#!/usr/bin/env python3
"""
Circular Mapper - 环状 DNA 精确比对工具
专为 plasmid-assembly-1.0.2 流程设计

功能：
1. 处理环状质粒参考序列（双拷贝策略）
2. 比对 long reads 到环状参考
3. 自动处理跨越连接点的 reads
4. 输出标准 BAM 格式

用法：
  # 单样本模式
  python circular_mapper.py -r ref.fasta -i reads.fastq -o output.bam
  
  # 批量处理 mode 0.2 结果
  python circular_mapper.py --batch /path/to/output_0.1.1/barcode_error_rate_0.2

作者: Senior Bioinformatics Engineer
日期: 2026-03-18
"""

import os
import sys
import argparse
import subprocess
import tempfile
import pathlib
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
import glob

try:
    from Bio import SeqIO
except ImportError:
    print("错误: 需要 Biopython。请安装: pip install biopython")
    sys.exit(1)

try:
    import pysam
except ImportError:
    print("错误: 需要 pysam。请安装: pip install pysam")
    sys.exit(1)


@dataclass
class CircularMappingResult:
    """比对结果统计"""
    sample_id: str
    total_reads: int
    mapped_reads: int
    circular_cross_reads: int
    mean_depth: float
    coverage_pct: float


class CircularMapper:
    """
    环状 DNA 比对器
    
    核心策略：双拷贝参考序列
    1. 将环状序列复制两份（头尾相接）
    2. 比对 reads 到双拷贝参考
    3. 坐标转换回原始参考
    4. 标记跨越连接点的 reads
    """
    
    def __init__(self, ref_fasta: str, threads: int = 4, 
                 min_identity: float = 0.85, preset: str = "map-ont",
                 enable_tr_realign: bool = False):
        """
        初始化环状比对器
        
        Args:
            ref_fasta: 参考序列 FASTA 文件
            threads: 比对线程数
            min_identity: 最小比对 identity
            preset: minimap2 preset (map-ont/map-pb/map-hifi)
            enable_tr_realign: 是否启用 TR 区域本地重比对 (v1.0.2.2)
        """
        self.ref_fasta = pathlib.Path(ref_fasta)
        self.threads = threads
        self.min_identity = min_identity
        self.preset = preset
        self.enable_tr_realign = enable_tr_realign
        
        # 加载参考序列信息
        self.ref_sequences = {}
        self.ref_lengths = {}
        self.doubled_fasta = None
        
        self._load_reference()
        # NOTE: rotation disabled in v1.0.2 - testing showed it moves SNP artifacts
        # to middle positions without reducing total count
        # self._rotate_reference()
        self._create_doubled_reference()
    
    def _load_reference(self):
        """加载原始参考序列"""
        if not self.ref_fasta.exists():
            print(f"[WARNING] 参考文件不存在: {self.ref_fasta}")
            return
        
        for record in SeqIO.parse(self.ref_fasta, "fasta"):
            self.ref_sequences[record.id] = str(record.seq)
            self.ref_lengths[record.id] = len(record.seq)
        
        if not self.ref_sequences:
            print(f"[WARNING] 参考文件为空或格式错误: {self.ref_fasta}")
            return
        
        print(f"[CircularMapper] 加载了 {len(self.ref_sequences)} 个参考序列")
        for ref_id, length in self.ref_lengths.items():
            print(f"  - {ref_id}: {length} bp")
    
    def _find_optimal_rotation(self, seq: str, window: int = 200, k: int = 15) -> int:
        """
        找到使 head/tail 重复度最低的旋转位点
        
        Args:
            seq: 原始序列
            window: 检查的 head/tail 长度
            k: k-mer 长度
            
        Returns:
            最优旋转偏移量 (0-based)
        """
        plen = len(seq)
        if plen <= window * 2:
            # 短序列不需要旋转
            return 0
        
        def rotate(s, offset):
            return s[offset:] + s[:offset]
        
        best_score = float('inf')
        best_offset = 0
        
        # 预计算所有 k-mers 的位置
        for offset in range(plen):
            rot = rotate(seq, offset)
            head = rot[:window]
            tail = rot[-window:]
            
            # 快速计算共享 k-mer 数
            head_kmers = set()
            for i in range(len(head) - k + 1):
                head_kmers.add(head[i:i+k])
            
            shared = 0
            for i in range(len(tail) - k + 1):
                if tail[i:i+k] in head_kmers:
                    shared += 1
            
            # 计算最长精确匹配 (简化版：只检查 tail 开头)
            max_match = 0
            for j in range(min(50, len(tail))):
                match_len = 0
                while j + match_len < len(tail) and tail[j + match_len] == head[match_len]:
                    match_len += 1
                max_match = max(max_match, match_len)
            
            score = shared + max_match * 5
            if score < best_score:
                best_score = score
                best_offset = offset
        
        return best_offset
    
    def _rotate_reference(self):
        """旋转参考序列到最优起始位点"""
        self.rotation_offsets = {}
        
        for ref_id, seq in list(self.ref_sequences.items()):
            offset = self._find_optimal_rotation(seq)
            self.rotation_offsets[ref_id] = offset
            
            if offset != 0:
                rotated_seq = seq[offset:] + seq[:offset]
                self.ref_sequences[ref_id] = rotated_seq
                print(f"[CircularMapper] 旋转参考序列 {ref_id}: offset={offset}, "
                      f"新起点=原位置{offset+1}")
            else:
                print(f"[CircularMapper] 参考序列 {ref_id} 无需旋转")
    
    def _create_doubled_reference(self) -> pathlib.Path:
        """创建双拷贝参考序列"""
        temp_fd, temp_path = tempfile.mkstemp(suffix=".doubled.fasta")
        
        with os.fdopen(temp_fd, 'w') as f:
            for ref_id, seq in self.ref_sequences.items():
                doubled_seq = seq + seq  # 双拷贝
                f.write(f">{ref_id}\n{doubled_seq}\n")
        
        self.doubled_fasta = pathlib.Path(temp_path)
        print(f"[CircularMapper] 创建双拷贝参考")
        return self.doubled_fasta
    
    def map_reads(self, fastq_path: str, output_bam: str, 
                  keep_temp: bool = False) -> CircularMappingResult:
        """
        比对 reads 到环状参考
        
        Args:
            fastq_path: 输入 FASTQ 文件
            output_bam: 输出 BAM 文件路径
            keep_temp: 是否保留临时文件
            
        Returns:
            CircularMappingResult 统计结果
        """
        fastq_path = pathlib.Path(fastq_path)
        output_bam = pathlib.Path(output_bam)
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        
        if not fastq_path.exists():
            print(f"[WARNING] FASTQ 文件不存在: {fastq_path}")
            return
        
        print(f"[CircularMapper] 开始比对...")
        print(f"  参考: {self.ref_fasta}")
        print(f"  输入: {fastq_path}")
        print(f"  输出: {output_bam}")
        
        # 步骤1: minimap2 比对
        temp_sam = output_bam.with_suffix('.temp.sam')
        self._run_minimap2(fastq_path, temp_sam)
        
        # 步骤2: 转换坐标并处理环化 reads
        temp_bam = output_bam.with_suffix('.temp.bam')
        stats = self._convert_and_filter_sam(temp_sam, temp_bam)
        
        # 步骤2.5: TR 区域本地重比对 (v1.0.2.2)
        # 参考 TRGT 算法：对 TR 区域的 reads 用打断重复的本地参考重新比对，
        # 消除 minimap2 在 tandem repeat 上的系统性比对错误。
        if self.enable_tr_realign:
            tr_bam = output_bam.with_suffix('.temp.tr.bam')
            self._tr_realign_postprocess(temp_bam, tr_bam)
            temp_bam = tr_bam
        
        # 步骤3: 排序和索引
        self._sort_and_index(temp_bam, output_bam)
        
        # 清理
        if not keep_temp:
            for f in [temp_sam, temp_bam, self.doubled_fasta]:
                if f and f.exists():
                    f.unlink()
        
        # 计算覆盖度
        coverage_stats = self._calculate_coverage(output_bam)
        
        print(f"[CircularMapper] 完成: {output_bam}")
        
        return CircularMappingResult(
            sample_id=output_bam.stem,
            total_reads=stats['total'],
            mapped_reads=stats['mapped'],
            circular_cross_reads=stats['circular_cross'],
            mean_depth=coverage_stats.get('mean_depth', 0),
            coverage_pct=coverage_stats.get('coverage_pct', 0)
        )
    
    def _run_minimap2(self, fastq_path: pathlib.Path, output_sam: pathlib.Path):
        """运行 minimap2 比对（只保留比对上的 reads）"""
        # minimap2 输出到临时 SAM，然后用 samtools 过滤未比对的 reads
        import tempfile
        
        cmd = [
            "minimap2",
            "-ax", self.preset,
            "-Y",                          # 保留 soft clip
            "--eqx",                       # 使用 X/= 代替 M
            "--secondary=no",              # 去除次级比对
            f"-t", str(self.threads),
            str(self.doubled_fasta),
            str(fastq_path)
        ]
        
        print(f"[CircularMapper] 运行 minimap2...")
        
        # 先用临时文件存储 minimap2 输出
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as tmp_sam:
            tmp_sam_path = tmp_sam.name
            result = subprocess.run(cmd, stdout=tmp_sam, stderr=subprocess.PIPE)
            if result.returncode != 0:
                stderr = result.stderr.decode() if result.stderr else "Unknown error"
                print(f"[WARNING] minimap2 失败: {stderr}")
                return
        
        # 使用 samtools 过滤未比对的 reads (-F 4 表示排除 flag 为 4 的 reads，即未比对的)
        filter_cmd = [
            "samtools", "view", "-h", "-F", "4", 
            "-o", str(output_sam),
            tmp_sam_path
        ]
        
        print(f"[CircularMapper] 过滤未比对 reads...")
        result = subprocess.run(filter_cmd, stderr=subprocess.PIPE)
        
        # 删除临时文件
        os.unlink(tmp_sam_path)
        
        if result.returncode != 0:
            stderr = result.stderr.decode() if result.stderr else "Unknown error"
            print(f"[WARNING] samtools 过滤失败: {stderr}")
            return
    
    def _split_cigar_at_position(self, cigar: List[Tuple[int, int]], 
                                   ref_pos: int, target_pos: int,
                                   read_seq: str, read_qual: List[int]) -> Tuple:
        """
        在指定参考位置分割 CIGAR
        
        Args:
            cigar: CIGAR 元组列表 [(op, length), ...]
            ref_pos: read 在参考序列上的起始位置 (0-based)
            target_pos: 分割位置的参考坐标 (0-based，即 ref_len)
            read_seq: read 序列
            read_qual: read 质量值列表
            
        Returns:
            (cigar1, cigar2, seq1, seq2, qual1, qual2, ref_len1, ref_len2)
        """
        cigar1 = []
        cigar2 = []
        seq1_parts = []
        seq2_parts = []
        qual1_parts = []
        qual2_parts = []
        
        current_ref_pos = ref_pos
        current_read_pos = 0
        split_found = False
        
        for op, length in cigar:
            # 定义哪些操作消耗参考序列，哪些消耗 read 序列
            consumes_ref = op in [0, 2, 3, 7, 8]   # M, D, N, =, X
            consumes_read = op in [0, 1, 4, 7, 8]  # M, I, S, =, X
            
            if split_found:
                # 已经找到分割点，剩余部分放入 cigar2
                cigar2.append((op, length))
                if consumes_read:
                    seq2_parts.append(read_seq[current_read_pos:current_read_pos + length])
                    qual2_parts.extend(read_qual[current_read_pos:current_read_pos + length])
                    current_read_pos += length
                continue
            
            if consumes_ref:
                # 此操作消耗参考序列
                op_end_ref_pos = current_ref_pos + length
                
                if op_end_ref_pos > target_pos:
                    # 此操作跨越分割点
                    split_offset = target_pos - current_ref_pos  # 在参考序列上的偏移量
                    
                    if consumes_read:
                        # 操作同时消耗 read 序列 (M, =, X)
                        # 按比例分割
                        first_part_len = split_offset
                        second_part_len = length - split_offset
                        
                        if first_part_len > 0:
                            cigar1.append((op, first_part_len))
                            seq1_parts.append(read_seq[current_read_pos:current_read_pos + first_part_len])
                            qual1_parts.extend(read_qual[current_read_pos:current_read_pos + first_part_len])
                        
                        if second_part_len > 0:
                            cigar2.append((op, second_part_len))
                            seq2_parts.append(read_seq[current_read_pos + first_part_len:current_read_pos + length])
                            qual2_parts.extend(read_qual[current_read_pos + first_part_len:current_read_pos + length])
                        
                        current_read_pos += length
                    else:
                        # 操作只消耗参考序列 (D, N)
                        first_part_len = split_offset
                        second_part_len = length - split_offset
                        
                        if first_part_len > 0:
                            cigar1.append((op, first_part_len))
                        if second_part_len > 0:
                            cigar2.append((op, second_part_len))
                    
                    split_found = True
                else:
                    # 此操作完全在分割点之前
                    cigar1.append((op, length))
                    if consumes_read:
                        seq1_parts.append(read_seq[current_read_pos:current_read_pos + length])
                        qual1_parts.extend(read_qual[current_read_pos:current_read_pos + length])
                        current_read_pos += length
                    
                    current_ref_pos = op_end_ref_pos
            else:
                # 此操作不消耗参考序列 (I, S, H, P)
                # 根据是否已分割决定放入哪一部分
                if consumes_read:
                    seq_part = read_seq[current_read_pos:current_read_pos + length]
                    qual_part = read_qual[current_read_pos:current_read_pos + length]
                    
                    if split_found:
                        # Soft clip 通常在开头或结尾
                        # 如果已经分割，检查位置决定
                        cigar2.append((op, length))
                        seq2_parts.append(seq_part)
                        qual2_parts.extend(qual_part)
                    else:
                        cigar1.append((op, length))
                        seq1_parts.append(seq_part)
                        qual1_parts.extend(qual_part)
                    
                    current_read_pos += length
        
        # 合并序列和质量
        seq1 = "".join(seq1_parts)
        seq2 = "".join(seq2_parts)
        qual1 = qual1_parts
        qual2 = qual2_parts
        
        # 验证 cigar1 和 seq1 长度
        cigar1_read_len = sum(l for o, l in cigar1 if o in [0, 1, 4, 7, 8])
        cigar2_read_len = sum(l for o, l in cigar2 if o in [0, 1, 4, 7, 8])
        
        if cigar1_read_len != len(seq1):
            print(f"[WARNING] CIGAR1 read length mismatch: cigar={cigar1_read_len}, seq={len(seq1)}")
        if cigar2_read_len != len(seq2):
            print(f"[WARNING] CIGAR2 read length mismatch: cigar={cigar2_read_len}, seq={len(seq2)}")
        
        # 计算参考序列长度
        ref_len1 = sum(l for o, l in cigar1 if o in [0, 2, 3, 7, 8])
        ref_len2 = sum(l for o, l in cigar2 if o in [0, 2, 3, 7, 8])
        
        return cigar1, cigar2, seq1, seq2, qual1, qual2, ref_len1, ref_len2

    def _convert_and_filter_sam(self, input_sam: pathlib.Path, output_bam: pathlib.Path) -> Dict:
        """转换 SAM 坐标并处理环化 reads"""
        print(f"[CircularMapper] 处理环化坐标...")
        
        # 创建 BAM header
        header = {
            'HD': {'VN': '1.6', 'SO': 'unsorted'},
            'SQ': [{'LN': length, 'SN': name} for name, length in self.ref_lengths.items()],
            'PG': [
                {'ID': 'CircularMapper', 'PN': 'CircularMapper', 'VN': '1.0.2.2'},
                {'ID': 'minimap2', 'PN': 'minimap2', 'VN': '2.30'}
            ]
        }
        
        infile = pysam.AlignmentFile(input_sam, "r")
        
        # =====================================================================
        # v1.0.2 修复：按 query_name 去重，解决 tandem repeat 导致的重复比对
        # =====================================================================
        from collections import defaultdict
        records_by_name = defaultdict(list)
        total_input = 0
        for read in infile:
            total_input += 1
            if read.is_unmapped:
                continue
            records_by_name[read.query_name].append(read)
        infile.close()
        
        # 对每个 query_name 选择最佳 record
        # 优先级: circular_cross > linear > second_copy
        # 对于 second_copy 的 supplementary alignment，恢复其 MAPQ 并去除 flag
        filtered_reads = []
        dedup_filtered = 0
        mapq_restored = 0
        for qname, records in records_by_name.items():
            # v1.0.2.1 bug fix: pysam AlignedSegment doesn't support dynamic attrs,
            # use per-query dict instead of method-level dict to avoid id() collisions
            restore_mapq_dict = {}
            if len(records) == 1:
                filtered_reads.append((records[0], restore_mapq_dict))
                continue
            
            scored = []
            for r in records:
                ref_name = r.reference_name
                if ref_name not in self.ref_lengths:
                    continue
                ref_len = self.ref_lengths[ref_name]
                orig_start = r.reference_start
                orig_end = r.reference_end if r.reference_end else orig_start
                crosses = (orig_start < ref_len) and (orig_end > ref_len)
                is_second = orig_start >= ref_len
                # 优先级: circular_cross(2) > linear(1) > second_copy(0)
                priority = 2 if crosses else (1 if not is_second else 0)
                scored.append((r, priority, is_second, r.query_alignment_length))
            
            if not scored:
                continue
            
            # 按优先级和长度排序
            scored.sort(key=lambda x: (x[1], x[3]), reverse=True)
            best_record, best_priority, best_is_second, _ = scored[0]
            
            # 如果选中的是 second_copy（通常是 supplementary，MAPQ=0），
            # 尝试从同组的 linear/circular_cross 记录中恢复 MAPQ
            if best_is_second:
                restore_mapq = None
                # 优先找同 query_name 的 linear (非 second_copy) 记录的 MAPQ
                for r, priority, is_second, _ in scored:
                    if not is_second and r.mapping_quality > 0:
                        restore_mapq = r.mapping_quality
                        break
                # 如果没找到，找任何非 supplementary 的 MAPQ
                if restore_mapq is None:
                    for r, priority, is_second, _ in scored:
                        if not r.is_supplementary and r.mapping_quality > 0:
                            restore_mapq = r.mapping_quality
                            break
                if restore_mapq is not None:
                    restore_mapq_dict[id(best_record)] = restore_mapq
                    mapq_restored += 1
            
            filtered_reads.append((best_record, restore_mapq_dict))
            dedup_filtered += len(records) - 1
        
        print(f"[CircularMapper] 去重: {total_input} records, {len(records_by_name)} reads, "
              f"过滤重复 {dedup_filtered}, MAPQ恢复 {mapq_restored}")
        
        outfile = pysam.AlignmentFile(output_bam, "wb", header=header)
        stats = {'total': total_input, 'mapped': 0, 'circular_cross': 0, 'filtered': dedup_filtered}
        
        for read_item in filtered_reads:
            if isinstance(read_item, tuple):
                read, restore_mapq_dict = read_item
            else:
                read = read_item
                restore_mapq_dict = {}
            # read is already mapped and filtered
            
            ref_name = read.reference_name
            if ref_name not in self.ref_lengths:
                continue
            
            ref_len = self.ref_lengths[ref_name]
            orig_start = read.reference_start
            orig_end = read.reference_end if read.reference_end else orig_start
            
            # 判断位置
            is_in_second_copy = orig_start >= ref_len
            crosses_origin = (orig_start < ref_len) and (orig_end > ref_len)
            
            # 检查 identity
            identity = self._calculate_identity(read)
            if identity < self.min_identity:
                stats['filtered'] += 1
                continue
            
            # =====================================================================
            # 恢复 MAPQ：双拷贝参考导致 minimap2 给大量 read MAPQ=0
            # （因为两个拷贝完全相同，次佳比对和最佳比对一样好）
            # =====================================================================
            mapq = read.mapping_quality
            if mapq == 0:
                # 优先使用去重阶段从同组记录恢复的 MAPQ
                restored = restore_mapq_dict.get(id(read))
                if restored is not None and restored > 0:
                    mapq = restored
                else:
                    # 单记录 read 或没有可用恢复值：
                    # 根据 alignment identity 估算合理的 MAPQ
                    # （双拷贝参考造成的人为 MAPQ=0 需要恢复）
                    if identity >= 0.99:
                        mapq = 60
                    elif identity >= 0.95:
                        mapq = 50
                    elif identity >= 0.90:
                        mapq = 40
                    elif identity >= 0.85:
                        mapq = 30
            
            # =====================================================================
            # 清除 supplementary flag (0x800)：双拷贝参考导致的 supplementary
            # 是人为的（同一个 read 在两个拷贝都有完全相同的匹配），不是
            # 真实的 chimera。包括 circular_cross 的分割记录也不应保留此 flag。
            # =====================================================================
            flag = read.flag & ~0x800
            
            # 转换坐标
            if is_in_second_copy:
                # 完全在第二拷贝，转换回第一拷贝
                new_read = pysam.AlignedSegment(outfile.header)
                new_read.query_name = read.query_name
                new_read.query_sequence = read.query_sequence
                new_read.query_qualities = read.query_qualities
                new_read.flag = flag
                new_read.reference_name = ref_name
                new_read.mapping_quality = mapq
                new_read.cigar = read.cigar
                new_read.reference_start = orig_start - ref_len
                new_read.tags = read.tags.copy() if read.tags else []
                new_read.set_tag('XC', 'second_copy', 'Z')
                new_read.set_tag('XI', round(identity, 4), 'f')
                outfile.write(new_read)
                stats['mapped'] += 1
                
            elif crosses_origin:
                # 跨越连接点！需要分割成两条 reads
                stats['circular_cross'] += 1
                
                # 获取 read 序列和质量
                read_seq = read.query_sequence
                read_qual = list(read.query_qualities) if read.query_qualities else []
                
                if not read_seq:
                    continue
                
                # 在 ref_len 处分割
                cigar1, cigar2, seq1, seq2, qual1, qual2, ref_len1, ref_len2 = \
                    self._split_cigar_at_position(
                        read.cigar, orig_start, ref_len, read_seq, read_qual
                    )
                
                # 创建第一条 read (从 orig_start 到 ref_len)
                if cigar1 and seq1:
                    read1 = pysam.AlignedSegment(outfile.header)
                    read1.query_name = f"{read.query_name}/1"
                    read1.query_sequence = seq1
                    read1.query_qualities = qual1 if qual1 else None
                    read1.flag = flag
                    read1.reference_name = ref_name
                    read1.reference_start = orig_start
                    read1.mapping_quality = mapq
                    read1.cigar = cigar1
                    read1.tags = read.tags.copy() if read.tags else []
                    read1.set_tag('XC', 'circular_cross_part1', 'Z')
                    read1.set_tag('XI', round(identity, 4), 'f')
                    outfile.write(read1)
                
                # 创建第二条 read (从 1 到环绕部分)
                if cigar2 and seq2:
                    read2 = pysam.AlignedSegment(outfile.header)
                    read2.query_name = f"{read.query_name}/2"
                    read2.query_sequence = seq2
                    read2.query_qualities = qual2 if qual2 else None
                    read2.flag = flag
                    read2.reference_name = ref_name
                    read2.reference_start = 0  # 从序列起始位置开始
                    read2.mapping_quality = mapq
                    read2.cigar = cigar2
                    read2.tags = read.tags.copy() if read.tags else []
                    read2.set_tag('XC', 'circular_cross_part2', 'Z')
                    read2.set_tag('XI', round(identity, 4), 'f')
                    outfile.write(read2)
                
                stats['mapped'] += 1
            else:
                # 正常比对到第一拷贝
                new_read = pysam.AlignedSegment(outfile.header)
                new_read.query_name = read.query_name
                new_read.query_sequence = read.query_sequence
                new_read.query_qualities = read.query_qualities
                new_read.flag = flag
                new_read.reference_name = ref_name
                new_read.mapping_quality = mapq
                new_read.cigar = read.cigar
                new_read.reference_start = orig_start
                new_read.tags = read.tags.copy() if read.tags else []
                new_read.set_tag('XC', 'linear', 'Z')
                new_read.set_tag('XI', round(identity, 4), 'f')
                outfile.write(new_read)
                stats['mapped'] += 1
        
        outfile.close()
        
        print(f"[CircularMapper] 统计: 总={stats['total']}, "
              f"比对={stats['mapped']}, 环化={stats['circular_cross']}, "
              f"过滤={stats['filtered']}")
        
        return stats
    
    def _calculate_identity(self, read) -> float:
        """计算比对 identity"""
        if not read.cigar:
            return 0.0
        
        match_len = 0
        total_len = 0
        
        for op, length in read.cigar:
            if op == 7:  # =
                match_len += length
                total_len += length
            elif op == 8:  # X
                total_len += length
            elif op == 0:  # M
                total_len += length
                match_len += length * 0.5  # 保守估计
            elif op in [1, 2]:  # I, D
                total_len += length
        
        return match_len / total_len if total_len > 0 else 0.0
    
    def _sort_and_index(self, input_bam: pathlib.Path, output_bam: pathlib.Path):
        """排序和索引 BAM"""
        # 排序
        cmd_sort = [
            "samtools", "sort",
            "-@", str(self.threads),
            "-o", str(output_bam),
            str(input_bam)
        ]
        subprocess.run(cmd_sort, check=True, capture_output=True)
        
        # 索引
        cmd_index = ["samtools", "index", str(output_bam)]
        subprocess.run(cmd_index, check=True, capture_output=True)
    
    def _tr_realign_postprocess(self, input_bam: pathlib.Path, output_bam: pathlib.Path):
        """
        TR 区域本地重比对后处理 (v1.0.2.2)。
        
        参考 TRGT 算法：对 tandem repeat 区域的 reads 提取后，用打断重复的
        本地参考重新比对，消除 minimap2 在 TR 上的系统性比对错误。
        
        核心策略：
        1. 构建 TR 本地参考 = head_chunk + NNNNNNNNNN + tail_chunk
        2. 提取覆盖 TR 区域的 reads
        3. minimap2 本地重比对
        4. 只保留重比对后 CIGAR 更简单的 reads，映射回原始坐标
        """
        TR_WINDOW = 250
        SPACER_LEN = 10
        
        ref_name = list(self.ref_sequences.keys())[0]
        ref_seq = self.ref_sequences[ref_name]
        ref_len = len(ref_seq)
        
        # 1. 构建 TR 本地参考
        head_chunk = ref_seq[:TR_WINDOW]
        tail_chunk = ref_seq[-TR_WINDOW:]
        tr_ref_seq = head_chunk + 'N' * SPACER_LEN + tail_chunk
        
        tr_ref_fasta = output_bam.with_suffix('.tr_ref.fasta')
        with open(tr_ref_fasta, 'w') as f:
            f.write(f'>TR_local\n{tr_ref_seq}\n')
        
        # 2. 提取覆盖 TR 区域的 reads
        in_bam = pysam.AlignmentFile(input_bam, 'rb')
        tr_reads = {}
        for read in in_bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            start = read.reference_start
            end = read.reference_end or start
            in_head = (start < TR_WINDOW and end > 0)
            in_tail = (start < ref_len and end > ref_len - TR_WINDOW)
            if in_head or in_tail:
                tr_reads[read.query_name] = read
        in_bam.close()
        
        if not tr_reads:
            import shutil
            shutil.copy(input_bam, output_bam)
            if tr_ref_fasta.exists():
                tr_ref_fasta.unlink()
            return
        
        # 3. 保存提取的 reads 到临时 fastq
        tr_fastq = output_bam.with_suffix('.tr_reads.fastq')
        with open(tr_fastq, 'w') as fout:
            for name, read in tr_reads.items():
                seq = read.query_sequence
                qual = read.query_qualities
                if not seq:
                    continue
                qual_str = ''.join(chr(min(q + 33, 126)) for q in qual) if qual else 'I' * len(seq)
                fout.write(f'@{name}\n{seq}\n+\n{qual_str}\n')
        
        # 4. 本地重比对
        tr_sam = output_bam.with_suffix('.tr_realigned.sam')
        cmd = [
            "minimap2", "-ax", "map-ont", "-Y", "--eqx", "--secondary=no",
            "-k", "15",
            "-t", str(min(4, self.threads)),
            str(tr_ref_fasta), str(tr_fastq)
        ]
        with open(tr_sam, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
            if result.returncode != 0:
                stderr = result.stderr.decode() if result.stderr else "Unknown error"
                print(f"[WARNING] TR realignment minimap2 failed: {stderr}")
        
        # 5. 解析重比对结果
        realigned = {}
        if tr_sam.exists():
            with pysam.AlignmentFile(tr_sam, 'r') as tr_bam:
                for read in tr_bam.fetch():
                    if read.is_unmapped or not read.cigar:
                        continue
                    
                    local_start = read.reference_start
                    local_end = read.reference_end or local_start
                    head_end = TR_WINDOW
                    tail_start = TR_WINDOW + SPACER_LEN
                    
                    # 映射回原始坐标
                    if local_end <= head_end:
                        orig_start = local_start
                        new_region = 'head'
                    elif local_start >= tail_start:
                        orig_start = ref_len - TR_WINDOW + (local_start - tail_start)
                        new_region = 'tail'
                    else:
                        continue  # 跨越 spacer，过滤
                    
                    cigar_ops = len(read.cigar)
                    indel_ops = sum(1 for op, _ in read.cigar if op in [1, 2])
                    
                    realigned[read.query_name] = {
                        'orig_start': orig_start,
                        'cigar': read.cigar,
                        'mapq': read.mapping_quality,
                        'cigar_ops': cigar_ops,
                        'indel_ops': indel_ops,
                        'region': new_region,
                    }
        
        # 6. 更新 BAM：只替换重比对后 CIGAR 更简单的 reads
        n_simplified = 0
        with pysam.AlignmentFile(input_bam, 'rb') as in_bam, \
             pysam.AlignmentFile(output_bam, 'wb', header=in_bam.header) as out_bam:
            for read in in_bam.fetch(until_eof=True):
                if read.query_name in realigned and read.query_name in tr_reads:
                    info = realigned[read.query_name]
                    orig = tr_reads[read.query_name]
                    
                    # 区域一致性检查：只替换映射到相同区域的 reads
                    orig_start = orig.reference_start
                    orig_end = orig.reference_end or orig_start
                    orig_in_head = (orig_start < TR_WINDOW)
                    orig_in_tail = (orig_end > ref_len - TR_WINDOW)
                    if orig_in_head and not orig_in_tail:
                        orig_region = 'head'
                    elif orig_in_tail and not orig_in_head:
                        orig_region = 'tail'
                    else:
                        orig_region = 'both'
                    
                    # 如果重比对后区域不一致，保留原比对
                    if info['region'] != orig_region and orig_region != 'both':
                        out_bam.write(read)
                        continue
                    
                    orig_ops = len(orig.cigar) if orig.cigar else 999
                    orig_indel = sum(1 for op, _ in orig.cigar if op in [1, 2]) if orig.cigar else 999
                    
                    # 只替换 CIGAR 更简单的
                    if info['cigar_ops'] < orig_ops or info['indel_ops'] < orig_indel:
                        n_simplified += 1
                        
                        new_read = pysam.AlignedSegment(out_bam.header)
                        new_read.query_name = read.query_name
                        new_read.query_sequence = read.query_sequence
                        new_read.query_qualities = read.query_qualities
                        new_read.flag = read.flag & ~0x800
                        new_read.reference_name = read.reference_name
                        new_read.reference_start = info['orig_start']
                        new_read.mapping_quality = info['mapq']
                        new_read.cigar = info['cigar']
                        new_read.tags = read.tags.copy() if read.tags else []
                        new_read.set_tag('TR', 1, 'i')
                        out_bam.write(new_read)
                        continue
                
                out_bam.write(read)
        
        print(f"[CircularMapper] TR重比对: {len(tr_reads)} reads 提取, "
              f"{len(realigned)} 成功重比对, {n_simplified} 简化替换")
        
        # 清理临时文件
        for f in [tr_ref_fasta, tr_fastq, tr_sam]:
            if f.exists():
                f.unlink()
    
    def _calculate_coverage(self, bam_file: pathlib.Path) -> Dict:
        """计算覆盖度统计"""
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
            
            total_length = sum(self.ref_lengths.values())
            covered_bases = 0
            total_depth = 0
            
            for ref_name in bam.references:
                ref_len = bam.get_reference_length(ref_name)
                for pileup in bam.pileup(ref_name):
                    if pileup.n > 0:
                        covered_bases += 1
                        total_depth += pileup.n
            
            bam.close()
            
            mean_depth = total_depth / covered_bases if covered_bases > 0 else 0
            coverage_pct = (covered_bases / total_length) * 100 if total_length > 0 else 0
            
            return {'mean_depth': mean_depth, 'coverage_pct': coverage_pct}
        except:
            return {'mean_depth': 0, 'coverage_pct': 0}


def process_single_sample(ref_fasta: str, fastq_path: str, output_bam: str, 
                          threads: int = 4, min_identity: float = 0.85):
    """处理单个样本"""
    mapper = CircularMapper(
        ref_fasta=ref_fasta,
        threads=threads,
        min_identity=min_identity,
        preset="map-ont"
    )
    
    result = mapper.map_reads(fastq_path, output_bam)
    return result


def batch_process(base_dir: str, error_rate: str = "0.2", threads: int = 4):
    """
    批量处理所有样本
    
    期望目录结构:
    base_dir/barcode_error_rate_{er}/{sample}/{sample}_plas_asm/{sample}_plassembler_plasmids.fasta
    base_dir/barcode_error_rate_{er}/{sample}/{sample}_plas_asm/plasmid_fastqs/plasmid_long.fastq
    """
    base_path = Path(base_dir)
    error_rate_dir = base_path / f"barcode_error_rate_{error_rate}"
    
    if not error_rate_dir.exists():
        print(f"错误: 目录不存在: {error_rate_dir}")
        return
    
    # 查找所有样本
    sample_dirs = [d for d in error_rate_dir.iterdir() if d.is_dir() and "Single" in d.name]
    sample_dirs = sorted(sample_dirs, key=lambda x: x.name)
    
    print(f"=" * 70)
    print(f"批量处理: {len(sample_dirs)} 个样本")
    print(f"Error Rate: {error_rate}")
    print(f="=" * 70)
    print()
    
    results = []
    
    for i, sample_dir in enumerate(sample_dirs, 1):
        sample_name = sample_dir.name
        
        # 自动查找输入文件
        ref_fasta = sample_dir / f"{sample_name}_plas_asm" / f"{sample_name}_plassembler_plasmids.fasta"
        fastq_path = sample_dir / f"{sample_name}_plas_asm" / "plasmid_fastqs" / "plasmids_long.fastq"
        output_bam = sample_dir / f"{sample_name}_circular_mapped.bam"
        
        print(f"[{i}/{len(sample_dirs)}] 处理: {sample_name}")
        
        # 检查输入文件
        if not ref_fasta.exists():
            print(f"  ⚠️  跳过: 参考文件不存在: {ref_fasta}")
            continue
        
        if not fastq_path.exists():
            print(f"  ⚠️  跳过: FASTQ 不存在: {fastq_path}")
            continue
        
        # 检查是否已处理
        if output_bam.exists():
            print(f"  ⏭️  跳过: 已存在: {output_bam}")
            continue
        
        try:
            result = process_single_sample(
                ref_fasta=str(ref_fasta),
                fastq_path=str(fastq_path),
                output_bam=str(output_bam),
                threads=threads
            )
            results.append(result)
            print(f"  ✅ 完成: {result.mapped_reads} reads 比对, "
                  f"{result.circular_cross_reads} 跨越连接点")
        except Exception as e:
            print(f"  ❌ 错误: {e}")
    
    # 输出汇总
    print()
    print("=" * 70)
    print("批量处理完成")
    print("=" * 70)
    print(f"成功处理: {len(results)}/{len(sample_dirs)} 样本")
    
    if results:
        total_reads = sum(r.total_reads for r in results)
        total_mapped = sum(r.mapped_reads for r in results)
        total_circular = sum(r.circular_cross_reads for r in results)
        
        print(f"总 reads: {total_reads}")
        print(f"总比对: {total_mapped} ({total_mapped/total_reads*100:.1f}%)")
        print(f"环化 reads: {total_circular}")


def main():
    parser = argparse.ArgumentParser(
        description="Circular Mapper - 环状质粒 reads 比对工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 单样本模式
  python circular_mapper.py -r ref.fasta -i reads.fastq -o output.bam
  
  # 批量处理 (自动查找 plassembler 输出)
  python circular_mapper.py --batch /path/to/output_0.1.1/barcode_error_rate_0.2
  
  # 指定线程数和最小 identity
  python circular_mapper.py -r ref.fasta -i reads.fastq -o output.bam -t 8 --min-identity 0.90
        """
    )
    
    # 单样本模式参数
    parser.add_argument('-r', '--reference', 
                        help='参考序列 FASTA 文件 (单样本模式)')
    parser.add_argument('-i', '--input', 
                        help='输入 FASTQ 文件 (单样本模式)')
    parser.add_argument('-o', '--output', 
                        help='输出 BAM 文件 (单样本模式)')
    
    # 批量模式参数
    parser.add_argument('--batch', 
                        help='批量处理目录 (包含 barcode_error_rate_X.XX 子目录)')
    parser.add_argument('--error-rate', default='0.2',
                        help='批量处理时的 error rate (默认: 0.2)')
    
    # 通用参数
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='线程数 (默认: 4)')
    parser.add_argument('--min-identity', type=float, default=0.85,
                        help='最小比对 identity (默认: 0.85)')
    
    args = parser.parse_args()
    
    # 检查模式
    if args.batch:
        # 批量模式
        batch_process(args.batch, args.error_rate, args.threads)
    elif args.reference and args.input and args.output:
        # 单样本模式
        result = process_single_sample(
            ref_fasta=args.reference,
            fastq_path=args.input,
            output_bam=args.output,
            threads=args.threads,
            min_identity=args.min_identity
        )
        
        print()
        print("=" * 60)
        print("比对结果")
        print("=" * 60)
        print(f"样本: {result.sample_id}")
        print(f"总 reads: {result.total_reads}")
        print(f"比对 reads: {result.mapped_reads}")
        print(f"环化 reads: {result.circular_cross_reads}")
        print(f"平均覆盖度: {result.mean_depth:.2f}X")
        print(f"覆盖度: {result.coverage_pct:.1f}%")
    else:
        parser.print_help()
        print("\n错误: 请指定 --batch 进行批量处理，或指定 -r/-i/-o 进行单样本处理")
        sys.exit(1)


if __name__ == "__main__":
    main()
