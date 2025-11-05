import pysam
import sys
from tqdm import tqdm

def find_high_qual_reads(bam_file, min_base_qual=45, max_reads=100):
    """
    查找BAM文件中存在base quality > min_base_qual的reads
    
    参数:
    bam_file: BAM文件路径
    min_base_qual: 最小base quality阈值
    max_reads: 最大输出reads数量
    """
    
    # 打开BAM文件
    bam = pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=40)
    
    count = 0
    found_reads = 0
    
    print(f"搜索BAM文件中存在base quality > {min_base_qual}的reads...")
    print("=" * 80)
    
    try:
        for read in tqdm(bam.fetch(until_eof=True), desc=f"reading {bam_file}"):
            
            if read.get_tag("np") < 3:
                continue
            
            count += 1
            
            # 检查read是否有quality信息
            if read.query_qualities is None:
                continue
            
            # 检查是否有任何一个base的quality > 阈值
            has_high_qual = any(qual > min_base_qual for qual in read.query_qualities)
            
            if has_high_qual:
                found_reads += 1
                print(f"\nRead #{found_reads}:")
                print(f"Read Name: {read.query_name}")
                print(f"Flags: {read.flag} (0x{read.flag:04x})")
                
                # 打印quality值统计
                quals = read.query_qualities
                max_qual = max(quals) if quals else 0
                min_qual = min(quals) if quals else 0
                avg_qual = sum(quals) / len(quals) if quals else 0
                
                print(f"Quality Stats - Max: {max_qual}, Min: {min_qual}, Avg: {avg_qual:.2f}")
                
                # 打印前20个base的quality值（如果read很长）
                print("First 20 base qualities:")
                qual_str = " ".join(f"{q:2d}" for q in quals[:20])
                print(qual_str, end=" ")
                
                if len(quals) > 20:
                    print(f"... and {len(quals) - 20} more bases")
                    
                print(read)
                
                print("-" * 60)
                
                # 如果找到的reads数量达到限制，停止搜索
                if found_reads >= max_reads:
                    print(f"\n已达到最大显示数量 ({max_reads})，停止搜索。")
                    break
            
            # 每处理10000个reads打印进度
            if count % 10000 == 0:
                print(f"已处理 {count} 个reads, 找到 {found_reads} 个符合条件的reads...")
    
    except KeyboardInterrupt:
        print("\n用户中断搜索")
    
    finally:
        bam.close()
    
    print(f"\n统计信息:")
    print(f"总共处理的reads: {count}")
    print(f"找到的符合条件的reads: {found_reads}")
    print(f"比例: {found_reads/count*100:.2f}%")

def detailed_quality_analysis(bam_file, min_base_qual=45, sample_size=50):
    """
    对符合条件的reads进行更详细的质量分析
    """
    bam = pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=40)
    
    high_qual_reads = []
    count = 0
    
    print(f"\n进行详细质量分析 (采样 {sample_size} 个reads)...")
    print("=" * 80)
    
    for read in tqdm(bam.fetch(until_eof=True), desc=f"reading {bam_file}"):
        
        if read.get_tag("np") < 3:
                continue
        
        if read.query_qualities is None:
            continue
        
        if any(qual > min_base_qual for qual in read.query_qualities):
            high_qual_reads.append(read)
            count += 1
            
            if count >= sample_size:
                break
    
    bam.close()
    
    if not high_qual_reads:
        print("未找到符合条件的reads")
        return
    
    # 分析质量分布
    all_high_quals = []
    for read in high_qual_reads:
        high_quals = [q for q in read.query_qualities if q > min_base_qual]
        all_high_quals.extend(high_quals)
    
    if all_high_quals:
        print(f"高质量base (> {min_base_qual}) 的统计:")
        print(f"  总数: {len(all_high_quals)}")
        print(f"  平均质量: {sum(all_high_quals)/len(all_high_quals):.2f}")
        print(f"  最大质量: {max(all_high_quals)}")
        print(f"  最小质量: {min(all_high_quals)}")
        
        # 质量值分布
        qual_bins = {}
        for qual in all_high_quals:
            bin_key = f"{qual}"
            qual_bins[bin_key] = qual_bins.get(bin_key, 0) + 1
        
        print(f"\n质量值分布:")
        for qual in sorted(qual_bins.keys(), key=int, reverse=True):
            print(f"  Q{qual}: {qual_bins[qual]}个base")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python script.py <bam_file> [min_base_qual] [max_reads]")
        print("示例: python script.py sample.bam 45 50")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    min_base_qual = int(sys.argv[2]) if len(sys.argv) > 2 else 45
    max_reads = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    
    # 检查BAM文件是否存在
    try:
        # 首先检查BAM索引文件
        pysam.index(bam_file)
    except:
        print(f"注意: BAM文件可能需要索引, 请确保 {bam_file}.bai 文件存在")
    
    # 查找高质量reads
    find_high_qual_reads(bam_file, min_base_qual, max_reads)
    
    # 进行详细分析
    detailed_quality_analysis(bam_file, min_base_qual)