use std::{
    cmp,
    collections::HashMap,
    fmt::Display,
    fs,
    io::{BufWriter, Write},
};

use gskits::{
    file_reader::{bed_reader::BedInfo, vcf_reader::VcfInfo},
    gsbam::bam_record_ext::BamRecordExt,
    pbar,
};
use indicatif::ProgressBar;
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};

use crate::{cli::AlignedBamParams, set_spin_pb};

use super::FastaData;

fn base_cnt_map_2_str(base_cnt: &HashMap<u8, usize>) -> String {
    let mut items = base_cnt
        .iter()
        .map(|(base, cnt)| (*base, *cnt))
        .collect::<Vec<_>>();
    items.sort_by_key(|v| v.0);

    let items = items
        .into_iter()
        .map(|(base, cnt)| format!("{}:{}", base as char, cnt))
        .collect::<Vec<_>>();
    items.join(",")
}

struct LocusStat {
    pos: usize,
    eq: usize,
    diff: HashMap<u8, usize>,
    ins: HashMap<u8, usize>,
    del: usize,
    depth: usize, // how much records that aligned to this position
    cur_base: String,
    next_base: String,
    cur_is_homo: i32,
    next_is_homo: i32,
    around_bases: String,
}

impl LocusStat {
    fn new(
        pos: usize,
        cur_base: String,
        next_base: String,
        around_bases: String,
        cur_is_homo: i32,
        next_is_homo: i32,
    ) -> Self {
        Self {
            pos,
            eq: 0,
            diff: HashMap::new(),
            ins: HashMap::new(),
            del: 0,
            depth: 0,
            cur_base,
            next_base,
            cur_is_homo,
            next_is_homo,
            around_bases,
        }
    }

    fn csv_header() -> &'static str {
        "pos\teq\tdiff\tins\tdel\tdepth\tcurBase\tnextBase\tcurIsHomo\tnextIsHomo\taroundBases\tdiffDetail\tinsDetail"
    }

    fn diff_str(&self) -> String {
        base_cnt_map_2_str(&self.diff)
    }

    fn ins_str(&self) -> String {
        base_cnt_map_2_str(&self.ins)
    }

    fn diff_tot(&self) -> usize {
        self.diff
            .values()
            .copied()
            .reduce(|acc, v| acc + v)
            .unwrap_or(0)
    }

    fn ins_tot(&self) -> usize {
        self.ins
            .values()
            .copied()
            .reduce(|acc, v| acc + v)
            .unwrap_or(0)
    }
}

impl Display for LocusStat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}",
            self.pos,
            self.eq,
            self.diff_tot(),
            self.ins_tot(),
            self.del,
            self.depth,
            self.cur_base,
            self.next_base,
            self.cur_is_homo,
            self.next_is_homo,
            self.around_bases,
            self.diff_str(),
            self.ins_str()
        )
    }
}

pub fn fact_ref_locus_info(
    args: &AlignedBamParams,
    output_dir: &str,
    hc_regions: Option<&BedInfo>,
    hc_variants: Option<&VcfInfo>,
    fasta_data: &FastaData,
    pbar: ProgressBar,
) {
    let bam_file = &args.bam;

    let mut bam_h = bam::IndexedReader::from_path(bam_file).unwrap();
    bam_h.set_threads(4).unwrap();

    let o_filepath = format!("{}/fact_aligned_bam_ref_locus_info.csv", output_dir);
    let o_file = fs::File::create(o_filepath).unwrap();
    let mut o_file_buff_writer = BufWriter::new(o_file);
    writeln!(
        &mut o_file_buff_writer,
        "refname\t{}",
        LocusStat::csv_header()
    )
    .unwrap();

    let pb = set_spin_pb(pbar, format!("fact_ref_locus_info"), pbar::DEFAULT_INTERVAL);

    for (refname, refseq) in fasta_data.get_ref_name2seq() {
        bam_h
            .fetch((refname, 0, refseq.len() as u64))
            .expect(&format!("fetch {} error", refname));
        let ref_seq_len = refseq.len();
        let ref_seq_bytes = refseq.as_bytes();

        let mut ref_locus_stat = (0..refseq.len())
            .into_iter()
            .map(|pos| {
                let around_start = if pos < 10 { 0 } else { pos - 10 };
                let around_end = cmp::min(pos + 11, ref_seq_len);

                let cur_is_homo = is_homo_pos(ref_seq_bytes, pos);

                let next_is_homo = if pos + 1 >= ref_seq_len {
                    false
                } else {
                    is_homo_pos(ref_seq_bytes, pos + 1)
                };

                LocusStat::new(
                    pos,
                    refseq[pos..pos + 1].to_string(),
                    if pos + 1 >= ref_seq_len {
                        "".to_string()
                    } else {
                        refseq[pos + 1..pos + 2].to_string()
                    },
                    format!(
                        "{}[{}]{}",
                        &refseq[around_start..pos],
                        &refseq[pos..pos + 1],
                        &refseq[pos + 1..around_end]
                    ),
                    if cur_is_homo { 1 } else { 0 },
                    if next_is_homo { 1 } else { 0 },
                )
            })
            .collect::<Vec<_>>();

        let refseq = refseq.as_bytes();

        for record in bam_h.records() {
            pb.inc(1);
            let record = record.unwrap();
            if record.is_secondary() || record.is_unmapped() || record.is_supplementary() {
                continue;
            }

            let record_ext = BamRecordExt::new(&record);

            let ref_start = record_ext.reference_start() as i64;
            let ref_end = record_ext.reference_end() as i64;
            let query_end = record_ext.query_alignment_end() as i64;
            // let ref_start = record.reference_start();
            // let ref_end = record.reference_end();

            ref_locus_stat
                .iter_mut()
                .skip(ref_start as usize)
                .take((ref_end - ref_start) as usize)
                .for_each(|stat| stat.depth += 1);

            let mut rpos_cursor = None;
            let mut qpos_cursor = None;

            let query_seq = record.seq().as_bytes();

            for [qpos, rpos] in record.aligned_pairs_full() {
                if qpos.is_some() {
                    qpos_cursor = qpos;
                }
                if rpos.is_some() {
                    rpos_cursor = rpos;
                }

                if let Some(rpos_cursor_) = rpos_cursor {
                    if rpos_cursor_ < ref_start {
                        continue;
                    }
                    if rpos_cursor_ >= ref_end {
                        break;
                    }
                } else {
                    continue;
                }

                if let Some(qpos_cursor_) = qpos_cursor {
                    if qpos_cursor_ >= query_end {
                        break;
                    }
                } else {
                    continue;
                }

                let rpos_cur_or_pre = rpos_cursor.unwrap() as usize;
                if let Some(hc_regions_) = hc_regions {
                    if !hc_regions_.point_within_region(refname, rpos_cur_or_pre) {
                        continue;
                    }
                }

                if let Some(vcf_regions_) = hc_variants {
                    if vcf_regions_.point_hit(refname, rpos_cur_or_pre) {
                        continue;
                    }
                }

                if qpos.is_none() {
                    // deletion
                    unsafe {
                        ref_locus_stat.get_unchecked_mut(rpos_cur_or_pre).del += 1;
                    }
                    continue;
                }

                if rpos.is_none() {
                    // insertion
                    unsafe {
                        let ins_base = *query_seq.get_unchecked(qpos.unwrap() as usize);
                        *ref_locus_stat
                            .get_unchecked_mut(rpos_cur_or_pre)
                            .ins
                            .entry(ins_base)
                            .or_insert(0) += 1;
                    }
                    continue;
                }

                unsafe {
                    if *refseq.get_unchecked(rpos.unwrap() as usize)
                        == *query_seq.get_unchecked(qpos.unwrap() as usize)
                    {
                        ref_locus_stat.get_unchecked_mut(rpos_cur_or_pre).eq += 1;
                    } else {
                        let diff_base = *query_seq.get_unchecked(qpos.unwrap() as usize);
                        *ref_locus_stat
                            .get_unchecked_mut(rpos_cur_or_pre)
                            .diff
                            .entry(diff_base)
                            .or_insert(0) += 1;
                    }
                }
            }
        }

        ref_locus_stat.iter().for_each(|locus_stat| {
            writeln!(&mut o_file_buff_writer, "{}\t{}", refname, locus_stat).unwrap()
        });
    }
    pb.finish();
}

fn is_homo_pos(bases: &[u8], pos: usize) -> bool {
    let mut res = false;
    let cur_base = bases[pos];
    if pos > 0 {
        res |= cur_base == bases[pos - 1];
    }

    if (pos + 1) < bases.len() {
        res |= cur_base == bases[pos + 1];
    }

    res
}
