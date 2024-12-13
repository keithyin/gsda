use std::{
    collections::HashMap,
    fmt::Display,
    fs,
    io::{BufWriter, Write},
};

use gskits::{
    file_reader::{bed_reader::BedInfo, vcf_reader::VcfInfo}, gsbam::bam_record_ext::BamRecordExt, pbar
};
use indicatif::ProgressBar;
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};

use crate::{cli::AlignedBamParams, set_spin_pb};

use super::FastaData;

struct BaseQStat {
    qvalue: u8,
    eq: usize,
    diff: usize,
    ins: usize,
    del: usize,
    depth: usize,
}

impl BaseQStat {
    fn new(qvalue: u8) -> Self {
        Self {
            qvalue,
            eq: 0,
            diff: 0,
            ins: 0,
            del: 0,
            depth: 0,
        }
    }

    fn csv_header() -> &'static str {
        "baseq\teq\tdiff\tins\tdel\tdepth"
    }
}

impl Display for BaseQStat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
            self.qvalue, self.eq, self.diff, self.ins, self.del, self.depth
        )
    }
}

/// qname refname qpos_start qpos_end rpos_start, rpos_end query ref
pub fn fact_baseq_stat(
    args: &AlignedBamParams,
    output_dir: &str,
    hc_regions: Option<&BedInfo>,
    hc_variants: Option<&VcfInfo>,
    fasta_data: &FastaData,
    pbar: ProgressBar
) {
    let bam_file = &args.bam;

    let mut bam_h = bam::IndexedReader::from_path(bam_file).unwrap();
    bam_h.set_threads(4).unwrap();

    let o_filepath = format!("{}/fact_baseq_stat.csv", output_dir);
    let o_file = fs::File::create(o_filepath).unwrap();
    let mut o_file_buff_writer = BufWriter::new(o_file);
    writeln!(
        &mut o_file_buff_writer,
        "refname\t{}",
        BaseQStat::csv_header()
    )
    .unwrap();

    let pb = set_spin_pb(pbar, format!("fact_baseq_stat"), pbar::DEFAULT_INTERVAL);

    for (refname, refseq) in fasta_data.get_ref_name2seq() {
        bam_h
            .fetch((refname, 0, refseq.len() as u64))
            .expect(&format!("fetch {} error", refname));

        let mut baseq2stat = HashMap::new();

        let refseq_bytes = refseq.as_bytes();

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

            let mut rpos_cursor = None;
            let mut qpos_cursor = None;

            let query_seq = record.seq().as_bytes();
            let qual = record.qual();

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

                if let Some(hc_regions_) = hc_regions {
                    if !hc_regions_.point_within_region(refname, rpos_cursor.unwrap() as usize) {
                        continue;
                    }
                }

                if let Some(hc_variants_) = hc_variants {
                    if hc_variants_.point_hit(refname, rpos_cursor.unwrap() as usize) {
                        continue;
                    }
                }

                let baseq = qual[qpos_cursor.unwrap() as usize];
                let stat = baseq2stat.entry(baseq).or_insert(BaseQStat::new(baseq));
                if qpos.is_some() {
                    stat.depth += 1;
                }
                if qpos.is_none() {
                    stat.del += 1;
                    continue;
                }
                if rpos.is_none() { // ins 需要加到 下一个 baseq 上
                    let next_qpos = qpos_cursor.unwrap() + 1;
                    if next_qpos >= query_end {
                        continue;
                    }

                    let next_baseq = qual[next_qpos as usize];
                    let stat = baseq2stat.entry(next_baseq).or_insert(BaseQStat::new(next_baseq));
                    stat.ins += 1;
                    continue;
                }

                unsafe {
                    if *refseq_bytes.get_unchecked(rpos.unwrap() as usize)
                        == *query_seq.get_unchecked(qpos.unwrap() as usize)
                    {
                        // eq
                        stat.eq += 1;
                    } else {
                        // diff
                        stat.diff += 1;
                    }
                }
                
            }
        }
        let mut baseq2stat = baseq2stat.into_iter().collect::<Vec<_>>();
        baseq2stat.sort_by_key(|v| v.0);
        baseq2stat.iter().for_each(|v| {
            writeln!(&mut o_file_buff_writer, "{}\t{}", refname, v.1).unwrap();
        });
    }
    pb.finish();
}
