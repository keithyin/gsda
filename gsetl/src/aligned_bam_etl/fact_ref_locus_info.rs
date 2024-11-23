use std::{
    cmp,
    fmt::Display,
    fs,
    io::{BufWriter, Write},
};

use gskits::{file_reader::{bed_reader::BedInfo, vcf_reader::VcfInfo}, pbar};
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};

use crate::cli::AlignedBamParams;

use super::FastaData;

struct LocusStat {
    pos: usize,
    eq: usize,
    diff: usize,
    ins: usize,
    del: usize,
    depth: usize, // how much records that aligned to this position
    around_bases: String,
}

impl LocusStat {
    fn new(pos: usize, around_bases: String) -> Self {
        Self {
            pos: pos,
            eq: 0,
            diff: 0,
            ins: 0,
            del: 0,
            depth: 0,
            around_bases: around_bases,
        }
    }

    fn csv_header() -> &'static str {
        "pos\teq\tdiff\tins\tdel\tdepth\taround_bases"
    }
}

impl Display for LocusStat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}",
            self.pos, self.eq, self.diff, self.ins, self.del, self.depth, self.around_bases
        )
    }
}

pub fn fact_ref_locus_info(
    args: &AlignedBamParams,
    output_dir: &str,
    hc_regions: Option<&BedInfo>,
    hc_variants: Option<&VcfInfo>,
    fasta_data: &FastaData,
) {
    let bam_file = &args.bam;

    let mut bam_h = bam::IndexedReader::from_path(bam_file).unwrap();
    bam_h.set_threads(4).unwrap();

    let o_filepath = format!("{}/fact_ref_locus_info.csv", output_dir);
    let o_file = fs::File::create(o_filepath).unwrap();
    let mut o_file_buff_writer = BufWriter::new(o_file);
    writeln!(
        &mut o_file_buff_writer,
        "refname\t{}",
        LocusStat::csv_header()
    )
    .unwrap();

    let pb = pbar::get_spin_pb(format!("fact_ref_locus_info"), pbar::DEFAULT_INTERVAL);

    for (refname, refseq) in fasta_data.get_ref_name2seq() {
        bam_h
            .fetch((refname, 0, refseq.len() as u64))
            .expect(&format!("fetch {} error", refname));
        let ref_seq_len = refseq.len();

        let mut ref_locus_stat = (0..refseq.len())
            .into_iter()
            .map(|pos| {
                let around_start = if pos < 10 { 0 } else { pos - 10 };
                let around_end = cmp::min(pos + 11, ref_seq_len);
                LocusStat::new(
                    pos,
                    format!(
                        "{} {} {}",
                        &refseq[around_start..pos],
                        &refseq[pos..pos + 1],
                        &refseq[pos + 1..around_end]
                    ),
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

            let ref_start = record.reference_start();
            let ref_end = record.reference_end();

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

                if qpos_cursor.is_none() {
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
                        ref_locus_stat.get_unchecked_mut(rpos_cur_or_pre).ins += 1;
                    }
                    continue;
                }

                unsafe {
                    if *refseq.get_unchecked(rpos.unwrap() as usize)
                        == *query_seq.get_unchecked(qpos.unwrap() as usize)
                    {
                        ref_locus_stat.get_unchecked_mut(rpos_cur_or_pre).eq += 1;
                    } else {
                        ref_locus_stat.get_unchecked_mut(rpos_cur_or_pre).diff += 1;
                    }
                }

                if rpos_cur_or_pre == (ref_end as usize - 1) {
                    break;
                }
            }
        }

        ref_locus_stat.iter().for_each(|locus_stat| {
            writeln!(&mut o_file_buff_writer, "{}\t{}", refname, locus_stat).unwrap()
        });
    }
    pb.finish();
}
