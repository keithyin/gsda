use std::{
    cmp,
    fmt::Display,
    fs,
    io::{BufWriter, Write},
};

use gskits::{
    file_reader::{bed_reader::BedInfo, vcf_reader::VcfInfo},
    pbar,
};
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};

use crate::cli::AlignedBamParams;

use super::FastaData;

struct ErrorLocusInfo {
    qstart: usize,
    qend: usize,
    rstart: usize,
    rend: usize,

    qseq: String,
    rseq: String,
}

impl ErrorLocusInfo {
    fn new(
        qstart: usize,
        qend: usize,
        rstart: usize,
        rend: usize,
        qseq: String,
        rseq: String,
    ) -> Self {
        Self {
            qstart,
            qend,
            rstart,
            rend,
            qseq,
            rseq,
        }
    }

    fn csv_header() -> &'static str {
        "qstart\rqend\rstart\trend\tqseq\tdepth\rseq"
    }
}

impl Display for ErrorLocusInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
            self.qstart, self.qend, self.rstart, self.rend, self.qseq, self.rseq
        )
    }
}

/// qname refname qpos_start qpos_end rpos_start, rpos_end query ref
pub fn fact_error_locus_info(
    args: &AlignedBamParams,
    output_dir: &str,
    hc_regions: Option<&BedInfo>,
    hc_variants: Option<&VcfInfo>,
    fasta_data: &FastaData,
) {
    let bam_file = &args.bam;

    let mut bam_h = bam::IndexedReader::from_path(bam_file).unwrap();
    bam_h.set_threads(4).unwrap();

    let o_filepath = format!("{}/fact_error_locus_info.csv", output_dir);
    let o_file = fs::File::create(o_filepath).unwrap();
    let mut o_file_buff_writer = BufWriter::new(o_file);
    writeln!(
        &mut o_file_buff_writer,
        "refname\t{}",
        ErrorLocusInfo::csv_header()
    )
    .unwrap();

    let pb = pbar::get_spin_pb(format!("fact_error_locus_info"), pbar::DEFAULT_INTERVAL);

    for (refname, refseq) in fasta_data.get_ref_name2seq() {
        bam_h
            .fetch((refname, 0, refseq.len() as u64))
            .expect(&format!("fetch {} error", refname));
        let ref_seq_len = refseq.len();

        let refseq_bytes = refseq.as_bytes();

        for record in bam_h.records() {
            pb.inc(1);
            let record = record.unwrap();
            if record.is_secondary() || record.is_unmapped() || record.is_supplementary() {
                continue;
            }

            let ref_start = record.reference_start();
            let ref_end = record.reference_end();

            let mut rpos_cursor = None;
            let mut qpos_cursor = None;

            let mut qdiff_start = None;
            let mut rdiff_start = None;
            let mut start_idx = None;

            let query_seq = record.seq().as_bytes();

            let mut err_locus_info = None;

            let aligned_pairs_full = record.aligned_pairs_full().into_iter().collect::<Vec<_>>();

            for (idx, [qpos, rpos]) in aligned_pairs_full.iter().enumerate() {
                let qpos = *qpos;
                let rpos = *rpos;

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

                if qpos.is_none() || rpos.is_none() {
                    // deletion
                    if qdiff_start.is_none() || rdiff_start.is_none() {
                        qdiff_start = qpos_cursor;
                        rdiff_start = rpos_cursor;
                        start_idx = Some(idx);
                    }
                    continue;
                }

                unsafe {
                    if *refseq_bytes.get_unchecked(rpos.unwrap() as usize)
                        == *query_seq.get_unchecked(qpos.unwrap() as usize)
                    {
                        // eq
                        if qdiff_start.is_some() {

                            

                            ErrorLocusInfo::new(
                                qdiff_start.unwrap() as usize,
                                qpos.unwrap() as usize,
                                rdiff_start.unwrap() as usize,
                                rpos.unwrap() as usize,
                                qseq,
                                rseq,
                            );
                        }

                        qdiff_start = None;
                        rdiff_start = None;
                        start_idx = None;
                    } else {
                        // diff
                        if qdiff_start.is_none() || rdiff_start.is_none() {
                            qdiff_start = qpos_cursor;
                            rdiff_start = rpos_cursor;
                        }
                    }
                }

                if rpos_cur_or_pre == (ref_end as usize - 1) {
                    break;
                }
            }
        }
    }
    pb.finish();
}
