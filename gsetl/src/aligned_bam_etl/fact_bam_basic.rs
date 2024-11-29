use std::{
    fmt::Display,
    fs,
    io::{BufWriter, Write},
};

use gskits::{
    gsbam::bam_record_ext::BamRecordExt,
    pbar,
};
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};

use crate::cli::AlignedBamParams;

use super::FastaData;

struct BasicInfo<'a> {
    qname: String,
    refname: &'a str,
    channel: u32,
    np: u32,
    rq: f32,
    iy: f32,
    ec: f32,
    rstart: usize,
    rend: usize,
    qstart: usize,
    qend: usize,
    qlen: usize,
    is_forward: bool,
    ori_start: usize,
    ori_end: usize
}

impl<'a> BasicInfo<'a> {
    pub fn new(
        qname: String,
        refname: &'a str,
        channel: u32,
        np: u32,
        rq: f32,
        iy: f32,
        ec: f32,
        rstart: usize,
        rend: usize,
        qstart: usize,
        qend: usize,
        qlen: usize,
        is_forward: bool,
        ori_start: usize,
        ori_end: usize
    ) -> Self {
        Self {
            qname,
            refname,
            channel,
            np,
            rq,
            iy,
            ec,
            rstart,
            rend,
            qstart,
            qend,
            qlen,
            is_forward,
            ori_start,
            ori_end
        }
    }

    pub fn csv_header() -> &'static str {
        "qname\trefname\tchannel\tnp\trq\tiy\tec\trstart\trend\tqstart\tqend\tqlen\tfwd\tori_start\tori_end"
    }
}

impl<'a> Display for BasicInfo<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}",
            self.qname,
            self.refname,
            self.channel,
            self.np,
            self.rq,
            self.iy,
            self.ec,
            self.rstart,
            self.rend,
            self.qstart,
            self.qend,
            self.qlen,
            self.is_forward,
            self.ori_start,
            self.ori_end
        )
    }
}

/// qname refname qpos_start qpos_end rpos_start, rpos_end query ref
pub fn fact_bam_basic(
    args: &AlignedBamParams,
    output_dir: &str,
    fasta_data: &FastaData,
) {
    let bam_file = &args.bam;

    let mut bam_h = bam::IndexedReader::from_path(bam_file).unwrap();
    bam_h.set_threads(10).unwrap();

    let o_filepath = format!("{}/fact_aligned_bam_bam_basic.csv", output_dir);
    let o_file = fs::File::create(o_filepath).unwrap();
    let mut o_file_buff_writer = BufWriter::new(o_file);
    writeln!(&mut o_file_buff_writer, "{}", BasicInfo::csv_header()).unwrap();

    let pb = pbar::get_spin_pb(format!("fact_bam_basic"), pbar::DEFAULT_INTERVAL);

    for (refname, refseq) in fasta_data.get_ref_name2seq() {
        bam_h
            .fetch((refname, 0, refseq.len() as u64))
            .expect(&format!("fetch {} error", refname));

        for record in bam_h.records() {
            pb.inc(1);
            let record = record.unwrap();
            if record.is_secondary() || record.is_unmapped() || record.is_supplementary() {
                continue;
            }

            let record_ext = BamRecordExt::new(&record);

            let rstart = record_ext.reference_start() as usize;
            let rend = record_ext.reference_end() as usize;
            let qstart = record_ext.query_alignment_start() as usize;

            // dump
            let record_ext = BamRecordExt::new(&record);
            let qname = record_ext.get_qname();
            let qend = record_ext.query_alignment_end();
            let qlen = record.seq_len_from_cigar(true);

            let be = record_ext.get_be().unwrap_or(vec![0, 0]);

            let basic_info = BasicInfo::new(
                qname,
                refname,
                record_ext.get_ch().unwrap_or(0),
                record_ext.get_np().unwrap_or(0),
                record_ext.get_rq().unwrap_or(-1.),
                record_ext.get_identity().unwrap_or(-1.),
                record_ext.get_coverage().unwrap_or(-1.),
                rstart,
                rend,
                qstart,
                qend,
                qlen,
                !record.is_reverse(),
                be[0] as usize,
                be[1] as usize
            );

            writeln!(&mut o_file_buff_writer, "{}", basic_info).unwrap();
        }
    }
    pb.finish();
}
