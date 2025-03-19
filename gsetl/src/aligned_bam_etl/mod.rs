use std::{collections::HashMap, fs, io::BufReader};

use gskits::{
    fastx_reader::{fasta_reader::FastaFileReader, fastq_reader::FastqReader},
    file_reader::{bed_reader::BedInfo, vcf_reader::VcfInfo},
    gsbam::bam_record_ext::{BamReader, BamRecord, BamRecordExt},
};
use rust_htslib::bam::Read;

use crate::cli::AlignedBamParams;

pub mod fact_record_stat;
pub mod fact_ref_locus_info;
pub mod fact_error_query_locus_info;
pub mod fact_bam_basic;
pub mod fact_baseq_stat;
pub mod fact_poly_info;

pub struct FastaData {
    ref_name2seq: HashMap<String, String>,
}

impl FastaData {
    pub fn new(filepath: &str) -> Self {
        let ref_name2seq =
            if filepath.ends_with("fa") || filepath.ends_with("fasta") || filepath.ends_with("fna")
            {
                let reader = FastaFileReader::new(filepath.to_string());
                reader
                    .into_iter()
                    .map(|record| (record.name, record.seq))
                    .collect::<HashMap<String, String>>()
            } else if filepath.ends_with("fq") || filepath.ends_with("fastq") {
                let reader = FastqReader::new(filepath.to_string());
                reader
                    .into_iter()
                    .map(|record| (record.name, record.seq))
                    .collect::<HashMap<String, String>>()
            } else if filepath.ends_with("bam") {
                let mut reader = BamReader::from_path(filepath).unwrap();
                reader.set_threads(4).unwrap();
                reader
                    .records()
                    .into_iter()
                    .map(|record| {
                        let record = record.unwrap();
                        let record_ext = BamRecordExt::new(&record);
                        (record_ext.get_qname(), record_ext.get_seq())
                    })
                    .collect::<HashMap<String, String>>()
            } else {
                panic!(
                    "invalid file format. fa/fasta/fna/fq/fastq/bam supported. but got {}",
                    filepath
                );
            };

        Self { ref_name2seq }
    }

    #[allow(unused)]
    pub fn get_ref_seq(&self, refname: &str) -> Option<&str> {
        self.ref_name2seq
            .get(refname)
            .and_then(|v| Some(v.as_str()))
    }

    pub fn get_ref_name2seq(&self) -> &HashMap<String, String> {
        &self.ref_name2seq
    }
}

pub fn get_hc_regions(region_file: Option<&str>) -> Option<BedInfo> {
    if let Some(region_file) = region_file {
        let reader = fs::File::open(region_file).unwrap();
        let mut buf_reader = BufReader::new(reader);
        Some(BedInfo::new(&mut buf_reader))
    } else {
        None
    }
}

pub fn get_hcvariants(vcf_file: Option<&str>) -> Option<VcfInfo> {
    if let Some(vcf_file) = vcf_file {
        let reader = fs::File::open(vcf_file).unwrap();
        let mut buf_reader = BufReader::new(reader);
        Some(VcfInfo::new(&mut buf_reader))
    } else {
        None
    }
}


pub fn audit(record: &BamRecord, param: &AlignedBamParams) -> bool {
    if record.is_unmapped() {
        return false;
    }

    if !param.use_seco && record.is_secondary() {
        return false;
    }

    if !param.use_supp && record.is_supplementary() {
        return false;
    }

    return true;
}