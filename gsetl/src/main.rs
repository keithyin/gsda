use std::{fs, io::Write, thread};

use aligned_bam_etl::{
    fact_bam_basic::fact_bam_basic, fact_baseq_stat::fact_baseq_stat,
    fact_error_query_locus_info::fact_error_query_locus_info, fact_poly_info::fact_poly_info,
    fact_record_stat::fact_record_stat, fact_ref_locus_info::fact_ref_locus_info, get_hc_regions,
    get_hcvariants, FastaData,
};
use clap::Parser;
use gskits::utils::command_line_str;

mod aligned_bam_etl;
mod cli;
mod non_aligned_bam_etl;
mod poly_n;

fn main() {
    let args = cli::Cli::parse();
    args.build_output_dir();

    let mut meta_file = fs::File::create(&format!("{}/meta.txt", args.output_dir)).unwrap();
    meta_file
        .write_all(format!("version: {}\n", env!("CARGO_PKG_VERSION")).as_bytes())
        .unwrap();
    meta_file
        .write_all(format!("cmd_line: {}\n", command_line_str()).as_bytes())
        .unwrap();

    match &args.commands {
        cli::Subcommands::AlignedBam(param) => {
            let bed_filepath = param.hc_regions_file.clone();
            let bed_thread = thread::spawn(move || {
                get_hc_regions(bed_filepath.as_ref().and_then(|v| Some(v.as_str())))
            });

            let vcf_filepath = param.vcf_file.clone();
            let vcf_thread = thread::spawn(move || {
                get_hcvariants(vcf_filepath.as_ref().and_then(|v| Some(v.as_str())))
            });

            let hc_regions = bed_thread.join().unwrap();
            let hc_variants = vcf_thread.join().unwrap();
            let reffasta = FastaData::new(&param.ref_file);

            fact_record_stat(
                param,
                &args.output_dir,
                hc_regions.as_ref(),
                hc_variants.as_ref(),
                &reffasta,
            );
            fact_ref_locus_info(
                param,
                &args.output_dir,
                hc_regions.as_ref(),
                hc_variants.as_ref(),
                &reffasta,
            );
            fact_error_query_locus_info(
                param,
                &args.output_dir,
                hc_regions.as_ref(),
                hc_variants.as_ref(),
                &reffasta,
            );
            fact_bam_basic(param, &args.output_dir, &reffasta);
            fact_baseq_stat(
                param,
                &args.output_dir,
                hc_regions.as_ref(),
                hc_variants.as_ref(),
                &reffasta,
            );

            fact_poly_info(
                param,
                &args.output_dir,
                hc_regions.as_ref(),
                hc_variants.as_ref(),
                &reffasta,
            );
        }
        cli::Subcommands::NonAlignedBam(_param) => {}
    }
}
