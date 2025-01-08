use std::{fs, io::Write, thread, time::Duration};

use aligned_bam_etl::{
    fact_bam_basic::fact_bam_basic, fact_baseq_stat::fact_baseq_stat,
    fact_error_query_locus_info::fact_error_query_locus_info, fact_poly_info::fact_poly_info,
    fact_record_stat::fact_record_stat, fact_ref_locus_info::fact_ref_locus_info, get_hc_regions,
    get_hcvariants, FastaData,
};
use clap::Parser;
use gskits::utils::command_line_str;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

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
            let multi_pb = MultiProgress::new();

            thread::scope(|s| {
                let args = &args;
                let hc_regions = hc_regions.as_ref();
                let hc_variants = hc_variants.as_ref();
                let reffasta = &reffasta;

                if param.fact_record_stat == 1 {
                    let pb = ProgressBar::new_spinner();
                    let pb = multi_pb.add(pb);
                    s.spawn(move || {
                        fact_record_stat(
                            param,
                            &args.output_dir,
                            hc_regions,
                            hc_variants,
                            reffasta,
                            pb,
                        )
                    });
                }

                if param.fact_ref_locus_info == 1 {
                    let pb = ProgressBar::new_spinner();
                    let pb = multi_pb.add(pb);
                    s.spawn(move || {
                        fact_ref_locus_info(
                            param,
                            &args.output_dir,
                            hc_regions,
                            hc_variants,
                            reffasta,
                            pb,
                        )
                    });
                }

                if param.fact_error_query_locus_info == 1 {
                    let pb = ProgressBar::new_spinner();
                    let pb = multi_pb.add(pb);
                    s.spawn(move || {
                        fact_error_query_locus_info(
                            param,
                            &args.output_dir,
                            hc_regions,
                            hc_variants,
                            reffasta,
                            pb,
                        )
                    });
                }

                if param.fact_bam_basic == 1 {
                    let pb = ProgressBar::new_spinner();
                    let pb = multi_pb.add(pb);
                    s.spawn(move || fact_bam_basic(param, &args.output_dir, reffasta, pb));
                }

                if param.fact_baseq_stat == 1 {
                    let pb = ProgressBar::new_spinner();
                    let pb = multi_pb.add(pb);
                    s.spawn(move || {
                        fact_baseq_stat(
                            param,
                            &args.output_dir,
                            hc_regions,
                            hc_variants,
                            reffasta,
                            pb,
                        )
                    });
                }

                if param.fact_poly_info == 1 {
                    let pb = ProgressBar::new_spinner();
                    let pb = multi_pb.add(pb);
                    s.spawn(move || {
                        fact_poly_info(
                            param,
                            &args.output_dir,
                            hc_regions,
                            hc_variants,
                            &reffasta,
                            pb,
                        )
                    });
                }
            });
        }
        cli::Subcommands::NonAlignedBam(_param) => {}
    }
}

pub fn set_spin_pb(pb: ProgressBar, message: String, interval: Duration) -> ProgressBar {
    pb.set_style(
        ProgressStyle::with_template("{msg} {spinner} {human_pos}  {per_sec} {elapsed}").unwrap(),
    );
    pb.enable_steady_tick(interval);
    pb.set_message(message);
    pb
}
