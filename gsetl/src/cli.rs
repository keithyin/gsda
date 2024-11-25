use std::{fs, path};

use clap::{self, Args, Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(version, about, long_about = None)]
pub struct Cli {
    #[arg(long = "outdir")]
    pub output_dir: String,
    #[arg(short = 'f', help = "remove everything in the ${outdir}")]
    pub force: bool,

    #[command(subcommand)]
    pub commands: Subcommands,
}

impl Cli {
    pub fn build_output_dir(&self) {
        if path::Path::new(&self.output_dir).exists() && !self.force {
            panic!(
                "output_dir: {} exists, use -f or change the output dir",
                self.output_dir
            );
        }

        if path::Path::new(&self.output_dir).exists() {
            fs::remove_dir_all(&self.output_dir)
                .expect(&format!("remove dir error {}", self.output_dir));
        }

        fs::create_dir_all(&self.output_dir)
            .expect(&format!("create dir error. {}", self.output_dir));
    }
}

#[derive(Debug, Subcommand, Clone)]
pub enum Subcommands {
    AlignedBam(AlignedBamParams),
    NonAlignedBam(NonAlignedBamParams),
}

#[derive(Debug, Args, Clone)]
pub struct AlignedBamParams {
    #[arg(long = "bam")]
    pub bam: String,

    #[arg(long = "ref-file", help="fasta/fastq/bam")]
    pub ref_file: String,

    #[arg(
        long = "hcregion",
        help = "bed file that contains the confidence regions"
    )]
    pub hc_regions_file: Option<String>,

    #[arg(long = "vcf", help = "vcf file that contains the variant locus")]
    pub vcf_file: Option<String>,

    #[arg(
        long = "rnames",
        help = "only the part of reference are considered, default: all"
    )]
    pub ref_names: Vec<String>,
}

#[derive(Debug, Args, Clone)]
pub struct NonAlignedBamParams {
    #[arg(long = "bam")]
    pub bam: String,
}
