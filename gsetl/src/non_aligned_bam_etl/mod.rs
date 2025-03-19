use std::{
    fmt::Display,
    fs::File,
    io::{BufWriter, Write},
    path,
};

use gskits::{
    dna::SEQ_NT4_TABLE,
    gsbam::bam_record_ext::BamRecordExt,
    pbar::{get_spin_pb, DEFAULT_INTERVAL},
};
use polars::{df, prelude::*};
use rust_htslib::bam::{Read, Reader, Record};

use crate::{cli::NonAlignedBamParams, utils};

/// 最终是要生成一个 Run 的聚合数据。下面的是每个 base 的结果作为中间结果
/// dw
/// dw-N
/// ar
/// ar-N
/// cq 中位数
/// cr 无脑均值
/// oe 中位数
pub fn fact_bam_basic2(bam_file: &str) {
    let mut reader = Reader::from_path(bam_file).unwrap();
    reader.set_threads(40).unwrap();

    let mut all_bases = vec![];
    let mut all_dwell_times = vec![];
    let mut all_arrival_times = vec![];
    let mut all_capture_rates = vec![];

    let mut all_channel_cqs = vec![];
    let mut all_channel_oes = vec![];
    let pb = get_spin_pb(format!("reading {bam_file}"), DEFAULT_INTERVAL);

    for record in reader.records() {
        let record = record.unwrap();
        pb.inc(1);
        let record_ext = BamRecordExt::new(&record);

        let bases = record_ext.get_seq();
        all_bases.extend_from_slice(
            &bases
                .as_bytes()
                .into_iter()
                .map(|v| (*v as char).to_string())
                .collect::<Vec<_>>(),
        );

        let dwell_times = record_ext
            .get_dw()
            .map(|v| v.into_iter().map(|v| v as i32).collect::<Vec<i32>>())
            .unwrap_or(vec![-1; bases.len()]);
        all_dwell_times.extend_from_slice(&dwell_times);

        let arrival_time = record_ext
            .get_ar()
            .map(|v| v.into_iter().map(|v| v as i32).collect::<Vec<i32>>())
            .unwrap_or(vec![-1; bases.len()]);
        all_arrival_times.extend_from_slice(&arrival_time);

        let capture_rates = record_ext.get_float_cr().unwrap_or(vec![-1.0; bases.len()]);
        all_capture_rates.extend_from_slice(&capture_rates);

        let cq = record_ext.get_cq().unwrap_or(-1.0);
        all_channel_cqs.push(cq);

        let oe = record_ext.get_oe().unwrap_or(-1.0);
        all_channel_oes.push(oe);
    }
    pb.finish();

    let base_level_df = df!(
        "base" => &all_bases,
        "dw" => &all_dwell_times,
        "ar" => &all_arrival_times,
        "cr" => &all_capture_rates
    )
    .unwrap();

    base_level_overall_aggr(base_level_df.clone().lazy());
    base_level_base_aggr(base_level_df.clone().lazy());
    speed(base_level_df.clone().lazy());

    let channel_level_df = df!(
        "oe" => &all_channel_oes,
        "cq" => &all_channel_cqs,
    )
    .unwrap();
    channel_level_overall_aggr(channel_level_df.clone().lazy());
}

fn base_level_overall_aggr(mut df: LazyFrame) {
    df = df.select([
        col("dw").mean().alias("dw-mean"),
        col("ar").mean().alias("ar-mean"),
        col("cr").mean().alias("cr-mean"),
    ]);

    println!("{:?}", df.collect().unwrap());
}

fn base_level_base_aggr(mut df: LazyFrame) {
    df = df.group_by([col("base")]).agg([
        col("dw").median().alias("dw-median"),
        col("ar").median().alias("ar-median"),
    ]);

    println!("{:?}", df.collect().unwrap());
}

fn channel_level_overall_aggr(mut df: LazyFrame) {
    df = df.select([
        col("oe").median().alias("oe-median"),
        col("cq").median().alias("cq-median"),
    ]);

    println!("{:?}", df.collect().unwrap());
}

fn speed(mut df: LazyFrame) {
    df = df
        .select([
            len().cast(DataType::Float64).alias("base_cnt"),
            ((col("dw").sum() + col("ar").sum()).cast(DataType::Float64) * lit(0.002_f64))
                .alias("secs"),
        ])
        .select([(col("base_cnt") / col("secs")).alias("speed")]);

    println!("{:?}", df.collect().unwrap());
}

#[derive(Debug, Default)]
struct Stat {
    dw_mean: f64,
    dw_acgt_mean: [f64; 4],

    ar_mean: f64,
    ar_acgt_mean: [f64; 4],

    cr_mean: f64,

    cq_mean: f64,
    all_oe: Vec<f64>,
    oe_median: f64,
    speed: f64,

    dw_add_ar: f64,

    acgt_cnt: [usize; 4],
    base_cnt: usize,
    channel_cnt: usize,
}

impl Stat {
    pub fn update(&mut self, rec: &Record) {
        let record_ext = BamRecordExt::new(rec);
        let bases = record_ext.get_seq();
        let dwell_times = record_ext
            .get_dw()
            .map(|v| v.into_iter().map(|v| v as f64).collect::<Vec<_>>())
            .unwrap_or(vec![-1.0; bases.len()]);
        let arrival_time = record_ext
            .get_ar()
            .map(|v| v.into_iter().map(|v| v as f64).collect::<Vec<_>>())
            .unwrap_or(vec![-1.0; bases.len()]);

        let capture_rates = record_ext
            .get_float_cr()
            .map(|v| v.into_iter().map(|v| v as f64).collect::<Vec<_>>())
            .unwrap_or(vec![-1.0; bases.len()]);
        let cq = record_ext.get_cq().unwrap_or(-1.0) as f64;
        let oe = record_ext.get_oe().unwrap_or(-1.0) as f64;

        self.update_base_level_info(
            bases.as_bytes(),
            &dwell_times,
            &arrival_time,
            &capture_rates,
        );
        self.update_channel_info(cq, oe);
    }

    fn update_base_level_info(&mut self, bases: &[u8], dws: &[f64], ars: &[f64], crs: &[f64]) {
        let mut base_cnt_cursor = self.base_cnt;
        crs.iter().copied().for_each(|cr| {
            self.cr_mean =
                (base_cnt_cursor as f64 * self.cr_mean + cr) / (base_cnt_cursor as f64 + 1.0);
            base_cnt_cursor += 1;
        });

        let mut base_cnt_cursor = self.base_cnt;
        let mut acgt_cnt_cursor = self.acgt_cnt.clone();
        bases.iter().zip(dws.iter()).for_each(|(&base, &dw)| {
            self.dw_mean =
                (base_cnt_cursor as f64 * self.dw_mean + dw) / (base_cnt_cursor as f64 + 1.0);
            base_cnt_cursor += 1;

            let base_idx = SEQ_NT4_TABLE[base as usize] as usize;
            self.dw_acgt_mean[base_idx] =
                (acgt_cnt_cursor[base_idx] as f64 * self.dw_acgt_mean[base_idx] + dw)
                    / (acgt_cnt_cursor[base_idx] as f64 + 1.0);
            acgt_cnt_cursor[base_idx] += 1;

            self.dw_add_ar += dw;
        });

        let mut base_cnt_cursor = self.base_cnt;
        let mut acgt_cnt_cursor = self.acgt_cnt.clone();
        bases.iter().zip(ars.iter()).for_each(|(&base, &ar)| {
            self.ar_mean =
                (base_cnt_cursor as f64 * self.ar_mean + ar) / (base_cnt_cursor as f64 + 1.0);
            base_cnt_cursor += 1;

            let base_idx = SEQ_NT4_TABLE[base as usize] as usize;
            self.ar_acgt_mean[base_idx] =
                (acgt_cnt_cursor[base_idx] as f64 * self.ar_acgt_mean[base_idx] + ar)
                    / (acgt_cnt_cursor[base_idx] as f64 + 1.0);
            acgt_cnt_cursor[base_idx] += 1;

            self.dw_add_ar += ar;
        });

        self.base_cnt = base_cnt_cursor;
        self.acgt_cnt = acgt_cnt_cursor;
    }

    fn update_channel_info(&mut self, cq: f64, oe: f64) {
        self.cq_mean =
            (self.channel_cnt as f64 * self.cq_mean + cq) / (self.channel_cnt as f64 + 1.0);
        self.all_oe.push(oe);

        self.channel_cnt += 1;
    }

    pub fn finish(&mut self) {
        let med = utils::median(&mut self.all_oe).unwrap_or(-1.0);
        self.oe_median = med;
        self.all_oe.clear();
        if self.dw_add_ar > 0.0 {
            self.speed = self.base_cnt as f64 / (self.dw_add_ar as f64 * 0.002);
        }
    }

    pub fn csv(&self) -> String {
        let mut o_str = String::new();
        o_str.push_str(&format!("dw-mean\t{:.2}\n", self.dw_mean));
        o_str.push_str(&format!("ar-mean\t{:.2}\n", self.ar_mean));
        o_str.push_str(&format!("cq-mean\t{:.2}\n", self.cq_mean));
        o_str.push_str(&format!("cr-mean\t{:.2}\n", self.cr_mean));
        o_str.push_str(&format!("oe-median\t{:.2}\n", self.oe_median));
        o_str.push_str(&format!("speed\t{:.2}\n", self.speed));

        ["A", "C", "G", "T"]
            .iter()
            .enumerate()
            .for_each(|(idx, &base)| {
                o_str.push_str(&format!(
                    "dw-mean-{}\t{:.2}\n",
                    base, self.dw_acgt_mean[idx]
                ));
            });

        ["A", "C", "G", "T"]
            .iter()
            .enumerate()
            .for_each(|(idx, &base)| {
                o_str.push_str(&format!(
                    "ar-mean-{}\t{:.2}\n",
                    base, self.ar_acgt_mean[idx]
                ));
            });

        o_str
    }
}

pub fn fact_bam_basic(param: &NonAlignedBamParams, output_dir: &str) {
    let bam_file = &param.bam;
    let filestem = path::Path::new(bam_file)
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();

    let oup_filepath = param
        .o_filepath
        .clone()
        .unwrap_or(format!("{output_dir}/{filestem}.non_aligned_aggr.csv"));

    let bam_threads = param.bam_threads.unwrap_or(num_cpus::get_physical());
    let mut reader = Reader::from_path(bam_file).unwrap();
    reader.set_threads(bam_threads).unwrap();

    let pb = get_spin_pb(format!("reading {bam_file}"), DEFAULT_INTERVAL);
    let mut stat = Stat::default();
    for record in reader.records() {
        let record = record.unwrap();
        pb.inc(1);
        stat.update(&record);
    }
    stat.finish();
    pb.finish();

    let out_file = File::create(&oup_filepath).unwrap();
    let mut out_writer = BufWriter::new(out_file);
    println!("{:?}", stat);

    writeln!(&mut out_writer, "name\tvalue").unwrap();
    write!(&mut out_writer, "{}", stat.csv()).unwrap();
}

#[cfg(test)]
mod test {
    use super::Stat;

    #[test]
    fn test_stat() {
        let mut stat = Stat::default();
        let bases = b"ACGTACGT";
        let dws = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let ars = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let crs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        stat.update_base_level_info(bases, &dws, &ars, &crs);

        stat.update_base_level_info(bases, &dws, &ars, &crs);

        println!("{:?}", stat);
    }
}
