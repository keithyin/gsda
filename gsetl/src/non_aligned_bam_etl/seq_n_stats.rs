use std::{
    fmt::Display,
    io::{BufWriter, Write},
};

use crossbeam::channel::{Receiver, Sender};
use gskits::{
    gsbam::bam_record_ext::BamRecordExt,
    pbar::{get_spin_pb, DEFAULT_INTERVAL},
};
use rust_htslib::bam::{Read, Record};

use crate::cli::NonAlignedBamSeqNStatsParams;

pub fn seq_n_stats_main(param: &NonAlignedBamSeqNStatsParams, output_dir: &str) {
    if param.length_percentile_thr.is_none() && param.length_thr.is_none() {
        panic!("--length-thr and --length-percentile-thr can't be all None");
    }

    let length_thr = param
        .length_thr
        .unwrap_or_else(|| length_percentile(&param.bam, param.length_percentile_thr.unwrap()));

    println!("length-thr={}", length_thr);

    let bam_path = &param.bam;
    let n = param.n;

    let out_dir = std::path::Path::new(output_dir);
    if !out_dir.exists() {
        std::fs::create_dir_all(out_dir).expect("create dir error");
    }
    let file_stem = std::path::Path::new(bam_path)
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();
    let output_filepath = out_dir.join(format!("{}.seq_n_stats.csv", file_stem));

    let first_n_out_fasta = out_dir.join(format!("{}.first-n.fasta", file_stem));
    let last_n_out_fasta = out_dir.join(format!("{}.last-n.fasta", file_stem));

    let mut writer = BufWriter::new(std::fs::File::create(output_filepath.clone()).unwrap());
    writeln!(
        &mut writer,
        "qname\tdw-first-n-median\tdw-last-n-median\tar-first-n-median\tar-last-n-median\tdw-first-n-mean\tdw-last-n-mean\tar-first-n-mean\tar-last-n-mean"
    )
    .unwrap();

    let mut first_n_writer = BufWriter::new(std::fs::File::create(first_n_out_fasta).unwrap());
    let mut last_n_writer = BufWriter::new(std::fs::File::create(last_n_out_fasta).unwrap());

    std::thread::scope(|thread_scope| {
        let (record_sender, record_recv) = crossbeam::channel::bounded(1000);
        thread_scope.spawn(move || {
            bam_reader(bam_path, Some(length_thr), record_sender);
        });

        let (stat_sender, stat_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..(num_cpus::get() * 4 / 5).max(10) {
            thread_scope.spawn({
                let record_recv = record_recv.clone();
                let stat_sender = stat_sender.clone();
                move || {
                    stats_worker(record_recv, stat_sender, n);
                }
            });
        }
        drop(record_recv);
        drop(stat_sender);

        let pb = get_spin_pb(
            format!("seq_n_stas: writing to {:?}", output_filepath),
            DEFAULT_INTERVAL,
        );
        for stat in stat_recv {
            pb.inc(1);
            writeln!(&mut writer, "{}", stat).unwrap();

            writeln!(&mut first_n_writer, ">{}\n{}", stat.qname, stat.first_n).unwrap();
            writeln!(&mut last_n_writer, ">{}\n{}", stat.qname, stat.last_n).unwrap();
        }
        pb.finish();
    });
}

#[derive(Debug, Clone)]
struct Stats {
    qname: String,
    ar_first_n_median: f64,
    ar_last_n_median: f64,
    dw_first_n_median: f64,
    dw_last_n_median: f64,

    ar_first_n_mean: f64,
    ar_last_n_mean: f64,
    dw_first_n_mean: f64,
    dw_last_n_mean: f64,
    first_n: String,
    last_n: String,
}

impl Display for Stats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname,
            self.dw_first_n_median,
            self.dw_last_n_median,
            self.ar_first_n_median,
            self.ar_last_n_median,
            self.dw_first_n_mean,
            self.dw_last_n_mean,
            self.ar_first_n_mean,
            self.ar_last_n_mean,
        )
    }
}

fn length_percentile(bam: &str, perentile: usize) -> usize {
    assert!(
        perentile <= 100,
        "expected perentile>=0 && perentile<=100. but got {}",
        perentile
    );

    let mut bam_reader = rust_htslib::bam::Reader::from_path(bam).unwrap();
    bam_reader.set_threads(num_cpus::get().max(4)).unwrap();

    let mut lengths = bam_reader
        .records()
        .into_iter()
        .map(|rec| rec.unwrap().seq_len())
        .collect::<Vec<usize>>();
    lengths.sort();
    if perentile == 0 {
        return 0;
    }
    if perentile == 100 {
        return usize::MAX;
    }

    let percentile_pos = ((lengths.len() as f32) * (perentile as f32 / 100.)).round() as usize;
    return lengths[percentile_pos];
}

fn bam_reader(bam: &str, length_thr: Option<usize>, sender: Sender<Record>) {
    let mut bam_reader = rust_htslib::bam::Reader::from_path(bam).unwrap();
    bam_reader
        .set_threads((num_cpus::get_physical() / 5).max(10))
        .unwrap();
    let mut record = rust_htslib::bam::Record::new();
    loop {
        if let Some(Ok(_)) = bam_reader.read(&mut record) {
            if record.seq_len() >= length_thr.unwrap_or(0) {
                sender.send(record).unwrap();
                record = rust_htslib::bam::Record::new();
            }
        } else {
            break;
        }
    }
}

fn stats_worker(recv: Receiver<Record>, sender: Sender<Stats>, n: usize) {
    for record in recv {
        let record_ext = BamRecordExt::new(&record);

        let seq = record_ext.get_seq();
        let n = n.min(seq.len());
        let first_n_seq = &seq[0..n];
        let last_n_seq = &seq[(seq.len() - n)..];

        let mut ar = record_ext.get_ar().unwrap();
        let mut dw = record_ext.get_dw().unwrap();
        let (ar_first_n_median, ar_last_n_median) = calculate_first_n_and_last_n_median(&mut ar, n);
        let (dw_first_n_median, dw_last_n_median) = calculate_first_n_and_last_n_median(&mut dw, n);

        let (ar_first_n_mean, ar_last_n_mean) = calculate_first_n_and_last_n_mean(&mut ar, n);
        let (dw_first_n_mean, dw_last_n_mean) = calculate_first_n_and_last_n_mean(&mut dw, n);

        sender
            .send(Stats {
                qname: record_ext.get_qname(),
                ar_first_n_median: ar_first_n_median,
                ar_last_n_median: ar_last_n_median,
                dw_first_n_median: dw_first_n_median,
                dw_last_n_median: dw_last_n_median,

                ar_first_n_mean: ar_first_n_mean,
                ar_last_n_mean: ar_last_n_mean,
                dw_first_n_mean: dw_first_n_mean,
                dw_last_n_mean: dw_last_n_mean,
                first_n: first_n_seq.to_string(),
                last_n: last_n_seq.to_string(),
            })
            .unwrap();
    }
}

fn calculate_first_n_and_last_n_median(values: &mut [u32], n: usize) -> (f64, f64) {
    let len = values.len();
    let n = n.min(len);

    if len <= 1 {
        return (0.0, 0.0);
    }

    let first_n_median = calculate_median_usize(&mut values[1..n]).unwrap_or(0.0);
    let last_n_median = calculate_median_usize(&mut values[(len - n).max(1)..]).unwrap_or(0.0);
    (first_n_median, last_n_median)
}

fn calculate_first_n_and_last_n_mean(values: &mut [u32], n: usize) -> (f64, f64) {
    let len = values.len();
    let n = n.min(len);

    if len <= 1 {
        return (0.0, 0.0);
    }

    let first_n_mean = calculate_mean_usize(&mut values[1..n]).unwrap_or(0.0);
    let last_n_mean = calculate_mean_usize(&mut values[(len - n).max(1)..]).unwrap_or(0.0);
    (first_n_mean, last_n_mean)
}

fn calculate_median_usize(numbers: &mut [u32]) -> Option<f64> {
    // 如果向量为空，则没有中位数。

    if numbers.is_empty() {
        return None;
    }

    // let v = numbers.iter().map(|v| *v as usize).sum::<usize>();
    // Some(v as f64 / numbers.len() as f64)

    // 对向量进行排序。
    numbers.sort();

    let mid = numbers.len() / 2;

    if numbers.len() % 2 == 1 {
        // 奇数个元素：中位数是中间的元素。
        Some(numbers[mid] as f64)
    } else {
        // 偶数个元素：中位数是中间两个元素的平均值。
        Some((numbers[mid - 1] + numbers[mid]) as f64 / 2.0)
    }
}

fn calculate_mean_usize(numbers: &mut [u32]) -> Option<f64> {
    // 如果向量为空，则没有中位数。

    if numbers.is_empty() {
        return None;
    }

    if numbers.is_empty() {
        return None;
    }

    let v = numbers.iter().map(|v| *v as usize).sum::<usize>();
    Some(v as f64 / numbers.len() as f64)
}
