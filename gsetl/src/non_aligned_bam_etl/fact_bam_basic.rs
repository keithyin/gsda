use std::{
    fs::{self},
    io::{BufWriter, Write},
    path,
};

use gskits::{
    dna::SEQ_NT4_TABLE,
    gsbam::bam_record_ext::BamRecordExt,
    pbar::{get_spin_pb, DEFAULT_INTERVAL},
};
use rust_htslib::bam::{Read, Reader, Record};

use crate::cli::NonAlignedBamParams;

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
        .unwrap_or(format!("{output_dir}/{filestem}.non_aligned_fact.csv"));

    let bam_threads = param.bam_threads.unwrap_or(num_cpus::get_physical());

    let mut reader = Reader::from_path(bam_file).unwrap();
    reader.set_threads(bam_threads).unwrap();

    let pb = get_spin_pb(format!("reading {bam_file}"), DEFAULT_INTERVAL);
    // ch base base_cnt dw_sum ar_sum cr_mean cq oe
    let o_f = fs::File::create(&oup_filepath).unwrap();
    let mut buf_wirter = BufWriter::new(o_f);
    writeln!(
        &mut buf_wirter,
        "qname\tbase\tbase_cnt\tdw_sum\tar_sum\tcr_mean\tcq\toe"
    )
    .unwrap();
    for record in reader.records() {
        let record = record.unwrap();
        dump_one_record(&record, &mut buf_wirter);
        pb.inc(1);
    }
    pb.finish();
}

/// ch base base_cnt dw_sum ar_sum cr_mean cq oe
fn dump_one_record(record: &Record, writer: &mut dyn Write) {
    let record_ext = BamRecordExt::new(record);
    let bases = record_ext.get_seq();
    let dwell_times = record_ext
        .get_dw()
        .map(|v| v.into_iter().map(|v| v as f64).collect::<Vec<_>>())
        .unwrap_or(vec![-1.0; bases.len()]);
    let arrival_time = record_ext
        .get_ar()
        .map(|v| v.into_iter().map(|v| v as f64).collect::<Vec<_>>())
        .unwrap_or(vec![-1.0; bases.len()]);
    // println!("{:?}", arrival_time.iter().take(20).collect::<Vec<_>>());

    let capture_rates = record_ext
        .get_float_cr()
        .map(|v| v.into_iter().map(|v| v as f64).collect::<Vec<_>>())
        .unwrap_or(vec![-1.0; bases.len()]);

    let cq = record_ext.get_cq().unwrap_or(-1.0) as f64;
    let oe = record_ext.get_oe().unwrap_or(-1.0) as f64;
    let qname = record_ext.get_qname();

    let mut acgt_dw_sum = [0.0; 4];
    let mut acgt_ar_sum = [0.0; 4];
    let mut acgt_cr_mean = [0.0; 4];
    let mut acgt_cnt = [0; 4];

    bases
        .as_bytes()
        .iter()
        .enumerate()
        .for_each(|(pos, &base)| {
            let idx = SEQ_NT4_TABLE[base as usize] as usize;
            // println!("base:{}, idx:{}", base as char, idx);
            let dw = dwell_times[pos];
            let ar = arrival_time[pos];
            let cr = capture_rates[pos];
            acgt_dw_sum[idx] += dw;
            acgt_ar_sum[idx] += ar;
            acgt_cnt[idx] += 1;
            acgt_cr_mean[idx] += cr;
        });

    acgt_cr_mean
        .iter_mut()
        .zip(acgt_cnt.iter())
        .for_each(|(cr, &cnt)| *cr = if cnt == 0 { 0.0 } else { *cr / cnt as f64 });
    // println!("{:?}", acgt_ar_sum);

    // ch base base_cnt dw_sum ar_sum cr_mean cq oe
    ["A", "C", "G", "T"]
        .iter()
        .enumerate()
        .for_each(|(idx, &base)| {
            let base_cnt = acgt_cnt[idx];
            let dw_sum = acgt_dw_sum[idx];
            let ar_sum = acgt_ar_sum[idx];
            let cr_mean = acgt_cr_mean[idx];
            writeln!(
                writer,
                "{qname}\t{base}\t{base_cnt}\t{dw_sum}\t{ar_sum}\t{cr_mean}\t{cq}\t{oe}"
            )
            .unwrap();
        });
}

#[cfg(test)]
mod test {
    use std::io::BufWriter;

    use gskits::gsbam::bam_record_ext::BamRecordExt;
    use rust_htslib::bam::{Read, Reader};

    use crate::non_aligned_bam_etl::fact_bam_basic::dump_one_record;

    #[test]
    fn test_stat() {
        let bam_file =
            "/data/ccs_data/data2025Q1/S_aureus_1h/20250207_Sync_Y0002_02_H01_Run0001_called.bam";
        let mut reader = Reader::from_path(bam_file).unwrap();
        reader.set_threads(10).unwrap();
        for record in reader.records() {
            let record = record.unwrap();
            let record_ext = BamRecordExt::new(&record);
            // println!("{}", record_ext.get_seq());

            let mut buf = vec![];
            let mut buf_writer = BufWriter::new(&mut buf);
            dump_one_record(&record, &mut buf_writer);
            drop(buf_writer);
            println!("{}", String::from_utf8(buf).unwrap());
            break;
        }
    }
}
