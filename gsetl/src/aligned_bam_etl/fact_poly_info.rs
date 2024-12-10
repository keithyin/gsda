use std::{
    fmt::Display,
    fs,
    io::{BufWriter, Write},
};

use gskits::{
    file_reader::{bed_reader::BedInfo, vcf_reader::VcfInfo},
    gsbam::bam_record_ext::BamRecordExt,
    pbar,
};
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read, Record};

use crate::{
    cli::AlignedBamParams,
    poly_n::{find_homopolymer_regions, position_relation, PosRelation},
};

use super::FastaData;

#[derive(Debug, PartialEq, Eq)]
struct RefPolyLocusInfo {
    rstart: usize,
    rend: usize,
    qstart: usize,
    qend: usize,
    ref_base: char,
    ref_repeats: usize,
    query_repeats: usize,
    qseq: String,
    query_clean: bool,
}

impl RefPolyLocusInfo {
    fn new(
        rstart: usize,
        rend: usize,
        qstart: usize,
        qend: usize,
        ref_base: char,
        ref_repeats: usize,
        qseq: String,
    ) -> Self {
        let query_repeats = qseq
            .as_bytes()
            .iter()
            .filter(|query_base| **query_base == ref_base as u8)
            .count();
        let query_clean = query_repeats == qseq.len();
        Self {
            rstart,
            rend,
            qstart,
            qend,
            ref_base,
            ref_repeats,
            query_repeats,
            qseq,
            query_clean,
        }
    }

    fn csv_header() -> &'static str {
        "rStart\trEnd\tqStart\tqEnd\trBase\trRepeats\tqRepeats\tqSeq\tqClean"
    }
}

impl Display for RefPolyLocusInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}",
            self.rstart, self.rend, self.qstart, self.qend, self.ref_base, self.ref_repeats, self.query_repeats, self.qseq, self.query_clean
        )
    }
}

/// qname refname qpos_start qpos_end rpos_start, rpos_end query ref
pub fn fact_poly_info(
    args: &AlignedBamParams,
    output_dir: &str,
    _hc_regions: Option<&BedInfo>,
    _hc_variants: Option<&VcfInfo>,
    fasta_data: &FastaData,
) {
    let bam_file = &args.bam;

    let mut bam_h = bam::IndexedReader::from_path(bam_file).unwrap();
    bam_h.set_threads(4).unwrap();

    let o_filepath = format!("{}/fact_poly_info.csv", output_dir);
    let o_file = fs::File::create(o_filepath).unwrap();
    let mut o_file_buff_writer = BufWriter::new(o_file);
    writeln!(
        &mut o_file_buff_writer,
        "rName\tqName\t{}",
        RefPolyLocusInfo::csv_header()
    )
    .unwrap();

    let pb = pbar::get_spin_pb(
        format!("fact_error_query_locus_info"),
        pbar::DEFAULT_INTERVAL,
    );

    for (refname, refseq) in fasta_data.get_ref_name2seq() {
        bam_h
            .fetch((refname, 0, refseq.len() as u64))
            .expect(&format!("fetch {} error", refname));

        let refseq_bytes = refseq.as_bytes();

        let ref_poly_regsion = find_homopolymer_regions(refseq_bytes);

        for record in bam_h.records() {
            pb.inc(1);
            let record = record.unwrap();
            let qname = BamRecordExt::new(&record).get_qname();
            let poly_infos = poly_info_extractor(&record, &ref_poly_regsion);
            if let Some(poly_infos) = poly_infos {
                poly_infos.into_iter()
                .for_each(|poly_info| {
                    writeln!(&mut o_file_buff_writer, "{}\t{}\t{}", refname, qname, poly_info).unwrap();
                });
            }
                
        }
    }
    pb.finish();
}

fn poly_info_extractor(
    record: &Record,
    ref_poly_region: &Vec<(usize, usize, u8)>,
) -> Option<Vec<RefPolyLocusInfo>> {
    // if record.is_secondary() || record.is_unmapped() || record.is_supplementary() {
    //     continue;
    // }
    let record_ext = BamRecordExt::new(&record);

    let ref_start = record_ext.reference_start() as i64;
    let ref_end = record_ext.reference_end() as i64;

    let query_end = record_ext.query_alignment_end() as i64;

    let mut rpos_cursor = None;
    let mut qpos_cursor = None;

    let query_seq = record.seq().as_bytes();
    let query_str = unsafe { String::from_utf8_unchecked(query_seq.clone()) };

    let mut poly_info = vec![];

    let aligned_pairs_full: Vec<[Option<i64>; 2]> =
        record.aligned_pairs_full().into_iter().collect::<Vec<_>>();

    let mut cur_poly_region_idx =
        if let Some(poly_idx) = move_poly_idx(ref_poly_region, 0, ref_start as usize) {
            if (ref_start as usize) > ref_poly_region[poly_idx].0 {
                poly_idx + 1
            } else {
                poly_idx
            }
        } else {
            return None;
        };

    if cur_poly_region_idx >= ref_poly_region.len() {
        return None;
    }

    let mut query_poly_seq = String::new();

    for [qpos, rpos] in aligned_pairs_full.iter() {
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

        if let Some(qpos_cursor_) = qpos_cursor {
            if qpos_cursor_ >= query_end {
                break;
            }
        } else {
            continue;
        }

        // 收尾
        if let Some(rpos_) = rpos.map(|v| v as usize) {
            let cur_poly_region = unsafe { ref_poly_region.get_unchecked(cur_poly_region_idx) };
            if rpos_ == (cur_poly_region.1 - 1) {
                if let Some(qpos_) = qpos.map(|v| v as usize) {
                    query_poly_seq.push_str(&query_str[qpos_..qpos_ + 1]);
                }

                // poly region end
                let qend = qpos_cursor.unwrap() as usize + 1;
                let qstart = qend - query_poly_seq.len();
                poly_info.push(RefPolyLocusInfo::new(
                    cur_poly_region.0,
                    cur_poly_region.1,
                    qstart,
                    qend,
                    cur_poly_region.2 as char,
                    cur_poly_region.1 - cur_poly_region.0,
                    query_poly_seq,
                ));

                // move
                cur_poly_region_idx += 1;
                if cur_poly_region_idx >= ref_poly_region.len() {
                    break;
                }

                query_poly_seq = String::new();
            }
        }

        let cur_region = unsafe { ref_poly_region.get_unchecked(cur_poly_region_idx) };
        // region 开始 & 持续
        let rpos_cur_or_pre = rpos_cursor.unwrap() as usize;

        if let Some(rpos_) = rpos.map(|v| v as usize) {
            if rpos_ < cur_region.0 {
                continue;
            }
            // rpos_ in region
            if let Some(qpos_) = qpos.map(|v| v as usize) {
                query_poly_seq.push_str(&query_str[qpos_..qpos_ + 1]);
            }
        } else {
            if (rpos_cur_or_pre + 1) < cur_region.0 {
                continue;
            }
            if let Some(qpos_) = qpos.map(|v| v as usize) {
                query_poly_seq.push_str(&query_str[qpos_..qpos_ + 1]);
            }
        }
    }

    Some(poly_info)
}

fn move_poly_idx(
    ref_poly_region: &Vec<(usize, usize, u8)>,
    cur_idx: usize,
    cur_pos: usize,
) -> Option<usize> {
    for idx in cur_idx..ref_poly_region.len() {
        match position_relation(&ref_poly_region[cur_idx], cur_pos) {
            PosRelation::Left => return Some(idx),
            PosRelation::Middle => return Some(idx),
            PosRelation::Right => (),
        }
    }

    return None;
}

#[cfg(test)]
mod test {
    use crate::{aligned_bam_etl::fact_poly_info::RefPolyLocusInfo, poly_n::find_homopolymer_regions};

    use super::poly_info_extractor;


    #[test]
    fn test_poly_info_extractor() {
        /*
           ref: ACCGGT
           qry: A-C--T
        */
        let refseq = b"ACCGGT";
        let mut record = rust_htslib::bam::Record::new();
        let seq = "ACT";
        record.set_pos(0);

        record.set(
            b"qname0",
            Some(&gskits::gsbam::cigar_ext::parse_cigar_string("1=1D1=2D1=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        let ref_poly_region = find_homopolymer_regions(refseq);
        let res = poly_info_extractor(&record, &ref_poly_region).unwrap();
        assert_eq!(res[0], RefPolyLocusInfo::new(1, 3, 1, 2, 'C', 2, "C".to_string()));
        assert_eq!(res[1], RefPolyLocusInfo::new(3, 5, 2, 2, 'G', 2, "".to_string()));

        
        /*
           ref: ACCGGT
           qry: A-C-GT
        */
        let refseq = b"ACCGGT";
        let mut record = rust_htslib::bam::Record::new();
        let seq = "ACGT";
        record.set_pos(0);

        record.set(
            b"qname0",
            Some(&gskits::gsbam::cigar_ext::parse_cigar_string("1=1D1=1D2=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        let ref_poly_region = find_homopolymer_regions(refseq);
        let res = poly_info_extractor(&record, &ref_poly_region).unwrap();
        assert_eq!(res[0], RefPolyLocusInfo::new(1, 3, 1, 2, 'C', 2, "C".to_string()));
        assert_eq!(res[1], RefPolyLocusInfo::new(3, 5, 2, 3, 'G', 2, "G".to_string()));


        /*
           ref: A--CCGGT
           qry: ACCCC-GT
        */
        let refseq = b"ACCGGT";
        let mut record = rust_htslib::bam::Record::new();
        let seq = "ACCCCGT";
        record.set_pos(0);

        record.set(
            b"qname0",
            Some(&gskits::gsbam::cigar_ext::parse_cigar_string("1=2I2=1D2=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        let ref_poly_region = find_homopolymer_regions(refseq);
        let res = poly_info_extractor(&record, &ref_poly_region).unwrap();
        assert_eq!(res[0], RefPolyLocusInfo::new(1, 3, 1, 5, 'C', 2, "CCCC".to_string()));
        assert_eq!(res[1], RefPolyLocusInfo::new(3, 5, 5, 6, 'G', 2, "G".to_string()));


        /*
           ref: A-TT
           qry: ATTT
        */
        let refseq = b"ATT";
        let mut record = rust_htslib::bam::Record::new();
        let seq = "ATTT";
        record.set_pos(0);

        record.set(
            b"qname0",
            Some(&gskits::gsbam::cigar_ext::parse_cigar_string("1=1I2=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        let ref_poly_region = find_homopolymer_regions(refseq);
        let res = poly_info_extractor(&record, &ref_poly_region).unwrap();
        assert_eq!(res[0], RefPolyLocusInfo::new(1, 3, 1, 4, 'T', 2, "TTT".to_string()));

        /*
           ref: AAA-TT
           qry:   ATTT
        */
        let refseq = b"AAATT";
        let mut record = rust_htslib::bam::Record::new();
        let seq = "ATTT";
        record.set_pos(2);

        record.set(
            b"qname0",
            Some(&gskits::gsbam::cigar_ext::parse_cigar_string("1=1I2=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        let ref_poly_region = find_homopolymer_regions(refseq);
        let res = poly_info_extractor(&record, &ref_poly_region).unwrap();
        assert_eq!(res[0], RefPolyLocusInfo::new(3, 5, 1, 4, 'T', 2, "TTT".to_string()));


    }
}
