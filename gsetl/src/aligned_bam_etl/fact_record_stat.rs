use std::{collections::HashMap, fs, io::{ BufWriter, Write}, sync::Arc, thread};

use gskits::{ file_reader::{bed_reader::BedInfo, vcf_reader::VcfInfo}, gsbam::bam_record_ext::BamRecordExt, pbar};
use indicatif::ProgressBar;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Cigar, Read};

use crate::{cli::AlignedBamParams, set_spin_pb};

use super::{audit, FastaData};



struct RecordReplica {
    qname: String,
    refname: Arc<String>,
    ch: u32,
    q_len: usize,
    passes: u32,
    rq: Option<f32>,
    cigars: Vec<bam::record::Cigar>,
    seq: String,
    align_ref_start: usize,
    align_ref_end: usize
}

impl RecordReplica {
    pub fn from_record(record: &bam::Record, refname: Arc<String>) -> Self {
        let cigars = record.cigar().take().0;
        let record_ext = BamRecordExt::new(record);

        Self {
            qname: record_ext.get_qname(),
            refname: refname,
            ch: record_ext.get_ch().unwrap(), 
            q_len: record_ext.get_seq().len(), 
            passes: record_ext.get_np().unwrap_or(1), 
            rq: record_ext.get_rq(), 
            cigars: cigars,
            seq: record_ext.get_seq(),
            align_ref_start: record.reference_start() as usize,
            align_ref_end: record.reference_end() as usize
        }
    }

    pub fn get_ch(&self) -> u32 {
        self.ch
    }

    pub fn get_seq_len(&self) -> usize {
        self.q_len
    }

    pub fn get_np(&self) -> u32 {
        self.passes
    }

    pub fn get_rq(&self) -> Option<f32> {
        self.rq
    }

    pub fn get_cigars(&self) -> &Vec<bam::record::Cigar> {
        &self.cigars
    }

    pub fn get_seq(&self) -> &str {
        &self.seq
    }

    pub fn get_aligned_ref_start(&self) -> usize {
        self.align_ref_start
    }

    pub fn get_aligned_ref_end(&self) -> usize {
        self.align_ref_end
    }
}


#[allow(unused)]
#[derive(Debug)]
struct Stat {
    qname: String,
    refname: Arc<String>,
    ch: u32,
    passes: u32,
    q_len: usize,
    m: usize,
    mm: usize,
    non_h_del: usize,
    h_del: usize,
    non_h_ins: usize,
    h_ins: usize,
    rq: Option<f32>,
    ignore_bps: usize,
    ref_start: usize,
    ref_end: usize,
    mm_ref_positions: Vec<usize>,
    ins_ref_positions: Vec<usize>,
    del_ref_positions: Vec<usize>
}

impl Stat {
    fn new(qname: String, refname: Arc<String>, ch: u32, q_len: usize, passes: u32, ref_start: usize, ref_end: usize, rq: Option<f32>) -> Self{
        Self { 
            qname: qname,
            refname: refname,
            ch: ch, 
            passes: passes,
            q_len: q_len, 
            m: 0, 
            mm: 0, 
            non_h_del: 0, 
            h_del: 0, 
            non_h_ins: 0, 
            h_ins: 0,
            rq: rq,
            ignore_bps: 0,
            ref_start: ref_start,
            ref_end: ref_end,
            mm_ref_positions: vec![],
            ins_ref_positions: vec![],
            del_ref_positions: vec![]
        }
    }
    
    fn ins_bp(&self) -> usize {
        self.non_h_ins + self.h_ins
    }

    fn del_bp(&self) -> usize {
        self.h_del + self.non_h_del
    }

    fn align_span(&self) -> usize {
        self.ins_bp() + self.del_bp() + self.m + self.mm
    }

    fn concordance(&self) -> f32 {
        (self.m as f32) / self.align_span() as f32
    }

    // qv is capped by 60
    fn concordance_qv(&self) -> f32 {
        let mut err_rate = 1.0 - self.concordance();
        err_rate = if err_rate < 1e-6 {1e-6} else {err_rate};
        -10. * err_rate.log10()
    }

    // qv coverage
    fn query_converage(&self) -> f32 {
        (self.ins_bp() + self.m + self.mm) as f32 / self.q_len as f32
    }
}

fn do_m_mm_stat(counter: &mut usize, start_pos: usize, n: usize, ref_name: &str, hc_regions: Option<&BedInfo>, hc_variants: Option<&VcfInfo>) -> (usize, Vec<usize>){
    let mut positions = vec![];
    let mut ignore_bps = 0_usize;
    if hc_regions.is_none() && hc_variants.is_none(){
        *counter += n as usize;
    } else {
        for shift in 0..n {
            let cur_pos = start_pos + shift;
            if valid_point(ref_name, cur_pos, hc_regions, hc_variants) {
                *counter += 1;
                positions.push(cur_pos);
            } else {
                ignore_bps += 1;
            }
        }
    }
    (ignore_bps, positions)
}


fn valid_insertion(ref_pos: usize, ref_name: &str, hc_regions: Option<&BedInfo>, hc_variants: Option<&VcfInfo>) -> bool {
    let mut valid = true;
    if let Some(hc_regions) = hc_regions {
        valid &= hc_regions.point_within_region(ref_name, ref_pos);
    }

    if let Some(hc_variants) = hc_variants {
        valid &= !hc_variants.point_hit(ref_name, ref_pos);
    }

    if ref_pos > 0 {
        if let Some(hc_regions) = hc_regions {
            valid &= hc_regions.point_within_region(ref_name, ref_pos - 1);
        }
    
        if let Some(hc_variants) = hc_variants {
            valid &= !hc_variants.point_hit(ref_name, ref_pos - 1);
        }
    }

    valid
}


fn valid_point(ref_name: &str, position: usize, hc_regions: Option<&BedInfo>, hc_variants: Option<&VcfInfo>) -> bool {

    let mut valid = true;
    if let Some(hc_regions) = hc_regions {
        valid &= hc_regions.point_within_region(ref_name, position);
    }

    if let Some(hc_variants) = hc_variants {
        valid &= !hc_variants.point_hit(ref_name, position);
    }
    
    valid  

}


fn stat_record_core(record: RecordReplica, references: &HashMap<String, String>, hc_regions: Option<&BedInfo>, hc_variants: Option<&VcfInfo>) -> Stat {
    let ref_name = record.refname.as_str();
    let ref_seq = references.get(ref_name).expect(&format!("refname:'{}' not found. valid refnames are:'{:?}'", ref_name, references.keys()));
    let query_seq = record.get_seq();

    let mut rpos_cursor = record.get_aligned_ref_start();
    let mut qpos_cursor = 0_usize;
    let mut stat = Stat::new(
            record.qname.clone(),
            record.refname.clone(),
        record.get_ch(), 
        record.get_seq_len(), 
        record.get_np(),
        record.get_aligned_ref_start(), 
        record.get_aligned_ref_end(), 
        record.get_rq());
    // println!("cigar:{:?}", record.get_cigars());
    for cig in record.get_cigars(){
        match cig {
            &Cigar::Equal(n) => {
                let (ignore_bps, _) = do_m_mm_stat(&mut stat.m, rpos_cursor, n as usize, ref_name, hc_regions, hc_variants);
                stat.ignore_bps += ignore_bps;
            },
            &Cigar::Diff(n) => {
                let (ignore_bps, mm_ref_pos) = do_m_mm_stat(&mut stat.mm, rpos_cursor, n as usize, ref_name, hc_regions, hc_variants);
                stat.ignore_bps += ignore_bps;
                stat.mm_ref_positions.extend(mm_ref_pos.into_iter());
            },
            &Cigar::Match(_) => panic!("Match not valid"),
            _ => (),
        }

        match cig {
            &Cigar::Ins(n) => {
                if valid_insertion(rpos_cursor, ref_name, hc_regions, hc_variants) {
                    let n = n as usize;
                    let mut is_hp = false;
                    let cur_bp_seq = &query_seq[qpos_cursor..qpos_cursor+n];
                    if rpos_cursor > 0 {
                        let cur_rseq = vec![&ref_seq.as_str()[rpos_cursor-1..rpos_cursor]; n].join("");
                        assert_eq!(cur_bp_seq.len(), cur_rseq.len());
                        is_hp |= cur_rseq.as_str() == cur_bp_seq;
                    }
                    if rpos_cursor < ref_seq.len() {
                        let cur_rseq = vec![&ref_seq.as_str()[rpos_cursor..rpos_cursor+1]; n].join("");
                        assert_eq!(cur_bp_seq.len(), cur_rseq.len());
                        is_hp |= cur_rseq.as_str() == cur_bp_seq;
                    }

                    if is_hp {
                        stat.h_ins += n;
                    } else {
                        stat.non_h_ins += n;
                    }
                    stat.ins_ref_positions.push(rpos_cursor);

                } else {
                    stat.ignore_bps += n as usize;
                }
            },
            &Cigar::Del(n) => {
                let n = n as usize;
                for shift in 0..n {
                    let cur_r_pos = rpos_cursor + shift;
                    if valid_point(ref_name, cur_r_pos, hc_regions, hc_variants) {
                        let mut is_hp = false;
                        let cur_ref_bp = ref_seq.as_bytes()[cur_r_pos];
                        if cur_r_pos > 0 {
                            is_hp |= cur_ref_bp == ref_seq.as_bytes()[cur_r_pos-1];
                        }
                        if cur_r_pos < ref_seq.len() {
                            is_hp |= cur_ref_bp == ref_seq.as_bytes()[cur_r_pos + 1];
                        }
                        if is_hp {
                            stat.h_del += 1;
                        } else {
                            stat.non_h_del += 1;
                        }
                        stat.del_ref_positions.push(cur_r_pos);
                    } else {
                        stat.ignore_bps += 1;
                    }
                }
            },

            _ => (),
        }

        match cig {
            &Cigar::Equal(n) | &Cigar::Diff(n)=> {
                rpos_cursor += n as usize;
                qpos_cursor += n as usize;
            },
            &Cigar::Del(n) | &Cigar::RefSkip(n)=> rpos_cursor += n as usize,
            &Cigar::Ins(n) | &Cigar::SoftClip(n) => qpos_cursor += n as usize,
            _ => panic!("? cigar:{:?}. {}", cig, cig),
        }
    }

    stat
}


pub fn fact_record_stat(args: &AlignedBamParams, output_dir: &str, hc_regions: Option<&BedInfo>, hc_variants: Option<&VcfInfo>, fasta_data: &FastaData, pbar: ProgressBar) {
    thread::scope(|thread_scope| {

        let (record_sender, record_receiver) = crossbeam::channel::bounded(200);
        thread_scope.spawn(move|| {
            let mut tid2refname = HashMap::new();
            let mut reader = bam::Reader::from_path(&args.bam).unwrap();
            reader.set_threads(10).unwrap();
            let header = bam::Header::from_template(reader.header());
            let header_view = bam::HeaderView::from_header(&header);
            for record in reader.records() {
                let record = record.unwrap();
                if !audit(&record, args) {
                    continue;
                }
                
                let tid = record.tid() as u32;
                if !tid2refname.contains_key(&tid) {
                    tid2refname.insert(tid, Arc::new(String::from_utf8(header_view.tid2name(tid).to_vec()).unwrap()));
                }

                record_sender.send(
                    RecordReplica::from_record(&record, tid2refname.get(&tid).unwrap().clone())
                    ).unwrap();
            }
        });

        let stat_threads = 8 as usize;
        let (stat_sender, stat_receiver) = crossbeam::channel::bounded(200);
        for _ in 0..stat_threads {
            let record_receiver_ = record_receiver.clone();
            let stat_sender_ = stat_sender.clone();
            thread_scope.spawn(move|| {
                for record in record_receiver_ {
                    let stat = stat_record_core(record, fasta_data.get_ref_name2seq(), hc_regions, hc_variants);
                    stat_sender_.send(stat).unwrap();
                }
            });
        }
        drop(record_receiver);
        drop(stat_sender);

        // write result
        let csv_filepath = format!("{}/fact_aligned_bam_record_stat.csv", output_dir);
        let csv_writer = fs::File::create(&csv_filepath).unwrap();
        let mut buf_csv_writer = BufWriter::new(csv_writer);
        let pbar = set_spin_pb(pbar, "fact_record_stat".to_string(), pbar::DEFAULT_INTERVAL);

        // writeln!(&mut buf_csv_writer, 
        //     "channel_id\treadLengthBp\tsubreadPasses\t\
        //     queryConverage\tpredictedConcordance\t\
        //     concordance\tconcordanceQv\tmatchBp\tmismatchBp\tnonHpInsertionBp\t\
        //     nonHpDeletionBp\thpInsertionBp\thpDeletionBp\tignoreBp\t\
        //     refStart\trefEnd\tmmRefPositions\tinsRefPositions\tdelRefPositions\tqname\trefname").unwrap();
        // for stat in stat_receiver {
        //     writeln!(&mut buf_csv_writer, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t\"{16}\"\t\"{17}\"\t\"{18}\"\t{19}\t{20}", 
        //         stat.ch, stat.q_len, stat.passes, 
        //         stat.query_converage(), stat.rq.unwrap_or(0.), stat.concordance(), 
        //         stat.concordance_qv(), stat.m, stat.mm, stat.non_h_ins, stat.non_h_del, 
        //         stat.h_ins, stat.h_del, stat.ignore_bps,
        //         stat.ref_start, stat.ref_end, stat.mm_ref_positions.iter().map(|v| v.to_string()).collect::<Vec<String>>().join(","),
        //         stat.ins_ref_positions.iter().map(|v| v.to_string()).collect::<Vec<String>>().join(","),
        //         stat.del_ref_positions.iter().map(|v| v.to_string()).collect::<Vec<String>>().join(","),
        //         stat.qname,
        //         stat.refname.as_str()
        //     ).unwrap();
        //     pbar.inc(1);
        // }

        writeln!(&mut buf_csv_writer, 
            "qname\tqueryCoverage\tconcordance\tconcordanceQv\tmatchBp\tmismatchBp\tnonHpInsertionBp\t\
            nonHpDeletionBp\thpInsertionBp\thpDeletionBp\tignoreBp").unwrap();
        for stat in stat_receiver {
            writeln!(&mut buf_csv_writer, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}", 
                stat.qname,
                stat.query_converage(), stat.concordance(), 
                stat.concordance_qv(), stat.m, stat.mm, stat.non_h_ins, stat.non_h_del, 
                stat.h_ins, stat.h_del, stat.ignore_bps,                
            ).unwrap();
            pbar.inc(1);
        }
        pbar.finish();

    });
}

