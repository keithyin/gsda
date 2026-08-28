# gsda


```
gsetl -f --outdir ${outdir} aligned-bam --bam ../asts/aligned.bam --ref-file ../asts/test_data/smc.bam
```

## fact_aligned_bam_bam_basic

* qname:
* refname
* channel:
* np:  num passes
* rq: predicted quality [0, 1], extracted from bam file
* iy: aligned identity, extracted from bam file
* ec: query effective coverage, extracted from bam file
* rstart: aligned reference start
* rend: aligned reference end
* qstart: aligned query start (if reverse alignment, this will be the start of the reversed seqeunce)
* qend: aligned query end (if reverse alignment, this will be the end of the reversed seqeunce)
* qlen: query length
* fwd: forward strand?
* ori_start: start position in the basecaller reads (extracting from sbr.bam)
* ori_end:  end position in the basecaller reads (extracting from sbr.bam)

## fact_aligned_bam_record_stat

* qname	
* queryCoverage: query effective coverage, compute from the cigar
* concordance: identity, compute from the cigar
* concordanceQv	: -10 * log10(1-concordance)
* matchBp
* mismatchBp	
* nonHpInsertionBp	
* nonHpDeletionBp	
* hpInsertionBp	
* hpDeletionBp	
* ignoreBp


## fact_aligned_bam_ref_locus_info

* refname	
* pos	
* eq	
* diff	
* ins	
* del	
* depth: how many records that aligned to this position
* curBase	
* nextBase
* curIsHomo
* nextIsHomo	
* aroundBases



## fact_error_query_locus_info

* qname	
* qstart	
* qend	
* rstart	
* rend	
* qseq	
* rseq


## 错误位点分析

```
gsmm2 align -q 20260728_240601Y0004_Run0004_demuxed.boost50.even-odd.smc_all_reads.post.bam -t ../ref_1k.fasta -p boost2ref-1k-rna --np-range 5:100



 asrtc --ref-fa ../ref_1k.fasta    -q 20260728_240601Y0004_Run0004_demuxed.boost50.even-odd.bam    -t 20260728_240601Y0004_Run0004_demuxed.boost50.even-odd.smc_all_reads.post.bam  -p boost50  --np-range 7:100 --oupIyT 0.8 --oupCovT 0.3



gsetl --outdir boost50-gsetl  aligned-bam --ref-file ../ref_1k.fasta --bam boost2ref-1k-rna.bam



python /root/projects/gsda/gseda/src/gseda/align_ana/gsetl_post_process/ref_locus_info_post_process.py  boost50-gsetl/fact_aligned_bam_ref_locus_info.csv


```