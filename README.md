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

