PY=/root/miniconda3/envs/py38/bin/python     # polars 1.8.2 + seaborn；系统 python3 也可
SCRIPT=/root/projects/gsda/third_party/gseda/src/gseda/ppl/rq_iy_analysis.py
export MPLBACKEND=Agg                        # 否则 matplotlib 找不到 display（只出图不弹窗，但会报错）
D=/data1/ccs_data/str-optimization/fourth-batch-of-data/20260831_250302Y0001_Run0001/withdi-barcode-partitioned
M='/data1/ccs_data/str-optimization/fourth-batch-of-data/STR第四批一代测序/merged_output'

$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.Barcode15.bam           --ref $M/260827STR-93-1.intersect.fa
$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.with_di.q20.Barcode15.bam --ref $M/260827STR-93-1.intersect.fa

$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.Barcode16.bam           --ref $M/260827STR-97-5.intersect.fa
$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.with_di.q20.Barcode16.bam --ref $M/260827STR-97-5.intersect.fa

$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.Barcode17.bam           --ref $M/260827STR-88-9.intersect.fa
$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.with_di.q20.Barcode17.bam --ref $M/260827STR-88-9.intersect.fa

$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.Barcode19.bam           --ref $M/260827STR-95-3.intersect.fa
$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.with_di.q20.Barcode19.bam --ref $M/260827STR-95-3.intersect.fa

$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.Barcode20.bam           --ref $M/260827STR-91-1.intersect.fa
$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.with_di.q20.Barcode20.bam --ref $M/260827STR-91-1.intersect.fa

$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.Barcode23.bam           --ref $M/260827STR-83-3.intersect.fa
$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.with_di.q20.Barcode23.bam --ref $M/260827STR-83-3.intersect.fa

$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.Barcode24.bam           --ref $M/260827STR-82-3.intersect.fa
$PY $SCRIPT --smc-bam $D/20260831_250302Y0001_Run0001_called-demuxed.withdi.smc_all_reads.with_di.q20.Barcode24.bam --ref $M/260827STR-82-3.intersect.fa
