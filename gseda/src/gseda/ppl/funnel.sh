# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa  \
#     --s_bam /data/four/run1/test_bam/20240929_240601Y0009_Run0001_called_ori.adapter.bam \
#     --smc_bams \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called_ori.smc_all_reads.bam \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called_ori3522.smc_all_reads.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa  \
#     --s_bam /data/four/run1/test_bam/20240929_240601Y0009_Run0001_called.adapter.bam \
#     --smc_bams \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called.smc_all_reads.bam \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called3522.smc_all_reads.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa  \
#     --s_bam /data/four/run1/test_bam/20240929_240601Y0009_Run0001_called_dbscan.adapter.bam \
#     --smc_bams \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called_dbscan.smc_all_reads.bam \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called_dbscan3522.smc_all_reads.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa  \
#     --s_bam /data/four/run1/test_bam/20240929_240601Y0009_Run0001_called_ori.adapter.bam \
#     --smc_bams \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called_ori3522.smc_all_reads.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa  \
#     --s_bam /data/four/run1/test_bam/20240929_240601Y0009_Run0001_called.adapter.bam \
#     --smc_bams \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called3522.smc_all_reads.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa  \
#     --s_bam /data/four/run1/test_bam/20240929_240601Y0009_Run0001_called_dbscan.adapter.bam \
#     --smc_bams \
#         /data/four/run1/test_bam/exp/20240929_240601Y0009_Run0001_called_dbscan3522.smc_all_reads.bam \
#     --sbr_qv_thr 10



# python subreads_smc_funnel_stat.py /data1/ccs_data/202603-rna-data/rna_data/RNA1K/RNA1k-nopoly-ccaa.fa  \
#     --s_bam /data1/ccs_data/202603-rna-data/rna_data/RNA1K/20260327_240601Y0004_Run0003_demuxed.bam \
#     --smc_bams \
#         /data1/ccs_data/202603-rna-data/rna_data/RNA1K/20260327_240601Y0004_Run0003_demuxed.0521.smc_all_reads.bam \
#     --sbr_qv_thr 10


# python subreads_smc_funnel_stat.py /data1/ccs_data/202603-rna-data/rna_data/RNA2K/RNA2K-nopoly-ccaa.fa  \
#     --s_bam /data1/ccs_data/202603-rna-data/rna_data/RNA2K/20260327_240601Y0004_Run0005_demuxed.bam \
#     --smc_bams \
#         /data1/ccs_data/202603-rna-data/rna_data/RNA2K/20260327_240601Y0004_Run0005_demuxed.0521.smc_all_reads.bam \
#     --sbr_qv_thr 10




# jinpu
# python subreads_smc_funnel_stat.py /data1/REF_GENOMES/ref_Saureus_ATCC25923.m.new.corrected.fasta  \
#     --s_bam /data1/ccs_data/2025Q1/eval-data/20250207_Sync_Y0003_08_H01_Run0001_called.adapter.bam \
#     --smc_bams \
#         /data1/ccs_data/2025Q1/eval-data/20250207_Sync_Y0003_08_H01_Run0001_called.adapter.baseline.smc_all_reads.bam \
#         /data1/ccs_data/2025Q1/eval-data/20250207_Sync_Y0003_08_H01_Run0001_called.adapter.gap-left-align.smc_all_reads.bam \
#     --sbr_qv_thr 10


# ecoli /data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter.gap-left-align
python subreads_smc_funnel_stat.py /data1/REF_GENOMES/MG1655.fa  \
    --s_bam /data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter.bam \
    --smc_bams \
        /data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter.baseline.smc_all_reads.bam \
        /data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter.gap-left-align.smc_all_reads.bam \
    --sbr_qv_thr 10