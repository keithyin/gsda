# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa \
#     --s_bam /data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q4/align-param-search/20240607_h84.icing.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta \
#     --s_bam /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q3/jinpu/smc-align-param-search/20240711_Sync_Y0006_02_H01_Run0001_called_m3_mm5_ins2_del2.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/align-param-search/20240711sbr2align-modif-smc-epo3-new.icing.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/ConsensusDcDw/20240711-smc3522.icing.bam \
#     --sbr_qv_thr 10


# /data/ccs_data/ccs_eval2024q3/jinpu/align-param-search/20240711sbr2align-modif-smc-epo3-new.icing.bam \
# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta \
#     --s_bam /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q3/jinpu/smc-align-param-search/20240711_Sync_Y0006_02_H01_Run0001_called_m3_mm5_ins2_del2.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/hg-ecoli-jinpu-np3-15-smc3522--ConsensusModelBiGRUDcDwArCr/20240711-smc3522.icing.bam \
#     --sbr_qv_thr 10
    
# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta \
#     --s_bam /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q3/jinpu/smc-align-param-search/20240711_Sync_Y0006_02_H01_Run0001_called_m3_mm5_ins2_del2.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/jinpu-test/20240711-baseline.icing.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/jinpu-test/20240711-virtual-softmax.icing.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/jinpu-test/20240711-dwarcr.icing.bam \
#     --sbr_qv_thr 10
#      /data/ccs_data/ccs_eval2024q3/jinpu/jinpu-test/20240711-ghm.icing.bam \

python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta \
    --s_bam /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.bam \
    --smc_bams \
        /data/ccs_data/ccs_eval2024q3/jinpu/smc-align-param-search/20240711_Sync_Y0006_02_H01_Run0001_called_m3_mm5_ins2_del2.smc_all_reads.bam \
        /data/ccs_data/smc-poa-3-5-2-2-noPolish/noPolish-icing/20240711-nopolish.icing.bam \
    --sbr_qv_thr 10

# /data/ccs_data/ccs_eval2024q3/jinpu/hg-ecoli-jinpu-np3-15-ConsensusModelBiGRU/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.icing.bam \

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.fasta  \
#     --s_bam /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.040icing.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta  \
#     --s_bam /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.align-param.smc.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/align-param-search/20240711sbr2align-modif-smc-epo3-new.icing.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/align-param-search/20240711sbr2align-modif-smc-epo3-new2.icing.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta  \
#     --s_bam /data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q3/jinpu/smc-align-param-search/20240711_Sync_Y0006_02_H01_Run0001_called_m3_mm5_ins2_del2.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q3/jinpu/model-v5/20240711_Sync_Y0006_02_H01_Run0001_called.smc_all_reads.bam \
#     --sbr_qv_thr 10

#        /data/ccs_data/ccs_eval2024q3/jinpu/align-param-search/20240711sbr2align-modif-smc-epo3.icing.bam \

# /data/ccs_data/ccs_eval2024q3/jinpu/align-param-search/20240711sbr2align-modif-smc-epo3.icing.bam \

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/kaipu/20241015_new/GRCh38_11_16.fa  \
#     --s_bam /data/ccs_data/kaipu/20241015_new/k2_5_subreads.bam \
#     --smc_bams \
#         /data/ccs_data/kaipu/20241015_new/k2_5_subreads.smc_all_reads.bam \
#         /data/ccs_data/kaipu/20241015_new/k2_5.icing.bam \
#     --sbr_qv_thr 10


# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/case-study/tr-at-error-kaipu/chr11_6000_ref-3002-G2A-2989InsAT.fa  \
#     --s_bam /data/ccs_data/case-study/tr-at-error-kaipu/k2_5_subreads.subset.bam \
#     --smc_bams \
#         /data/ccs_data/case-study/tr-at-error-kaipu/k2_5_subreads.subset.smc_all_reads.bam \
#         /data/ccs_data/case-study/tr-at-error-kaipu/k2_5_sbr.subset.icing.bam \
#     --sbr_qv_thr 10


# /data/ccs_data/ccs_eval2024q4/baseline-reexport-epo9/20240607_Sync_H84_01_H01_Run0002_adapter.icing.bam \
# /data/ccs_data/ccs_eval2024q4/baseline-reexport-epo9/20240607_Sync_H84_01_H01_Run0002_adapter.icing.bam \

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa \
#     --s_bam /data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q4/icing041/result.icing.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa \
#     --s_bam /data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0003_adapter.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0003_adapter.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q4/baseline-all-sbrs/20240607_Sync_H84_01_H01_Run0003_adapter.icing.bam \
#     --sbr_qv_thr 10


# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa \
#     --s_bam /data/ccs_data/ccs_eval2024q4/20240608_Sync_H84_01_H01_Run0001_adapter.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q4/20240608_Sync_H84_01_H01_Run0001_adapter.smc_all_reads.bam \
#         /data/ccs_data/ccs_eval2024q4/baseline-all-sbrs/20240608_Sync_H84_01_H01_Run0001_adapter.icing.bam \
#     --sbr_qv_thr 10



# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa  \
#     --s_bam /data/ccs_data/smc-upgrade/20240930_240601Y0009_Run0001_called.adapter.bam \
#     --smc_bams \
#         /data/ccs_data/smc-upgrade/exp/20240930_240601Y0009_Run0001_called.baseline.smc_all_reads.bam \
#         /data/ccs_data/smc-upgrade/exp/20240930_240601Y0009_Run0001_called.smc_all_reads.bam \
#     --sbr_qv_thr 10


# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa  \
#     --s_bam /data/ccs_data/smc-upgrade/20240930_240601Y0009_Run0001_called_dbscan.adapter.bam \
#     --smc_bams \
#         /data/ccs_data/smc-upgrade/exp/20240930_240601Y0009_Run0001_called_dbscan.baseline.smc_all_reads.bam \
#         /data/ccs_data/smc-upgrade/exp/20240930_240601Y0009_Run0001_called_dbscan.smc_all_reads.bam \
#     --sbr_qv_thr 10

# python subreads_smc_funnel_stat_speed_up.py /data/ccs_data/MG1655.fa \
#     --s_bam /data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter.bam \
#     --smc_bams \
#         /data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter.smc_all_reads.bam \
#     --sbr_qv_thr 10

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