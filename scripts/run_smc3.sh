#!/bin/bash


        
# 执行 docker 命令
docker run -it -v `pwd`:/data --gpus all \
    192.168.3.38:5000/algo/adacus:smc5.5.0_adapter_demux0.0.4_barcode_remover1.0.3_smicing0.4.4_bmi_0.1.4 \
        smc \
        /data/20260114_250701Y0005_Run0002/20260114_250701Y0005_Run0002_adapter.bam \
        /data/20260114_250701Y0005_Run0002/20260114_250701Y0005_Run0002_adapter.hmm.ByStrand \
        --maxLength 10000000 \
        --minLength 10 \
        --hmmModelSpec Kit__500_Chem__1_BC__1_PW3_v4 \
        --byStrand 1 


# 执行 docker 命令
docker run -it -v `pwd`:/data --gpus all \
    192.168.3.38:5000/algo/adacus:smc5.5.0_adapter_demux0.0.4_barcode_remover1.0.3_smicing0.4.4_bmi_0.1.4 \
        smc \
        /data/20260114_250701Y0005_Run0003/20260114_250701Y0005_Run0003_adapter.bam \
        /data/20260114_250701Y0005_Run0003/20260114_250701Y0005_Run0003_adapter.hmm.ByStrand \
        --maxLength 10000000 \
        --minLength 10 \
        --hmmModelSpec Kit__500_Chem__1_BC__1_PW3_v4 \
        --byStrand 1 



docker run -it -v `pwd`:/data --gpus all \
    192.168.3.38:5000/algo/adacus:smc5.5.0_adapter_demux0.0.5_barcode_remover1.0.3_smicing0.5.2_bmi_0.1.5 \
        smc \
        /data/20260122_240901Y0005_Run0001_adapter.bam \
        /data/20260122_240901Y0005_Run0001_adapter.hmm.byStrand \
        --maxLength 10000000 \
        --minLength 10 \
        --hmmModelSpec Kit__500_Chem__1_BC__1_PW3_v4 \
        --byStrand 1 