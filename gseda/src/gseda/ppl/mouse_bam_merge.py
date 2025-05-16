
"""
将多个 adapter.bam 和 smc_all_reads.bam 文件合并成一个. dw 和 ar 特征需要从 called bam 中进行透传
"""

def main():
    adapter_bams = [
        "/data/ccs_data/little-mouse/20250515_240601Y0005_Run0001_called.bam",
        "/data/ccs_data/little-mouse/20250515_240601Y0005_Run0002_called.bam"
    ]
    
    sbr_bams = [
        "/data/ccs_data/little-mouse/20250515_240601Y0005_Run0001_adapter.bam",
        "/data/ccs_data/little-mouse/20250515_240601Y0005_Run0002_adapter.bam"
    ]
    
    smc_bams = [
        "/data/ccs_data/little-mouse/20250515_240601Y0005_Run0001-rerun.smc_all_reads.bam",
        "/data/ccs_data/little-mouse/20250515_240601Y0005_Run0002-rerun.smc_all_reads.bam"
    ]
    
    
    pass

if __name__ == "__main__":
    pass
