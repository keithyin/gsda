"""Configuration for GSEDA Server"""

import os
from typing import List


class Settings:
    """Application settings"""

    # Application name
    APP_NAME = "GSEDA - Genomic Data Analysis Tools"
    APP_VERSION = "1.0.0"

    # API configuration
    API_PREFIX = "/api"
    DEBUG = os.getenv("DEBUG", "False").lower() == "true"

    # CLI tools configuration
    GSADA_PACKAGE = "gseda"
    CLI_TOOL_TIMEOUT = 3600  # 1 hour in seconds

    # Paths
    SERVER_ROOT = os.path.dirname(__file__)
    PROJECT_ROOT = os.path.dirname(os.path.dirname(SERVER_ROOT))

    # Static files
    STATIC_DIR = os.path.join(SERVER_ROOT, "static")
    TEMPLATE_DIR = os.path.join(SERVER_ROOT, "templates")

    # Frontend directory
    FRONTEND_BUILD_DIR = os.path.join(STATIC_DIR, "dist")

    # Tools configuration - list of all CLI tools from pyproject.toml
    # Format: (command_name, module_path)
    TOOLS_CONFIG: List[tuple] = [
        ("reads-quality-stats", "gseda.ppl.reads_quality_stats_v2:main_cli"),
        ("reads-quality-stats-v3", "gseda.ppl.reads_quality_stats_v3:main_cli"),
        ("reads-time-err", "gseda.ppl.reads_time_err_ana:main_cli"),
        ("reads-quality-hp-tr", "gseda.ppl.reads_quality_stats_hp_tr:main_cli"),
        ("reads-quality-hp", "gseda.ppl.reads_quality_stats_hp:main_cli"),
        ("sequencing-report", "gseda.ppl.sequencing_report:main_cli"),
        ("sequencing-report-v2", "gseda.ppl.sequencing_report_v2:main_cli"),
        ("sequencing-report-overall", "gseda.ppl.sequencing_report_overall:main_cli"),
        ("seq-n-stats", "gseda.ppl.seq_n_stats:main_cli"),
        ("seq-n-stats-v2", "gseda.ppl.seq_n_stats_v2:main_cli"),
        ("low-high-q-quality-stats", "gseda.ppl.low_high_q_quality_stats:main_cli"),
        ("bam-basic-stat", "gseda.bam_ana.bam_basic_stat:main_cli"),
        ("fastx-basic-stat", "gseda.fastx_ana.fastx_basic_stat:main_cli"),
        ("msa-view", "gseda.msa_view.msa_view_using_pileup:main_cli"),
        ("phreq-ana", "gseda.ppl.phreq_analysis:main"),
        ("wga", "gseda.ppl.whole_genome_alignment:main_cli"),
        ("bam2fx", "gseda.bam_filter.bam2fx:main_cli"),
        ("smc-mem-est", "gseda.bam_ana.mem_est:main_cli"),
        ("low-q-analysis", "gseda.ppl.low_q_analysis:main_cli"),
        ("homo-and-str-ratio", "gseda.ppl.homo_and_str_region_coverage:main_cli"),
        ("macebell-ratio", "gseda.ppl.macebell_ratio:main_cli"),
        ("dump-specified-barcode-sbr", "gseda.bam_filter.dump_sbr_according_barcode:main_cli"),
    ]


settings = Settings()
