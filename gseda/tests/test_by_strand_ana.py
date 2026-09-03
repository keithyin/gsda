import pysam
import os
import subprocess
import json
import pathlib

def create_mock_bam(filename, reads_data):
    header = { 'HD': {'VN': '1.0'},
               'SQ': [{'SN': 'ref', 'LN': 1000}] }

    with pysam.AlignmentFile(filename, "wb", header=header) as out_bam:
        for seq, ch in reads_data:
            record = pysam.AlignedSegment()
            record.query_sequence = seq
            record.set_tag("ch", ch)
            record.is_unmapped = True
            out_bam.write(record)

def test_analysis():
    true_data = [
        ("ATGCATGC", 1),
        ("ATGCATGC", 2),
        ("ATGCATGC", 3),
        ("ATGCATGC", 4),
    ]
    false_data = [
        ("ATGCATGC", 1),
        ("ATGCATGT", 2),
        ("ATGCATC", 3),
        ("ATGCATGC", 5),
    ]
    bam_true = "test_true.bam"
    bam_false = "test_false.bam"
    out_prefix = "test_out"

    try:
        create_mock_bam(bam_true, true_data)
        create_mock_bam(bam_false, false_data)

        cmd = [
            "python3", "src/gseda/ppl/ByStrandTrueAna.py",
            "--bam-true", bam_true,
            "--bam-false", bam_false,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)

        summary_path = f"{out_prefix}.summary.json"
        if not os.path.exists(summary_path):
            print(f"Error: {summary_path} not created")
            return

        with open(summary_path, 'r') as f:
            summary = json.load(f)

        print(f"Summary: {summary}")
        assert summary["total_molecules"] == 3
        assert summary["total_comparisons"] == 3
        assert summary["perfect_matches_count"] == 1
        print("Verification successful!")

    finally:
        for f in [bam_true, bam_false, f"{out_prefix}.details.tsv", f"{out_prefix}.summary.json"]:
            if os.path.exists(f):
                os.remove(f)

if __name__ == "__main__":
    test_analysis()
