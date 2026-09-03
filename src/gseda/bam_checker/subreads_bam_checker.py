import pysam
from tqdm import tqdm

def check_continuous(subraeds_bam_path: str):
    channel_positions = {}
    with pysam.AlignmentFile(
        filename=subraeds_bam_path, check_sq=False, mode="rb", threads=40
    ) as in_bam:
        
        for idx, record in tqdm(enumerate(in_bam.fetch(until_eof=True)), desc=f"reading {subraeds_bam_path}"):
            ch = record.get_tag("ch")
            channel_positions.setdefault(ch, [])
            channel_positions[ch].append(idx)

    for ch, positions in tqdm(channel_positions.items(), desc="checking"):
        positions = sorted(positions)
        if len(positions) > 1:
            succ = True
            pre_pos = positions[0]
            for cur_pos in positions[1:]:
                if (cur_pos - pre_pos) > 1:
                    succ = False
                    break
                pre_pos = cur_pos
            if not succ:
                print(f"not continuous. ch:{ch}, positions:{positions}")
                
def main():
    p = "/data/ccs_data/bad-adapter-bam/20250523_240601Y0011_Run0002_adapter.bam"
    check_continuous(p)
    pass
if __name__ == "__main__":
    main()
