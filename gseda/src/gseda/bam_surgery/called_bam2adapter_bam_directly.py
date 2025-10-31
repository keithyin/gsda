import pysam
from tqdm import tqdm
import array


def convert(in_bam_path: str, out_bam_path: str, out_header):
    # interested_channel_and_start_end = """
    # 131235  0-36
    # 32400  0-33
    # 32496  0-32
    # 49264  0-32
    # 41168  0-36
    # 33024 0-37
    # 55264  0-33
    # 63136  0-35
    # 73264  0-53
    # 74224  0-54
    # 81664 0-31
    # 81504  0-58
    # 81904  0-31
    # 90048  0-30
    # 90176 0-31
    # 97632 0-35
    # 67488  0-36
    # 67248 0-32
    # 67072 0-32
    # 72480 0-30 
    # """
    
    # items = interested_channel_and_start_end.split("\n")
    # items = [item.strip() for item in items if item.strip() != ""]
    # print(items)
    
    # items = {int(item.split(" ")[0]): (int(item.split(" ")[-1].split("-")[0]),  int(item.split(" ")[-1].split("-")[1]) + 1) for item in items}
    # print(items)

    with pysam.AlignmentFile(in_bam_path, mode="rb", check_sq=False, threads=40) as in_bam, pysam.AlignmentFile(out_bam_path, check_sq=False, mode="wb", threads=40, header=out_header) as out_bam:
        cnt = -1
        for record in tqdm(in_bam.fetch(until_eof=True), desc=f"processing {in_bam_path}"):
            qname = record.query_name
            ch = int(qname.split("_")[1])
            # if ch not in items:
            #     continue
            
            # start, end = items[ch]
            start, end = 0, len(record.query_sequence)
            
            cnt += 1
            record.query_sequence = record.query_sequence[start:end]
            
            record.query_name = f"read_0/0/subread/{cnt}"
            if record.query_qualities is None:
                record.query_qualities = array.array('H', [50] * (end-start))
            
            record.set_tag("np", 1)
            record.set_tag("cx", 3)
            record.set_tag("ch", 0)
            record.set_tag("rq", 0.8)
            record.set_tag("cr", array.array(
                'B', [0] * (end - start)))
            
            ar = record.get_tag("ar")
            record.set_tag("ar", ar[start:end])
            
            dw = record.get_tag('dw')
            record.set_tag("dw", dw[start:end])
            
            record.set_tag("sn", array.array('f', [20, 20, 20, 20]))
            record.set_tag("be", array.array(
                'I', [0, end-start]))

            record.set_tag("RG", "0425")

            out_bam.write(record)

            pass

    pass


def main():
    in_bam = "/data1/ccs_data/rna/try2_call.bam"
    out_bam = "/data1/ccs_data/rna/try2_call.adapter.bam"

    header = {
        'HD': {'VN': '1.5', 'SO': 'unknown', "bn": "0.01"},
        'RG': [{
            'ID': '0425',
            'PL': 'Gseq500',
            'rn': '20251015_240601Y0004_Run0003',
            'if': 'RECORDTYPE=SUBREAD;KITV=0.0.500;CHV=dna_l958_pa54;BCV=2.3.3;MODEL=L958_PA54_lstm_10_256_v1.27.trt;BCCLUSTER=dbscan'
        }]
    }

    convert(in_bam_path=in_bam, out_bam_path=out_bam, out_header=header)

    pass


if __name__ == "__main__":
    """
    docker run -t --rm -v`pwd`:/data 192.168.3.38:5000/algo/adacus:smc5.5.0_adapter_demux0.0.3_barcode_remover1.0.2_smicing0.4.2_bmi_0.1.4 smc /data/20251015_240601Y0004_Run0003_adapter.bam --hmmModelSpec Kit__500_Chem__1_BC__1_PW3_v4
    """
    main()


