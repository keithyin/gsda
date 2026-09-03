from datetime import datetime
import re
import matplotlib.pyplot as plt
import glob
from typing import List, Tuple
import polars as pl
import os


def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "30000"
    os.environ["POLARS_FMT_STR_LEN"] = "100"


def get_files(root_dir: str) -> List[str]:
    files = list(glob.glob(f"{root_dir}/*called-metric/*gsmm2-hp-aggr.csv"))
    return files


def read_datas(files: List[str]) -> List[Tuple[str, pl.DataFrame]]:
    result = []
    for fname in files:
        data = pl.read_csv(fname, separator="\t")
        result.append([fname, data])
    return result


def filter_data(data: pl.DataFrame) -> pl.DataFrame:
    return data.filter(pl.col("called") < (pl.col("true_cnt") + pl.lit(3))
                       ).filter(pl.col("called") > (pl.col("true_cnt") - pl.lit(3)))


def extract_date(filename: str) -> datetime.date:
    """从文件名提取日期(yyyymmdd格式)"""
    match = re.search(r"/(\d{8})_", filename)
    if match:
        return datetime.strptime(match.group(1), "%Y%m%d").date()
    raise ValueError(f"无法从文件名提取日期: {filename}")


def process_and_plot(data_list: List[Tuple[str, pl.DataFrame]], output_dir: str = "plots"):
    """主处理函数：合并数据并生成可视化图表"""
    # 数据预处理
    processed = []
    for filename, df in data_list:
        try:
            date = filename
            processed.append(df.with_columns(
                pl.lit(date).alias("date"),
                pl.col("called").cast(str)
            ))
        except ValueError as e:
            print(f"警告：{str(e)}")

    if not processed:
        print("无有效数据可处理")
        return

    combined = pl.concat(processed).sort("date")
    os.makedirs(output_dir, exist_ok=True)

    # 按true_cnt分组处理
    for cnt, cnt_df in combined.group_by("true_cnt"):
        fig, axes = plt.subplots(4, 2, figsize=(20, 15))
        fig.suptitle(f"True Count = {cnt}", fontsize=16)

        bases = ["A", "C", "G", "T"]
        tags = ["pure", "mixed"]

        # 绘制4x2子图矩阵
        for i, base in enumerate(bases):
            for j, tag in enumerate(tags):
                ax = axes[i][j]
                subset = cnt_df.filter(
                    (pl.col("true_base") == base) &
                    (pl.col("tag") == tag)
                )

                # 按called值分组绘制曲线
                for called, group in subset.group_by("called"):
                    sorted_df = group.sort("date")
                    ax.plot(
                        sorted_df["date"],
                        sorted_df["ratio_within_motif"],
                        marker="o",
                        label=f"called={called}"
                    )

                ax.set_title(f"Base {base} | Tag {tag}")
                ax.set_xlabel("Date")
                ax.set_ylabel("Ratio Within Motif")
                ax.legend()
                ax.grid(True)

        plt.tight_layout()
        output_path = os.path.join(output_dir, f"true_cnt_{cnt}.png")
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"已生成图表: {output_path}")


def main():
    files = files = list(
        glob.glob("/data/ccs_data/sa/*smc_all_reads-metric/*gsmm2-hp-aggr.csv"))
    file_and_datas = read_datas(files)
    file_and_datas = [[v[0], filter_data(v[1])] for v in file_and_datas]
    file_and_datas = sorted(file_and_datas, key=lambda x: x[0])
    for i in range(len(file_and_datas)):
        file_and_datas[i][0] = i

    process_and_plot(file_and_datas, output_dir="plotdir-sa-smc")


if __name__ == "__main__":
    polars_env_init()
    main()
