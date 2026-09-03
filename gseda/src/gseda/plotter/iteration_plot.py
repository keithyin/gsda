import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def do_plot(file_path: str, position=3.5):
    
    # # 设置全局字体样式
    # plt.rcParams.update({
    #     'xtick.labelsize': 12,
    #     'ytick.labelsize': 12,
    #     'xtick.labelweight': 'bold',
    #     'ytick.labelweight': 'bold',
    #     'axes.labelweight': 'bold'
    # })
    
    fig = plt.figure(figsize=(20, 10))
    axes = fig.add_subplot(1, 1, 1)
    axes.tick_params(axis='both', width=4, labelsize=16)
    axes.grid(visible=True)
    data = pl.read_csv(file_path, separator=",")
    data = data.to_pandas()
    sns.lineplot(data, x="idx", y="value", hue="tag", ax=axes, marker="o", markersize=10,   markeredgewidth=1.5,  linewidth=2 )
    xmin = data["idx"].min()
    xmax = data["idx"].max()
    
    axes.set_ylabel("PhreQ", fontdict={"fontsize": 14})
    axes.set_xlabel("Iter", fontdict={"fontsize": 14})

    axes.set_xticks(np.arange(xmin-1, xmax + 2, 1))  # +1确保包含最大值
    axes.set_xlim(xmin=xmin-1, xmax=xmax+1)  # 设置x轴范围

    axes.axvline(x=position, color='red', linestyle='--', linewidth=1)

    y_min, y_max = axes.get_ylim()  # 获取y轴范围元组
    axes.text(position, y_max*0.95, 'DatasetChange', fontdict={"fontsize": 14, "fontweight": "bold"},
              rotation=0, verticalalignment='top',
              bbox=dict(facecolor='white', alpha=0.8))

    plt.savefig("{}.png".format(file_path))
    


def plot_time(file_path: str):
    fig = plt.figure(figsize=(20, 10))
    axes = fig.add_subplot(1, 1, 1)
    axes.tick_params(axis='both', width=4, labelsize=16)
    
    axes.grid(visible=True)
    data = pl.read_csv(file_path, separator=",")
    data = data.to_pandas()
    xmin = data["idx"].min()
    xmax = data["idx"].max()

    axes.set_xticks(np.arange(xmin-1, xmax + 2, 1))  # +1确保包含最大值
    axes.set_xlim(xmin=xmin-1, xmax=xmax+1)  # 设置x轴范围
    axes.set_ylabel("Time(min)")
    axes.set_xlabel("Iter")
    sns.lineplot(data, x="idx", y="time", ax=axes, marker="o", markersize=10,   markeredgewidth=1.5,  linewidth=2)
    fig.savefig("{}.png".format(file_path))
    pass


def main():
    do_plot("/root/projects/gsda/gseda/src/gseda/plotter/smc_iter_data.csv", position=3.5)
    do_plot("gseda/src/gseda/plotter/icing_iter_data.csv", position=6.5)
    plot_time("/root/projects/gsda/gseda/src/gseda/plotter/speed_data.csv")
    pass


if __name__ == "__main__":
    main()
