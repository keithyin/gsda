import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import argparse


def plot_rq_iy_scatter(df: pl.DataFrame):
    figure = plt.figure(figsize=(10, 10))
    axs = figure.add_subplot(1, 1, 1)
    plt.sca(axs)
    plt.grid(True, linestyle=":", linewidth=0.5, color="gray")

    df = df.to_pandas()
    sns.scatterplot(df, x="phreq_rq", y="phreq_iy", ax=axs)
    sns.kdeplot(
        df,
        x="phreq_rq",
        y="phreq_iy",
        ax=axs,
        cmap="RdYlGn",
        fill=True,
        alpha=0.5,
        label="Density",
        cbar=True,
    )

    axs.set_xticks(list(range(0, 30, 1)))
    axs.set_yticks(list(range(0, 30, 1)))
    axs.set_xlim((0, 30))
    axs.set_ylim((0, 30))
    axs.set_xlabel("PredictedChannelQ", fontdict={"size": 16})
    axs.set_ylabel("SbrDiscQ", fontdict={"size": 16})
    perfect_line = pl.DataFrame(
        {
            "x": list(range(0, 30)),
            "y": list(range(0, 30)),
        }
    )

    sns.lineplot(
        perfect_line.to_pandas(), x="x", y="y", ax=axs, color="blue", linestyle="--"
    )
    figure.savefig(fname="discrenpy_rq_scatter.png")


def stat(fname: str):
    df = pl.read_csv(fname, separator="\t")
    df = df.filter(pl.col("np").eq(3))
    plot_rq_iy_scatter(df)
    pass


def main(args):
    stat(args.data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("data")
    main(parser.parse_args())
