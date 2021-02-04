#!/usr/bin/env python
# coding: utf-8


import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from collections import OrderedDict
from sklearn.cluster import KMeans

import plotly.graph_objects as go
import plotly
from plotly.subplots import make_subplots

# def visualize(analysis_dir, ratio_scores_path, z_scores_path, annotated_path, htmp_size):
def visualize(args):
    # argparsing
    analysis_dir = args.analysis_directory
    ratio_scores_path = args.ratio_scores
    z_scores_path = args.z_scores
    annotated_path = args.annotated
    htmp_size = args.size

    # load data
    annotated = pd.read_csv(annotated_path, sep="\t", low_memory=False, header=None)
    annotated.columns = ["chr", "start", "end", "gene", "size", "mapp", "gc"]
    annotated.index = [
        str(annotated["chr"][i])
        + "_"
        + str(annotated["start"][i])
        + "_"
        + str(annotated["end"][i])
        + "_"
        + str(annotated["gene"][i])
        for i in range(len(annotated.index))
    ]
    ratio_scores_df = pd.read_csv(ratio_scores_path, sep="\t", low_memory=False, index_col=0)
    z_scores_scaled_df = pd.read_csv(z_scores_path, sep="\t", low_memory=False, index_col=0)
    RS = ratio_scores_df.fillna(1.0).copy()

    def splits(index, htmp_size):
        indexchrs = [s.split("_")[0] for s in index]
        chrs = list(OrderedDict.fromkeys(indexchrs))
        chrstarts = [indexchrs.index(s) for s in chrs]
        chrstarts.append(len(indexchrs))
        arr = np.abs(np.array(chrstarts) - (len(index) / 2))
        split_pos = list(arr).index(min(arr)) + 1

        if len(chrs) <= 1 or len(index) <= htmp_size:
            return [index]
        else:
            return [
                *splits(index[: chrstarts[:split_pos][-1]], htmp_size),
                *splits(index[chrstarts[:split_pos][-1] :], htmp_size),
            ]

    def get_plots(df, annotated):
        n_clusters = math.ceil(len(df.columns) / 40)
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(
            df[df.std(axis=1) < np.nanquantile(df.std(axis=1), 0.95)].T
        )

        clustered_samples = pd.DataFrame(list(zip(df.columns, kmeans.labels_))).sort_values(by=1)[0]

        df = df[clustered_samples].copy()

        htmp = go.Heatmap(
            z=[np.clip(df.iloc[i, :], 0, 2) for i in range(len(df.index))],
            x=df.columns,
            y=df.index,
            hoverongaps=False,
        )

        mgc = go.Scatter(y=df.index, x=annotated.loc[df.index, "gc"])
        mmapp = go.Scatter(y=df.index, x=annotated.loc[df.index, "mapp"])
        msize = go.Scatter(y=df.index, x=annotated.loc[df.index, "size"])

        sp = make_subplots(
            rows=1,
            cols=4,
            subplot_titles=("ratios", "mapp", "gc", "log_size"),
            column_widths=[0.85, 0.05, 0.05, 0.05],
            shared_yaxes=True,
        )

        sp.add_trace(htmp, row=1, col=1)
        sp.add_trace(mmapp, row=1, col=2)
        sp.add_trace(mgc, row=1, col=3)
        sp.add_trace(msize, row=1, col=4)

        sp.update_layout(
            width=len(df.columns) * 3 + 300,
            height=len(df.index) * 3 + 300,
        )
        return sp

    print("plotting heatmaps...")
    # parallelize?
    for i, split in enumerate(splits(RS.index, htmp_size)):
        plotly.offline.plot(
            get_plots(RS.loc[split, :], annotated),
            filename=str(pathlib.Path(analysis_dir) / f"ratio_scores_extended_{i}.html"),
        )
        print(f"ratio_scores_extended_{i}.html finished")

    # analyze relation GC content and st. dev. of z-scores after grouping
    print("analyze relation GC content and st. dev. of z-scores after grouping.")

    z_scores_scaled_df["std"] = [
        z_scores_scaled_df.loc[i, :].std() for i in z_scores_scaled_df.index
    ]
    annotated["std"] = z_scores_scaled_df["std"]
    z_scores_scaled_df = z_scores_scaled_df.drop(columns="std")

    ratio_scores_df["std"] = [ratio_scores_df.loc[i, :].std() for i in ratio_scores_df.index]
    annotated["ratio"] = ratio_scores_df["std"]
    ratio_scores_df = ratio_scores_df.drop(columns="std")

    plt.figure(figsize=(8, 5))
    plt.scatter(annotated["std"], annotated["gc"], marker=".", alpha=0.1, label="sd of z-scores")
    plt.ylabel("GC content")
    plt.xlabel("sd of z-scores")
    plt.title("GC content vs z-scores after grouping")
    plt.legend()
    plt.savefig(pathlib.Path(analysis_dir) / "ANALYSIS_GC_vs_SD.pdf", format="pdf")

    plt.figure(figsize=(8, 5))
    plt.scatter(
        annotated["ratio"], annotated["gc"], marker=".", alpha=0.1, label="sd of ratio scores"
    )
    plt.ylabel("GC content")
    plt.xlabel("sd of ratios")
    plt.title("GC content vs sd of ratio scores after grouping")
    plt.legend()
    plt.savefig(str(pathlib.Path(analysis_dir) / "ANALYSIS_ratio_vs_gc.pdf"), format="pdf")
