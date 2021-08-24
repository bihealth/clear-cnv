#!/usr/bin/env python
# coding: utf-8

import pathlib

import pandas as pd
import numpy as np

import multiprocessing as mp
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

from clearCNV import util
from clearCNV import cnv_arithmetics as ca


def cnv_calling(args):
    # argparsing
    intsv_path = args.coverages
    analysis_dir = args.analysis_directory
    matchscores_path = args.matchscores
    calls_path = args.cnv_calls
    ratio_scores_path = args.ratio_scores
    z_scores_path = args.z_scores
    EXPECTED_CNV_RATE = args.expected_artefacts
    SAMPLE_SCORE_FACTOR = args.sample_score_factor
    MINIMUM_SAMPLE_GROUP = args.minimum_group_sizes
    ZSCALE = args.zscale if args.zscale >= 0 and args.zscale <= 2 else 0.65
    CORES = min([args.cores, mp.cpu_count()]) if args.cores else mp.cpu_count()
    DUP_CUTOFF = args.dup_cutoff  # 1.35
    DEL_CUTOFF = args.del_cutoff  # 0.75
    TRANSPROB = args.trans_prob  # 0.001
    PLOT_REGIONS = args.plot_regions

    # load data
    D0 = util.load_dataframe(intsv_path)
    Matchscores = pd.read_csv(matchscores_path, sep="\t", low_memory=False, header=0, index_col=0)

    # adaptive threshold for group sizes
    minmatchscore = np.median(Matchscores.median()) * SAMPLE_SCORE_FACTOR
    Matchscores_bools = Matchscores < minmatchscore
    selected_samples = Matchscores_bools[Matchscores_bools.sum() >= MINIMUM_SAMPLE_GROUP].index
    failed_samples = Matchscores_bools.index.difference(selected_samples)
    np.fill_diagonal(Matchscores_bools.values, False)
    # optional for neater code downstream
    Matchscores_bools_selected = Matchscores_bools.loc[selected_samples, selected_samples]
    Matchscores_selected = Matchscores.loc[selected_samples, selected_samples]

    if len(selected_samples) == 0:
        print(
            "ERROR: NO SAMPLES IN ANALYSIS. There were no samples selected to perform CNV calling. Try to pick a lower value for MINIMUM_GROUP_SIZES or a higher value for SAMPLE_SCORE_FACTOR. Aborting."
        )
        raise Exception(
            "No samples selected to smaple group. Try to pick a lower value for MINIMUM_GROUP_SIZES or a higher value for SAMPLE_SCORE_FACTOR."
        )

    for p in [calls_path, z_scores_path, ratio_scores_path]:
        pathlib.Path(p).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(analysis_dir).mkdir(parents=True, exist_ok=True)

    # evaluation
    # if len(failed_samples) > 0:
    plt.figure(figsize=(6, 4))
    plt.title("Threshold finding of sample groups")
    x = plt.hist(Matchscores.median().sort_values(), bins=30)
    plt.vlines(minmatchscore, 0, max(x[0]), color="orange", label="threshold")
    plt.xlim(0, 1)
    plt.ylabel("number of samples")
    plt.xlabel("median matchscore")
    plt.legend()
    plt.savefig(pathlib.Path(analysis_dir) / "ANALYSIS_group_sizes_threshold.pdf", format="pdf")

    plt.figure(figsize=(6, 4))
    plt.title("Sample group sizes and cutoff")
    x = plt.hist(
        Matchscores_bools.sum().sort_values(),
        bins=int(len(selected_samples) / 3) + 1,
        range=(0, len(selected_samples)),
    )
    plt.vlines(MINIMUM_SAMPLE_GROUP, 0, max(x[0]), color="darkred", label="cutoff")
    plt.xlim(0, len(selected_samples))
    plt.ylabel("number of samples")
    plt.xlabel("group size")
    plt.legend()
    plt.savefig(pathlib.Path(analysis_dir) / "ANALYSIS_group_sizes_cutoff.pdf", format="pdf")

    util.print_clustermap(
        Matchscores_selected,
        path=pathlib.Path(analysis_dir) / "ANALYSIS_clustermap_filtered.pdf",
        title="clustered heat map of selected sample distances",
    )

    # separate chr subsets
    DS0 = util.filter_by_coverage(D0[selected_samples])
    INDEX = DS0.index.copy()
    SAMPLES = DS0.columns.copy()
    DS1 = util.normalize_within_sample(DS0)
    DS1.index = range(len(INDEX))
    I = np.array([(not (s.startswith("X") or s.startswith("Y"))) for s in INDEX])
    DA = (DS1.to_numpy()[I], I)
    I = np.array([s.startswith("X") for s in INDEX])
    DX = (DS1.to_numpy()[I], I)
    I = np.array([s.startswith("Y") for s in INDEX])
    DY = (DS1.to_numpy()[I], I)

    # define data containers
    n, m = len(DS1.columns), len(DS1.index)
    z_scores_scaled = np.zeros((n, m))
    # e_scores        = np.zeros((n,m))
    ratio_scores = 1 + np.zeros((n, m))
    sample_scores = np.zeros(n)

    # parallelized execution
    def collect_result(result):
        sample_id = result[0]
        index = result[1]
        nonlocal sample_scores
        nonlocal ratio_scores
        nonlocal z_scores_scaled
        # E_scores_df.loc[result[1].index,sample]        = result[1]
        sample_scores[sample_id] = result[2] * (sum(index) / len(INDEX))
        ratio_scores[sample_id, index] = result[3]
        z_scores_scaled[sample_id, index] = result[4]

    #  parallelized autosomal calling
    pool = mp.Pool(CORES)
    print("multi threaded calling autosomal targets")
    for i in range(len(selected_samples)):
        pool.apply_async(
            util.calling_cnv,
            args=(
                i,
                DA[0],
                Matchscores_bools_selected,
                DA[1],
                EXPECTED_CNV_RATE,
                ZSCALE,
            ),
            callback=collect_result,
        )
    pool.close()
    pool.join()

    #  calculate expected X-cov ratio and prepare match scores
    if sum(DX[1]) > 0 and np.max(DX[0].mean(axis=0)) / np.min(DX[0].mean(axis=0)) > 1.25:
        kmeans = KMeans(n_clusters=2).fit(DX[0].mean(axis=0).reshape(-1, 1))
        Mx = Matchscores_bools_selected.copy()
        for i, x in enumerate(kmeans.labels_):
            sample = Matchscores_bools_selected.index[i]
            if x:
                Mx[sample] &= kmeans.labels_ == x
            else:
                Mx[sample] &= (1 - kmeans.labels_) == x
        Mx &= Mx.T
        sample_scores_buffer = sample_scores.copy()
        #  --- parallelized X-chromosome calling ---  #
        print("multi threaded calling on X-chromosome")
        pool = mp.Pool(CORES)
        for i in range(len(selected_samples)):
            if sum(Mx.iloc[i]) >= MINIMUM_SAMPLE_GROUP:
                pool.apply_async(
                    util.calling_cnv,
                    args=(
                        i,
                        DX[0],
                        Mx,
                        DX[1],
                        EXPECTED_CNV_RATE,
                        ZSCALE,
                    ),
                    callback=collect_result,
                )
        pool.close()
        pool.join()
        sample_scores += sample_scores_buffer

    #  calculate expected Y-cov ratio and prepare match scores
    if sum(DY[1]) > 0 and np.max(DY[0].mean(axis=0)) / np.min(DY[0].mean(axis=0)) > 1.25:
        kmeans = KMeans(n_clusters=2).fit(DY[0].mean(axis=0).reshape(-1, 1))
        My = Matchscores_bools_selected.copy()
        for i, y in enumerate(kmeans.labels_):
            sample = Matchscores_bools_selected.index[i]
            My[sample] &= kmeans.labels_ == y
        My &= My.T
        sample_scores_buffer = sample_scores.copy()
        #  --- parallelized Y-chromosome calling ---  #
        print("multi threaded calling on Y-chromosome")
        pool = mp.Pool(CORES)
        for i in range(len(selected_samples)):
            if sum(My.iloc[i]) >= MINIMUM_SAMPLE_GROUP:
                pool.apply_async(
                    util.calling_cnv,
                    args=(
                        i,
                        DY[0],
                        My,
                        DY[1],
                        EXPECTED_CNV_RATE,
                        ZSCALE,
                    ),
                    callback=collect_result,
                )
        pool.close()
        pool.join()
        sample_scores += sample_scores_buffer

    # extract single exon CNVs
    print("compute HMM on z-scores")
    # HMM guided CNV candidates
    # =========================================================================
    # HMM PARAMETRIZATION
    mean_del = np.std(z_scores_scaled.flatten(), ddof=1) * (-3.0)
    mean_dup = np.std(z_scores_scaled.flatten(), ddof=1) * 4.0
    mean_median = np.median(z_scores_scaled.flatten())
    # HMM guided CNV candidates
    # =========================================================================
    means = np.array([[mean_del], [mean_median], [mean_dup]])
    p_trans = TRANSPROB
    p_stay = 1.0 - (2.0 * TRANSPROB)
    transitionprobs = np.array(
        [[p_stay, p_trans, p_trans], [p_trans, p_stay, p_trans], [p_trans, p_trans, p_stay]]
    )
    # HMM = np.array([util.hmm_scores(sample, probs, transitionprobs) for sample in z_scores_scaled])
    # =========================================================================
    pool = mp.Pool(CORES)
    HMM = np.array(
        [
            pool.apply(util.hmm_scores, args=([sample, means, transitionprobs]))
            for sample in z_scores_scaled
        ]
    )
    pool.close()
    # =========================================================================

    print("iterate HMM-labels and condense CNVs")
    # extract CNVs from HMM
    cnv_carrying_indexes = np.array(range(z_scores_scaled.shape[0]))[
        [sum((s != 1)) > 0 for s in HMM]
    ]
    pool = mp.Pool(CORES)
    CNVs = list(
        zip(
            cnv_carrying_indexes,
            [
                pool.apply(
                    util.merge_score_cnvs,
                    args=(
                        [
                            HMM[i, :],
                            z_scores_scaled[i, :],
                            ratio_scores[i, :],
                            INDEX,
                            DEL_CUTOFF,
                            DUP_CUTOFF,
                        ]
                    ),
                )
                for i in cnv_carrying_indexes
            ],
        )
    )
    pool.close()

    # format and save results
    RD = pd.DataFrame(
        [
            [SAMPLES[sample_cnv[0]], *cnv.to_list(), sample_scores[CNVs[i][0]]]
            for i, sample_cnv in enumerate(CNVs)
            for cnv in sample_cnv[1]
        ],
        columns=[
            "sample",
            "chr",
            "start",
            "end",
            "gene",
            "aberration",
            "size",
            "score",
            "sample_score",
        ],
    ).sort_values(by="score", ascending=False)

    # extract single exon CNVs
    RZ = pd.DataFrame(z_scores_scaled, columns=INDEX, index=SAMPLES).T
    RR = pd.DataFrame(ratio_scores, columns=INDEX, index=SAMPLES).T
    S = pd.DataFrame(sample_scores, index=SAMPLES).T
    SINGLE_EXONS = pd.concat(
        [
            pd.DataFrame(
                [
                    [
                        RZ.columns[c],
                        *RZ.index[i].split("_"),
                        "DEL",
                        1,
                        RZ.iloc[i, c],
                        float(S[RZ.columns[c]]),
                    ]
                    for i, c in zip(*np.where((RZ < -3.5) & (RR < DEL_CUTOFF)))
                ],
                columns=[
                    "sample",
                    "chr",
                    "start",
                    "end",
                    "gene",
                    "aberration",
                    "size",
                    "score",
                    "sample_score",
                ],
            ),
            pd.DataFrame(
                [
                    [
                        RZ.columns[c],
                        *RZ.index[i].split("_"),
                        "DUP",
                        1,
                        RZ.iloc[i, c],
                        float(S[RZ.columns[c]]),
                    ]
                    for i, c in zip(*np.where((RZ > 4.5) & (RR > DUP_CUTOFF)))
                ],
                columns=[
                    "sample",
                    "chr",
                    "start",
                    "end",
                    "gene",
                    "aberration",
                    "size",
                    "score",
                    "sample_score",
                ],
            ),
        ]
    )
    print("extract single exon CNVs")
    BIG_HITS = ca.flatten_list(
        [
            l.get_hits()
            for l in ca.load_cnvs_from_df(
                RD[["sample", "chr", "start", "end", "gene", "aberration", "score", "size"]]
            )
        ]
    )
    SINGLE_HITS = ca.flatten_list(
        [
            l.get_hits()
            for l in ca.load_cnvs_from_df(
                SINGLE_EXONS[
                    ["sample", "chr", "start", "end", "gene", "aberration", "score", "size"]
                ]
            )
        ]
    )
    FS = pd.DataFrame(Matchscores_bools.loc[failed_samples].sum(axis=1), columns=["group_size"])
    FS["median_sample_score"] = Matchscores.loc[failed_samples].median(axis=1)
    X = [
        [*h.to_list()[:6], h.to_list()[7], h.to_list()[6], *S[h.to_list()[0]]]
        for h in ca.hitsA_not_in_hitsB(SINGLE_HITS, BIG_HITS)
    ]
    FINAL = pd.concat([RD, pd.DataFrame(X, columns=RD.columns)])
    FINAL["score"] = FINAL["score"].astype(float).apply(lambda x: float(abs(x)))
    FINAL = FINAL.sort_values(by="score", ascending=False)
    FINAL.index = list(range(FINAL.shape[0]))

    print("saving results...")
    FINAL.to_csv(calls_path, sep="\t", index=False)
    RZ.to_csv(z_scores_path, sep="\t")
    RR.to_csv(ratio_scores_path, sep="\t")
    S.T.to_csv(pathlib.Path(analysis_dir) / "sample_scores.tsv", sep="\t")
    Matchscores_bools.to_csv(pathlib.Path(analysis_dir) / "samplegroups.tsv", sep="\t")
    FS.to_csv(pathlib.Path(analysis_dir) / "failed_samples.tsv", sep="\t")

    if PLOT_REGIONS:
        print("plotting all called CNVs with sample groups...")
        factor = 0.08
        final_regions = ["-".join(FINAL.loc[i, ["chr", "start", "end"]]) for i in FINAL.index]
        for i, region in enumerate(final_regions):
            print(f"plotting {i} of %d" % len(final_regions))
            sample = FINAL["sample"].iloc[i]
            r = util.select_region(region, RR, sample, Matchscores_bools_selected, buffer=2).clip(
                0, 2
            )
            plt.figure(figsize=(factor * len(r.columns), factor * len(r.index)))
            plt.pcolor(r)
            plt.title(f"{sample} - {region}")
            plt.ylabel("targets")
            plt.xlabel("samples")
            plt.savefig(
                pathlib.Path(analysis_dir) / str(sample + "_" + region + ".pdf"), format="pdf"
            )
            plt.close()

    print("done!")
    return RD, RZ, RR
