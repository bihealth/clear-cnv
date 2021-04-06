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
    SENSITIVITY = args.sensitivity if args.sensitivity >= 0 and args.sensitivity <= 1 else 0.7
    CORES = min([args.cores, mp.cpu_count()]) if args.cores else mp.cpu_count()
    DUP_CUTOFF = args.dup_cutoff #1.35
    DEL_CUTOFF = args.del_cutoff #0.75

    # load data
    D0 = util.load_dataframe(intsv_path)
    Matchscores = pd.read_csv(matchscores_path, sep="\t", low_memory=False, header=0, index_col=0)

    # adaptive threshold for group sizes
    minmatchscore = np.median(Matchscores.median())*SAMPLE_SCORE_FACTOR
    Matchscores_bools = Matchscores < minmatchscore
    selected_samples = Matchscores_bools[Matchscores_bools.sum() >= MINIMUM_SAMPLE_GROUP].index
    failed_samples = Matchscores_bools.index.difference(selected_samples)
    # optional for neater code downstream
    Matchscores_bools_selected = Matchscores_bools.loc[selected_samples,selected_samples]
    Matchscores_selected = Matchscores.loc[selected_samples,selected_samples]

    # evaluation
    #if len(failed_samples) > 0:
    with open(pathlib.Path(analysis_dir) / "failed_samples.txt", 'w') as f:
        print(
            "%.d samples failed to have a sufficient group size. They are ignored in cnv calling."
            % len(failed_samples),
            file=f
        )
        for s in failed_samples:
            print(s,file=f)

    plt.figure(figsize=(6,4))
    plt.title("Threshold finding of sample groups")
    x = plt.hist(Matchscores.median().sort_values(),bins=30)
    plt.vlines(minmatchscore,0,max(x[0]),color='orange',label="threshold")
    plt.xlim(0,1)
    plt.ylabel("number of samples")
    plt.xlabel("median matchscore")
    plt.legend()
    plt.savefig(pathlib.Path(analysis_dir) / "ANALYSIS_group_sizes_threshold.pdf", format="pdf")


    plt.figure(figsize=(6,4))
    plt.title("Sample group sizes and cutoff")
    x = plt.hist(Matchscores_bools.sum().sort_values(),bins=round(len(selected_samples)/3),range=(0,len(selected_samples)))
    plt.vlines(MINIMUM_SAMPLE_GROUP,0,max(x[0]),color='darkred',label="cutoff")
    plt.xlim(0,len(selected_samples))
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
                SENSITIVITY,
            ),
            callback=collect_result,
        )
    pool.close()
    pool.join()

    #  calculate expected X-cov ratio and prepare match scores
    if sum(DX[1]) > 0 and np.max(DX[0].mean(axis=0)) / np.min(DX[0].mean(axis=0)) > 1.5:
        kmeans = KMeans(n_clusters=2).fit(DX[0].mean(axis=0).reshape(-1, 1))
        Mx = Matchscores_bools_selected.copy()
        for i, x in enumerate(kmeans.labels_):
            sample = Matchscores_bools_selected.index[i]
            if x:
                Mx[sample] &= kmeans.labels_ == x
            else:
                Mx[sample] &= (1 - kmeans.labels_) == x

        sample_scores_buffer = sample_scores.copy()
        #  --- parallelized X-chromosome calling ---  #
        print("multi threaded calling on X-chromosome")
        pool = mp.Pool(CORES)
        for i in range(len(selected_samples)):
            pool.apply_async(
                util.calling_cnv,
                args=(
                    i,
                    DX[0],
                    Matchscores_bools_selected,
                    DX[1],
                    EXPECTED_CNV_RATE,
                    SENSITIVITY,
                ),
                callback=collect_result,
            )
        pool.close()
        pool.join()
        sample_scores += sample_scores_buffer

    #  calculate expected Y-cov ratio and prepare match scores
    if sum(DY[1]) > 0 and np.max(DY[0].mean(axis=0)) / np.min(DY[0].mean(axis=0)) > 1.5:
        kmeans = KMeans(n_clusters=2).fit(DY[0].mean(axis=0).reshape(-1, 1))
        My = Matchscores_bools_selected.copy()
        for i, y in enumerate(kmeans.labels_):
            sample = Matchscores_bools_selected.index[i]
            My[sample] &= kmeans.labels_ == y

        sample_scores_buffer = sample_scores.copy()
        #  --- parallelized Y-chromosome calling ---  #
        print("multi threaded calling on Y-chromosome")
        pool = mp.Pool(CORES)
        for i in range(len(selected_samples)):
            pool.apply_async(
                util.calling_cnv,
                args=(
                    i,
                    DY[0],
                    My,
                    DY[1],
                    EXPECTED_CNV_RATE,
                    SENSITIVITY,
                ),
                callback=collect_result,
            )
        pool.close()
        pool.join()
        sample_scores += sample_scores_buffer

    # extract single exon CNVs
    S = pd.DataFrame(sample_scores, index=SAMPLES).T
    DFZ = pd.DataFrame(z_scores_scaled, columns=INDEX, index=SAMPLES).T
    DFR = pd.DataFrame(ratio_scores, columns=INDEX, index=SAMPLES).T
    SINGLE_EXONS = pd.concat(
        [
            pd.DataFrame(
                [
                    [
                        DFZ.columns[c],
                        *DFZ.index[i].split("_"),
                        "DEL",
                        1,
                        DFZ.iloc[i, c],
                        float(S[DFZ.columns[c]]),
                    ]
                    for i, c in zip(*np.where((DFZ < -5) & (DFR < DEL_CUTOFF)))
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
                        DFZ.columns[c],
                        *DFZ.index[i].split("_"),
                        "DUP",
                        1,
                        DFZ.iloc[i, c],
                        float(S[DFZ.columns[c]]),
                    ]
                    for i, c in zip(*np.where((DFZ > 6) & (DFR > DUP_CUTOFF)))
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

    print("compute HMM on z-scores")
    # HMM guided CNV candidates
    # =========================================================================
    # --- fix lockdown bug --- #
    # pool = mp.Pool(CORES)
    probs = np.array([[-2.5], [0.0], [3.5]])
    transitionprobs = np.array(
        [[0.9999, 0.0001, 0.0], [0.0001, 0.9998, 0.0001], [0.0, 0.0001, 0.9999]]
    )
    # HMM = np.array([pool.apply(util.hmm_scores, args=([sample,probs,transitionprobs]))
    #                for sample in z_scores_scaled])
    # pool.close()
    HMM = np.array([util.hmm_scores(sample, probs, transitionprobs) for sample in z_scores_scaled])
    # --- fix lockdown bug --- #
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
                            DUP_CUTOFF
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
                    for i, c in zip(*np.where((RZ < -5) & (RR < DEL_CUTOFF)))
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
                    for i, c in zip(*np.where((RZ > 6) & (RR > DUP_CUTOFF)))
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
    print("extract single exon CNVs.")
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
    X = [
        [*h.to_list()[:6], h.to_list()[7], h.to_list()[6], *S[h.to_list()[0]]]
        for h in ca.hitsA_not_in_hitsB(SINGLE_HITS, BIG_HITS)
    ]
    FINAL = pd.concat([RD, pd.DataFrame(X, columns=RD.columns)])
    FINAL["score"] = FINAL["score"].astype(float)#.apply(x)
    FINAL["score"] = FINAL["score"].apply(lambda x: float(abs(x)))
    FINAL = FINAL.sort_values(by="score", ascending=False)

    print("saving results...")
    FINAL.to_csv(calls_path, sep="\t", index=False)
    RZ.to_csv(z_scores_path, sep="\t")
    RR.to_csv(ratio_scores_path, sep="\t")
    S.T.to_csv(pathlib.Path(analysis_dir) / "sample_scores.tsv", sep= '\t')
    # analysis_dir
    print("done!")
    return RD, RZ, RR
