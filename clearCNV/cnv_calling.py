#!/usr/bin/env python
# coding: utf-8

# In[14]:


import argparse

import pandas as pd
import numpy as np

import multiprocessing as mp
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import math

from clearCNV import util
from clearCNV import cnv_arithmetics as ca


# In[2]:


r"""
parser = argparse.ArgumentParser(description="CNV calling script. Output is a single cnv.calls file in tsv format. Some quality control plots are added to analysis directory in the process.")
parser.add_argument("-p", "--panel",               help="Name of the data set(or panel)",                                                      required=True,  type=str)
parser.add_argument("-c", "--coverages",           help="Coverages file in tsv format",                                                        required=True,  type=str)
parser.add_argument("-a", "--analysis_directory",  help="Path to the directory, where analysis files are stored",                              required=True,  type=str)
parser.add_argument("-m", "--matchscores",         help="matchscores.tsv file generated with matchscores.py",                                  required=True,  type=str)
parser.add_argument("-C", "--cnv_calls",           help="Output cnv.calls file formatted in tsv format",                                       required=True,  type=str)
parser.add_argument("-r", "--ratio_scores",        help="Output ratio scores file in tsv format. Best kept together with cnv.calls",           required=True,  type=str)
parser.add_argument("-z", "--z_scores",            help="Output z-scores file in tsv format. Best kept together with cnv.calls",               required=True,  type=str)
parser.add_argument("-x", "--expected_artifacts",  help="Expected ratio of CNVs or artifacs in target fragment counts",                        required=False, type=float, default=0.02)
parser.add_argument("-u", "--minimum_sample_score",help="A lower threshold results in better fitting, but smaller calling groups",             required=False, type=float, default=0.15)
parser.add_argument("-g", "--minimum_group_sizes", help="Minimum group size per CNV calling group per match scores",                           required=False, type=int,   default=33)
parser.add_argument("-s", "--sensitivity",         help="A higher sensitivity results in more CNV calls. Can only be 0.0 <= sens <= 1.0",      required=False, type=float, default=0.7)
parser.add_argument("--cores",                     help="Number of cpu cores used in parallel processing. Default: determined automatically.", required=False, type=int,   default=0)

args = parser.parse_args()

panel                = args.panel
intsv_path           = args.coverages
analysis_dir         = args.analysis_directory
matchscores_path     = args.matchscores
calls_path           = args.cnv_calls
ratio_scores_path    = args.ratio_scores
z_scores_path        = args.z_scores
EXPECTED_CNV_RATE    = args.expected_artifacts
MINIMUM_SAMPLE_SCORE = args.minimum_sample_score
MINIMUM_SAMPLE_GROUP = args.minimum_group_sizes
SENSITIVITY          = args.sensitivity if args.sensitivity >= 0 and args.sensitivity <= 1 else 0.0
CORES                = min([args.cores, mp.cpu_count()]) if args.cores else mp.cpu_count()
"""


# In[15]:


r"""
panel = "TAADv2"

intsv_path       = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/CALLING/%s/cov/coverages.tsv"%panel

# the resulting matchscore matrix is saved here
calls_path       = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/CALLING/%s/results/cnv.calls"%panel

# all analysis results are saved here
analysis_dir     = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/CALLING/%s/analysis/"%panel

matchscores_path = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/CALLING/%s/analysis/matchscores.tsv"%panel

ratio_scores_path = analysis_dir+'ratio_scores.tsv'
z_scores_path = analysis_dir+'z_scores.tsv'

EXPECTED_CNV_RATE    = 0.02
MINIMUM_SAMPLE_GROUP = 33
SENSITIVITY = 0.7
MINIMUM_SAMPLE_SCORE = 0.15
CORES = 4
"""


# In[20]:


# D,Z,R = cnv_calling(panel, intsv_path, analysis_dir, matchscores_path,
#        calls_path, ratio_scores_path, z_scores_path, EXPECTED_CNV_RATE,
#        MINIMUM_SAMPLE_GROUP, MINIMUM_SAMPLE_SCORE, SENSITIVITY, CORES)


# In[21]:


# def cnv_calling(panel, intsv_path, analysis_dir, matchscores_path,
#        calls_path, ratio_scores_path, z_scores_path, EXPECTED_CNV_RATE,
#        MINIMUM_SAMPLE_GROUP, MINIMUM_SAMPLE_SCORE, SENSITIVITY, CORES):
#
def cnv_calling(args):
    # argparsing
    panel = args.panel
    intsv_path = args.coverages
    analysis_dir = args.analysis_directory
    matchscores_path = args.matchscores
    calls_path = args.cnv_calls
    ratio_scores_path = args.ratio_scores
    z_scores_path = args.z_scores
    EXPECTED_CNV_RATE = args.expected_artifacts
    MINIMUM_SAMPLE_SCORE = args.minimum_sample_score
    MINIMUM_SAMPLE_GROUP = args.minimum_group_sizes
    SENSITIVITY = args.sensitivity if args.sensitivity >= 0 and args.sensitivity <= 1 else 0.0
    CORES = min([args.cores, mp.cpu_count()]) if args.cores else mp.cpu_count()

    # load data
    D0 = util.load_dataframe(intsv_path)
    Matchscores = pd.read_csv(matchscores_path, sep="\t", low_memory=False, header=0, index_col=0)

    # adaptive threshold for group sizes
    MINIMUM_SAMPLE_SCORE_O = MINIMUM_SAMPLE_SCORE
    suggested_MINIMUM_SAMPLE_SCORE = np.mean(
        Matchscores.sum(axis=0) / (len(Matchscores.columns) - 1)
    )
    if suggested_MINIMUM_SAMPLE_SCORE < MINIMUM_SAMPLE_SCORE:
        MINIMUM_SAMPLE_SCORE = (suggested_MINIMUM_SAMPLE_SCORE + MINIMUM_SAMPLE_SCORE * 2) / 3
    print("chose %.04f as MINIMUM_SAMPLE_SCORE" % MINIMUM_SAMPLE_SCORE)

    # define calling groups
    matchmatrix = Matchscores < MINIMUM_SAMPLE_SCORE

    # select all underfitting samples
    sample_group_sizes = matchmatrix.sum(axis=0) - 1
    selected_samples_matrix = matchmatrix[sample_group_sizes >= MINIMUM_SAMPLE_GROUP]
    failed_samples = matchmatrix.index.difference(selected_samples_matrix.index)
    Matchscores_selected = Matchscores[selected_samples_matrix.index]
    Matchscores = Matchscores_selected.T[selected_samples_matrix.index].T
    Matchscores_bools = Matchscores <= MINIMUM_SAMPLE_SCORE
    Matchscores_bools = Matchscores_bools  # * (1-np.identity(len(Matchscores_bools)))).astype(bool)

    # evaluation
    plt.figure(figsize=(8, 4), dpi=200)
    n, bins, patches = plt.hist(
        sample_group_sizes,
        int(len(sample_group_sizes) / 2),
        label=str(len(failed_samples)) + " failed samples",
    )
    plt.vlines(
        MINIMUM_SAMPLE_GROUP,
        0,
        max(n) * 0.75,
        label="cutoff at %.d" % MINIMUM_SAMPLE_GROUP,
        color="darkred",
    )
    plt.xlabel("sample group size")
    plt.ylabel("samples")
    plt.legend()
    plt.savefig(analysis_dir + "ANALYSIS_group_sizes_cutoff.pdf", format="pdf")

    print(
        "%.d samples failed to have a sufficient group size. They are ignored in cnv calling."
        % len(failed_samples)
    )
    for f in failed_samples:
        print(f)

    util.print_clustermap(
        Matchscores,
        path=(analysis_dir + "ANALYSIS_clustermap_filtered.pdf"),
        title="clustered heat map of selected sample distances",
    )

    # separate chr subsets
    DS0 = util.filter_by_coverage(D0[Matchscores.columns])
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
    for i in range(len(Matchscores_bools.columns)):
        pool.apply_async(
            util.calling_cnv,
            args=(i, DA[0], Matchscores_bools, DA[1], EXPECTED_CNV_RATE, SENSITIVITY,),
            callback=collect_result,
        )
    pool.close()
    pool.join()

    #  calculate expected X-cov ratio and prepare match scores
    if sum(DX[1]) > 0 and np.max(DX[0].mean(axis=0)) / np.min(DX[0].mean(axis=0)) > 1.5:
        kmeans = KMeans(n_clusters=2).fit(DX[0].mean(axis=0).reshape(-1, 1))
        Mx = Matchscores_bools.copy()
        for i, x in enumerate(kmeans.labels_):
            sample = Matchscores_bools.index[i]
            if x:
                Mx[sample] &= kmeans.labels_ == x
            else:
                Mx[sample] &= (1 - kmeans.labels_) == x

        sample_scores_buffer = sample_scores.copy()
        #  --- parallelized X-chromosome calling ---  #
        print("multi threaded calling on X-chromosome")
        pool = mp.Pool(CORES)
        for i in range(len(Matchscores_bools.columns)):
            pool.apply_async(
                util.calling_cnv,
                args=(i, DX[0], Matchscores_bools, DX[1], EXPECTED_CNV_RATE, SENSITIVITY,),
                callback=collect_result,
            )
        pool.close()
        pool.join()
        sample_scores += sample_scores_buffer

    #  calculate expected Y-cov ratio and prepare match scores
    if sum(DY[1]) > 0 and np.max(DY[0].mean(axis=0)) / np.min(DY[0].mean(axis=0)) > 1.5:
        kmeans = KMeans(n_clusters=2).fit(DY[0].mean(axis=0).reshape(-1, 1))
        My = Matchscores_bools.copy()
        for i, y in enumerate(kmeans.labels_):
            sample = Matchscores_bools.index[i]
            My[sample] &= kmeans.labels_ == y

        sample_scores_buffer = sample_scores.copy()
        #  --- parallelized Y-chromosome calling ---  #
        print("multi threaded calling on Y-chromosome")
        pool = mp.Pool(CORES)
        for i in range(len(Matchscores_bools.columns)):
            pool.apply_async(
                util.calling_cnv,
                args=(i, DY[0], My, DY[1], EXPECTED_CNV_RATE, SENSITIVITY,),
                callback=collect_result,
            )
        pool.close()
        pool.join()
        sample_scores += sample_scores_buffer

    print("compute HMM on z-scores")
    # HMM guided CNV candidates
    # =========================================================================
    # --- fix lockdown bug --- #
    # pool = mp.Pool(CORES)
    probs = np.array([[-3.0], [0.0], [5.0]])
    transitionprobs = np.array([[0.99, 0.01, 0.0], [0.001, 0.998, 0.001], [0.0, 0.001, 0.999]])
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
                    args=([HMM[i, :], z_scores_scaled[i, :], ratio_scores[i:,], INDEX]),
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
    print("saving results...")
    RD.to_csv(calls_path, sep="\t", index=False)
    RZ.to_csv(z_scores_path, sep="\t")
    RR.to_csv(ratio_scores_path, sep="\t")
    print("done!")
    return RD, RZ, RR
