#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import regex as re
from sklearn.neighbors import KernelDensity
from sklearn.linear_model import LinearRegression
import seaborn as sns
import scipy.cluster.hierarchy as hc
from hmmlearn import hmm
from clearCNV import cnv_arithmetics as ca


def load_dataframe(path, sep="\t", annotated=None):
    D0 = pd.read_csv(path, sep="\t", low_memory=False)
    D0["index"] = [
        str(D0.iloc[i, 0])
        + "_"
        + str(D0.iloc[i, 1])
        + "_"
        + str(D0.iloc[i, 2])
        + "_"
        + str(D0.iloc[i, 3])
        for i in range(len(D0.index))
    ]
    if not (
        D0.columns[0] in ["chr", "Chr", "chromosome", "Chromosome"]
        and D0.columns[1] in ["start", "Start"]
        and D0.columns[2] in ["end", "End"]
        and D0.columns[3] in ["gene", "Gene"]
    ):
        raise Exception(
            "File needs to be formatted like: chr\tstart\tend\tgene\tsample_0\tsample_1\tsample_n"
        )
    D0.set_index("index", inplace=True)
    D0 = D0.iloc[:, 4:]
    D0.index = [s[3:] if s.startswith("chr") else s for s in D0.index]
    if annotated:
        D0 = D0.mul(annotated["mapp"] >= 1.0, axis=0)
    return D0


def load_annotated(path, header=["chr", "start", "end", "gene", "size", "mapp", "gc"]):
    annotated = pd.read_csv(path, sep="\t", low_memory=False, header=None)
    annotated.columns = header
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
    annotated.index = [s[3:] if s.startswith("chr") else s for s in annotated.index]
    return annotated


def trimmed_std(df, axis=0, alpha=0.05):
    if axis not in (0, 1):
        raise ValueError("Invalid axis %d" % axis)
    if axis == 0:
        n = len(df.index)
    else:
        n = len(df.columns)
    df = df.copy()
    df.values.sort(axis=axis)
    lower = math.ceil(n * alpha)
    upper = math.floor(n * (1.0 - alpha))
    if axis == 0:
        return np.asarray(np.std(df.iloc[lower:upper, :], axis=axis, ddof=1))
    else:
        return np.asarray(np.std(df.iloc[:, lower:upper], axis=axis, ddof=1))


def trimmed_std_np(np_matrix, axis=1, alpha=0.1):
    D = np.sort(np_matrix, axis=axis)
    if axis:
        S = D[:, math.ceil(D.shape[axis] * alpha / 2) : -math.ceil(D.shape[axis] * alpha / 2)]
    else:
        S = D[math.ceil(D.shape[axis] * alpha / 2) : -math.ceil(D.shape[axis] * alpha / 2), :]
    return np.std(np.asarray(S), ddof=1, axis=axis)


def filter_by_coverage(D0_, threshold=5):
    # samples
    D0_ = D0_.transpose()[D0_.median(axis=0) > threshold].transpose()
    # exons
    D0_ = D0_[D0_.median(axis=1) > threshold]
    return D0_


def trim(D, axis=0, alpha=0.1):
    if axis:
        return D[
            trimmed_std_np(D, axis=axis) < np.quantile(trimmed_std_np(D, axis=axis), 1 - alpha), :
        ]
    else:
        return D[
            :, trimmed_std_np(D, axis=axis) < np.quantile(trimmed_std_np(D, axis=axis), 1 - alpha)
        ]


# within sample
def normalize(D, axis=1):
    return np.apply_along_axis(lambda v: v / np.median(v[np.nonzero(v)]), axis, np.asarray(D))


def normalize_within_exon(D0_):
    c, i = D0_.columns, D0_.index
    x = np.apply_along_axis(lambda v: v / np.median(v[np.nonzero(v)]), 1, np.asarray(D0_))
    return pd.DataFrame(x, columns=c, index=i)


def normalize_within_sample(D0_):
    c, i = D0_.columns, D0_.index
    x = np.apply_along_axis(lambda v: v / np.median(v[np.nonzero(v)]), 0, np.asarray(D0_))
    return pd.DataFrame(x, columns=c, index=i)


def normalize_chromosomewise(DA, INDEX):
    chroms = np.array([s.split("_")[0] for s in INDEX[DA[1]]])
    chroms_set = sorted([int(s) for s in set(chroms)])
    DD = np.zeros(DA[0].shape)
    chroms_mat = np.zeros((len(chroms), len(set(chroms)))).astype(bool)
    for i, c in enumerate(chroms_set):
        chroms_mat[:, i] = [v == str(c) for v in chroms]
    for i, c in enumerate(chroms_set):
        x = normalize(DA[0][chroms_mat[:, i]], axis=1)
        DD[chroms_mat[:, i], :] = x + 1 - np.median(x, axis=0)
    return DD


def matchscore(i, df, expected_cnv_rate):
    scores = np.zeros(len(df.columns))
    for j in range(len(df.columns)):
        if i > j:
            scores[j] = np.sum(
                sorted(abs(df.iloc[:, i] - df.iloc[:, j]))[
                    0 : math.ceil((1 - expected_cnv_rate) * len(df.index))
                ]
            ) / ((1 - expected_cnv_rate) * len(df.index))
    return (i, scores)


def print_clustermap(df, path, title, fileformat="pdf"):
    linkage = hc.linkage(df, method="average")
    clumap = sns.clustermap(df, row_linkage=linkage, col_linkage=linkage, figsize=(12, 12))
    ax = clumap.ax_heatmap
    ax.set_yticks(ticks=[])
    ax.set_xticks(ticks=[])
    plt.title(title)
    plt.savefig(path, format=fileformat)


# def hmm_predict(D):
#    df = D.copy()
#    model = hmm.GaussianHMM(n_components=3, covariance_type="full")
#    model.startprob_ = np.array([0.0001, 0.9998, 0.0001])
#    model.transmat_ = np.array(
#        [[0.9999, 0.0001, 0.0], [0.0001, 0.9998, 0.0001], [0.0, 0.0001, 0.9999]]
#    )
#    model.means_ = np.array([[-2], [0.0], [2]])
#    model.covars_ = np.tile(np.identity(1), (3, 1, 1))
#    for col in D.columns:
#        df[col] = model.predict([[s] for s in D[col]])
#    return df


def escore(target, rs):
    k = math.ceil(len(target) * 0.025)
    a = np.array(sorted(target)[k:-k]).reshape(-1, 1)
    kde = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(a)
    s = np.linspace(0, 2, 100)
    e = kde.score_samples(s.reshape(-1, 1))
    n = np.where(e == max(e))[0][0]
    # right side
    if rs >= 1:
        x = s[n:].reshape((-1, 1))
        y = e[n:]
        Lmodel = LinearRegression().fit(x, y)
        return (np.exp(-Lmodel.predict(np.array(rs).reshape(-1, 1))) - 1)[0]
    # left side
    else:
        x = s[:n].reshape((-1, 1))
        y = e[:n]
        Lmodel = LinearRegression().fit(x, y)
        return (-np.exp(-Lmodel.predict(np.array(rs).reshape(-1, 1))) + 1)[0]


# takes np.array as input objects
def escore_sample(D, sample):
    rs = np.zeros(len(sample))
    for i, t in enumerate(sample):
        rs[i] = escore(D[:, i], t)
    return rs


def pcolor(D, factor=0.08):
    plt.figure(figsize=(factor * len(D.columns), factor * len(D.index)))
    plt.pcolor(D)
    plt.show()


def npcolor(D, factor=0.08):
    pcolor(pd.DataFrame(D), factor)


def hmm_scores(
    X,
    means=np.array([[-3.5], [0.0], [4.5]]),
    transitionprobs=np.array([[0.999, 0.001, 0.0], [0.001, 0.998, 0.001], [0.0, 0.001, 0.999]]),
):
    model = hmm.GaussianHMM(n_components=3, covariance_type="full")
    model.startprob_ = np.array([0.001, 0.998, 0.001])
    model.transmat_ = np.array(transitionprobs)
    model.means_ = np.array(means)
    model.covars_ = np.tile(np.identity(1), (3, 1, 1))
    return model.predict([np.array([s]) for s in X])


def score_cnv(start, end, z_scores):
    x = np.abs(z_scores[start:end])
    return np.mean(x)  # sum(np.log(x)/np.log(sum(x)))


def rscore_cnv(start, end, r_scores):
    x = r_scores[start:end]
    return np.mean(x)


# hmm_preds,z_scores,r_scores,index = HMM[9],z_scores_scaled[9],ratio_scores[9],INDEX
def merge_score_cnvs(hmm_preds, z_scores, r_scores, index, del_cutoff, dup_cutoff):
    CNVs = []
    length = 0 if hmm_preds[0] == 1 else 1
    start = 0
    for i in range(1, len(hmm_preds)):
        if hmm_preds[i] != hmm_preds[i - 1]:
            # w -> cnv
            if hmm_preds[i - 1] == 1:
                length = 1
                start = i
            # cnv -> w or cnv
            if hmm_preds[i - 1] != 1:
                rs = rscore_cnv(start, i, r_scores)
                # ensure that the ratio score is sufficiently high
                if rs <= del_cutoff or rs >= dup_cutoff:
                    CNVs.append(
                        ca.CNV(
                            *index[start].split("_")[:2],
                            *index[i - 1].split("_")[2:4],
                            "DUP" if hmm_preds[i - 1] == 2 else "DEL",
                            length - 1,
                            score_cnv(start, i, z_scores)
                        )
                    )
                # cnv -> cnv
                if abs(hmm_preds[i] - hmm_preds[i - 1]) > 1:
                    length = 1
                    start = i
        if hmm_preds[i] != 1:
            length += 1
    return CNVs


def center_samples(D, limit=50):
    INDEX, SAMPLES = D.index, D.columns
    if D.shape[0] > limit:
        D = normalize(D, axis=0)
    D = pd.DataFrame(normalize(D, axis=1) + 1 - np.median(normalize(D, axis=1), axis=0))
    D.index, D.columns = INDEX, SAMPLES
    return D


def center_samples_np(D, limit=50):
    if D[0].shape[0] > limit:
        D = normalize(D, axis=0)
    return normalize(D, axis=1) + 1 - np.median(normalize(D, axis=1), axis=0)


def calling_cnv(index_sample, D, MS, index, EXPECTED_CNV_RATE, SENSITIVITY=0.7):
    sample_group = MS.iloc[index_sample]
    # per target median normalized sample
    # select sample group
    D1 = D[:, sample_group]
    # normalize sample per target
    ratio_values = D[:, index_sample] / np.median(D1, axis=1)
    # normalize sample group per target
    D2 = normalize(D1, axis=1)
    np.nan_to_num(D2, copy=False, nan=1.0)
    D3 = trim(D2, axis=0, alpha=0.1)

    # new metric instead of z-scores
    zscores = (ratio_values - 1.0) / (
        trimmed_std_np(D3, axis=1, alpha=EXPECTED_CNV_RATE) ** (2 - SENSITIVITY)
    )

    sample_score = 1.0 / (
        trimmed_std_np(zscores[:, None], axis=0, alpha=EXPECTED_CNV_RATE)[0] ** (2 - SENSITIVITY)
    )
    # print(index_sample,D1.shape,sample_score)
    z_scores_scaled = zscores * sample_score
    return (index_sample, index, sample_score, ratio_values, z_scores_scaled)


def turntransform(v):
    v -= min(v)
    v /= max(v)
    v += np.arange(1, 0, -1 / len(v))
    return v


def select_region(region, df, sample, matchscores_bools_selected, buffer=2):
    cols = df.columns
    if sample:
        # cols = np.append(np.array(df.loc[:,matchscores_bools_selected.loc[:,sample]].columns),sample)
        cols = np.array(df.loc[:, matchscores_bools_selected.loc[:, sample]].columns)
    target = re.split("[^0-9]", region.replace(" ", ""))
    index = np.array(
        [
            (
                c[0] == target[0]
                and int(c[1]) >= int(target[1]) - 5
                and int(c[2]) <= int(target[2]) + 5
            )
            for c in [s.split("_") for s in df.index]
        ]
    )
    n = index.shape[0]
    lower = np.argmax(index == 1)
    upper = np.argmax(index[lower:] == 0) + lower
    index[max(0, lower - buffer) : min(upper + buffer, n)] = True
    # df.loc[index,cols]
    r = df.loc[index, cols].copy()
    if not sample:
        return r
    spacer = pd.DataFrame(
        [float(round(np.mean(r.mean())))] * r.shape[0], index=r.index, columns=["spacer"]
    )
    rsample = pd.DataFrame(df.loc[index, sample])
    return pd.concat([r.T, spacer.T, rsample.T]).T
