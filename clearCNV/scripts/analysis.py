#!/usr/bin/env python
# coding: utf-8

# In[43]:

import argparse
import math
import pandas as pd
import numpy  as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from scipy.spatial.distance import cdist
import os.path
import sys


# In[46]:


r"""panel = "TAADv1"
intsv_path       = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/%s/cov/coverages.tsv"%panel
# all analysis results are saved here
analysis_dir     = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/%s/analysis/"%panel
matchscores_path = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/%s/analysis/matchscores.tsv"%panel
expected_cnv_rate= 0.02"""


# In[42]:

parser = argparse.ArgumentParser(description="Computes quality control plots into analysis directory")
parser.add_argument("-p", "--panel", help="Name of the data set (or panel)",required=True, type=str)
parser.add_argument("-c", "--coverages", help="Coverages file in tsv format",required=True, type=str)
parser.add_argument("-a", "--analysis_directory", help="Output path to the directory, where analysis files are stored", required=True, type=str)
parser.add_argument("-m", "--matchscores", help="matchscores.tsv file generated with matchscores.py", required=True, type=str)
parser.add_argument("-x", "--expected_artifacts", help="Expected ratio of CNVs or artifacs in target fragment counts", required=False, type=float, default=0.02)

args = parser.parse_args()

panel             = args.panel
intsv_path        = args.coverages
analysis_dir      = args.analysis_directory
matchscores_path  = args.matchscores
expected_cnv_rate = args.expected_artifacts

# In[45]:


r"""intsv_path       = snakemake.input[0]
matchscores_path = snakemake.input[1]

# all analysis results are saved here
analysis_dir     = snakemake.params[0]
expected_cnv_rate= float(snakemake.params[1])"""


# In[47]:


D0    = pd.read_csv(intsv_path, sep='\t', low_memory=False)


# In[48]:


D0['index']=[str(D0.iloc[i,0])+'_'+str(D0.iloc[i,1])+'_'+str(D0.iloc[i,2])+'_'+str(D0.iloc[i,3]) for i in range(len(D0.index))]
D0.set_index('index', inplace=True)
D0 = D0.iloc[:,5:]


# In[50]:


D0a = D0[[(not (s.startswith('chrX') or s.startswith('X') or s.startswith('chrY') or s.startswith('Y'))) for s in D0.index]]
D0s = D0[[(s.startswith('chrX') or s.startswith('X')) for s in D0.index]]


# In[51]:


def filter_by_coverage(D0_):
    # samples
    D0_ = D0_.transpose()[D0_.median(axis = 0) > 1].transpose()
    # exons
    D0_ = D0_[D0_.median(axis = 1) > 1]
    return D0_

def normalize_within_exon(D0_):
    return D0_.div(D0_.median(axis=1), axis=0)

def normalize_within_sample(D0_):
    return D0_.div(D0_.median(axis=0), axis=1)

def trimmed_std(df, axis=0, alpha=0.05):
    """Trimmed standard deviation along axis with alpha-trimming."""
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
        return np.asarray(np.std(df.iloc[lower:upper,:], axis=axis, ddof=1))
    else:
        return np.asarray(np.std(df.iloc[:,lower:upper], axis=axis, ddof=1))

def matchscore(i, df, expected_cnv_rate):
    scores = np.zeros(len(df.columns))
    for j in range(len(df.columns)):
        if i > j:
            scores[j] = np.sum(sorted(abs(df.iloc[:,i] - df.iloc[:,j]))[0:math.ceil((1-expected_cnv_rate)*len(df.index))]) / ((1-expected_cnv_rate) * len(df.index))
    return (i,scores)

def collect_result(result):
    global preMatchscores
    preMatchscores.iloc[:,result[0]] = result[1]


# In[52]:


D1a = filter_by_coverage(D0a)
D2a = normalize_within_sample(D1a)
D3a = normalize_within_exon(D2a)


# In[53]:


D1s = filter_by_coverage(D0s)
D2s = normalize_within_sample(D1s)
D3s = normalize_within_exon(D2s)


# In[54]:


D1 = filter_by_coverage(D0)
D2 = normalize_within_sample(D1)
D3 = normalize_within_exon(D2)


# In[11]:


r"""tstep = len(D3a.index) / 20
cstep = len(D3a.columns) / 20

plt.figure(figsize=(cstep,tstep))
_ = plt.pcolor(D3a)
_ = plt.yticks(range(0,len(D3a.index),math.ceil(tstep/10)),[s.split('_')[3] for s in list(D3a.index)][::math.ceil(tstep/10)])
_ = plt.xticks([list(D3a.columns).index(s) for s in D3a.std(axis=0).sort_values()[-10:].index],D3a.std(axis=0).sort_values()[-10:].index,rotation=90)
_ = plt.title("Whole data set, showing top 10 std samples")
plt.savefig(analysis_dir+"whole_data.pdf", format='pdf')
"""


# In[55]:


if os.path.isfile(matchscores_path):
    Matchscores = pd.read_csv(matchscores_path, sep='\t', index_col = 0, low_memory=False)
    print("found matchscores file.")
else:
    preMatchscores = pd.DataFrame(np.zeros((len(D3a.columns),len(D3a.columns))))
    pool = mp.Pool(mp.cpu_count())
    print("calculate matchscores")
    for i in range(len(D3a.columns)):
        pool.apply_async(matchscore, args=(i, D3a, expected_cnv_rate), callback=collect_result)
    pool.close()
    pool.join()
    Matchscores = preMatchscores.add(preMatchscores.T)
    Matchscores.index,Matchscores.columns = D3a.columns, D3a.columns
    # save match scores
    Matchscores.to_csv(matchscores_path, sep='\t', header=True, index=True)
    print("matchscores saved to "+matchscores_path)


# In[13]:


# Analyze


# In[59]:


x_ratios = D2.loc[D3s.index,:].mean()/D2.loc[D3a.index,:].mean()


# In[64]:


if x_ratios.max()/x_ratios.min() > 1.5:
    print("ChrX is unequally distributed.")
    # assume there are unequal X-chromosome numbers
    N = [2,3,4,5]
    X = []
    scores = np.zeros(len(N))
    for i,n in enumerate(N):
        kmeans = KMeans(n_clusters=n).fit(np.array(x_ratios).reshape(-1, 1))
        X.append(kmeans)
        scores[i] = sum([sorted(kmeans.cluster_centers_)[i+1][0]-sorted(kmeans.cluster_centers_)[i][0] for i in range(len(kmeans.cluster_centers_)-1)])/len(kmeans.cluster_centers_)
    kmeans = X[list(scores).index(max(list(scores)))]
    n_clusters = N[list(scores).index(max(list(scores)))]
    print("Chose %d clusters for X-copy numbers"%(n_clusters))
    X=[]

    mean_x_ratios = [x_ratios[kmeans.labels_ == i].mean() for i in range(n_clusters)]
    H = [x_ratios[kmeans.labels_ == i] / min(mean_x_ratios) for i in range(n_clusters)]

    plt.figure(figsize=(6,4), dpi=200)
    plt.title("Clustered distributions of ChrX.")
    plt.xlabel("Average coverage ratio of ChrX-targets")
    for h in H:
        plt.hist(h, bins=40,range=(min([h.min() for h in H]),max([h.max() for h in H])))
    plt.savefig(analysis_dir+"ANALYSIS_chrX_clustering.pdf", format='pdf')
else:
    print("No targtes on chromosome X found, skipping X-ratio calculations.")


# In[65]:


std_trimmed_within_exon = trimmed_std(D3a, axis=1, alpha=0.05)
limit = np.quantile(std_trimmed_within_exon, 0.95)

plt.figure(figsize=(6,4), dpi=200)
n, bins, patches = plt.hist(std_trimmed_within_exon, 50, range=(0,1), label="target deviations")
plt.vlines(limit, 0, max(n)*0.75, label="%.3f = %.0fth percentile."%(limit, 100.0* 0.95),color='darkred')

for i in range(sum(sum([bins < limit])), len(patches)):
    plt.setp(patches[i], 'facecolor', 'grey')

plt.xlabel("deviations")
plt.title("Deviations of normalised fragments per exon counts")
leg = plt.legend()
plt.savefig(analysis_dir+"ANALYSIS_exon_deviations.pdf", format='pdf')


# In[66]:


std_trimmed_within_sample = trimmed_std(D3, axis=0, alpha=0.05)
limit = np.quantile(std_trimmed_within_sample, 0.95)

plt.figure(figsize=(6,4), dpi=200)
n, bins, patches = plt.hist(std_trimmed_within_sample, 20, label="sample deviations")
plt.vlines(limit, 0, max(n)*0.75, label="%.3f = %.0fth percentile"%(limit,100.0* 0.95),color='darkred')

for i in range(sum(sum([bins < limit])), len(patches)):
    plt.setp(patches[i], 'facecolor', 'grey')
plt.xlabel("deviations")
plt.ylabel("samples")
plt.title("Deviations of normalised fragments per sample counts")
leg = plt.legend()
plt.savefig(analysis_dir+"ANALYSIS_sample_deviations.pdf", format='pdf')


# In[68]:


plt.figure(figsize=(6,4), dpi=200)
#plt.hist(np.sqrt(Matchscores.sum(axis=1)), max(10, math.ceil(len(D3.columns)/10)),label='sqrt(matchscores)')
plt.hist(sorted(Matchscores.sum(axis=1) / (len(D3s.columns)-1)), max(10, math.ceil(len(D3s.columns)/10)),label='sqrt(matchscores)')
plt.xlabel("normalized average scores")
plt.ylabel("samples")
plt.title("Match Scores Distribution.")
plt.legend()
plt.savefig(analysis_dir+"ANALYSIS_match_scores_distribution.pdf", format='pdf')


# In[28]:


linkage = hc.linkage(Matchscores, method='average')
#sns.set(color_codes=True)
clumap = sns.clustermap(Matchscores, row_linkage=linkage, col_linkage=linkage, figsize=(12, 12))
ax = clumap.ax_heatmap
ax.set_yticks(ticks=[])
ax.set_xticks(ticks=[])
plt.title("clustered heat map of sample distances")
plt.savefig(analysis_dir+"ANALYSIS_clustermap_unfiltered.pdf", format='pdf') # used in snakemake!
