#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
import os.path
import math
from clearCNV import util


# In[ ]:


r"""
parser = argparse.ArgumentParser(description="Matchscore calculation script.")
parser.add_argument("-p", "--panel", help="Name of the data set (or panel)",required=True, type=str)
parser.add_argument("-c", "--coverages", help="Coverages file in tsv format",required=True, type=str)
parser.add_argument("-m", "--matchscores", help="Output matchscores.tsv file", required=True, type=str)
parser.add_argument("-x", "--expected_artifacts", help="Expected ratio of CNVs or artifacs in target fragment counts", required=False, type=float, default=0.02)
parser.add_argument("--cores", help="Number of cpu cores used in parallel processing. Default: determined automatically.", required=False, type=int,   default=0)

args = parser.parse_args()

panel             = args.panel
intsv_path        = args.coverages
matchscores_path  = args.matchscores
expected_cnv_rate = args.expected_artifacts
CORES             = min([args.cores, mp.cpu_count()]) if args.cores else mp.cpu_count()
"""


# In[2]:


r"""
panel = "TAADv2"

intsv_path       = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/CALLING/%s/cov/coverages.tsv"%panel
#mappability_path = "/vol/sshfs/vmay/bih_cluster/fast/work/projects/cubit/18.12/static_data/db/goldenpath/variable/GRCh37/wgEncodeCrgMapabilityAlign100mer.eq_1.bed"

matchscores_path = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/CALLING/%s/analysis/matchscores.tsv"%panel

expected_cnv_rate= 0.02
"""


# In[7]:


# def matchscores(panel, intsv_path, matchscores_path, expected_cnv_rate, fast=True):
def matchscores(args):
    # parse args
    panel = args.panel
    intsv_path = args.coverages
    matchscores_path = args.matchscores
    expected_cnv_rate = args.expected_artifacts
    CORES = min([args.cores, mp.cpu_count()]) if args.cores else mp.cpu_count()
    fast = True

    D0 = util.load_dataframe(intsv_path)
    D3a = util.normalize_within_exon(
        util.normalize_within_sample(
            util.filter_by_coverage(
                D0[[not (s.startswith("X") or s.startswith("Y")) for s in D0.index]]
            )
        )
    )
    # center biased samples
    D3a = util.center_samples(D3a)

    # only take every i-th target into account to speed up
    if fast:
        X = D3a.T[D3a.index[:: min([len(D3a.index), math.floor(len(D3a.index) / 1000)])]].T
    else:
        X = D3a

    preMatchscores = pd.DataFrame(np.zeros((len(X.columns), len(X.columns))))

    def collect_result(result):
        nonlocal preMatchscores
        preMatchscores.iloc[:, result[0]] = result[1]

    pool = mp.Pool(CORES)

    print(
        f"calc matchscores with {CORES} threads at an expected artifact rate "
        f"of {expected_cnv_rate}"
    )
    for i in range(len(X.columns)):
        pool.apply_async(util.matchscore, args=(i, X, expected_cnv_rate), callback=collect_result)
    pool.close()
    pool.join()

    Matchscores = preMatchscores.add(preMatchscores.T)
    Matchscores.index, Matchscores.columns = X.columns, X.columns
    Matchscores.to_csv(matchscores_path, sep="\t", header=True, index=True)
    return Matchscores
