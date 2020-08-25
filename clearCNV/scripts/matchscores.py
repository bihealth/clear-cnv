#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import pandas as pd
import numpy  as np
import multiprocessing as mp
import os.path
import math
import sys
import util


# In[ ]:


parser = argparse.ArgumentParser(description="Matchscore calculation script.")
parser.add_argument("-p", "--panel", help="Name of the data set (or panel)",required=True, type=str)
parser.add_argument("-c", "--coverages", help="Coverages file in tsv format",required=True, type=str)
parser.add_argument("-m", "--matchscores", help="Output matchscores.tsv file", required=True, type=str)
parser.add_argument("-x", "--expected_artifacts", help="Expected ratio of CNVs or artifacs in target fragment counts", required=False, type=float, default=0.02)

args = parser.parse_args()

panel             = args.panel
intsv_path        = args.coverages
matchscores_path  = args.matchscores
expected_cnv_rate = args.expected_artifacts


# In[70]:


r"""
panel = "SDAG1"

intsv_path       = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/CALLING/%s/cov/coverages.tsv"%panel
#mappability_path = "/vol/sshfs/vmay/bih_cluster/fast/work/projects/cubit/18.12/static_data/db/goldenpath/variable/GRCh37/wgEncodeCrgMapabilityAlign100mer.eq_1.bed"

matchscores_path = "/vol/sshfs/vmay/bih_cluster/fast/work/users/vmay_m/workflow/FCC/CALLING/%s/analysis/matchscores.tsv"%panel

expected_cnv_rate= 0.02
"""


# In[ ]:


"""intsv_path       = snakemake.input[0]

matchscores_path = snakemake.output[0]

# all analysis results are saved here
expected_cnv_rate= float(snakemake.params[0])"""


# In[71]:


D0 = util.load_dataframe(intsv_path)


# In[72]:


D0a = D0[[(not (s.startswith('chrX') or s.startswith('X') or s.startswith('chrY') or s.startswith('Y') or s.startswith('MT') or s.startswith('M'))) for s in D0.index]]
D0s = D0[[(s.startswith('chrX') or s.startswith('X')) for s in D0.index]]


# In[73]:


def collect_result(result):
    global preMatchscores
    preMatchscores.iloc[:,result[0]] = result[1]


# In[74]:


D1a = util.filter_by_coverage(D0a)
D2a = util.normalize_within_sample(D1a)
D3a = util.normalize_within_exon(D2a)


# In[ ]:


#preMatchscores = pd.DataFrame(np.zeros((len(D3a.columns),len(D3a.columns))))
#pool = mp.Pool(mp.cpu_count())
#print("calculate matchscores")
#for i in range(len(D3a.columns)):
#    pool.apply_async(util.matchscore, args=(i, D3a, expected_cnv_rate), callback=collect_result)
#pool.close()
#pool.join()
#Matchscores = preMatchscores.add(preMatchscores.T)
#Matchscores.index,Matchscores.columns = D3a.columns, D3a.columns
## save match scores
#Matchscores.to_csv(matchscores_path, sep='\t', header=True, index=True)
#print("matchscores saved")
#M0 = Matchscores


# In[149]:


X = D3a.T[D3a.index[::min([len(D3a.index),math.ceil(len(D3a.index)/1000)])]].T


# In[151]:


preMatchscores = pd.DataFrame(np.zeros((len(X.columns),len(X.columns))))
pool = mp.Pool(mp.cpu_count())
print("calculate matchscores")
for i in range(len(X.columns)):
    pool.apply_async(util.matchscore, args=(i, X, expected_cnv_rate), callback=collect_result)
pool.close()
pool.join()
Matchscores = preMatchscores.add(preMatchscores.T)
Matchscores.index,Matchscores.columns = X.columns, X.columns
# save match scores
Matchscores.to_csv(matchscores_path, sep='\t', header=True, index=True)
print("matchscores saved")
M1 = Matchscores
