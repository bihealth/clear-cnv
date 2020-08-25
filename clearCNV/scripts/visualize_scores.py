#!/usr/bin/env python
# coding: utf-8

# In[16]:

import argparse
import sys
import pandas as pd
import numpy  as np
import matplotlib.pyplot as plt
import math
from sklearn.cluster import KMeans

import plotly.graph_objects as go
import plotly
from plotly.subplots import make_subplots


# In[ ]:


"""analysis_dir = "/vol/sshfs/vmay/bih_cluster/fast/users/vmay_m/work/workflow/FCC/CALLING/BM/analysis/"
z_scores_path = "/vol/sshfs/vmay/bih_cluster/fast/users/vmay_m/work/workflow/FCC/CALLING/BM/analysis/z_scores.tsv"
ratio_scores_path = "/vol/sshfs/vmay/bih_cluster/fast/users/vmay_m/work/workflow/FCC/CALLING/BM/analysis/ratio_scores.tsv"
annotated_path = "/vol/sshfs/vmay/bih_cluster/fast/users/vmay_m/work/workflow/FCC/CALLING/BM/results/annotations.bed"
"""


# In[ ]:

parser = argparse.ArgumentParser(description="CNV calling script. Output is a single cnv.calls file in tsv format. Some quality control plots are added to analysis directory in the process.")
parser.add_argument("-a", "--analysis_directory", help="Path to the directory, where analysis files are stored", required=True, type=str)
parser.add_argument("-r", "--ratio_scores", help="Ratio scores file in tsv format, generated in cnv_calling.py", required=True, type=str)
parser.add_argument("-z", "--z_scores", help="Z-scores file in tsv format, generated in cnv_calling.py", required=True, type=str)
parser.add_argument("-n", "--annotated", help="CALL_GROUPS.py analysis_directory z_scores.tsv ratio_scores.tsv annotations.bed", required=True, type=str)

args = parser.parse_args()

analysis_dir      = args.analysis_directory
ratio_scores_path = args.ratio_scores
z_scores_path     = args.z_scores
annotated_path    = args.annotated

# In[5]:


annotated = pd.read_csv(annotated_path, sep='\t', low_memory=False, header=None)
annotated.columns=["chr","start","end","gene","size","mapp","gc"]
annotated.index = [str(annotated["chr"][i])+'_'+str(annotated["start"][i])+'_'+str(annotated["end"][i])+'_'+str(annotated["gene"][i]) for i in range(len(annotated.index))]


# In[9]:


ratio_scores_df    = pd.read_csv(ratio_scores_path, sep='\t', low_memory=False,index_col = 0)


# In[12]:


z_scores_scaled_df = pd.read_csv(z_scores_path, sep='\t', low_memory=False,index_col = 0)


# In[14]:


#=========================================================================================
# Deep analysis with plotly
#=========================================================================================


# In[19]:


RS = ratio_scores_df.fillna(1.0).copy()
n_clusters=math.ceil(len(RS.columns)/40)
kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(RS[RS.std(axis=1) < np.nanquantile(RS.std(axis=1),0.95)].T)

clustered_samples = pd.DataFrame(list(zip(RS.columns,kmeans.labels_))).sort_values(by=1)[0]

df = RS[clustered_samples].copy()

limit = len(df.index)

htmp = go.Heatmap(
                   z=[np.clip(df.iloc[i,:],0,2) for i in range(0,limit)],
                   x=df.columns,
                   y=df.index[:limit],
                   hoverongaps = False)

mgc   = go.Scatter(y=df.index[:limit], x=annotated['gc'])
mmapp = go.Scatter(y=df.index[:limit], x=annotated['mapp'])
msize = go.Scatter(y=df.index[:limit], x=annotated['size'])

sp = make_subplots(
    rows=1, cols=4,
    subplot_titles=("ratios", "mapp", "gc","log_size"),
    column_widths=[0.85, 0.05, 0.05, 0.05],
    shared_yaxes=True)

sp.add_trace(htmp,
              row=1, col=1)
sp.add_trace(mmapp,
              row=1, col=2)
sp.add_trace(mgc,
              row=1, col=3)
sp.add_trace(msize,
              row=1, col=4)

sp.update_layout(
    width=len(df.columns)*7.5,
    height=limit*5,
)
plotly.offline.plot(sp, filename=analysis_dir+"ratio_scores_extended.html")


# In[21]:


#=========================================================================================
# analyze relation GC content and st. dev. of z-scores after grouping
#=========================================================================================

print("analyze relation GC content and st. dev. of z-scores after grouping.")


# In[22]:


z_scores_scaled_df["std"] = [z_scores_scaled_df.loc[i,:].std() for i in z_scores_scaled_df.index]
annotated["std"] = z_scores_scaled_df["std"]
z_scores_scaled_df = z_scores_scaled_df.drop(columns="std")


# In[23]:


ratio_scores_df["std"] = [ratio_scores_df.loc[i,:].std() for i in ratio_scores_df.index]
annotated["ratio"] = ratio_scores_df["std"]
ratio_scores_df = ratio_scores_df.drop(columns="std")


# In[24]:


plt.figure(figsize=(8,5))
plt.scatter(annotated["std"],annotated["gc"],marker=".",alpha=0.1,label="sd of z-scores")
plt.ylabel("GC content")
plt.xlabel("sd of z-scores")
plt.title("GC content vs z-scores after grouping")
plt.legend()
plt.savefig(analysis_dir+"ANALYSIS_GC_vs_SD.pdf",format='pdf')


# In[25]:


plt.figure(figsize=(8,5))
plt.scatter(annotated["ratio"],annotated["gc"],marker=".",alpha=0.1,label="sd of ratio scores")
plt.ylabel("GC content")
plt.xlabel("sd of ratios")
plt.title("GC content vs sd of ratio scores after grouping")
plt.legend()
plt.savefig(analysis_dir+"ANALYSIS_ratio_vs_gc.pdf",format='pdf')
