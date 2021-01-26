"""Handle access to the data."""

import typing

import os
import attr
import cattr
import pickle

import pathlib
import pandas as pd
import numpy as np
from math import sqrt

from collections import Counter
import subprocess

from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import AgglomerativeClustering

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from logzero import logger

from .. import util
from .. import cnv_arithmetics as ca
from .cache import cache
from . import settings
from . import ui_plots


#: Whether or not to pickle data (for development only).
PICKLE_DATA = False


@attr.s(auto_attribs=True, frozen=True)
class Data:
    BED: typing.Any
    D: typing.Any
    META: typing.Any
    X: typing.Any
    allsamples: typing.Any
    dropoutsamples: typing.Any
    samples: typing.Any
    panels: typing.Any
    bedsmatrix: typing.Any


def _load_coverage(path, sep="\t"):
    D0 = pd.read_csv(path, sep="\t", low_memory=False)
    if not (
        D0.columns[0] in ["chr", "Chr", "chromosome", "Chromosome"]
        and D0.columns[1] in ["start", "Start"]
        and D0.columns[2] in ["end", "End"]
    ):
        raise Exception(
            "File needs to be formatted like: chr\tstart\tend\tsample_0\tsample_1\tsample_n"
        )
    D0 = D0.iloc[:, 3:]
    return D0


def _find_batches(xd, savepath, factor=0.99, _n=1, _bic=None, _df=None):
    xds = xd[["X", "Y"]].to_numpy()
    xindex = xd.index
    GM = GaussianMixture(n_components=_n, n_init=10).fit(xds)
    # print(GM.bic(xds))
    xds_df = pd.DataFrame(xds, index=xd.index, columns=["X", "Y"])
    xds_df["clustering"] = GM.predict(xds)
    BIC = GM.bic(xds)

    if _bic and factor * _bic <= BIC:
        # render figures
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=_df, x="X", y="Y", hue="clustering")
        plt.savefig(savepath)

        _df.to_csv(str(savepath) + ".tsv", sep="\t", header=True, index=True)

        return (_n - 1, _bic, _df)
    else:
        return _find_batches(xd, savepath, factor, _n + 1, BIC, xds_df)


@cache.memoize()
def load_all_data(us: settings.UntangleSettings):
    if PICKLE_DATA and os.path.exists("all_data.bin"):
        with open("all_data.bin", "rb") as f:
            return pickle.load(f)
    logger.info("Starting to load data ...")
    logger.info("    - %s", settings.PATH_BEDFILE)
    BED = pd.read_csv(settings.PATH_BEDFILE, sep="\t", header=None).astype(str)
    BED = BED[[b.isdigit() for b in BED.iloc[:, 0]]]
    BED = BED.astype(int)

    logger.info("    - %s", settings.PATH_COVERAGES)
    D = _load_coverage(settings.PATH_COVERAGES)
    D = D.T[BED.index].T

    META = pd.read_csv(settings.PATH_METAFILE, sep="\t")
    META = META.sort_values(by="panel")

    logger.info("    - allsamples")
    allsamples = pd.DataFrame()
    for i in META.index:
        s = META.loc[i, :]
        bamfiles = [line.rstrip("\n") for line in open(s["bamsfiles"])]
        allsamples = pd.concat(
            [
                allsamples,
                pd.DataFrame(
                    list(
                        zip(
                            [s["panel"]] * len(bamfiles),
                            [pathlib.Path(bf).stem for bf in bamfiles],
                            [bf for bf in bamfiles],
                        )
                    )
                ),
            ]
        )
    allsamples = allsamples.set_index(1)
    allsamples.columns = ["panel", "path"]

    D = D[allsamples.index]

    d = D.columns[((D > 1).sum() == 0) | (D.sum() <= us.threshold)]
    dropoutsamples = pd.DataFrame(allsamples.T[d].T["path"])

    D1 = D.drop(columns=dropoutsamples.index)
    D.index = BED.apply(lambda row: "chr%s:%s" % (row[0], row[1]), axis=1)

    # X = D1 > us.threshold
    X = D1.ge((BED.iloc[:, 2] - BED.iloc[:, 1]) / us.threshold, axis=0)

    logger.info("    - bedsmatrix")
    samples = allsamples.drop(dropoutsamples.index)
    panels = sorted(set(META["panel"]))
    bedsmatrix = np.zeros((len(panels), BED.shape[0]))

    for j, panel in enumerate(panels):
        bed_panel_path = META.iloc[j, 2]
        # def merge_bedfiles(beds,bedfile):
        command = [
            "bedtools",
            "intersect",
            "-wa",
            "-a",
            settings.PATH_BEDFILE,
            "-b",
            bed_panel_path,
        ]
        p = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
        )
        output, error = p.communicate()

        btargets = [ca.Target(*BED.iloc[i, :], "gene") for i in BED.index]
        ptargets = [ca.Target(*line.split("\t"), "gene") for line in output.split("\n")[:-1]]
        for i, t in enumerate(btargets):
            if t in ptargets:
                bedsmatrix[j][i] = 1
    bedsmatrix = pd.DataFrame(bedsmatrix, index=panels)

    logger.info("... done loading data.")

    result = Data(
        BED=BED,
        D=D,
        X=X,
        META=META,
        allsamples=allsamples,
        dropoutsamples=dropoutsamples,
        samples=samples,
        panels=panels,
        bedsmatrix=bedsmatrix,
    )
    if PICKLE_DATA:
        with open("all_data.bin", "wb") as f:
            pickle.dump(result, f)
    return result


@cache.memoize()
def compute_pca(us: settings.UntangleSettings):
    logger.info("Computing PCA ...")
    data = load_all_data(us)
    pca_simple = PCA(n_components=2, random_state=us.pca_seed)
    pca2df = pd.DataFrame(pca_simple.fit_transform(data.X.dropna(0).T))
    pca2df.columns = ["X", "Y"]
    pca2df["panel"] = list(data.samples["panel"])
    logger.info("... done computing PCA.")
    return pca2df


@cache.memoize()
def compute_tsne(us: settings.UntangleSettings):
    logger.info("Computing tSNE ...")
    data = load_all_data(us)
    pca = PCA(n_components=us.pca_components, random_state=us.pca_seed,)
    XT = pca.fit_transform(data.X.dropna(0).T)

    X_embedded = TSNE(n_components=2, random_state=us.pca_seed).fit_transform(XT)

    XD = pd.DataFrame(X_embedded)
    XD.columns = ["X", "Y"]
    XD["panel"] = list(data.samples["panel"])
    # for paper demo only
    # XD["panel"] = [paneldict[p] for p in XD["panel"]]
    logger.info("... done computing tSNE.")
    return XD

@cache.memoize()
def compute_acluster(us: settings.UntangleSettings):
    logger.info("Computing Agglomerative Clustering ...")
    data = load_all_data(us)
    XD = compute_tsne(us)
    clustering = AgglomerativeClustering(n_clusters=len(data.panels)).fit(XD[["X", "Y"]].to_numpy())
    # TODO: yikes, inplace update with cached data...
    XD["clustering"] = clustering.labels_
    # majority vote for panel assignments
    cluster_panel_dict = {i:Counter(data.samples[clustering.labels_ == i]["panel"]).most_common()[0][0] for i in range(len(data.panels))}
    XD["new_assignments"] = list(map(lambda x: cluster_panel_dict[x],clustering.labels_))
    logger.info("... done computing Agglomerative Clustering.")
    return XD


@cache.memoize()
def compute_panelcoldict(us: settings.UntangleSettings):
    data = load_all_data(us)
    clustercolors = cm.get_cmap("viridis", len(data.panels))
    panelcoldict = {
        data.panels[i]: c for i, c in enumerate(clustercolors(np.linspace(0, 1, len(data.panels))))
    }
    return panelcoldict


@cache.memoize()
def compute_clustercoldict(us: settings.UntangleSettings):
    data = load_all_data(us)
    clustercolors = cm.get_cmap("viridis", len(data.panels))
    clustercoldict = {
        i: c for i, c in enumerate(clustercolors(np.linspace(0, 1, len(data.panels))))
    }
    return clustercoldict
