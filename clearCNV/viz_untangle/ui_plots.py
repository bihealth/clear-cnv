import base64
import io
import tempfile

import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from logzero import logger

import seaborn as sns

sns.set_theme(color_codes=True)


from .cache import cache

# matplotlib.use('Agg')


@cache.memoize()
def plot_clustermap_panels_as_base64(data, panelcoldict):
    logger.info("Creating panels clustermap ...")
    n = min([data.X.shape[0], 1000])
    m = min([data.X.shape[1], 1000])
    plt.figure(figsize=(8, 5))
    sns.clustermap(
        data.X.iloc[:: math.ceil(data.X.shape[0] / n), :: math.ceil(data.X.shape[1] / m)],
        col_colors=list(map(lambda x: panelcoldict[x], data.samples["panel"])),
    )
    logger.info("... done creating panels clustermap.")
    with tempfile.NamedTemporaryFile(suffix=".png") as f:
        plt.savefig(f, format="png")
        f.flush()
        f.seek(0)
        return base64.b64encode(f.read()).decode("ascii")


@cache.memoize()
def plot_clustermap_clustering_as_base64(data, XD, clustercoldict):
    logger.info("Creating clustering clustermap ...")
    n = min([data.X.shape[0], 1000])
    m = min([data.X.shape[1], 1000])
    plt.figure(figsize=(8, 5))
    sns.clustermap(
        data.X.iloc[:: math.ceil(data.X.shape[0] / n), :: math.ceil(data.X.shape[1] / m)],
        col_colors=list(map(lambda x: clustercoldict[x], XD["clustering"])),
    )

    logger.info("... done creating clustering clustermap.")
    with tempfile.NamedTemporaryFile(suffix=".png") as f:
        plt.savefig(f, format="png")
        f.flush()
        f.seek(0)
        return base64.b64encode(f.read()).decode("ascii")

@cache.memoize()
def plot_clustermap_batches_as_base64(data, XD, clustercoldict):
    logger.info("Creating batches clustermap ...")
    n = min([data.X.shape[0], 1000])
    m = min([data.X.shape[1], 1000])
    plt.figure(figsize=(8, 5))
    data = load_all_data(us)
    D1 = data.D1
    D2 = util.normalize_within_sample(D1)
    D3 = util.normalize_within_exon(D2)
    sns.clustermap(
        D3.iloc[:: math.ceil(D3.shape[0] / n), :: math.ceil(D3.shape[1] / m)],
        col_colors=list(map(lambda x: clustercoldict[x], XD["clustering"])),
    )
    logger.info("... done creating batches clustermap.")
    with tempfile.NamedTemporaryFile(suffix=".png") as f:
        plt.savefig(f, format="png")
        f.flush()
        f.seek(0)
        return base64.b64encode(f.read()).decode("ascii")
