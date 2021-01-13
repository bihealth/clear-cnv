import base64
import io
import tempfile

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
    plt.figure(figsize=(8, 5))
    sns.clustermap(data.X, col_colors=list(map(lambda x: panelcoldict[x], data.samples["panel"])))
    logger.info("... done creating panels clustermap.")
    with tempfile.NamedTemporaryFile(suffix=".png") as f:
        plt.savefig(f, format="png")
        f.flush()
        f.seek(0)
        return base64.b64encode(f.read()).decode("ascii")


@cache.memoize()
def plot_clustermap_clustering_as_base64(data, XD, clustercoldict):
    logger.info("Creating clustering clustermap ...")
    plt.figure(figsize=(8, 5))
    sns.clustermap(data.X, col_colors=list(map(lambda x: clustercoldict[x], XD["clustering"])))
    logger.info("... done creating clustering clustermap.")
    with tempfile.NamedTemporaryFile(suffix=".png") as f:
        plt.savefig(f, format="png")
        f.flush()
        f.seek(0)
        return base64.b64encode(f.read()).decode("ascii")
