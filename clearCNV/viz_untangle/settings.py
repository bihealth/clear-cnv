"""Settings for the Dash-based visual untangling.

While it is general best practice not to have global state such as this settings module, the
Flask-based architecture of Dash forces us to use it.
"""

import pathlib

import attr
import cattr

#: Text to display on top left.
HOME_BRAND = "clear-CNV"

#: Base input path.
PATH_BASE = None

#: Path to metafile.
PATH_METAFILE = None

#: Path to coverages.
PATH_COVERAGES = None

#: Path to targets BED.
PATH_BEDFILE = None

#: Path to batch definition.
PATH_BATCHES = None

#: Helper function to setup the paths based on an input path.
def setup_paths(path):
    global PATH_BASE, PATH_METAFILE, PATH_COVERAGES, PATH_BEDFILE, PATH_BATCHES
    PATH_BASE = pathlib.Path(path)
    PATH_METAFILE = PATH_BASE / "meta.tsv"
    PATH_COVERAGES = PATH_BASE / "coverages.tsv"
    PATH_BEDFILE = PATH_BASE / "union.bed"
    PATH_BATCHES = PATH_BASE / "batches"


#: Whether or not to preload data.
CACHE_PRELOAD_DATA = True


@attr.s(frozen=True, auto_attribs=True)
class UntangleSettings:
    threshold: int = 50
    min_group_size: int = 20
    pca_components: int = 20
    batch_factor: float = 0.985
    pca_seed: int = 100


#: Configuration for untangling.
UNTANGLE_SETTINGS = UntangleSettings()

#: Whether or not to enable debugging.
DEBUG = False

#: The type of the cache to use from {"filesystem", "redis"}
CACHE_TYPE = "filesystem"
#: Default cache timeout.
CACHE_DEFAULT_TIMEOUT = 300
#: For "filesystem" cache: directory to store the cache in.
CACHE_DIR = None
#: For "redis" cache: the URL to use for connecting to the cache.
CACHE_REDIS_URL = None
#: whether or not to preload data on startup
CACHE_PRELOAD_DATA = True
