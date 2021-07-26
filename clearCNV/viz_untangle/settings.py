"""Settings for the Dash-based visual reassignment.

While it is general best practice not to have global state such as this settings module, the
Flask-based architecture of Dash forces us to use it.
"""

import pathlib

import attr
import cattr

#: Text to display on top left.
HOME_BRAND = "clear-CNV"

#: Base input path.
# PATH_BASE = None

#: Path to metafile.
PATH_METAFILE = None

#: Path to coverages.
PATH_COVERAGES = None

#: Path to targets BED.
PATH_BEDFILE = None

#: Path to batch definition.
BATCH_OUTPUT_PATH = None

#: Helper function to setup the paths based on an input path.
def setup_paths(args):
    global PATH_METAFILE, PATH_COVERAGES, PATH_BEDFILE, BATCH_OUTPUT_PATH
    # PATH_BASE = pathlib.Path(path)
    PATH_METAFILE = args.metafile
    PATH_COVERAGES = args.coverages
    PATH_BEDFILE = args.bedfile
    BATCH_OUTPUT_PATH = args.new_panel_assignments_directory


#: Whether or not to preload data.
CACHE_PRELOAD_DATA = True


@attr.s(frozen=True, auto_attribs=True)
class reassignSettings:
    threshold: int = 50
    pca_components: int = 20
    batch_num: str = "2"
    pca_seed: int = 100
    colormap = "tab20"


#: Configuration for reassignment.
reassign_SETTINGS = reassignSettings()

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

#: Output path for batch separation.
# BATCH_OUTPUT_PATH = None
