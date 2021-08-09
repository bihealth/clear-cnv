"""Setup the SCelVis Dash application.

When importing this module, the app is built and configurd.  Thus, it is important that before
this module is imported, the values in ``.settings`` must already have been setup.
"""

import os

import dash
import flask
from flask import request, helpers
from logzero import logger


from . import callbacks, settings, store, cache
from ..__init__ import __version__
from .ui import build_layout

#: Path to assets.
ASSETS_FOLDER = os.path.join(os.path.dirname(__file__), "assets")

#: The Flask application to use.
app_flask = flask.Flask(__name__)

#: The Dash application to run.
app = dash.Dash(
    __name__,
    # Use our specific Flask app
    server=app_flask,
    # All files from "assets/" will be served as "/assets/*"
    assets_folder=ASSETS_FOLDER,
)

# Setup the cache.
cache.setup_cache(app)

# Set app title
app.title = "clear-CNV v%s" % __version__

# Serve assets locally
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

app.config.suppress_callback_exceptions = True

# Setup the application's main layout.
app.layout = build_layout()

# Register the callbacks with the app.
#
callbacks.register_control_to_buffer(app)
callbacks.register_buffer_to_graphs(app)
callbacks.register_save_button(app)
