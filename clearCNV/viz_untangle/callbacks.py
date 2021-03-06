"""Registration of callbacks for the UI."""

import base64
import os.path
import uuid
import json

import dash
import dash_html_components as html
from dash.dependencies import Input, Output, State
import pandas as pd
from logzero import logger
from werkzeug.utils import secure_filename

from . import ui, settings, store
from .settings import UntangleSettings


def register_control_to_buffer(app):
    """Register callsbacks for controls."""

    @app.callback(
        Output("buffer-cluster-params", "children"),
        [
            Input("input-thresh", "value"),
            Input("input-min-group-size", "value"),
            Input("input-pca-components", "value"),
            Input("input-batch-factor", "value"),
        ],
    )
    def parameters_to_json(threshold, min_group_size, pca_components, batch_factor):
        return json.dumps(
            {
                "threshold": threshold,
                "min_group_size": min_group_size,
                "pca_components": pca_components,
                "batch_factor": batch_factor,
            }
        )


def register_buffer_to_graphs(app):
    """Register callsbacks for controls."""

    @app.callback(
        [
            Output("graph-pca", "figure"),
            Output("graph-tsne", "figure"),
            Output("graph-agg-clust", "figure"),
        ],
        [Input("buffer-cluster-params", "children")],
    )
    def parameters_to_json(params_str):
        us = UntangleSettings(**json.loads(params_str))
        return [
            ui.render_pca(us),
            ui.render_tsne(us),
            ui.render_agg_clust(us),
        ]
