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
from .settings import reassignSettings, BATCH_OUTPUT_PATH


def register_control_to_buffer(app):
    """Register callsbacks for controls."""

    @app.callback(
        Output("buffer-cluster-params", "children"),
        [
            Input("input-thresh", "value"),
            Input("input-pca-components", "value"),
            Input("input-batch-factor", "value"),
        ],
    )
    def parameters_to_json(threshold, pca_components, batch_num):
        return json.dumps(
            {
                "threshold": threshold,
                "pca_components": pca_components,
                "batch_num": batch_num,
            }
        )


def register_buffer_to_graphs(app):
    """Register callsbacks for controls."""

    @app.callback(
        [
            Output("graph-pca", "figure"),
            Output("graph-tsne", "figure"),
            Output("graph-agg-clust", "figure"),
            Output("container-image-cluster-panels", "children"),
            Output("container-image-cluster-clustering", "children"),
            # Output("container-image-batch-clustering", "children"),
            Output("container-image-batch-separation", "children"),
        ],
        [Input("buffer-cluster-params", "children")],
    )
    def parameters_to_json(params_str):
        us = reassignSettings(**json.loads(params_str))
        return [
            ui.render_pca(us),
            ui.render_tsne(us),
            ui.render_agg_clust(us),
            ui.render_image_cluster_panels(us),
            ui.render_image_cluster_clustering(us),
            # ui.render_image_batches_clustering(us),
            ui.render_images_batch_separation(us),
        ]


def register_save_button(app):
    """Register button for the save button."""

    @app.callback(
        Output("buffer-save-batches-text", "children"),
        [
            Input("input-save-batches-button", "n_clicks"),
            Input("buffer-cluster-params", "children"),
        ],
    )
    def on_click_button(n_clicks, params_str):
        if n_clicks is not None:
            us = reassignSettings(**json.loads(params_str))
            print("clicked save", n_clicks, us)
            # how to get us ????
            store.save_results(us, n_clicks)
            return "%d: saving to %s" % (n_clicks or 0, BATCH_OUTPUT_PATH)
