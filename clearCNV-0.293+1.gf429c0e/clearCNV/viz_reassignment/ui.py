"""Definitions of Dash layout."""

import os.path

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
from matplotlib.colors import rgb2hex

import seaborn as sns

from .cache import cache
from . import settings
from . import store
from . import ui_plots
from ..__init__ import __version__
from .settings import reassignSettings


def render_navbar():
    """Render the site navbar"""
    return dbc.Navbar(
        dbc.Container(
            children=[
                # Use row and col to control vertical alignment of logo / brand
                dbc.NavbarBrand(
                    [html.I(className="fas fa-home mr-1"), settings.HOME_BRAND], id="page-brand"
                ),
            ]
        ),
        # className="mb-5",
        color="primary",
        dark=True,
        id="page-navbar",
    )


@cache.memoize()
def render_pca(us: reassignSettings):
    return px.scatter(
        store.compute_pca(us).sort_values(by="panel"),
        x="X",
        y="Y",
        color="panel",
        color_discrete_sequence=[rgb2hex(v) for v in store.compute_clustercoldict(us).values()],
    )


@cache.memoize()
def render_tsne(us: reassignSettings):
    return px.scatter(
        store.compute_tsne(us).sort_values(by="panel"),
        x="X",
        y="Y",
        color="panel",
        color_discrete_sequence=[rgb2hex(v) for v in store.compute_clustercoldict(us).values()],
    )


@cache.memoize()
def render_agg_clust(us: reassignSettings):
    return px.scatter(
        store.compute_acluster(us)[0].sort_values(by="new_assignments"),
        x="X",
        y="Y",
        color="new_assignments",
        color_discrete_sequence=[rgb2hex(v) for v in store.compute_clustercoldict(us).values()],
    )


@cache.memoize()
def render_image_cluster_panels(us: reassignSettings):
    data = store.load_all_data(us)
    img_b64 = ui_plots.plot_clustermap_panels_as_base64(data, store.compute_panelcoldict(us))
    return html.Img(src="data:image/png;base64,%s" % img_b64, className="img-responsive")


@cache.memoize()
def render_image_cluster_clustering(us: reassignSettings):
    data = store.load_all_data(us)
    img_b64 = ui_plots.plot_clustermap_clustering_as_base64(
        data, *store.compute_acluster(us), store.compute_clustercoldict(us)
    )
    return html.Img(src="data:image/png;base64,%s" % img_b64, className="img-responsive")


# @cache.memoize()
# def render_image_batches_clustering(us: reassignSettings):
#    data = store.load_all_data(us)
#    img_b64 = ui_plots.plot_clustermap_batches_as_base64(
#        data, store.compute_acluster(us), store.compute_clustercoldict(us)
#    )
#    return html.Img(src="data:image/png;base64,%s" % img_b64, className="img-responsive")


@cache.memoize()
def render_images_batch_separation(us: reassignSettings):
    return [
        dcc.Graph(figure=px.scatter(p, x="X", y="Y", color="batch"), id="batch-separation-%i" % i)
        for (i, p) in enumerate(store.compute_batches(us))
    ]


def render_main_content():
    """Render page main content"""

    # thresh: int = 50
    # pca_components: int = 20
    # batch_num: float = 0.985
    # pca_seed: int = 100

    controls = dbc.Card(
        children=[
            html.P("Clustering Parameters", className="font-weight-bold"),
            dbc.FormGroup(
                children=[
                    dbc.Label("Threshold"),
                    dbc.Input(
                        id="input-thresh",
                        type="number",
                        value=settings.reassign_SETTINGS.threshold,
                        debounce=True,
                    ),
                ],
            ),
            dbc.FormGroup(
                children=[
                    dbc.Label("PCA Components"),
                    dbc.Input(
                        id="input-pca-components",
                        type="number",
                        value=settings.reassign_SETTINGS.pca_components,
                        debounce=True,
                    ),
                ],
            ),
            dbc.FormGroup(
                children=[
                    dbc.Label("number of batches (clusters)"),
                    dbc.Input(
                        id="input-batch-factor",
                        type="text",
                        value=settings.reassign_SETTINGS.batch_num,
                        debounce=True,
                    ),
                ],
            ),
            html.Div(id="buffer-cluster-params", style={"display": "none"}),
            dbc.Button("Save", id="input-save-batches-button"),
            html.P(id="buffer-save-batches-text"),
        ],
        body=True,
    )

    def make_card(id_, type_="graph"):
        if type_ == "graph":
            return dbc.Card(
                dbc.CardBody(
                    dbc.Spinner([dcc.Graph(id=id_)], color="primary"),
                    className="mt-3 text-center",
                ),
                className="border-top-0 rounded-0",
            )
        else:
            return dbc.Card(
                dbc.CardBody(
                    dbc.Spinner(color="primary", id=id_),
                    className="mt-3 text-center",
                ),
                className="border-top-0 rounded-0",
            )

    tabs_content = [
        dbc.Tab(make_card("graph-pca"), label="PCA"),
        dbc.Tab(make_card("graph-tsne"), label="tSNE"),
        dbc.Tab(make_card("graph-agg-clust"), label="Agg. Clust."),
        dbc.Tab(
            make_card("container-image-cluster-panels", type_="html"), label="Original assignment"
        ),
        dbc.Tab(
            make_card("container-image-cluster-clustering", type_="html"), label="New assignment"
        ),
        # dbc.Tab(
        #    make_card("container-image-batch-clustering", type_="html"), label="Batches clustering"
        # ),
        dbc.Tab(make_card("container-image-batch-separation", type_="html"), label="Batches"),
    ]

    return html.Div(
        children=[
            dbc.Row(
                children=[
                    dbc.Col(
                        children=[controls],
                        id="controls",
                        md=3,
                    ),
                    dbc.Col(children=[dbc.Tabs(tabs_content)], id="plots", md=9),
                ]
            )
        ],
        className="container pt-3",
    )


def render_footer():
    """Render page footer"""
    return html.Footer(
        html.Div(
            children=[
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.Span(
                                    "clear-CNV v%s by BIH CUBI" % __version__,
                                    className="text-muted",
                                )
                            ],
                            className="col-5",
                        ),
                        html.Div(
                            children=[
                                html.A(
                                    children=[
                                        html.I(className="fas fa-book mr-1"),
                                        "Manual & Tutorial",
                                    ],
                                    href="https://clear-cnv.readthedocs.io",
                                    className="text-muted mr-3",
                                ),
                                html.A(
                                    children=[
                                        html.I(className="fas fa-globe-europe mr-1"),
                                        "CUBI Homepage",
                                    ],
                                    href="https://www.cubi.bihealth.org",
                                    className="text-muted mr-3",
                                ),
                                html.A(
                                    children=[
                                        html.I(className="fab fa-github mr-1"),
                                        "GitHub Project",
                                    ],
                                    href="https://github.com/bihealth/clear-cnv",
                                    className="text-muted",
                                ),
                            ],
                            className="col-7 text-right",
                        ),
                    ],
                    className="row",
                )
            ],
            className="container",
        ),
        className="footer",
    )


def build_layout():
    """Build the overall Dash app layout"""
    return html.Div(
        children=[
            # Represents the URL bar, doesn't render anything.
            dcc.Location(id="url", refresh=False),
            # Navbar, content, footer.
            render_navbar(),
            render_main_content(),
            render_footer(),
        ],
        id="_dash-app-content",  # make pytest-dash happy
    )
