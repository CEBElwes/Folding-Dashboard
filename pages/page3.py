from dash import Dash, html, dcc, Output, Input
import dash, dash_table
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns

app = Dash(__name__)

## Read in data
ddg_info = pd.read_csv("ddg_info.csv")


### ----------------------
# Layout
layout = html.Div(children=[
    html.Br(),
    html.H1(children='Folding Energies for Genes and Variants'),

    html.Div(children='''
        Here you can explore the (mis)folding energies.
    '''),


    html.Div([

        # Graph container
        html.Div([
            ## Dropdown for gene
            dcc.Dropdown(id="gene",
                         value="-1",
                         searchable=True,
                         placeholder="Choose gene...",
                         clearable=True),
            ## Gene graph
            dcc.Graph(id="ddg_by_gene"),
        ], style={'width': '49%', 'display': 'inline-block', 'padding': 10}),

        html.Div([
            ## Dropdown for variant
            dcc.Dropdown(id="variant",
                         value="XXX",
                         searchable=True,
                         placeholder="Choose variant...",
                         clearable=True),
            ## Variant graph
            dcc.Graph(id="ddg_by_variant"),
        ], style={'width': '52%', 'display': 'inline-block', 'padding': 10})
    ], style={'display': 'flex'}),
])

#-------------------------------------------
## Callbacks

## Callback for Dropdown for gene
@app.callback(
    Output(component_id="gene", component_property="options"),
    Input(component_id="gene", component_property="value")
)
def set_dropdown_options_page3_1(_):
    gene_dropdown = [{'label': i, 'value': i} for i in ddg_info.index.tolist()]
    return gene_dropdown

@app.callback(
    Output(component_id="variant", component_property="options"),
    Input(component_id="variant", component_property="value")
)
def set_dropdown_options_page3_2(_):
    variant_dropdown = [{'label': i, 'value': i} for i in ddg_info.columns]
    return variant_dropdown


## Callback for mutation rankings per metabolite
@app.callback(
    Output(component_id="drug_sensitivity_by_pathway", component_property="figure"),
    Input(component_id="pathway", component_property="value")
)
def drug_sensitivity_by_pathway_plot(pathway_name):
    pathway_dataframe = drugsensitivity_shorthouse
    pathway = pd.DataFrame(pathway_dataframe.loc[pathway_name].sort_values())
    pathway.columns = ["log10(Pvalue) * correlation direction"]
    pathway["Drug"] = pathway.index
    pathway = pathway.reset_index(drop=True)
    pathway["Drug Rank"] = pathway.index + 1

    pathway["Association"] = np.where((pathway["log10(Pvalue) * correlation direction"] >= 0), "Resistance",
                                      "Sensitivity")

    scatterplot = px.scatter(pathway, x="Drug Rank", y="log10(Pvalue) * correlation direction", hover_name="Drug",
                             title="Association of " + pathway_name + " activity with drug resistance/sensitivity",
                             template="simple_white", color="Association", color_discrete_sequence=["blue", "red"])
    return scatterplot


@app.callback(
    Output(component_id="pathway_ranking_by_drug", component_property="figure"),
    Input(component_id="drug", component_property="value")
)
def pathway_ranking_by_drug_plot(drug_name):
    pathway_dataframe = drugsensitivity_shorthouse.T
    pathway = pd.DataFrame(pathway_dataframe.loc[drug_name].sort_values())
    pathway.columns = ["log10(Pvalue) * correlation direction"]
    pathway["SMPDB Pathway"] = pathway.index
    pathway = pathway.reset_index(drop=True)
    pathway["Pathway Rank"] = pathway.index + 1

    pathway["Association"] = np.where((pathway["log10(Pvalue) * correlation direction"] >= 0), "Resistance",
                                      "Sensitivity")

    scatterplot = px.scatter(pathway, x="Pathway Rank", y="log10(Pvalue) * correlation direction",
                             hover_name="SMPDB Pathway",
                             title="Association of resistance to " + drug_name + " with SMPDB pathway activity levels",
                             template="simple_white", color="Association", color_discrete_sequence=["blue", "red"])
    return scatterplot

if __name__ == '__main__':
    app.run_server(debug=True, port=8052)
