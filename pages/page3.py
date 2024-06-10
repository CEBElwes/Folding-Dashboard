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
drugsensitivity_shorthouse = pd.read_csv("./Data/Pathway_direction_pvalue_shorthouse.csv", index_col = "Pathway")
drugsensitivity_cherkaoui = pd.read_csv("./Data/Pathway_direction_pvalue_cherkaoui.csv", index_col = "Pathway")

### ----------------------
# Layout
layout = html.Div(children=[
    html.Br() ,
    html.H1(children='Relationship between Metabolic Pathways and Drug Sensitivity'),

    html.Div(children='''
        Here you can explore the relationships between metabolic pathways and drug sensitivity.
        A positive value indicates that activity of the pathway is associated with an increased resistance to a drug.
        Use dropdown menus to explore specific pathways and drugs.
    '''),
    html.Br(),
    html.Div(children='''
        This data has been normalised in two differing ways - please see the relevant publications for details, but toggle between them below - default is Shorthouse et al.
        ''',style={'textAlign': 'center'}),

    dcc.RadioItems(
    value = 'shorthouse',
    id = 'dataset_type'),
    html.Div([

        # Graph container
        html.Div([
    ## Dropdown for pathways
            dcc.Dropdown(id="pathway2"
                 , value="Citric Acid Cycle"
                 , searchable=True
                 , placeholder="Choose pathway..."
                 , clearable=True),
                 ## Pathway ranking graph
            dcc.Graph(id="drug_sensitivity_by_pathway"),
            ],style={'width': '49%', 'display': 'inline-block', 'padding': 10}),

        html.Div([
            ## Dropdown for TFs
            dcc.Dropdown(id="drug"
                 , value="Cisplatin"
                 , searchable=True
                 , placeholder="Choose Drug..."
                 , clearable=True),
    ## TF ranking graph
    dcc.Graph(id="pathway_ranking_by_drug"),
    ],style={'width': '52%', 'display': 'inline-block', 'padding': 10})
    ], style={'display': 'flex'}),
])

#-------------------------------------------
## Callbacks

## Callback for Dropdown for pathway
@app.callback(
    Output(component_id = "pathway", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)
def set_dropdown_options_page3_1(dataset):
    if dataset == "shorthouse":
        pathway_dropdown = [{'label': i, 'value': i} for i in drugsensitivity_shorthouse.index.tolist()[1:]]
    return pathway_dropdown

@app.callback(
    Output(component_id = "drug", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)
def set_dropdown_options_page3_2(dataset):
    if dataset == "shorthouse":
        drug_dropdown = [{'label': i, 'value': i} for i in drugsensitivity_shorthouse.columns]
    return drug_dropdown


## Callback for mutation rankings per metabolite
@app.callback(
    Output(component_id = "drug_sensitivity_by_pathway", component_property = "figure"),
    [Input(component_id = "pathway", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)

def drug_sensitivity_by_pathway_plot(pathway_name, dataset):
    if dataset == "shorthouse":
        pathway_dataframe = drugsensitivity_shorthouse
    pathway = pd.DataFrame(pathway_dataframe.loc[pathway_name].sort_values())
    pathway.columns = ["log10(Pvalue) * correlation direction"]
    pathway["Drug"] = pathway.index
    pathway = pathway.reset_index()
    pathway["Drug Rank"] = pathway.index+1

    pathway["Association"] = np.where((pathway["log10(Pvalue) * correlation direction"] >= 0), "Resistance",
                                                 "Sensitivity")

    scatterplot = px.scatter(pathway, x = "Drug Rank", y = "log10(Pvalue) * correlation direction", hover_name = "Drug"
                             ,title = "Association of " + pathway_name + " activity with drug resistance/sensitivity", template = "simple_white",
                              color = "Association", color_discrete_sequence = ["blue", "red"])
    return scatterplot


@app.callback(
    Output(component_id = "pathway_ranking_by_drug", component_property = "figure"),
    [Input(component_id = "drug", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)

def pathway_ranking_by_drug_plot(drug_name, dataset):
    if dataset == "shorthouse":
        pathway_dataframe = drugsensitivity_shorthouse.T
    pathway = pd.DataFrame(pathway_dataframe.loc[drug_name].sort_values())
    pathway.columns = ["log10(Pvalue) * correlation direction"]
    pathway["SMPDB Pathway"] = pathway.index
    pathway = pathway.reset_index()
    pathway["Pathway Rank"] = pathway.index+1

    pathway["Association"] = np.where((pathway["log10(Pvalue) * correlation direction"] >= 0), "Resistance",
                                                 "Sensitivity")

    scatterplot = px.scatter(pathway, x = "Pathway Rank", y = "log10(Pvalue) * correlation direction", hover_name = "SMPDB Pathway"
                             ,title = "Association of resistance to " + drug_name + " with SMPDB pathway activity levels", template = "simple_white",
                             color = "Association", color_discrete_sequence = ["blue", "red"])
    return scatterplot

if __name__ == '__main__':
    app.run_server(debug=True, port = 8052)