# -*- coding: utf-8 -*-
import pandas as pd
import io, sys, os, datetime
import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import base64
import pickle
import plotly.express as px
import matplotlib
matplotlib.use('Agg')
from dash_extensions import Download
from dash_extensions.snippets import send_file
import json
import time
import subprocess
from pathlib import Path

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BaseFolder = "./"
FULL_PATH = os.path.abspath(BaseFolder)+"/"

sys.path.append(FULL_PATH)
from netcom.netcom import pathwayEnrichment, EdgeR_to_seeds, simulation

#Delete all results older then 5 days
deleteOlFiles_command = "find "+FULL_PATH+"Results/* -type d -ctime +5 -exec rm -rf {} \;"
os.system(deleteOlFiles_command)

try:
     with open(BaseFolder+"data/DB/DB.pickle", 'rb') as handle:
         DB = pickle.load(handle)
except:
    DB = pd.read_pickle(BaseFolder+"data/DB/DB.pickle")      
  
df_el_ = DB['full_enzymes_labels_jun.txt']
df_ecMapping_ = DB['ec_reac_mapping_jun.txt']
df_reactions_ = DB['reactions_3_balanced.txt']
df_ec_to_compoundIndex_ = DB['compound_labels_jun.txt']


def read_edgeR(df):
    try:
        df_edgeR_grouped = df.groupby("association")['X'].apply(list).to_frame()
    except:
        df_edgeR_grouped = df.groupby("association")['enzyme'].apply(list).to_frame()
    return df_edgeR_grouped


def createParametersDict(folder):    
    start = time.time()
    #write default parameters pickle file. this will be changed by the user later.
    parametersDict = {}
    parametersDict["drop_fragment_with_size"] = 1
    parametersDict["filter_hubness"] = 50
    parametersDict["soft_color_A"] = "green"
    parametersDict["dark_color_A"] = "lime"
    parametersDict["corrected_p-val"] = 0.05
    parametersDict["enrichment_results_slice"] = [0, 100]
    parametersDict["figure_name"] = "Figure"
    parametersDict["network_layout_iter"] = 75
    parametersDict["treatment_col"] = ""
    parametersDict["comparison_col"] = ""
    parametersDict["Not_associated_col"] = ""
    parametersDict["Min_entities_Enrichment"] = 3
    parametersDict["Max_entities_Enrichment"] = 25
    parametersDict["Enriched_pathways"] = []    
    parametersDict["Final_folder"] = "" 
    
    folder = str(folder).strip("\"").strip("\'")
    
    f = open(folder+"parametersDict.json", "w", encoding="utf8")
    json.dump([parametersDict], f)
    f.close()   
    try:
        os.system("rm "+folder+"main_proccess_results_html.pkl")
    except:
        print()
        
    end = time.time()
    print("createParametersDict time (sec):")
    print(end - start)

    

def loadParametersDict(folder):
    start = time.time()
    folder = str(folder).strip("\"").strip("\'")
    f = open(folder+"parametersDict.json", "r")
    output = json.load(f)
    f.close()
    end = time.time()
    print("loadParametersDict time (sec):")
    print(end - start)    
    return output[0]
  

def update_parameters(val, col, folder):
    start = time.time()   
    try:
        folder = str(folder).strip("\"").strip("\'")
        
        f = open(folder+"parametersDict.json", "r")
        parametersDict = json.load(f)[0]
        f.close()
    except:
        folder = str(folder).strip("\"").strip("\'")
        createParametersDict(folder)
        f = open(folder+"parametersDict.json", "r")
        parametersDict = json.load(f)[0]
        f.close()        
    
    parametersDict[col] = val
    f = open(folder+"parametersDict.json", "w", encoding="utf8")
    json.dump([parametersDict], f)
    f.close()
    
    end = time.time()
    print("update_parameters time (sec):")
    print(end - start)        
    
    
    
def presentDatasetStatistics(folder):
    start = time.time()
    print("loading edger")
    df = pd.read_csv(folder+"raw_input_edger.csv")
    
    colorsDict={}
    colorsDict["treatment_col"] = "blue"
    colorsDict["comparison_col"] = "red"
    colorsDict["Not_associated"] = "Not_associated"
    print("prep colors")
    df["Treatment color"] = df["association"].replace(colorsDict)
    
    #VOLCANO PLOT    
    try:
        volcano = px.scatter(df, x="logFC", y="PValue",color="Treatment color",
                             hover_name="enzyme", log_y=True)
        try:
            labels = df[["association"]].value_counts().index
            values = df[["association"]].value_counts().values
        except:
            labels = df["association"].value_counts().index
            values = df["association"].value_counts().values
            
        pieChart = go.Figure(data=[go.Pie(labels=labels, values=values)])

        pvalHist = px.histogram(df, x="PValue")


        descriptionGraphs = html.Div(dbc.Row([
                                    dbc.Col(dcc.Graph(
                                        id='volcano-scatter',
                                        #style={'display': 'inline-block'},
                                        figure=volcano
                                        )),
                                    
                                    dbc.Col(dcc.Graph(
                                        id='pie-chart',
                                        #style={'display': 'inline-block'},
                                        figure=pieChart
                                        )),
                                    
                                    dbc.Col(dcc.Graph(
                                        id='pval-hist',
                                        #style={'display': 'inline-block'},
                                        figure=pvalHist
                                        )),                       
                                ])
                        )
                
        end = time.time()
        print("presentDatasetStatistics time (sec):")
        print(end - start)
        #calculate enrichment for keep_pathway.txt file creation
        parametersDict = loadParametersDict(folder)
        #print(parametersDict)
        folder = str(folder).strip("\"").strip("\'")   
        FinalFolder_ = folder   
        Seeds_A_input, T1_seeds_tag, ECs_A_input, Seeds_B_input, T2_seeds_tag, ECs_B_input, Seeds_All_input, ECs_All_input=EdgeR_to_seeds(edgeR_row_location=FinalFolder_+"raw_input_edger.csv",
                                                                                                              col_treatment_1=parametersDict["treatment_col"],
                                                                                                              col_treatment_2=parametersDict["comparison_col"],
                                                                                                              outputFolder=FinalFolder_,
                                                                                                                  input_sep=",")
     

        with open(FinalFolder_+"keep_pathways.txt", 'w') as f:
            f.write("\n")    
                
        pathways_enzymes_A = pathwayEnrichment(BaseFolder=BaseFolder,
                                        FinalFolder=folder,
                                        DE_ecs_list_=ECs_A_input,
                                        All_B=ECs_B_input,
                                        All_ecs_list_=ECs_A_input+ECs_B_input,
                                        input_type="enzymes", # "compound" or "enzyme"
                                        outputfilename=parametersDict["treatment_col"],
                                        minEntitiesInPathway=parametersDict["Min_entities_Enrichment"], 
                                        maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                        drop_pathways=False,
                                        keep_pathways=FinalFolder_+"keep_pathways.txt"
                                        )
        
        
        pathways_enzymes_B = pathwayEnrichment(BaseFolder=BaseFolder,
                                        FinalFolder=folder,
                                        DE_ecs_list_=ECs_B_input,
                                        All_B=ECs_A_input,
                                        All_ecs_list_=ECs_A_input+ECs_B_input,
                                        input_type="enzymes", # "compound" or "enzyme"
                                        outputfilename=parametersDict["comparison_col"],
                                        minEntitiesInPathway=parametersDict["Min_entities_Enrichment"],
                                        maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                        drop_pathways=False,
                                        keep_pathways=FinalFolder_+"keep_pathways.txt"                                    
                                        )
        
        return descriptionGraphs   
    except Exception as e:
        print(e)    
        end = time.time()
        print("presentDatasetStatistics time (sec):")
        print(end - start)        
    

def CreateBarPlot(folder):
    start = time.time()    
    parametersDict = loadParametersDict(folder)
    folder = str(folder).strip("\"").strip("\'")
    FinalFolder_ = folder
    df = pd.read_csv(FinalFolder_+"raw_input_edger.csv")
    df=read_edgeR(df)   
    Seeds_A_input, T1_seeds_tag, ECs_A_input, Seeds_B_input, T2_seeds_tag, ECs_B_input, Seeds_All_input, ECs_All_input=EdgeR_to_seeds(edgeR_row_location=FinalFolder_+"raw_input_edger.csv",
                                                                                                          col_treatment_1=parametersDict["treatment_col"],
                                                                                                          col_treatment_2=parametersDict["comparison_col"],
                                                                                                          outputFolder=FinalFolder_,
                                                                                                          input_sep=",")

    with open(FinalFolder_+"keep_pathways.txt", 'w') as f:
        f.write("\n")    
            
    pathways_enzymes_A = pathwayEnrichment(BaseFolder=BaseFolder,
                                    FinalFolder=folder,
                                    DE_ecs_list_=ECs_A_input,
                                    All_B=ECs_B_input,
                                    All_ecs_list_=ECs_A_input+ECs_B_input,
                                    input_type="enzymes", # "compound" or "enzyme"
                                    outputfilename=parametersDict["treatment_col"],
                                    minEntitiesInPathway=parametersDict["Min_entities_Enrichment"], 
                                    maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                    drop_pathways=False,
                                    keep_pathways=FinalFolder_+"keep_pathways.txt"
                                    )
        
    pathways_enzymes_B = pathwayEnrichment(BaseFolder=BaseFolder,
                                    FinalFolder=folder,
                                    DE_ecs_list_=ECs_B_input,
                                    All_B=ECs_A_input,
                                    All_ecs_list_=ECs_A_input+ECs_B_input,
                                    input_type="enzymes", # "compound" or "enzyme"
                                    outputfilename=parametersDict["comparison_col"],
                                    minEntitiesInPathway=parametersDict["Min_entities_Enrichment"],
                                    maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                    drop_pathways=False,
                                    keep_pathways=FinalFolder_+"keep_pathways.txt"                                    
                                    )
    
    df_enzymes_A = pd.read_csv(FinalFolder_+parametersDict["treatment_col"]+"_pathway.csv")
    df_enzymes_B = pd.read_csv(FinalFolder_+parametersDict["comparison_col"]+"_pathway.csv")
    
    pathways_seed_rich = pathwayEnrichment(BaseFolder=BaseFolder,
                                    FinalFolder=folder,
                                    DE_ecs_list_=T1_seeds_tag.copy(),
                                    All_B=T2_seeds_tag.copy(),
                                    All_ecs_list_=T1_seeds_tag.copy()+T2_seeds_tag.copy(),
                                    input_type="metabolites", # "metabolites" or "enzyme"
                                    outputfilename=parametersDict["treatment_col"],
                                    minEntitiesInPathway=parametersDict["Min_entities_Enrichment"],
                                    maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                    drop_pathways=False,
                                    keep_pathways=FinalFolder_+"keep_pathways.txt"
                                    )

    pathways_seed_poor = pathwayEnrichment(BaseFolder=BaseFolder,
                                    FinalFolder=FinalFolder_,
                                    DE_ecs_list_=T2_seeds_tag.copy(),
                                    All_B=T1_seeds_tag.copy(),
                                    All_ecs_list_=T1_seeds_tag.copy()+T2_seeds_tag.copy(),
                                    input_type="metabolites", # "metabolites" or "enzyme"
                                    outputfilename=parametersDict["comparison_col"],
                                    minEntitiesInPathway=parametersDict["Min_entities_Enrichment"],
                                    maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                    drop_pathways=False,
                                    keep_pathways=FinalFolder_+"keep_pathways.txt"
                                    )
    
    df_compounds_A = pd.read_csv(FinalFolder_+parametersDict["treatment_col"]+"_pathway.csv")
    df_compounds_B = pd.read_csv(FinalFolder_+parametersDict["comparison_col"]+"_pathway.csv")

    listOfAllPathways = df_enzymes_A["Pathway"].values.tolist()
    listOfAllPathways += df_enzymes_B["Pathway"].values.tolist()
    listOfAllPathways += df_compounds_A["Pathway"].values.tolist()
    listOfAllPathways += df_compounds_B["Pathway"].values.tolist()
    
    update_parameters(list(set(listOfAllPathways)), "Enriched_pathways", folder)
 
    BarPltEnzymes = dcc.Graph(id='g1',
                              figure={'data': [
                                            {'x': df_enzymes_A["Pathway"].values.tolist(), 'y': df_enzymes_A["Count"].values.tolist(), 'type': 'bar', 'name':parametersDict["treatment_col"]},
                                            {'x': df_enzymes_B["Pathway"].values.tolist(), 'y': df_enzymes_B["Count"].values.tolist(), 'type': 'bar', 'name':parametersDict["comparison_col"]}
                                            ],
                                    'layout':{
                                        'xaxis': {'autorange': True, 'title': 'X Axis', 'automargin': True}
                                        }},
                                    style={'height': '800px', 'width': '1000px'}
        )


    BarPltCompounds = dcc.Graph(id='g1',
                              figure={'data': [
                                            {'x': df_compounds_A["Pathway"].values.tolist(), 'y': df_compounds_A["Count"].values.tolist(), 'type': 'bar', 'name':parametersDict["treatment_col"]},
                                            {'x': df_compounds_B["Pathway"].values.tolist(), 'y': df_compounds_B["Count"].values.tolist(), 'type': 'bar', 'name':parametersDict["comparison_col"]}
                                            ],
                                    'layout':{
                                        'xaxis': {'autorange': True, 'title': 'X Axis', 'automargin': True}
                                        }},
                                    style={'height': '800px', 'width': '1000px'}
        )
    
   
    BarPlts = html.Div([
        html.Div([
            html.Div([
                html.H3('Enzymes - pathways histogram'),
                BarPltEnzymes
            ], className="six columns"),
            
            html.Div([
                html.H3('Compounds - pathways histogram'),
                BarPltCompounds
            ], className="six columns"),
        ],
            className="row")
    ])


    end = time.time()
    print("CreateBarPlot time (sec):")
    print(end - start)        
                     
    return BarPlts

def slider_enrch_min():   
    return html.Div([
                    dbc.FormGroup(
                    [
                        dbc.Label("Entities number in a pathway [range, used for enrichment analysis]",
                                  html_for="slider-enrch-min",
                                  style={'text-align': 'center', 'font-size': '100%', 'text-transform': 'uppercase'}),
                        dcc.RangeSlider(id="slider-enrch-min",
                               min=1,
                               max=100,
                               step=1,
                               marks={i: str(i) for i in range(100)},
                               value=[5, 25])
                    ]
                ),
                    html.Div(id='slider-output-enrch-min-container', style={'margin-bottom': 20})
                    ],
                    style={'padding': '10px 10px 10px 10px'}
                )

def slider_node_hubness():
    return html.Div([
                    dbc.FormGroup(
                    [
                        dbc.Label("Limit node hubness [default = 50]",
                                  html_for="slider-hubness",
                                  style={'text-align': 'center', 'font-size': '100%', 'text-transform': 'uppercase'}
                                  ),
                        dcc.Slider(id="slider-hubness",
                                   min=1,
                                   max=100,
                                   step=1,
                                   marks={i: str(i) for i in range(100)},
                                   value=50),
                    ]
                ),
                    html.Div(id='slider-output-hubness-container', style={'margin-bottom': 20})
                    ]
                )

def slider_network_iter():
    return html.Div([
                    dbc.FormGroup(
                    [
                        dbc.Label("Set network layout iterations [default = 75]",
                                  html_for="slider-iter",
                                  style={'text-align': 'center', 'font-size': '100%', 'text-transform': 'uppercase'}
                                  ),
                        dcc.Slider(id="slider-iter",
                                   min=1,
                                   max=200,
                                   step=5,
                                   marks={i: str(i) for i in range(0, 200, 5)},
                                   value=75),
                    ]
                ),
                    html.Div(id='slider-output-iter-container', style={'margin-bottom': 20})
                    ]
                )

def select_colors_seeds():
    return html.Div([
                    dbc.FormGroup(
                    [
                        dbc.Label("Environmental resource node color",
                                  html_for="seeds-color-dropdown", width=4),
                        dbc.Col(
                            dcc.Dropdown(
                                id="seeds-color-dropdown",
                                options=[
                                    {'label': 'Green', 'value': 'green'},
                                    {'label': 'Orange', 'value': 'orange'},
                                    {'label': 'Blue', 'value': 'blue'},
                                    {'label': 'Red', 'value': 'red'},
                                    {'label': 'Goldenrod', 'value': 'goldenrod'},
                                    {'label': 'Magenta', 'value': 'magenta'},
                                    {'label': 'Medium Purple', 'value': 'mediumpurple'},
                                    {'label': 'Chocolate', 'value': 'chocolate'},
                                    {'label': 'Khaki', 'value': 'khaki'},                    
                                ],
                                value='green',
                            ),
                            width=10,
                        ),
                    ],
                    row=True,
                ),
                    html.Div(id='output-seedcolor-container', style={'padding': 10})   
                ],
                    style={"width": "30%",
                                   #'display':'inline-block',
                                   'text-align':'center',
                                   #'padding-left':'35%', 
                                   #'padding-right':'35%'
                                   }
                    )

def select_colors_unique():
    return html.Div([
                    dbc.FormGroup(
                    [
                        dbc.Label("Unique node color", html_for="unique-color-dropdown", width=4),
                        dbc.Col(
                            dcc.Dropdown(
                                id="unique-color-dropdown",
                                options=[
                                    {'label': 'Lime', 'value': 'lime'},
                                    {'label': 'Salmon', 'value': 'salmon'},
                                    {'label': 'Cyan', 'value': 'cyan'},
                                    {'label': 'Pink', 'value': 'pink'},
                                    {'label': 'Yellow', 'value': 'yellow'},
                                    {'label': 'Gold', 'value': 'gold'},
                                    {'label': 'Light Purple', 'value': 'lightpurple'},
                                    {'label': 'Gray', 'value': 'gray'},
                                    {'label': 'Light Blue', 'value': 'lightblue'},                    
                                ],
                                value='lime',
                            ),
                            width=10, 
                        ),
                    ],
                    row=True,
                ),
                    html.Div(id='output-uniquecolor-container', style={'padding': 10})   
                ],
                    style={"width": "30%",
                                   #'display':'inline-block',
                                   'text-align':'center',
                                   #'padding-left':'35%', 
                                   #'padding-right':'35%'
                                   }
                    )

def execution_button():
    return html.Div([
                    html.Div(dcc.Input(id='input-on-submit', type='hidden')),
                    html.Button('Submit', id='submit-val', n_clicks=0),
                    html.Div(id='container-button-basic'),
                    #html.Div(id='container-button-basic', children=''),
                    dcc.Interval(id='final-results-listener', interval=5000, n_intervals=0, max_intervals=180),#fires 180 times * 5 sec = 15minutes
                    html.Div(id='container-final-results')
                ])

def pathways_dropout(folder):
    #Read enrichment csv files
    folder = str(folder).strip("\"").strip("\'")

    parametersDict = loadParametersDict(folder)
    df_compounds_A = pd.read_csv(folder+parametersDict["treatment_col"]+"_pathway.csv")
    df_compounds_B = pd.read_csv(folder+parametersDict["comparison_col"]+"_pathway.csv")

    df_compounds_A = df_compounds_A.set_index("Pathway")
    df_compounds_B = df_compounds_B.set_index("Pathway")
    
    #iterate over the options and add the counts to label
    df_pathwayCounts = pd.DataFrame(index=parametersDict["Enriched_pathways"],
                                    columns=[parametersDict["treatment_col"]+"_count",
                                             parametersDict["comparison_col"]+"_count"])    
    #populate treatment col index 
    for i in df_compounds_A.index:
        #print("A")
        #print(i)
        try:
            df_pathwayCounts.at[i, parametersDict["treatment_col"]+"_count"] = df_compounds_A.loc[i, "Count"]
        except Exception as e:
            print(e)
        
            
    for i in df_compounds_B.index:
        #print("B")        
        #print(i)
        try:
            df_pathwayCounts.at[i, parametersDict["comparison_col"]+"_count"] = df_compounds_B.loc[i, "Count"]
        except Exception as e:
            print(e)
    
    
    df_pathwayCounts = df_pathwayCounts.reset_index()
    df_pathwayCounts = df_pathwayCounts.fillna(0)
    #print(df_pathwayCounts)
    opts = [{'label' : str(i[0])+" "+str(i[1])+" "+str(i[2]) , 'value' : i[0]} for i in df_pathwayCounts.values.tolist()]    
    #opts = [{'label' : i, 'value' : i} for i in parametersDict["Enriched_pathways"]]
    #print("options:")
    #print(opts)
    #print("values:")
    #print(parametersDict["Enriched_pathways"])
    return dbc.Container(
                children=[
                    html.Details([
                    html.Summary('Select pathways to dropout - (Pathway / Number-of-compounds-'+parametersDict["treatment_col"]+' / Number-of-compounds-'+parametersDict["comparison_col"]+')'),
                    html.Br(),
                    dbc.Col([
                        dcc.Checklist(id='my-checklist1',
                            options=opts,
                            value=parametersDict["Enriched_pathways"],
                            labelStyle = {'display': 'block'},
                            )  
                        ])
                    ])
                ])

#Create dropdown
def Explore_enriched_pathways_treatment(folder):
    #REAd all subgraphs object
    parametersDict = loadParametersDict(folder)
    
    with open(BaseFolder+folder+parametersDict["treatment_col"]+"_allSubGraphs.pickle", 'rb') as handle:
        allSubGraphs = pickle.load(handle) 
    
    
    
    opts = [{'label' : i[0] , 'value' : i[0]} for i in allSubGraphs]    
    #opts = [{'label' : i, 'value' : i} for i in parametersDict["Enriched_pathways"]]

    return html.Div(
                children=[
                    html.H3('Select enriched pathway to explore its subgraph- '+parametersDict["treatment_col"]),
                    dbc.Col([
                        dcc.Dropdown(id='explore-enriched-pathways-treatment-dropdown',
                            options=opts,
                            #value=parametersDict["Enriched_pathways"],
                            )  
                        ])
                    ])

def Explore_enriched_pathways_control(folder):
    #REAd all subgraphs object
    parametersDict = loadParametersDict(folder)
    
    with open(BaseFolder+folder+parametersDict["comparison_col"]+"_allSubGraphs.pickle", 'rb') as handle:
        allSubGraphs = pickle.load(handle) 
    
    opts = [{'label' : i[0] , 'value' : i[0]} for i in allSubGraphs]    

    return html.Div(
                children=[
                    html.H3('Select enriched pathway to explore its subgraph- '+parametersDict["comparison_col"]),
                    dbc.Col([
                        dcc.Dropdown(id='explore-enriched-pathways-control-dropdown',
                            options=opts,
                            #value=parametersDict["Enriched_pathways"],
                            )  
                        ])
                    ])

image_dep_filename = './assets/dep_sign.png' 
image_volcani_filename = './assets/Minhal_logos_english.png' 
image_banner_options_filename = './assets/banner_option2.jpg' 
   
fixed_header = html.Div(
                    style={'background-image': "url(/assets/bg.jpg)",
                           "padding": "60px"
                           },                        
                    children=[
                            dbc.Row(
                                style={'display':'flex'},
                                children=[
                                    html.Img(style={'height':'197px',
                                                    'width':'170px'},
                                             src=image_dep_filename),
                                    html.Div(
                                             children=[html.Div(style={"border":"2px black solid",
                                                                       'margin': '10px 20px 10px 20px',
                                                                       'backgroundColor':'white'},
                                                                children=html.H1(children=["Web Tools"],
                                                                                 style={'marginBottom': 25,
                                                                                        'marginTop': 25,
                                                                                        'marginLeft': 300,
                                                                                        'marginRight': 300,
                                                                                        'backgroundColor':'white'})           
                                                                ),
                                                       html.H2(children="Freilich Lab",
                                                               style={'text-align':'center'})
                                                       ]),
                                    html.Img(style={'height':'187px',
                                                    'width':'189px'},
                                             src=image_volcani_filename),
                                    ],
                                justify="center",
                                align="center")
                ])

header_links = html.Div(children=[
                        dbc.Row([
                            dbc.Button("Home",
                                href='https://www.freilich-lab.com/',
                                style={"margin-right": "15px"}),
                            dbc.Button("Members",
                                href='https://www.freilich-lab.com/members/',
                                style={"margin-right": "15px"}),
                            dbc.Button("Projects",
                                href='https://www.freilich-lab.com/projects',
                                style={"margin-right": "15px"}),
                            dbc.Button("Gallery",
                                href='https://www.freilich-lab.com/gallery',
                                style={"margin-right": "15px"}),
                            dbc.Button("Contact",
                                href='https://www.freilich-lab.com/contact',
                                style={"margin-right": "15px"})                                         
                            ],
                            justify="center",
                            align="center"),
                        
                        html.Img(src=image_banner_options_filename,)                        
                        
                        ],
    style={'textAlign': 'center'})

main_header = html.H1(style={"margin-left": "60px"},
                      children='NetCom - A Tool For Functional Interpretation Of Metagenomic Data')

network_parameters_header = html.H2(children='Network parameters')
differential_abundance_header = html.H2(children='Differential abundance analysis')

text_f = open(BaseFolder+"data/introduction.txt", "r")
introduction_text = text_f.read()
text_f.close()

introduction_text = html.Div([
                        html.H4(introduction_text, style={'whiteSpace': 'pre-wrap',
                                                          "margin-right": "60px",
                                                          "margin-left": "60px"}),
                        html.Div([
                            html.A('For more information please check our publication Tal et al Microorganisms 2021', href='https://www.mdpi.com/2076-2607/8/6/840', style={'font-size':'20.8px'}),
                            html.Div(children='Data in the example file was taken from metagenomics sequencing of the root environment of wheat (Treatment 1 – \'root\') and the more distant soil not under direct effect of the plant (Treatment 2 – \'soil\').', style={'font-size':'20.8px'}),
                            html.A('Data was taken from Ofek-Lalzar et al, Nature Communication 2014', href= 'https://www.nature.com/articles/ncomms5950', style={'font-size':'20.8px'}),
                            html.H4(children="Please notice, a run takes 1-15 minutes (depends on your data). DO NOT refresh or close the tab until the compressed results file is downloaded.", style={'color': 'green'})
                            ],                        
                            )
            ])

def parameters_form(): 
    return dbc.Form([differential_abundance_header,
                        slider_enrch_min(),
                        network_parameters_header,
                        select_colors_seeds(),
                        select_colors_unique(),
                        slider_node_hubness(),
                        slider_network_iter(),                            
                        execution_button(),
                        html.Div(dcc.Store(id='cache', storage_type='local'))]) 


#Upload file
#############

def ul():    
    return html.Div([                 
                    dcc.Upload(
                        id='upload-data',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Files')
                        ]),
                        style={
                            'width': '100%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '2px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                        },
                        # Allow multiple files to be uploaded
                        multiple=False
                    ),
                    html.Div(id='output-data-upload', children=""),                    
                    ])
                
def parse_contents(contents, filename, date, folder):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
       FinalFolder_ = folder            
       df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep="\t")
       df.to_csv(FinalFolder_+"raw_input_edger.csv")
    except Exception as e:
        print(e)
        return print('There was an error processing this file.')
    return df

  
download_link = html.Div([html.Button("Download results", id="btn"), Download(id="download")])
download_instructions_link = html.Div([html.Button("Download instructions file", id="btn-inst"), Download(id="download-inst")])
download_example_link = html.Div([html.Button("Download example file", id="btn-examp"), Download(id="download-examp")])
download_link_placeholder = html.Div(id="download_btn_placeholder")
download_treatment_enrichment_link = html.Div([html.Button("Download pathway enrichment file", id="btn-enrichment_treatment"), Download(id="download-treatment-enrichment")])
download_comparison_enrichment_link = html.Div([html.Button("Download pathway enrichment file", id="btn-enrichment_comparison"), Download(id="download-comparison-enrichment")])

# =============================================================================

#THE APP
##############

external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, url_base_pathname='/netcom/')
server = app.server
app.title = 'NetCom'

app.layout = html.Div(
            style={"margin-right": "60px",
                   "margin-left": "60px"},
            children=[fixed_header,
                       header_links,
                       main_header,
                       introduction_text,
                       dbc.Row(children=[download_instructions_link, download_example_link]),
                       html.Div(ul()),                    
                       parameters_form(),
                       download_link_placeholder                      
                    ])

@app.callback(Output('cache', 'data'), Input('upload-data', 'contents'))
def store_data(value):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if 'upload-data' in changed_id:
        try:
            os.mkdir(BaseFolder + "Results")
        except:
            print()   
        d = datetime.datetime.today()
        try:
            FinalFolder_ = BaseFolder + "Results/" + str(d.year) + str(d.month) + str(d.day) + str(d.hour) + str(
                d.minute) + str(
                d.second) + "/"
        except:
            FinalFolder_ = BASE_DIR + "/Results/" + str(d.year) + str(d.month) + str(d.day) + str(d.hour) + str(d.minute) + str(
                d.second) + "/"       
        try:
            os.mkdir(FinalFolder_)    
        except:
            print()
        
        print("#### " + FinalFolder_)
        dumped_folder = json.dumps(FinalFolder_)
        return dumped_folder

@app.callback(Output('output-data-upload', 'children'),
              [Input('upload-data', 'contents'), Input('cache', 'data')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])
def update_output(list_of_contents, folder, list_of_names, list_of_dates):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if 'upload-data' in changed_id:    
        time.sleep(2)
        folder = json.loads(folder)
        folder = str(folder).strip("\"").strip("\'")
        print("#### #### " + folder)   
        if list_of_contents is not None:  
            createParametersDict(folder)
            print("folder: "+folder)   
            df_grouped = read_edgeR(parse_contents(list_of_contents, list_of_names, list_of_dates, folder=folder))
            df_grouped.to_csv(folder+"input_edger.csv")
            d = pd.read_csv((folder+"input_edger.csv"))                                
            labels = d["association"].values.tolist()                       
            treatment_options = []
            for i in labels:
                if (str(i) == "Not_associated") or (str(i) == "Not associated"):
                    update_parameters(i, "Not_associated_col", folder)
                else:
                    treatment_options.append(i)
                
            update_parameters(treatment_options[0], "treatment_col", folder)  
            update_parameters(treatment_options[1], "comparison_col", folder)
            
            #print(treatment_options)
            print("callback update_output")
            return presentDatasetStatistics(folder)

@app.callback(
    dash.dependencies.Output('slider-output-enrch-min-container', 'children'),
    [dash.dependencies.Input('slider-enrch-min', 'value'), Input('cache', 'data')])
def update_enrichment_min(value, folder):

    update_parameters(value[0], "Min_entities_Enrichment", folder)
    update_parameters(value[1], "Max_entities_Enrichment", folder)
    print("callback update_enrichment_min")   
    xx = 'Range selected = {}'.format(str(value))
    return [xx, CreateBarPlot(folder), pathways_dropout(folder)]

@app.callback(
    dash.dependencies.Output('slider-output-hubness-container', 'children'),
    [dash.dependencies.Input('slider-hubness', 'value'), Input('cache', 'data')])
def update_hubness(value, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    update_parameters(value, "filter_hubness", folder)
    print("callback update_hubness")     
    return 'Max hubness selected = {}'.format(value)

@app.callback(
    dash.dependencies.Output('slider-output-iter-container', 'children'),
    [dash.dependencies.Input('slider-iter', 'value'), Input('cache', 'data')])
def update_iter(value, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    update_parameters(value, "network_layout_iter", folder)
    print("callback update_iter")   
    return 'Number of iterations selected = {}'.format(value)

#color radiobuttons output callback
@app.callback(
    dash.dependencies.Output('output-uniquecolor-container', 'children'),
    [dash.dependencies.Input('unique-color-dropdown', 'value'), Input('cache', 'data')])
def update_unique_color(value, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    update_parameters(value, "dark_color_A", folder)
    print("callback update_unique_color")
    return ""
    
@app.callback(
    dash.dependencies.Output('output-seedcolor-container', 'children'),
    [dash.dependencies.Input('seeds-color-dropdown', 'value'), Input('cache', 'data')])
def update_seeds_color(value, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    update_parameters(value, "soft_color_A", folder)
    print("callback update_seeds_color")
    return ""
      
#Button
@app.callback(
    dash.dependencies.Output('container-button-basic', 'children'),
    [dash.dependencies.Input('submit-val', 'n_clicks'), Input('cache', 'data')],
    [dash.dependencies.State('input-on-submit', 'value')])
def update_button(n_clicks, folder, value):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    try:
        if n_clicks > 1:
            cmd_str="python3 "+BaseFolder+"main_proccess.py "+folder
            proc = subprocess.Popen([cmd_str], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
            print("callback update_button")  
            return dbc.Spinner(spinner_style={"width": "6rem", "height": "6rem"})
    except:
        return ""

@app.callback(Output("download", "data"), [Input("btn", "n_clicks"), Input('cache', 'data')])
def download_file(n_nlicks, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    print("callback download_file")
    return send_file(folder+"AllCompressed/results.zip")

##Add callbacks for download instructions and download example
@app.callback(Output("download-examp", "data"), [Input("btn-examp", "n_clicks")])
def download_example_file(n_nlicks):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0] 
    if 'btn-examp' in changed_id:
        print("callback download_file")
        return send_file(BaseFolder+"data/example_input.txt")
    else:
        print()

@app.callback(Output("download-inst", "data"), [Input("btn-inst", "n_clicks")])
def download_instructions_file(n_nlicks):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0] 
    if 'btn-inst' in changed_id:
        print("callback download_file")
        return send_file(BaseFolder+"data/instructions.txt")
    else:
        print()

@app.callback(Output("download-treatment-enrichment", "data"),
              [Input("btn-enrichment_treatment", "n_clicks"), Input('cache', 'data')])
def download_treatment_enrichment_file(n_nlicks, folder):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0] 
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    if 'btn-enrichment_treatment' in changed_id:
        ParametersDict = loadParametersDict(folder)
        print("callback download_file")
        return send_file(folder+ParametersDict["treatment_col"]+"_Metabolites_pathway.csv")
    else:
        print()

@app.callback(Output("download-comparison-enrichment", "data"),
              [Input("btn-enrichment_comparison", "n_clicks"), Input('cache', 'data')])
def download_comparison_instructions_file(n_nlicks, folder):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0] 
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")        
    if 'btn-enrichment_comparison' in changed_id:
        ParametersDict = loadParametersDict(folder)
        print("callback download_file")
        return send_file(folder+ParametersDict["comparison_col"]+"_Metabolites_pathway.csv")
    else:
        print()     


@app.callback(Output("download-comparison-highres", "data"),
              [Input("btn-highres_comparison", "n_clicks"), Input('cache', 'data')])
def download_comparison_highres_file(n_nlicks, folder):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0] 
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")        
    if 'btn-highres_comparison' in changed_id:
        ParametersDict = loadParametersDict(folder)
        print("callback download_file")
        return send_file(folder+ParametersDict["comparison_col"]+"_Network.png")
    else:
        print()     

@app.callback(Output("download-treatment-highres", "data"),
              [Input("btn-highres_treatment", "n_clicks"), Input('cache', 'data')])
def download_treatment_highres_file(n_nlicks, folder):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0] 
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")        
    if 'btn-highres_treatment' in changed_id:
        ParametersDict = loadParametersDict(folder)
        print("callback download_file")
        return send_file(folder+ParametersDict["treatment_col"]+"_Network.png")
    else:
        print()          
        
@app.callback(Output('my-checklist1', 'children'),
              [Input('my-checklist1', 'value'), Input('cache', 'data')])
def prepare_data(value, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    if value:
        update_parameters(value, "Enriched_pathways", folder)
        #parametersDict["Enriched_pathways"] = value
        #print(parametersDict["Enriched_pathways"])
        with open(folder+"keep_pathways.txt", 'w') as f:
            for item in value:
                f.write("%s\n" % item)          
        return loadParametersDict(folder)["Enriched_pathways"]

@app.callback(Output('container-final-results', "children"),
              [Input('final-results-listener', "n_intervals"),
               Input('cache', 'data')])
def check_results_creation(n_intervals,folder):
   folder = json.loads(folder)
   folder = str(folder).strip("\"").strip("\'")   
   my_file = Path(folder+"main_proccess_results_html.pkl")
   print("waiting for main proccess results")
   print(folder+"main_proccess_results_html.pkl")
   if my_file.is_file():
       print("File was created!")
       file = open(folder+"main_proccess_results_html.pkl", 'rb')
       data = pickle.load(file)
       file.close()
       exploreSubNetworksDropdowns = dbc.Row([
                                       dbc.Col([download_treatment_enrichment_link,
                                               Explore_enriched_pathways_treatment(folder), 
                                               html.Div(id="treatmet-subgraph")]),
                                       dbc.Col([download_comparison_enrichment_link,
                                               Explore_enriched_pathways_control(folder), 
                                               html.Div(id="control-subgraph")])
                                       ])
       return [data, exploreSubNetworksDropdowns]

@app.callback(
            Output("treatmet-subgraph", "children"),
            [Input('explore-enriched-pathways-treatment-dropdown', "value"),
             Input('cache', 'data')]
            )
def present_subgraph_treatment(value, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'") 
    ParametersDict = loadParametersDict(folder)
    #load the list object of all subgraphs
    file = open(folder+ParametersDict["treatment_col"]+"_allSubGraphs.pickle", 'rb')
    subgraphs = pickle.load(file)
    file.close()    
    #select the subgraph according to value
    for i in subgraphs:
        if i[0] == value:
            #return the html object
            return i[1]
    
@app.callback(
            Output("control-subgraph", "children"),
            [Input('explore-enriched-pathways-control-dropdown', "value"),
             Input('cache', 'data')]
            )
def present_subgraph_control(value, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'") 
    ParametersDict = loadParametersDict(folder)
    #load the list object of all subgraphs
    file = open(folder+ParametersDict["comparison_col"]+"_allSubGraphs.pickle", 'rb')
    subgraphs = pickle.load(file)
    file.close()    

    #select the subgraph according to value
    for i in subgraphs:
        if i[0] == value:
            #return the html object
            return i[1]

@app.callback(
            Output('container-button-basic', component_property='style'),
            [Input('final-results-listener', "n_intervals"),
             Input('cache', 'data')]
            )
def toggle_collapse(n_intervals, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    my_file = Path(folder+"main_proccess_results_html.pkl")
    if my_file.is_file():
       return {'display':'none'}
    else:
       return {'display':'inline'}

@app.callback(
            Output("download_btn_placeholder", "children"),
            [Input('final-results-listener', "n_intervals"),
             Input('cache', 'data')]
            )
def download_finished_results_btn(n_intervals, folder):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    my_file = Path(folder+"main_proccess_results_html.pkl")
    if my_file.is_file():
       return download_link

#Make similar callback to shutdown sniffer
@app.callback(
            Output('final-results-listener', "disabled"),
            [Input('final-results-listener', "n_intervals"),
             Input('cache', 'data')],
            [State('final-results-listener', 'disabled')]
            )
def stop_sniffer_when_results_finished(n_intervals, folder, disabled_state):
    folder = json.loads(folder)
    folder = str(folder).strip("\"").strip("\'")    
    my_file = Path(folder+"main_proccess_results_html.pkl")
    if my_file.is_file():
       return not disabled_state
   

if __name__ == '__main__':
    app.run_server(debug=False, use_reloader=False,host='0.0.0.0')
