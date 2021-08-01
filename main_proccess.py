#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 09:33:17 2021

@author: ofir
"""

import numpy as np
import pandas as pd
import sys, os, random
import plotly.graph_objects as go
import dash_core_components as dcc
import dash_html_components as html
import networkx as nx
import dash_bootstrap_components as dbc
import base64
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#import ast
import plotly as py
import json
import time
from dash_extensions import Download

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BaseFolder = "./"
FULL_PATH = os.path.abspath(BaseFolder)+"/"

sys.path.append(FULL_PATH)
from netcom.netcom import pathwayEnrichment, EdgeR_to_seeds, simulation


try:
     with open(BaseFolder+"data/DB/DB.pickle", 'rb') as handle:
         DB = pickle.load(handle)
except:
    DB = pd.read_pickle(BaseFolder+"data/DB/DB.pickle")      
  
df_el_ = DB['full_enzymes_labels_jun.txt']
df_ecMapping_ = DB['ec_reac_mapping_jun.txt']
df_reactions_ = DB['reactions_3_balanced.txt']
df_ec_to_compoundIndex_ = DB['compound_labels_jun.txt']




def loadParametersDict(folder):
    start = time.time()
    
    f = open(folder+"parametersDict.json", "r")
    output = json.load(f)
    f.close()
    end = time.time()
    print("loadParametersDict time (sec):")
    print(end - start)    
    return output[0]
   
def CreateCompoundsNetwork_3D(folder, prefix,
                            G,
                            patches,
                            node_colors_Dict,
                            EnrichedNodesColorsDict,
                            FinalFigureName,
                            idd,
                            network_iter,                           
                            df_el = df_el_.copy(),
                            df_ecMapping = df_ecMapping_.copy(),
                            df_reactions = df_reactions_.copy(),
                            df_ec_to_compoundIndex = df_ec_to_compoundIndex_.copy(),
                            ):
   
    start = time.time()
    
    parametersDict = loadParametersDict(folder)
    FinalFolder_ = folder
    
    
    edges,colors = zip(*nx.get_edge_attributes(G,'color').items())
    edges,widths = zip(*nx.get_edge_attributes(G,'width').items())
    edges,texts = zip(*nx.get_edge_attributes(G,'texts').items())
    
    #texts for hover
    
    node_color = []
    node_alpha = []
    node_size = []
    count_colored_edges = 0
    
    for i in colors:
        if i != "gray": count_colored_edges += 1
    
    for node in list(G.nodes()):
        node_color.append(node_colors_Dict[node])
        
        if node_colors_Dict[node] == "gray":
            node_alpha.append(0.5)
            node_size.append(25)
        else:
            node_alpha.append(1)
            node_size.append(35)
     
   
    nodes_of_largest_component  = max(nx.connected_components(G), key = len)
    largest_component = G.subgraph(nodes_of_largest_component)
    #For reproduceability
    seed = 123
    random.seed(seed)
    np.random.seed(seed)
   
    pos1 = nx.spring_layout(G,k=0.075,iterations=network_iter, dim=3, seed=seed)
    
    pos2 = nx.spring_layout(G, pos=pos1,fixed=nodes_of_largest_component,
                           k=0.001,iterations=0, dim=3, seed=seed)





    
    ###VISUALIZATION

    
# =============================================================================

    PatchNodeNames = []
    EnrichedNodesLabelsDict = {}
    
    for i, vali in enumerate(patches):
        for j in vali[1]:
                PatchNodeNames.append(j)           
        
        if len(PatchNodeNames) > 0:
            for cl in PatchNodeNames:
                EnrichedNodesLabelsDict[cl] = vali[0]
                             
        PatchNodeNames = []

    #add all other compounds (which are not in enriched pathways) to colors dictionary
    for i in G.nodes():
        if not i in EnrichedNodesColorsDict:
            EnrichedNodesColorsDict[i] = "gray"

    for i in G.nodes():
        if not i in EnrichedNodesLabelsDict:
            EnrichedNodesLabelsDict[i] = " "


    with open(BaseFolder+"data/DB/kegg_data/KeggNumToLabelsDict.pickle", 'rb') as handle:
        KeggNumToLabelsDict = pickle.load(handle)  

    
    #Create nodes color list for the 3d trace3 enrichment circles
    enrichment_node_colors = [EnrichedNodesColorsDict[i] for i in G.nodes()]
        
    enrichment_node_labels = []
    for i in list(G.nodes()):
        try:
            enrichment_node_labels.append(KeggNumToLabelsDict[i]+ " " + EnrichedNodesLabelsDict[i])
        except:
            print(i+" was not found in kegg2Label dict")
            enrichment_node_labels.append(i+ " " + EnrichedNodesLabelsDict[i])
      
    #Plotly visualization 3D
    N=len(G.nodes())
    Xn=[list(pos2.values())[k][0] for k in range(N)]# x-coordinates of nodes
    Yn=[list(pos2.values())[k][1] for k in range(N)]# y-coordinates
    Zn=[list(pos2.values())[k][2] for k in range(N)]# z-coordinates
    Xe=[]
    Ye=[]
    Ze=[]
    for e in G.edges():
        Xe+=[pos2[e[0]][0],pos2[e[1]][0], None]# x-coordinates of edge ends
        Ye+=[pos2[e[0]][1],pos2[e[1]][1], None]
        Ze+=[pos2[e[0]][2],pos2[e[1]][2], None]

    
    print("Building trace1")
    trace1=go.Scatter3d(x=Xe,
                   y=Ye,
                   z=Ze,
                   mode='lines',
                   line=dict(color=colors, width=3),
                   text=texts,
                   hoverinfo='text'
                   )
    
    print("Building trace2")
    trace2=go.Scatter3d(x=Xn,
                   y=Yn,
                   z=Zn,
                   mode='markers',
                   name='Compounds',
                   marker=dict(symbol='circle',
                                 size=10,
                                 color=node_color,
                                 line=dict(color=enrichment_node_colors, width=8)
                                 ),
                   #text=list(G.nodes()),
                   text=enrichment_node_labels,
                   hoverinfo='text'
                   )
    
    axis=dict(showbackground=False,
              showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              title=''
              )
    
    layout = go.Layout(
             #title="Network of Compounds",             
             width=1000,
             height=1000,
             showlegend=False,
             scene=dict(
                 xaxis=dict(axis),
                 yaxis=dict(axis),
                 zaxis=dict(axis),
            ),
         margin=dict(
            t=100
        ),
        hovermode='closest',
        annotations=[
               dict(
               showarrow=False,
                text="",
                xref='paper',
                yref='paper',
                x=0,
                y=0.1,
                xanchor='left',
                yanchor='bottom',
                font=dict(
                size=14
                )
                )
            ],    )
    
    print("Building 3d fig")
    fig_Compounds = go.Figure(data=[trace1, trace2],
                  layout=layout
                     ) 

    title = FinalFigureName

    fig_Compounds.update_layout(
        title=title,
        font=dict(
            #family="Courier New, monospace",
            size=18,
            #color="RebeccaPurple"
        )
    )


    graph_compounds_3D = html.Div([
                            #html.H2(title),
                            dcc.Graph(
                                style={'height': '80%', 'width': '80%'},
                                id=idd,
                                figure=fig_Compounds)
                                ], className="four columns")    

    try:
        print("Saving html")   
        py.offline.plot(fig_Compounds, filename=FinalFolder_+"3D_network_"+prefix+"_"+title.replace(" ","_")+".html", auto_open=False)  
    except Exception as e: print(e)
    
    
    
    end = time.time()
    print("CreateCompoundsNetwork time (sec):")
    print(end - start)        
                     
    
    return graph_compounds_3D


def CreateCompoundsNetwork(folder, prefix,
                            edgeR_grouped,
                            All_compounds_A,
                            All_compounds_B,
                            Seeds_A,
                            Seeds_B,
                            ECs_A,
                            ECs_B,
                            ECs_All,
                            soft_color,#seeds
                            dark_color,# unique,  
                            filter_hubness,
                            FinalFigureName,
                            idd,
                            drop_fragment_with_size,
                            network_iter,                           
                            df_el = df_el_.copy(),
                            df_ecMapping = df_ecMapping_.copy(),
                            df_reactions = df_reactions_.copy(),
                            df_ec_to_compoundIndex = df_ec_to_compoundIndex_.copy(),
                            ):
   
    start = time.time()
    
    parametersDict = loadParametersDict(folder)
    FinalFolder_ = folder
    
    All_compounds_B = All_compounds_B+Seeds_B
    All_compounds_A = All_compounds_A+Seeds_A
    
    IndexEnzymeList = df_el[["Index", "Enzyme_Code"]].values.tolist()
    DictEL = {}
    
    for i in IndexEnzymeList:
        DictEL[i[1]] = i[0]
        DictEL[i[0]] = i[1]
        
    EClist_tmp = []
    Errors = []
    
    for i in ECs_All:
        try:
            EClist_tmp.append(DictEL[i])
        except:
            Errors.append(i)

    df_ecMapping = df_ecMapping[df_ecMapping["Index"].isin(EClist_tmp)]
    ListOfRelevantReactions = df_ecMapping["Reactions"].values.tolist()
    
    df_ecMapping["Reactions"] = [ i.split(":")[1].strip().split(" ") for i in ListOfRelevantReactions]
    
    #This to be used as Reaction to Enzyme_index 
    df_ecMapping_expolded = df_ecMapping[["Index", "Reactions"]].explode("Reactions")
    
    flat_list = []
            
    for i in ListOfRelevantReactions:
        for j in i.split(":")[1].strip().split(" "):
            try:
                flat_list.append(int(j.strip()))
            except:
                _=""
    
    #Save enzyme information of each reaction.
    df_reactions = df_reactions[df_reactions["Reaction_index"].isin(flat_list)]
    
    #Fix reaction directions - flip =
    tempList = []
    
    for i in df_reactions.values.tolist():
        if i[1] == "=":
            x = i
            x[3], x[2] = x[2], x[3]
            tempList.append(x)
    
    df_reactions = pd.concat([df_reactions, pd.DataFrame(tempList, columns=list(df_reactions.columns))])
        
    l = df_reactions["Reaction_Left"].values.tolist()
    r = df_reactions["Reaction_Right"].values.tolist()

    #l = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Reaction_Left"].values.tolist()]
    #r = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Reaction_Right"].values.tolist()]

          
    DictComp = {}
    
    for i in df_ec_to_compoundIndex[["Index", "Compound_Code"]].values.tolist():
        DictComp[i[1]] = i[0]
        DictComp[i[0]] = i[1]
            
    df = pd.DataFrame([l, r, df_reactions["Reaction_index"].values.tolist()]).T           
    df = df.explode(0)         
    df = df.explode(1)
    cols = ["Left", "Right", "Reaction"]
    df.columns = cols
    
    #map reaction to enzyme     
    ReactionToEnzymeDict = { i[1]:i[0] for i in df_ecMapping_expolded.values.tolist() }

    def map_enz(x):
        return ReactionToEnzymeDict[str(x)]
    
    print(ReactionToEnzymeDict)
    
    df["Enzyme"] = df["Reaction"].apply(map_enz)
    
    df_grouped = df.groupby(["Left","Right"])["Enzyme"].apply(list)
    df_grouped = df_grouped.to_frame().reset_index()
    
    def enzInd_to_enzEC(l):
        return [DictEL[i] for i in l] 
    
    df_grouped["Enzyme"] = df_grouped["Enzyme"].apply(enzInd_to_enzEC)
    
    
    #Fix Compound IDs
    def CompInd_to_CompC(x):
        try:
            return DictComp[int(x)]
        except:
            return np.nan
    
    df_grouped["Left"] = df_grouped["Left"].apply(CompInd_to_CompC)
    df_grouped["Right"] = df_grouped["Right"].apply(CompInd_to_CompC)
    
    df_grouped = df_grouped[(df_grouped["Left"].isin(All_compounds_A)) &
                            (df_grouped["Right"].isin(All_compounds_A))]    
    
    #Set node color
    node_colors_Dict = {}
    
    all_nodes = list(set(df_grouped["Left"].values.tolist()+df_grouped["Right"].values.tolist()))
    
    for i in Seeds_A:
        node_colors_Dict[i] = soft_color
    
    for i in list(np.setdiff1d(All_compounds_A, All_compounds_B)):
        if not i in node_colors_Dict:
            node_colors_Dict[i] = dark_color
    
    for i in All_compounds_B:
        if not i in node_colors_Dict:
            node_colors_Dict[i] = "gray"
    
    for i in all_nodes:
        if not i in node_colors_Dict:
            node_colors_Dict[i] = "gray"
    
    #set edges color
    df_grouped["Color"] = "gray"
    
    def isEnzymeUnique(l):
        if (len(list(set(l).intersection(set(ECs_A)))) > 0) and (len(list(set(l).intersection(set(ECs_B)))) == 0):
            return dark_color
        else:
            return "gray"
    
    df_grouped["Color"] = df_grouped["Enzyme"].apply(isEnzymeUnique)
    
    
    def setEdgeWidthByColor(x):
        if x != "gray":
            return 1
        else:
            return 0.75
        
    df_grouped["Edge_width"] = df_grouped["Color"].apply(setEdgeWidthByColor)
    
    patches = pathwayEnrichment(BaseFolder=BaseFolder,
                                FinalFolder=FinalFolder_,
                                DE_ecs_list_=All_compounds_A,
                                All_B = All_compounds_B,
                                All_ecs_list_=All_compounds_A+All_compounds_B,
                                input_type="metabolites", # "metabolites" or "enzymes"
                                outputfilename=FinalFigureName+"_Compounds",
                                minEntitiesInPathway=parametersDict["Min_entities_Enrichment"],
                                maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                keep_pathways=FinalFolder_+"keep_pathways.txt"
                                )
    
    patches_seeds = pathwayEnrichment(BaseFolder=BaseFolder,
                                        FinalFolder=FinalFolder_,
                                        DE_ecs_list_=Seeds_A,
                                        All_B = All_compounds_B,
                                        All_ecs_list_=All_compounds_A+All_compounds_B,
                                        input_type="metabolites", # "metabolites" or "enzymes"
                                        outputfilename=FinalFigureName+"_Resources",
                                        minEntitiesInPathway=parametersDict["Min_entities_Enrichment"],
                                        maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                        keep_pathways=FinalFolder_+"keep_pathways.txt"
                                        )
    
    patches_enzymes = pathwayEnrichment(BaseFolder=BaseFolder,
                                        FinalFolder=FinalFolder_,
                                        DE_ecs_list_=ECs_A,
                                        All_B = ECs_B,
                                        All_ecs_list_=ECs_A+ECs_B,
                                        input_type="enzymes", # "metabolites" or "enzymes"
                                        outputfilename=FinalFigureName+"_Enzymes",
                                        minEntitiesInPathway=parametersDict["Min_entities_Enrichment"],
                                        maxEntitiesInPathway=parametersDict["Max_entities_Enrichment"],
                                        keep_pathways=FinalFolder_+"keep_pathways.txt"
                                        )    
    
    #Save all enrichment analyses to files
    print("saving enrichment")
# =============================================================================
#     try:
#         with open(FinalFolder_+FinalFigureName+"_enrichment_unique_metabolites.txt", 'w') as f:
#             for item in patches:
#                 f.write("%s\n" % item)   
#     except Exception as e:
#         print(e)
#         print("no data to save")
# =============================================================================
        
# =============================================================================
#     try:
#         print("saving enrichment")
#         with open(FinalFolder_+FinalFigureName+"_enrichment_resource_metabolites.txt", 'w') as f:
#             for item in patches_seeds:
#                 f.write("%s\n" % item)   
#     except:
#         print("no data to save")
# =============================================================================
        
# =============================================================================
#     try:
#         print("saving enrichment")
#         with open(FinalFolder_+FinalFigureName+"_enrichment_enzymes.txt", 'w') as f:
#             for item in patches_enzymes:
#                 f.write("%s\n" % item)   
#     except:
#         print("no data to save")
# =============================================================================
        
    #filter by user input
    #try:
    patches.sort(key=len, reverse=True)    
    lenPatches = len(patches)
    lower_frac = parametersDict["enrichment_results_slice"][0]
    higher_frac = parametersDict["enrichment_results_slice"][1]

    patches = patches[int((lenPatches/100)*lower_frac):int((lenPatches/100)*higher_frac)]
    #except:
    #    patches = []
    
    #Filter df by nodes appearance
    temp_df = df_grouped["Left"].value_counts().reset_index()
    keeplist = temp_df[temp_df["Left"] <= filter_hubness]["index"].values.tolist()
    temp_df = df_grouped["Right"].value_counts().reset_index()
    keeplist += temp_df[temp_df["Right"] <= filter_hubness]["index"].values.tolist()
    keeplist = list(set(keeplist))
    
    df_grouped = df_grouped[(df_grouped["Left"].isin(keeplist)) & (df_grouped["Right"].isin(keeplist))]
    
    
    #Network
    G = nx.Graph()
    
    for i in df_grouped[["Left", "Right", "Color","Edge_width", "Enzyme"]].values.tolist():
        #print("adding edges")
        G.add_edge(i[0], i[1], color=i[2], width=i[3], texts=" ".join(list(set(i[4]))))#Added Enzyme
        
    #Filter out isolates
    #G.remove_nodes_from(list(nx.isolates(G)))
    
    #Filter out fragments
    print("filtering out small components")
    for component in list(nx.connected_components(G)):
        if len(component)<drop_fragment_with_size:
            for node in component:  
                G.remove_node(node)
    
    edges,colors = zip(*nx.get_edge_attributes(G,'color').items())
    edges,widths = zip(*nx.get_edge_attributes(G,'width').items())
    edges,texts = zip(*nx.get_edge_attributes(G,'texts').items())
    
    #texts for hover
    node_color = []
    node_alpha = []
    node_size = []
    count_colored_edges = 0
    
    for i in colors:
        if i != "gray": count_colored_edges += 1
    
    for node in list(G.nodes()):
        node_color.append(node_colors_Dict[node])
        
        if node_colors_Dict[node] == "gray":
            node_alpha.append(0.5)
            node_size.append(25)
        else:
            node_alpha.append(1)
            node_size.append(35)
 
    nodes_of_largest_component  = max(nx.connected_components(G), key = len)
    largest_component = G.subgraph(nodes_of_largest_component)
    
    #For reproduceability
    seed = 123
    random.seed(seed)
    np.random.seed(seed)
   
    pos1 = nx.spring_layout(G,k=0.075,iterations=network_iter, dim=3, seed=seed)
    
    pos2 = nx.spring_layout(G, pos=pos1,fixed=nodes_of_largest_component,
                           k=0.001,iterations=0, dim=3, seed=seed)

    pos1_2d = nx.spring_layout(G,k=0.075,iterations=network_iter, dim=2, seed=seed)
    
    pos2_2d = nx.spring_layout(G, pos=pos1_2d,fixed=nodes_of_largest_component,
                           k=0.001,iterations=0, dim=2, seed=seed)
    

    ###VISUALIZATION
    #2d matplotlib visualization
    #Labels only seeds and uniques
    nodes_labeldict = {}
    for i in G.nodes():
        if (i in Seeds_A) or (i in list(np.setdiff1d(All_compounds_A, All_compounds_B))):
            nodes_labeldict[i] = i
        else: nodes_labeldict[i] = ""
  

    fig = plt.figure(figsize=(50,50))
    
    nx.draw_networkx_edges(G, pos2_2d, alpha=0.25, width=widths, edge_color=colors, arrows=False)
    nx.draw_networkx_nodes(G, pos2_2d, node_color=node_color, node_size=node_size,
                                   alpha=node_alpha)
    nx.draw_networkx_labels(G,pos2_2d,nodes_labeldict,font_size=8,font_color='black')
    
# =============================================================================

    ax=plt.gca()
    PatchNodeNames = []
    #all_patch_colors =  list(mcd.XKCD_COLORS.keys())[::7]
    all_patch_colors = ['orange', 'green', 'blue', 'red', 'aqua', 'magenta', 'yellow', 'lawngreen', 'gray', 'palevioletred']
    EnrichedNodesColorsDict = {}
    EnrichedNodesLabelsDict = {}
    handles = []
    
    for i, vali in enumerate(patches):
        #node to position
        for j in vali[1]:
                PatchNodeNames.append(j)
        
        #Create legend colors -
        if len(PatchNodeNames) > 0:
            handles.append(mpatches.Patch(color=all_patch_colors[i], label=vali[0]))
            #add colored circles as patches
            for cl in PatchNodeNames:
                circle = plt.Circle(pos2_2d[cl], 0.01, color=all_patch_colors[i], fill=True, alpha=0.3)
                ax.add_artist(circle)
                #build colors dict
                try:
                    EnrichedNodesColorsDict[cl] = all_patch_colors[i]                
                except:
                    EnrichedNodesColorsDict[cl] = "wheat"     
                            
                EnrichedNodesLabelsDict[cl] = vali[0]
                
        PatchNodeNames = []

    handles2 = []
    #legend 2
    handles2.append(mpatches.Patch(color=parametersDict["soft_color_A"], label="Environmental resource"))
    handles2.append(mpatches.Patch(color=parametersDict["dark_color_A"], label="Compound"))
    
    # plot the legend
    leg1 = plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=('upper left'), fontsize=45, title="Pathways")
    for lh in leg1.legendHandles: 
        lh.set_alpha(0.5)
    
    leg1.get_title().set_fontsize('45') #legend 'Title' fontsize
    ax.add_artist(leg1)    
    leg2 = plt.legend(handles=handles2, bbox_to_anchor=(1.05, 0.25), loc=('upper left'), fontsize=45, title="Entities")
    leg2.get_title().set_fontsize('45') #legend 'Title' fontsize
    ax.add_artist(leg2)    

    #add all other compounds (which are not in enriched pathways) to colors dictionary
    for i in G.nodes():
        if not i in EnrichedNodesColorsDict:
            EnrichedNodesColorsDict[i] = "gray"

    for i in G.nodes():
        if not i in EnrichedNodesLabelsDict:
            EnrichedNodesLabelsDict[i] = " "

    with open(BaseFolder+"data/DB/kegg_data/KeggNumToLabelsDict.pickle", 'rb') as handle:
        KeggNumToLabelsDict = pickle.load(handle)  

    #Create nodes color list for the 3d trace3 enrichment circles
    enrichment_node_colors = [EnrichedNodesColorsDict[i] for i in G.nodes()]
    
    enrichment_node_labels = []
    for i in list(G.nodes()):
        try:
            enrichment_node_labels.append(KeggNumToLabelsDict[i]+ " " + EnrichedNodesLabelsDict[i])
        except:
            #print(i+" was not found in kegg2Label dict")
            enrichment_node_labels.append(i+ " " + EnrichedNodesLabelsDict[i])
    
    plt.axis('equal')    

    #Constant legend
    txt = "    Edges represent enzymes; nodes represent metabolites. Colored edges represent differentially abundant enzymes;\n\
    colored nodes represent environmental resources and treatment specific compounds (colors are selected by the user).\n\
    Nodes' background colors (wider circles around the nodes) represent pathways that are enriched\n\
    (FDR adjusted P value <= 0.05) with network components (nodes) that are unique to the treated samples."
    
    fig.text(.13, .05, txt, ha='left', fontsize=45)

    #Save 2d matplotlib image to unique folder in Results    
    fig.savefig(FinalFolder_+prefix+'_Network.png', bbox_inches='tight')
        
    plt.close() 
    image2D_base64 = base64.b64encode(open(FinalFolder_+prefix+'_Network.png', 'rb').read()).decode('ascii')
    image2D = html.Div(html.Img(src='data:image/png;base64,{}'.format(image2D_base64), style={'height':'100%', 'width':'100%'}))

    #Plotly visualization 3D
    N=len(G.nodes())
    Xn=[list(pos2.values())[k][0] for k in range(N)]# x-coordinates of nodes
    Yn=[list(pos2.values())[k][1] for k in range(N)]# y-coordinates
    Zn=[list(pos2.values())[k][2] for k in range(N)]# z-coordinates
    Xe=[]
    Ye=[]
    Ze=[]
    for e in G.edges():
        Xe+=[pos2[e[0]][0],pos2[e[1]][0], None]# x-coordinates of edge ends
        Ye+=[pos2[e[0]][1],pos2[e[1]][1], None]
        Ze+=[pos2[e[0]][2],pos2[e[1]][2], None]

    print("Building trace1")
    trace1=go.Scatter3d(x=Xe,
                   y=Ye,
                   z=Ze,
                   mode='lines',
                   line=dict(color=colors, width=3),
                   text=texts,
                   hoverinfo='text'
                   )
    
    print("Building trace2")
    trace2=go.Scatter3d(x=Xn,
                   y=Yn,
                   z=Zn,
                   mode='markers',
                   name='Compounds',
                   marker=dict(symbol='circle',
                                 size=10,
                                 color=node_color,
                                 line=dict(color=enrichment_node_colors, width=8)
                                 ),
                   #text=list(G.nodes()),
                   text=enrichment_node_labels,
                   hoverinfo='text'
                   )
    
    axis=dict(showbackground=False,
              showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              title=''
              )
    
    layout = go.Layout(
             #title="Network of Compounds",             
             width=1000,
             height=1000,
             showlegend=False,
             scene=dict(
                 xaxis=dict(axis),
                 yaxis=dict(axis),
                 zaxis=dict(axis),
            ),
         margin=dict(
            t=100
        ),
        hovermode='closest',
        annotations=[
               dict(
               showarrow=False,
                text="",
                xref='paper',
                yref='paper',
                x=0,
                y=0.1,
                xanchor='left',
                yanchor='bottom',
                font=dict(
                size=14
                )
                )
            ],    )
    
    print("Building 3d fig")
    fig_Compounds = go.Figure(data=[trace1, trace2],
                  layout=layout
                     ) 

    title = FinalFigureName

    fig_Compounds.update_layout(
        title=title,
        font=dict(
            size=18,
        )
    )

    graph_compounds = html.Div([
                            #html.H2(title),
                            dcc.Graph(
                                style={'height': '80%', 'width': '80%'},
                                id=idd,
                                figure=fig_Compounds)
                                ], className="four columns")    
   
    try:
        print("Saving html")   
        py.offline.plot(fig_Compounds, filename=FinalFolder_+"3D_network_"+title.replace(" ","_")+".html", auto_open=False)  
    except Exception as e: print(e)
    
    #Here calculate and save all sub-graphs, later will be presented instantly in the app. 
    #Save all sub-graphs as a dict. keys=pathways, values=html object 
    allEnrichedNetworks = []
    for counter, patch in enumerate(patches):
        
        H = G.subgraph(patch[1])
        
        #Reduce H according to patch
        subgraph = CreateCompoundsNetwork_3D(folder=folder,
                                             prefix=prefix+"_"+str(counter),
                                            G=H,
                                            patches=[patch],
                                            node_colors_Dict=node_colors_Dict,
                                            EnrichedNodesColorsDict=EnrichedNodesColorsDict,
                                            FinalFigureName=patch[0],
                                            idd='subplot-pathwayEnriched',
                                            network_iter=20,                           
                                            df_el = df_el_.copy(),
                                            df_ecMapping = df_ecMapping_.copy(),
                                            df_reactions = df_reactions_.copy(),
                                            df_ec_to_compoundIndex = df_ec_to_compoundIndex_.copy(),
                                            )
        
        allEnrichedNetworks.append([patch[0], subgraph])
        print("Construction all subgraphs")
    with open(folder+prefix+'_allSubGraphs.pickle', 'wb') as handle:
        pickle.dump(allEnrichedNetworks, handle, protocol=pickle.DEFAULT_PROTOCOL)

  
    end = time.time()
    print("CreateCompoundsNetwork time (sec):")
    print(end - start)        
                     
    return [graph_compounds, image2D]

def main_proccess(folder, df_el_=df_el_.copy(), df_ecMapping_=df_ecMapping_.copy(), df_reactions_=df_reactions_.copy(), df_ec_to_compoundIndex_=df_ec_to_compoundIndex_.copy()):
    print("starting main process")
    start = time.time()

    try:
        os.mkdir(BaseFolder+"Results")
    except:
        print()
    
    df_ = pd.read_csv(folder+"input_edger.csv")
    
    parametersDict = loadParametersDict(folder)
    
    FinalFolder_ = folder
    ###A###
   
    df_edgeR_grouped=df_.set_index("association")
    Seeds_A_input, T1_seeds_tag, ECs_A_input, Seeds_B_input, T2_seeds_tag, ECs_B_input, Seeds_All_input, ECs_All_input=EdgeR_to_seeds(edgeR_row_location=FinalFolder_+"raw_input_edger.csv",
                                                                                                          col_treatment_1=parametersDict["treatment_col"],
                                                                                                          col_treatment_2=parametersDict["comparison_col"],
                                                                                                          outputFolder=FinalFolder_,
                                                                                                          input_sep=",")
  
    sim_df_T1, steps_df_T1 = simulation(input1=ECs_All_input.copy(), input2=[T1_seeds_tag.copy()], resfolder=FinalFolder_, prefix=parametersDict["treatment_col"])
    sim_df_T2, steps_df_T2 = simulation(input1=ECs_All_input.copy(), input2=[T2_seeds_tag.copy()], resfolder=FinalFolder_, prefix=parametersDict["comparison_col"])

    All_compounds_A_sim = sim_df_T1.values.tolist()[0][0]
    All_compounds_B_sim = sim_df_T2.values.tolist()[0][0]
    
    #save to files
    text_file = open(FinalFolder_+parametersDict["treatment_col"]+"_compounds.txt", "w")
    text_file.writelines("; ".join(All_compounds_A_sim))
    text_file.close()    

    #save to files
    text_file = open(FinalFolder_+parametersDict["comparison_col"]+"_compounds.txt", "w")
    text_file.writelines("; ".join(All_compounds_B_sim))
    text_file.close()    

    #try:
    graph_compounds_A = CreateCompoundsNetwork(folder=folder,
                        edgeR_grouped=df_edgeR_grouped.copy(),
                        idd='compound_A-network-graph',
                        All_compounds_A = All_compounds_A_sim,
                        All_compounds_B = All_compounds_B_sim,
                        Seeds_A = T1_seeds_tag.copy(),
                        Seeds_B = T2_seeds_tag.copy(),
                        ECs_A = ECs_A_input.copy(),
                        ECs_B = ECs_B_input.copy(),
                        ECs_All = ECs_All_input.copy(),
                        df_el = df_el_.copy(),
                        df_ecMapping = df_ecMapping_.copy(),
                        df_reactions = df_reactions_.copy(),
                        df_ec_to_compoundIndex = df_ec_to_compoundIndex_.copy(),
                        drop_fragment_with_size = parametersDict["drop_fragment_with_size"],#size of fragments to drop
                        filter_hubness = parametersDict["filter_hubness"],#filter by max number of node connections
                        soft_color = parametersDict["soft_color_A"],#seeds
                        dark_color = parametersDict["dark_color_A"],
                        FinalFigureName=parametersDict["treatment_col"],
                        network_iter=parametersDict["network_layout_iter"],
                        prefix=parametersDict["treatment_col"])
    
    graph_compounds_A_2D = graph_compounds_A[1]
    graph_compounds_A = graph_compounds_A[0]
    
    graph_compounds_B = CreateCompoundsNetwork(folder=folder,
                        edgeR_grouped=df_edgeR_grouped.copy(),
                        idd='compound_B-network-graph',
                        All_compounds_A = All_compounds_B_sim,
                        All_compounds_B = All_compounds_A_sim,
                        Seeds_A = T2_seeds_tag.copy(),
                        Seeds_B = T1_seeds_tag.copy(),
                        ECs_A = ECs_B_input.copy(),
                        ECs_B = ECs_A_input.copy(),
                        ECs_All = ECs_All_input.copy(),
                        df_el = df_el_.copy(),
                        df_ecMapping = df_ecMapping_.copy(),
                        df_reactions = df_reactions_.copy(),
                        df_ec_to_compoundIndex = df_ec_to_compoundIndex_.copy(),
                        drop_fragment_with_size = parametersDict["drop_fragment_with_size"],#size of fragments to drop
                        filter_hubness = parametersDict["filter_hubness"],#filter by max number of node connections
                        soft_color = parametersDict["soft_color_A"],#seeds
                        dark_color = parametersDict["dark_color_A"],
                        FinalFigureName=parametersDict["comparison_col"],
                        network_iter=parametersDict["network_layout_iter"],
                        prefix=parametersDict["comparison_col"])
       
    graph_compounds_B_2D = graph_compounds_B[1]
    graph_compounds_B = graph_compounds_B[0]
       
    #except Exception as e:
    #    print(e) 

    #zip all results
    from zipfile import ZipFile
    from os.path import basename
    # create a ZipFile object
    try:
        os.mkdir(FinalFolder_+"AllCompressed")
    except:
        print()    
        
    with ZipFile(FinalFolder_+"AllCompressed/results.zip", 'w') as zipObj:
       # Iterate over all the files in directory
       for folderName, subfolders, filenames in os.walk(FinalFolder_):
           for filename in filenames:
               #create complete filepath of file in directory
               filePath = os.path.join(folderName, filename)
               # Add file to zip
               print(filePath)
               if (not "results.zip" in filePath) and \
                   (not "parametersDict.json" in filePath) and \
                   (not "raw_input_edger.csv" in filePath) and \
                   (not "input_edger.csv" in filePath) and \
                   (not "keep_pathways.txt" in filePath) and \
                   (not "All_ECs.txt" in filePath) and \
                   (not "simulation_steps" in filePath) and \
                   (not "Final_results" in filePath) and \
                   (not "main_proccess_results_html" in filePath) and \
                   (not parametersDict["treatment_col"]+"_pathway.csv" in filePath) and \
                   (not parametersDict["comparison_col"]+"_pathway.csv" in filePath) and \
                   (not "allSubGraphs" in filePath):
                       zipObj.write(filePath, basename(filePath))
       zipObj.write("./data/Readme.txt", "Readme.txt")
    os.system("rm "+FinalFolder_+"*.txt")
    os.system("rm "+FinalFolder_+"*.csv")
    
    print("results are compressed")
    titles = dbc.Row([
                   dbc.Col(html.H3(parametersDict["treatment_col"], style={'whiteSpace': 'pre-line'})),
                   dbc.Col(html.H3(parametersDict["comparison_col"], style={'whiteSpace': 'pre-line'})) 
                    ])

    end = time.time()
    print("main_proccess time (sec):")
    print(end - start)
    
    download_treatment_enrichment_link = html.Div([html.Button("Download high resolution file", id="btn-highres_treatment"), Download(id="download-treatment-highres")])
    download_comparison_enrichment_link = html.Div([html.Button("Download high resolution file", id="btn-highres_comparison"), Download(id="download-comparison-highres")])

    #replace link with download button    
    res = html.Div([
                titles,
                dbc.Row([
                    dbc.Col([graph_compounds_A_2D,
                             download_treatment_enrichment_link]),
                    dbc.Col([graph_compounds_B_2D,
                             download_comparison_enrichment_link])
                        ]),        
                html.Div([
                    dbc.Row([
                        dbc.Col(graph_compounds_A),
                        dbc.Col(graph_compounds_B)
                        ])
                ])
            ])


    with open(folder+'main_proccess_results_html.pkl', 'wb') as output:
        pickle.dump(res, output, pickle.DEFAULT_PROTOCOL)

print(sys.argv[1])
main_proccess(folder=sys.argv[1])
