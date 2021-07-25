#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:26:56 2020

@author: ofir
"""

import pickle
import ast
import random
import networkx as nx
import numpy as np
from scipy.stats import fisher_exact
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import matplotlib.patches as mpatches



def pathwayEnrichment(BaseFolder,
                      FinalFolder,
                      DE_ecs_list_,
                      All_B,
                      All_ecs_list_,
                      input_type, # "metabolites" or "enzymes"
                      outputfilename,
                      minEntitiesInPathway,
                      maxEntitiesInPathway,
                      enrichmentMode="pathway", #"pathway" or "module"
                      correction_algorithm="fdr_bh",
                      produceOutputfiles=True,
                      drop_pathways=False,
                      keep_pathways=False
                      ): 
    try:
        print("running enrichment")
        drop_pathways_list = []
        if drop_pathways != False:
            with open(drop_pathways) as f:
                drop_pathways_list = f.readlines()
                drop_pathways_list = [x.strip().strip("\n") for x in drop_pathways_list] 
                #print(drop_pathways_list)
    
    
        keep_pathways_list = []
        if keep_pathways != False:
            with open(keep_pathways) as f:
                keep_pathways_list = f.readlines()
                keep_pathways_list = [x.strip().strip("\n") for x in keep_pathways_list] 
                #print(keep_pathways_list)
                
                
        if enrichmentMode == "pathway":
            pathwayToEntityDict = {}           
            if input_type == "enzymes":
                print("Enzyme - open pathwayToEnzDict")
                with open(BaseFolder+"data/DB/kegg_data/pathwayToEnzDict.pickle", 'rb') as handle:
                    pathwayToEntityDict = pickle.load(handle)
            if input_type == "metabolites":
                print("Compound - open pathwayToCompoundDict")            
                with open(BaseFolder+"data/DB/kegg_data/pathwayToCompoundDict.pickle", 'rb') as handle:
                    pathwayToEntityDict = pickle.load(handle)              
    
        elif enrichmentMode == "module":
            pathwayToEntityDict = {}    
            if input_type == "enzymes":
                print("Enzyme - open moduleToEnzDict")
                with open(BaseFolder+"data/DB/kegg_data/moduleToEnzDict.pickle", 'rb') as handle:
                    pathwayToEntityDict = pickle.load(handle)
            if input_type == "metabolites":
                print("Compound - open moduleToCompoundDict")            
                with open(BaseFolder+"data/DB/kegg_data/moduleToCompoundDict.pickle", 'rb') as handle:
                    pathwayToEntityDict = pickle.load(handle)              
        
        print("calculating Fisher")
    
        def calculate_fisher(pathway_ecs_list,
                             pathway_name):
            
            InPathwayDE = list(set(DE_ecs_list_).intersection(set(pathway_ecs_list)))
            InPathwayNotDE = list(set(All_B).intersection(set(pathway_ecs_list)))            
            NotInPathwayDE = list(np.setdiff1d(DE_ecs_list_,pathway_ecs_list))
            NotInPathwayNotDE = list(np.setdiff1d(All_B, InPathwayNotDE))
            
            oddsratio, pvalue = fisher_exact([[len(InPathwayDE), len(InPathwayNotDE)],
                                              [len(NotInPathwayDE), len(NotInPathwayNotDE)]],
                                             alternative="greater")
            
            return oddsratio, pvalue, InPathwayDE, len(InPathwayDE), len(InPathwayNotDE), len(NotInPathwayDE), len(NotInPathwayNotDE)
        
        results = []
        for i in list(pathwayToEntityDict.keys()):
            Odds, Pval, enzymesDE, l_InPathwayDE, l_InPathwayNotDE, l_NotInPathwayDE, l_NotInPathwayNotDE = calculate_fisher(pathway_name = i,
                                                                                                             pathway_ecs_list = pathwayToEntityDict[i])
                                                                                
            results.append([i, Odds, Pval, enzymesDE, l_InPathwayDE, l_InPathwayNotDE, l_NotInPathwayDE, l_NotInPathwayNotDE])
            
        df_pathway_enrichment = pd.DataFrame(results, columns=["Pathway", "Odds", "P-val", "entities", "l_InPathwayDE", "l_InPathwayNotDE", "l_NotInPathwayDE", "l_NotInPathwayNotDE"])
        #print(df_pathway_enrichment)
        if len(drop_pathways_list) > 1:
            df_pathway_enrichment = df_pathway_enrichment[~df_pathway_enrichment["Pathway"].isin(drop_pathways_list)]
        
        if len(keep_pathways_list) > 1:
            df_pathway_enrichment = df_pathway_enrichment[df_pathway_enrichment["Pathway"].isin(keep_pathways_list)]
        
        df_pathway_enrichment["Count"] = df_pathway_enrichment["entities"].apply(lambda x: len(x))
        #discard empty lists
        df_pathway_enrichment = df_pathway_enrichment[df_pathway_enrichment.astype(str)['entities'] != '[]']
        
        #multiple test correction
        from statsmodels.stats.multitest import multipletests
              
        print("correcting P vals for multiple tests")
        df_pathway_enrichment = df_pathway_enrichment[df_pathway_enrichment["Count"] >= minEntitiesInPathway]
        df_pathway_enrichment = df_pathway_enrichment[df_pathway_enrichment["Count"] <= maxEntitiesInPathway]
        corrected_pvals = multipletests(pvals=df_pathway_enrichment['P-val'],
                                       alpha=0.05,
                                       method=correction_algorithm,
                                       is_sorted=False,
                                       returnsorted=False)
        
        df_pathway_enrichment["P-val_corrected"] = corrected_pvals[1]
        df_pathway_enrichment['P-val_corrected'] = df_pathway_enrichment['P-val_corrected'].astype(float)
        df_pathway_enrichment = df_pathway_enrichment.sort_values(['P-val_corrected'])       
        
        #Add compound names (real names, not C numbers) columns to df_pathway_enrichment
        with open(BaseFolder+"data/DB/kegg_data/KeggNumToLabelsDict.pickle", 'rb') as handle:
            KeggNumToLabelsDict = pickle.load(handle) 
            
        def CnumberToCompoundName(l):
            new_l = []
            for i in l:
                try:
                    new_l.append(KeggNumToLabelsDict[i])
                except:
                    #print(i+" is not in Dict")
                    print()
            return new_l
        
        if input_type=="metabolites":
            df_pathway_enrichment["entities_names"] = df_pathway_enrichment["entities"].apply(CnumberToCompoundName)
        
        if produceOutputfiles==True:
            df_pathway_enrichment.to_csv(FinalFolder+outputfilename+"_"+enrichmentMode+".csv")
               
        df_pathway_enrichment = df_pathway_enrichment[df_pathway_enrichment["P-val_corrected"] <= 0.05]
        patches = df_pathway_enrichment[["Pathway", "entities"]].values.tolist()       
    except Exception as e: print(e)
        
    if not 'patches' in locals():
        patches = []
        return None
    else:
        return patches
 

def EdgeR_to_seeds(edgeR_row_location, col_treatment_1, col_treatment_2, outputFolder, input_sep=","):
    edgeR_row = pd.read_csv(edgeR_row_location, sep=input_sep)
    BaseFolder = "./" #where the netcom.py file is
   
    try:
        df_edgeR_grouped = edgeR_row.groupby("association")['X'].apply(list).to_frame()
    except:
        df_edgeR_grouped = edgeR_row.groupby("association")['enzyme'].apply(list).to_frame()    

    #fix list strings to list
    try:
        enzymes_fixed = []
        for i in df_edgeR_grouped["enzyme"].values.tolist():
            #print(i)
            enzymes_fixed.append(ast.literal_eval(i))
        
        df_edgeR_grouped["enzyme"] = enzymes_fixed
    except:
        _=""
    
    df_edgeR_grouped_T = df_edgeR_grouped.T
    df_edgeR_grouped_T_cols = list(df_edgeR_grouped_T.columns)
    associated_cols = [col_treatment_1, col_treatment_2]
    not_associated_cols = [col for col in df_edgeR_grouped_T_cols if not col in associated_cols]
    print("col_treatment_1: "+ col_treatment_1)
    T1 = df_edgeR_grouped_T[col_treatment_1].values.tolist()[0]
    T2 = df_edgeR_grouped_T[col_treatment_2].values.tolist()[0]
    NA = df_edgeR_grouped_T[not_associated_cols[0]].values.tolist()[0]
    ALL = T1+T2+NA
      
    def extract_seeds(EClist):
        
        try:
             with open(BaseFolder+"data/DB/DB.pickle", 'rb') as handle:
                 DB = pickle.load(handle)
        except:
            DB = pd.read_pickle(BaseFolder+"data/DB/DB.pickle")      
          
        df_el = DB['full_enzymes_labels_jun.txt']
        df_ecMapping = DB['ec_reac_mapping_jun.txt']
        df_reactions = DB['reactions_3_balanced.txt']
        df_ec_to_compoundIndex = DB['compound_labels_jun.txt']
        rlist = df_ecMapping["Reactions"].values.tolist()   
        df_ecMapping["Reactions"] = [i.split(":")[1].lstrip().rstrip().split(" ") for i in rlist]    
        IndexEnzymeList = df_el[["Index", "Enzyme_Code"]].values.tolist()
        DictEL = {}
        
        for i in IndexEnzymeList:
            DictEL[i[1]] = i[0]
        
        #fix list strings to list
        #edgeR_row = ast.literal_eval(edgeR_row)
        
        #EClist = edgeR_row[0]
        print("Extracting EClist")
    
        EClist_ = []
        
        for i in EClist:
            try:
                EClist_.append(DictEL[i])
            except:
                _=""
                
        df_ecMapping = df_ecMapping[df_ecMapping["Index"].isin(EClist_)]
        ListOfRelevantReactions = df_ecMapping["Reactions"].values.tolist()
        flat_list = []
                
        for i in ListOfRelevantReactions:
            for j in i:
                try:
                    flat_list.append(int(j))
                except:
                    _=""
                    
        df_reactions = df_reactions[df_reactions['Reaction_index'].isin(flat_list)]
        #l = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Reaction_Left"].values.tolist()]
        #r = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Reaction_Right"].values.tolist()]
        l = df_reactions["Reaction_Left"].values.tolist()
        r = df_reactions["Reaction_Right"].values.tolist()     
        DictComp = {}
        
        for i in df_ec_to_compoundIndex[["Index", "Compound_Code"]].values.tolist():
            DictComp[i[0]] = i[1]        
        
        df_tmp = pd.DataFrame()
        df_tmp["l"]=l
        df_tmp["r"]=r
        df_tmp = df_tmp.explode("l").explode("r")
        df_tmp_grouped = pd.DataFrame({'count' : df_tmp.groupby( ["l","r"] ).size()}).reset_index()
    
        G_list = [ (i[0], i[1], i[2]) for i in df_tmp_grouped.values.tolist() ]
        G_seeds = nx.DiGraph()
        G_seeds.add_weighted_edges_from(G_list) 
        seeds = [len(c) for c in sorted(nx.strongly_connected_components_recursive(G_seeds), key=len)]    
        seeds += list(set(list(np.setdiff1d(df_tmp["l"].values,df_tmp["r"].values))))    
        MySeeds = [str(DictComp[int(i)]) for i in seeds]
        CompList = list(set(df_tmp["l"].values.tolist()+df_tmp["r"].values.tolist()))
        CompList = [str(DictComp[int(i)]) for i in CompList]
        MySeeds = list(set(MySeeds))
        CompList = list(set(CompList))
        EClist = list(set(EClist))
        print("Seeds were calculated")
        return MySeeds, CompList, EClist
 
    T1_seeds, T1_compounds, T1_ECs = extract_seeds(T1)
    T2_seeds, T2_compounds, T2_ECs = extract_seeds(T2)
    ALL_seeds, ALL_compounds, ALL_ECs = extract_seeds(ALL)
    
    T1_seeds_tag = list(set(T1_seeds).intersection(set(ALL_seeds)))
    T2_seeds_tag = list(set(T2_seeds).intersection(set(ALL_seeds)))
    
    print("Seeds were calculated")

    text_file = open(outputFolder+col_treatment_1+"_ECs.txt", "w")
    text_file.writelines("; ".join(T1_ECs))
    text_file.close()    
    
    text_file = open(outputFolder+col_treatment_2+"_ECs.txt", "w")
    text_file.writelines("; ".join(T2_ECs))
    text_file.close()        
    
    text_file = open(outputFolder+col_treatment_1+"_resources.txt", "w")
    text_file.writelines("; ".join(T1_seeds_tag))
    text_file.close()
    
    text_file = open(outputFolder+col_treatment_2+"_resources.txt", "w")
    text_file.writelines("; ".join(T2_seeds_tag))
    text_file.close()    

    text_file = open(outputFolder+"All_ECs.txt", "w")
    text_file.writelines("; ".join(ALL_ECs))
    text_file.close()    

    print("Seeds were extracted")  
    return T1_seeds, T1_seeds_tag, T1_ECs, T2_seeds, T2_seeds_tag, T2_ECs, ALL_seeds, ALL_ECs 


def simulation(input1, input2, resfolder, prefix):  
    BaseFolder = "./" #where the netcom.py file is
    try:
        with open(BaseFolder+"data/DB/DB.pickle", 'rb') as handle:
            DB = pickle.load(handle)
    except:
        DB = pd.read_pickle(BaseFolder+"data/DB/DB.pickle")    
        
    def is_consist(LeftSides, RightSides, directions, BucketOfCompounds, Reaction_Direction):
        IsConsist = []
        new = BucketOfCompounds
        for row in range(len(LeftSides)):
            if ((set(LeftSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '>'))) or \
               ((set(LeftSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '='))):   
                IsConsist.append(1)
                new = new+RightSides[row]
            elif ((set(RightSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '='))):       
                IsConsist.append(1)
                new = new+LeftSides[row]         
            else:
                IsConsist.append(0)
        return IsConsist, new
    
    print("Running simulation")
    organisms = ["All_ECs"]
    enzymes = input1    
    
    input1_df = pd.DataFrame(columns=["Organism", "Enzymes_List"])
    input1_df["Organism"] = organisms
    input1_df["Enzymes_List"] = input1_df["Enzymes_List"].astype("object")
      
    ## PATCH - FIX INPUT2##
    D={}

    for i in DB['compound_labels_jun.txt'][["Index", "Compound_Code"]].values.tolist():
        D[str(i[1]).strip()] = str(i[0]).strip()
        
    for i, vali in enumerate(input2):
        for j, valj in enumerate(vali):
            input2[i][j] = D[str(valj).strip()]
    
    #####    
    #Dictionary - relevant enzymes of each organism
    OrgEnzDict = {}
    
    #Dictionary - Compound index to Compound Code
    CompoundDict = {}

    for i in DB['compound_labels_jun.txt'][["Index", "Compound_Code"]].values.tolist():
        CompoundDict[i[0]] = i[1]
     
    #enviornments:
    listOfinput2OfOrganisms = []
    newOrganismsCompounds = []
    
    SimSteps_df = pd.DataFrame(columns=["Cycle", "Organism", "InitialEnv", "numberOfCompounds", "Compounds"])
    simulation_steps = []
    
    for count_,k in enumerate(input2):     
        #For each Organism:
        Enzymes_List = []
        k = [int(kk) for kk in k]
        newOrganismsCompounds = []
        for i, vali in enumerate(input1_df["Organism"].values):
            Enzymes_List = enzymes      
                
            #Enzymes_list code to index:
            ll = DB["full_enzymes_labels_jun.txt"][["Index", "Enzyme_Code"]].values.tolist()
            d = {}
            for j in ll:
                d[j[1]] = j[0]
                
            Enzymes_List_ = []
            for j in Enzymes_List:
                try:
                    Enzymes_List_.append(d[j])    
                except:
                    _=""
            #Find reactions consist input Compounds: check if all compounds in envirnment (k) are in the same side of the equation
            #AND an enzyme (j) is in the list --> TAKE the other side of equation according to equilibrium > < =   
            df = DB["ec_reac_mapping_jun.txt"][DB["ec_reac_mapping_jun.txt"]['Index'].isin(Enzymes_List_)]

            RelevantReactions = df["Reactions"].values.tolist()
            RelevantReactions = [x.lstrip().rstrip().split(":")[1].lstrip().rstrip().split(" ") for x in RelevantReactions]
            RelevantReactions = [item for sublist in RelevantReactions for item in sublist]
    
            for x, valx in enumerate(RelevantReactions):
                try:
                    RelevantReactions[x] = int(valx)
                except:
                    _=""
                                      
            RelevantReactions = [x for x in RelevantReactions if x != '']
    
            OrgEnzDict[i] = RelevantReactions
            df = DB["reactions_3_balanced.txt"][DB["reactions_3_balanced.txt"]['Reaction_index'].isin(RelevantReactions)]
            
            RL = df["Reaction_Left"].values.tolist()
            RR = df["Reaction_Right"].values.tolist()
            Reaction_Direction_ = df["Reaction_Direction"].values.tolist()       
            OnlyNewCompounds = []
            newCompounds = []
            prevCompounds=k
            C = 0
                   
            IC, newCompounds = is_consist(RL, RR, Reaction_Direction_, prevCompounds, Reaction_Direction=Reaction_Direction_)
            simulation_steps.append([count_, C, vali, k, len(set(newCompounds)), list(set(newCompounds))])
            
            OnlyNewCompounds = list(set(newCompounds).intersection(set(prevCompounds)) )
            
            while set(newCompounds) != set(prevCompounds): 
                OnlyNewCompounds = list(set(newCompounds).intersection(set(prevCompounds)) )          
                simulation_steps.append([count_, C, vali, k, len(set(newCompounds)), list(set(newCompounds))])
                prevCompounds = list(set(newCompounds))
                IC, newCompounds = is_consist(RL, RR, Reaction_Direction_, prevCompounds,Reaction_Direction=Reaction_Direction_)
                C += 1
    
        #1 - Filter by relevant enzyme (done) 2 - filter by relevant reactions (only reactions which are subsets of the final env)
        #3 - flatten the list
        
            newCompoundsCodes = [CompoundDict[k] for k in newCompounds]    
        
            if C > 0 :
                newOrganismsCompounds.append(list(set(newCompoundsCodes)))
            else:
                newOrganismsCompounds.append([])
                
        listOfinput2OfOrganisms.append(newOrganismsCompounds)     
    SimSteps_df = pd.DataFrame(simulation_steps, columns=["EnvIndex", "Cycle", "Organism", "InitialEnv", "numberOfCompounds", "Compounds"])
    Final_df = pd.DataFrame(listOfinput2OfOrganisms, columns=organisms)
    
    SimSteps_df.to_csv(resfolder+prefix+"_simulation_steps.csv")
    Final_df.to_csv(resfolder+prefix+"_Final_results.csv")    
    return Final_df, SimSteps_df


def CreateCompoundsNetwork_2D(BaseFolder,
                              FinalFolder,
                              prefix,
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
                              drop_fragment_with_size,
                              network_iter,
                              drop_pathways=False,
                              minEntitiesInPathway=4
                              ):

    try:
         with open(BaseFolder+"data/DB/DB.pickle", 'rb') as handle:
             DB = pickle.load(handle)
    except:
        DB = pd.read_pickle(BaseFolder+"data/DB/DB.pickle")      
      
    df_el = DB['full_enzymes_labels_jun.txt']
    df_ecMapping = DB['ec_reac_mapping_jun.txt']
    df_reactions = DB['reactions_3_balanced.txt']
    df_ec_to_compoundIndex = DB['compounds_lables_jun_1.txt']   
    
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
            _=""
    
    df_ecMapping = df_ecMapping[df_ecMapping["Index"].isin(EClist_tmp)]
    ListOfRelevantReactions = df_ecMapping["Reactions"].values.tolist()
    
    df_ecMapping["Reactions_list"] = [ i.split(":")[1].strip().split(" ") for i in ListOfRelevantReactions]
    
    #This to be used as Reaction to Enzyme_index 
    df_ecMapping_expolded = df_ecMapping[["Index", "Reactions_list"]].explode("Reactions_list")
    
    flat_list = []
            
    for i in ListOfRelevantReactions:
        for j in i.split(":")[1].strip().split(" "):
            try:
                flat_list.append(int(j.strip()))
            except:
                _=""
    
    #Save enzyme information of each reaction.
    df_reactions = df_reactions[df_reactions["Reactions"].isin(flat_list)]
    
    #Fix reaction directions - flip =
    tempList = []
    for i in df_reactions.values.tolist():
        if i[1] == "=":
            x = i
            x[3], x[2] = x[2], x[3]
            tempList.append(x)
    
    df_reactions = pd.concat([df_reactions, pd.DataFrame(tempList, columns=list(df_reactions.columns))])
    
    l = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Left"].values.tolist()]
    r = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Right"].values.tolist()]
      
    DictComp = {}
    
    for i in df_ec_to_compoundIndex[["Index", "Compound_Code"]].values.tolist():
        DictComp[i[1]] = i[0]
        DictComp[i[0]] = i[1]
        
    df = pd.DataFrame([l, r, df_reactions["Reactions"].values.tolist()]).T           
    df = df.explode(0)         
    df = df.explode(1)
    cols = ["Left", "Right", "Reaction"]
    df.columns = cols
    
    #map reaction to enzyme (col 4)
    ReactionToEnzymeDict = { i[1]:i[0] for i in df_ecMapping_expolded.values.tolist() }

    def map_enz(x):
        return ReactionToEnzymeDict[str(x)]
    
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
    
    #Save only rows contain compounds from the all_compounds list 
    df_grouped = df_grouped[(df_grouped["Left"].isin(All_compounds_A)) &
                            (df_grouped["Right"].isin(All_compounds_A))]    
    
    #Set node color
    node_colors_Dict = {}    
    all_nodes = list(set(df_grouped["Left"].values.tolist()+df_grouped["Right"].values.tolist()))
    #all_nodes = All_compounds_B
    
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
        
    All_compounds_A_unique = list(np.setdiff1d(All_compounds_A, All_compounds_B))
    All_compounds_B_unique = list(np.setdiff1d(All_compounds_B, All_compounds_A))
    
    patches = pathwayEnrichment(BaseFolder=BaseFolder,
                                FinalFolder=FinalFolder,
                                DE_ecs_list_=All_compounds_A_unique,
                                All_B = All_compounds_B_unique,
                                All_ecs_list_=All_compounds_A_unique+All_compounds_B_unique,
                                input_type="compound", # "compound" or "enzyme"
                                outputfilename=FinalFigureName+"_Compounds",
                                minEntitiesInPathway=minEntitiesInPathway,
                                produceOutputfiles=False,
                                drop_pathways=drop_pathways,
                                )
    
    #Filter df by nodes appearance
    temp_df = df_grouped["Left"].value_counts().reset_index()
    keeplist = temp_df[temp_df["Left"] <= filter_hubness]["index"].values.tolist()
    temp_df = df_grouped["Right"].value_counts().reset_index()
    keeplist += temp_df[temp_df["Right"] <= filter_hubness]["index"].values.tolist()
    keeplist = list(set(keeplist))
    
    df_grouped = df_grouped[(df_grouped["Left"].isin(keeplist)) & (df_grouped["Right"].isin(keeplist))]

    #Network
    G = nx.Graph()
    
    print("adding edges")
    for i in df_grouped[["Left", "Right", "Color","Edge_width", "Enzyme"]].values.tolist():
        G.add_edge(i[0], i[1], color=i[2], width=i[3], texts=" ".join(list(set(i[4]))))#Added Enzyme
        
    #Filter out isolates
    print("filtering out small components")    
    #Filter out fragments
    for component in list(nx.connected_components(G)):
        if len(component)<drop_fragment_with_size:
            for node in component:
                G.remove_node(node)
    
    edges,colors = zip(*nx.get_edge_attributes(G,'color').items())
    edges,widths = zip(*nx.get_edge_attributes(G,'width').items())
    edges,texts = zip(*nx.get_edge_attributes(G,'texts').items())#texts for hover
    
    node_color = []
    node_alpha = []
    node_size = []
    count_colored_edges = 0
    
    for i in colors:
        if i != "gray": count_colored_edges += 1
    
    for node in list(G.nodes()):
        node_color.append(node_colors_Dict[node])
        
        if node_colors_Dict[node] == "gray":
            node_alpha.append(0.25)
            node_size.append(25)
        else:
            node_alpha.append(0.75)
            node_size.append(35)
 
    nodes_of_largest_component  = max(nx.connected_components(G), key = len)
    largest_component = G.subgraph(nodes_of_largest_component)
    
    #For reproduceability
    seed = 123
    random.seed(seed)
    np.random.seed(seed)
    
    pos1_2d = nx.spring_layout(G,k=0.075,iterations=network_iter, dim=2, seed=seed)
    
    pos2_2d = nx.spring_layout(G, pos=pos1_2d,fixed=nodes_of_largest_component,
                           k=0.001,iterations=0, dim=2, seed=seed)
    

    ###VISUALIZATION
    #2d matplotlib visualization
    #Labels only seeds and uniques
    nodes_labeldict = {}
    for i in list(G.nodes()):
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
    all_patch_colors =  list(mcd.XKCD_COLORS.keys())[::7]
    EnrichedNodesLabelsDict = {}
    handles = []
    handles_nodes = []
    UsedCompoundsList = []
    for i, vali in enumerate(patches):
        #node to position
        for j in vali[1]:
            PatchNodeNames.append(j)
        
        #Create legend colors -
        if len(PatchNodeNames) > 0:                    
            #add colored circles as patches
            switch = 0
            for cl in PatchNodeNames:
                if (cl in list(G.nodes())) and not (cl in UsedCompoundsList):
                    switch = 1
                    #print("in G")
                    # plot circles using the RGBA colors
                    circle = plt.Circle(pos2_2d[cl], 0.01, color=all_patch_colors[i], fill=True, alpha=0.3)
                    #print(cl+" "+str(all_patch_colors[i])+" "+str(vali[0]))
                    ax.add_artist(circle)                    
                    EnrichedNodesLabelsDict[cl] = vali[0]
                    #print(cl)
            
            if switch == 1:
                handles.append(mpatches.Patch(color=all_patch_colors[i], label=vali[0])) 
                switch = 0
        
        UsedCompoundsList += PatchNodeNames     
        PatchNodeNames = []
 
    # plot the legend
    l1 = plt.legend(handles=handles, bbox_to_anchor=(1.05,0.9), loc=2, title_fontsize=45, fontsize=40, title="Pathways" )
    #print("legend: "+soft_color)
    handles_nodes.append(mpatches.Patch(color=soft_color, label="Resource"))
    #print("legend: "+dark_color)    
    handles_nodes.append(mpatches.Patch(color=dark_color, label="Unique"))
    
    plt.legend(handles=handles_nodes, bbox_to_anchor=(1.05, 1), loc=2, title_fontsize=45, fontsize=40, title="Nodes" ) 
    
    plt.gca().add_artist(l1)
    
    Title = prefix.strip().strip("_")+" Network"
    plt.title(Title, fontdict=None, fontsize=50, ha='center')
    plt.axis('equal')    
    
    #Save 2d matplotlib image to unique folder in Results    
    fig.savefig(FinalFolder+prefix+'_Network.png', bbox_inches='tight')
    
    plt.show() 