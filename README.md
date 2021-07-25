# NetCom 
This code is a local version of the tool described in the paper:
NetCom: a network-based tool for predicting metabolic activities of
microbial communities based on interpretation of metagenomics data
Authors: Ofir Tal, Rotem Bartuv, Maria Vetcos, Shlomit Medina, Jiandong
Jiang and Shiri Freilich

and availible online: 
https://freilich-lab-tools.com/netcom/

# Instructions
Input file: EdgeR results file, tab seperated consist the columns [enzyme	logFC	logCPM	PValue	FDR	association].
After a file is uploaded, a new section appears, enables the selection of the control and the experiment treatments. 
The selection of control and experiment treatments drives the calculatin of basic statistics which are now shown.
Entities number in a pathway: select a range of minimum and maximum number of entities linked with a pathway in the enrichment analysis.  
The selection of the range above a drives the enrichment analysis calculation and a new option is then enabled - 
Select pathways to dropout: select which pathways to drop out of the network.
Environmental resource node color: select the color of the compounds which are essential for the network to develop (seeds).
Unique node color: select the color of the compounds which are differential abunded due to the treatment.
Limit node hubness: Select the maximum allowed number of connected edges to a node.
Set network layout iterations: NetworkX network layout itterations. Higher number would produce arranged network, but consume more time to calculate (up to few minutes).

# The tool is implemented in Python 3.8 (Linux OS, Anaconda) using the following packages:
pandas 1.0.5
plotly 4.14.1
dash 1.18.1
dash_core_components 1.14.1
dash_html_components 1.1.1
networkx 2.4
dash_bootstrap_components 0.11.1
matplotlib 3.2.2
dash_extensions 0.0.41
numpy 1.18.5
scipy 1.5.0
statsmodels 0.11.1

# Parameters
fisher exact test - scipy.stats.fisher_exact, alternative hypothesis 'greater'

Test results and p-value correction for multiple tests - statsmodels.stats.multitest.multipletests default parameters, alpha=0.05

# Output files:
A dual set of results would be produced, one for each treatment, consists the following files:
keep_pathways.txt
raw_input_edger.csv
input_edger.csv
All_ECs.txt
<treatment_name>_simulation_steps.csv
<treatment_name>_resources_pathway.csv
<treatment_name>_resources.txt
<treatment_name>_pathway.txt
<treatment_name>_Network.png
<treatment_name>_Final_results.csv
<treatment_name>_Enzymes_pathway.csv
<treatment_name>_enrichment_unique_metabolites.txt
<treatment_name>_enrichment_resource_metabolites.txt
<treatment_name>_enrichment_enzymes.txt
<treatment_name>_ECSs.txt
<treatment_name>_metabolites_pathway.csv
<treatment_name>_metabolites.txt
3D_network_<treatment_name>.html

# run locally - recommended approach
mkdir netcom
cd netcom
virtualenv netcom
source netcom/bin/activate 
pip install -r requirements.txt

execute by:
python app.py

open browser at http://0.0.0.0:8050/netcom/

upload file for analysis



