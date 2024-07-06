# UVM_CRISPR_DepMap_Analysis

## Description
Analyses to compare MAGeCK beta scores of our in-house Uveal Melanoma (UVM) CRISPR screen with DepMap data.

## Data
Raw data are stored in "data":
- "annotation.tsv": Downloaded annotations from DepMap
- "avana.tsv": MAGeCK beta scroes of AVANA data from DepMap
- "MAGeCK_gene_corrected_beta.tsv": MAGeCK beta scores from the UVM screen
- "pan_genes.csv": list of pan essential genes from ...?

## Analyses
The analysis scripts are in "scripts", which call functions in "src".  
(Ignore "converted_notebook.R"!)

### Data Processing
Processed data are stored in "processed_data".

"filter_depmap.R":
- Divides DepMap data into two groups: "SKCM" and "Pan-cancer"
- Filters DepMap and UVM data for only common genes

### Calculate Fold Changes and Plotting
Results tables are in "processed_data".
Output plots are in "plots".

"uvm_vs_depmap.R":
- Ranks the beta scores in UVM, SKCM and Pan-cancer groups
- Calculates fold changes for 'UVM vs SKCM' and 'UVM vs Pan-cancer'
- Plots boxplots of ranks for the top 10 significant genes
- Plots volcano plots of fold changes, highlighting significant genes

