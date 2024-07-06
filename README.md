# UVM_CRISPR_DepMap_Analysis

## Description
Analysis to compare MAGeCK beta scores of our in-house Uveal Melanoma CRISPR screen with DepMap

## Data
Raw data is stored in "data"
- annotation.tsv - downloaded annotations from DepMap
- avana.tsv - MAGeCK beta scroes of avana data from DepMap
- MAGeCK_gene_corrected_beta.tsv - MAGeCK beta scores of UVM screen
- pan_genes.csv - list of pan essential genes from ...?

## Analyses
The scripts for analyses are in "scripts", which call functions in "src"
### Data Processing
Processed data is stored in "processed_data"
- "scripts/filter_depmap.R"

### Plotting
Output plots are in "plots"
- "src/boxplot.R" 
- "src/volcano_plot.R"

