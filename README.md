# UVM_CRISPR_DepMap_Analysis

## Description
Analyses to compare MAGeCK beta scores of Uveal Melanoma (UVM) CRISPR screens with the Avana DepMap dataset.

## Data
Raw data are stored in `data`:
- `annotation.tsv`: Downloaded annotations from DepMap
- `avana.tsv`: MAGeCK beta scroes of AVANA data from DepMap
- `MAGeCK_gene_corrected_beta.tsv`: MAGeCK beta scores from the UVM screen
- `pan_genes.csv`: list of pan essential genes from ...?

## Analyses
The analysis scripts are in `scripts`, which call functions in `src`.

### Data Processing
Processed data are stored in `processed_data`.

`filter_depmap.R`:
- Divides DepMap data into two groups: "SKCM" and "Pan-cancer"
- Filters DepMap and UVM data for only common genes

### Calculate Fold Changes and Plotting
Results tables are in `results`.  
Output plots are in `plots`.

`uvm_vs_depmap.R`:
- Ranks the beta scores in UVM and Pan-cancer groups
- Calculates fold changes for 'UVM vs Pan-cancer'
- Plots boxplots of ranks for the top 10 UVM specific significant genes
- Plots volcano plots of fold changes, highlighting significantly upregulated or downregulated genes

