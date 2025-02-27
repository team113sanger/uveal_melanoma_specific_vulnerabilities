# UVM_VS_Pan_Cancer_MAGeCK_Analysis

## Description
Analyses to compare MAGeCK beta scores of 10 Uveal Melanoma (UVM) CRISPR screens with Pan-Cancer screens from the Avana DepMap dataset.

## Raw Data
Raw MAgeCK beta scores are stored in `data/raw`:
- `avana.tsv`: MAGeCK beta scroes of DepMap screens
- `MAGeCK_gene_corrected_beta.tsv`: MAGeCK beta scores of UVM screens  

Metadata are stored in `data/metadata`
- `annotation.tsv`: Downloaded annotations from DepMap
- `pan_genes.csv`: List of pan essential genes identified using ProdeTool (https://github.com/cantorethomas/prodeTool)
- `avana_removed_cell_lines.tsv`: List of cell lines removed from avana.tsv for Pan-Cancer group in analysis  

## Analyses

### Data Pre-processing
Processed MAGeCK beta scores are in `data/processed`  
`scripts/01.preprocessing/prepare_beta_scores.R` was used to:
- Filter for only common genes between the UVM and Avana libraries.
- Filter Avana data to exclude melanoma and uveal melanoma cell lines, creating a 'pan-cancer' group for analysis.

### Calculate Fold Changes
Analysis scripts are in `scripts/02.analysis`  
Output results tables are in `results/tables`.  

`scripts/02.analysis/uvm_vs_pan_cancer_analysis.R`:
- Ranks the beta scores in UVM and Pan-cancer groups
- Calculates fold changes for 'UVM vs Pan-cancer'
- Finds top 10 UVM specific genes

### Plotting 
Plotting scripts are in `scripts/03.visualisation`  
Output plots are in `results/figures`.

`scripts/03.visualisation/uvm_vs_pan_cancer_plotting.R`:
- Plots boxplots of ranks for the top 10 UVM specific significant genes
- Plots volcano plots of fold changes, highlighting significantly upregulated or downregulated genes