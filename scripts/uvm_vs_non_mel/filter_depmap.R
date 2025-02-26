library(dplyr)
library(readr)
library(vroom)

# Load depmap data
anno <- read.table(
  "data/annotation.tsv",
  sep = " ", header = TRUE
) |>
  select(DepMap_ID, primary_disease, Subtype)

avana <- vroom("data/avana.tsv")
avana <- na.omit(avana)

# Load UVM data
uvm_beta_scores <- read_tsv("data/MAGeCK_gene_corrected_beta.tsv")

# Divide avana data into pan-cancer group
non_melanoma_lines <- anno |>
  filter(!primary_disease %in% c("Unknown", "Non-Cancerous")) |>  
  filter(!grepl("melanoma", Subtype, ignore.case = TRUE)) |>  
  pull(DepMap_ID)

anno |>
  filter(DepMap_ID %in% melanoma_lines)

avana_nMel <- avana |>
    select(
      genes,
      any_of(non_melanoma_lines)
    )

# Save file with cell lines removed from avana
removed_lines <- anno |>
  filter(!DepMap_ID %in% non_melanoma_lines) |>
  filter(DepMap_ID %in% colnames(avana))
write_tsv(removed_lines, file.path("processed_data", "removed_lines.tsv"))

# Select only common genes in UVM data and depmap
common_genes <- intersect(uvm_beta_scores[["genes"]], avana[["genes"]])

filter_genes <- function(df, genes_to_keep) {
  filtered_df <- df |>
    filter(genes %in% genes_to_keep) |>
    arrange(genes)

  return(filtered_df)
}

process_and_write_files <- function(df, filename) {
  scores_df <- filter_genes(df, common_genes)
  write_tsv(scores_df, file.path("processed_data", filename))
}

datasets <- list(avana_nMel, uvm_beta_scores)
filenames <- list(
  "pan_cancer_scores.tsv",
  "uvm_scores.tsv"
)

purrr::walk2(datasets, filenames, process_and_write_files)
