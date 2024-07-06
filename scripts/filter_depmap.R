library(dplyr)
library(readr)

# Load depmap data
anno <- read.table("data/annotation.tsv", sep = " ", header = TRUE) %>%
  select("DepMap_ID", "primary_disease")

avana <- read.table("data/avana.tsv",
  sep = " ",
  header = TRUE
)
colnames(avana) <- gsub("\\.", "-", colnames(avana))
avana <- na.omit(avana)

# Load UVM data
uvm_beta_scores <- read_tsv("data/MAGeCK_gene_corrected_beta.tsv")

# Tally depmap cancer types
primary_disease_counts <- table(anno$primary_disease)
primary_disease_counts

# Divide avana data into SKCM and Pan-cancer
filter_avana_by_cancer_type <- function(avana_df, annotations_df, cancer_types) {
  desired_cell_lines <- annotations_df %>%
    filter(primary_disease %in% cancer_types) %>%
    pull(DepMap_ID)

  filtered_avana <- avana_df %>%
    select(genes, any_of(desired_cell_lines))

  return(filtered_avana)
}

avana_sk <- filter_avana_by_cancer_type(avana, anno, "Skin Cancer")

non_skin_cancer_lines <- anno %>%
  filter(!primary_disease %in% c("Skin Cancer", "Unknown", "Non-Cancerous")) %>%
  distinct(primary_disease) %>%
  pull(primary_disease)
avana_Nsk <- filter_avana_by_cancer_type(avana, anno, non_skin_cancer_lines)

# Select only common genes in UVM data and depmap
common_genes <- intersect(uvm_beta_scores$genes, avana$genes)

avana_sk_scores <- avana_sk %>%
  filter(genes %in% common_genes)
write_tsv(avana_sk_scores, "processed_data/avana_sk_scores.tsv")

avana_Nsk_scores <- avana_Nsk %>%
  filter(genes %in% common_genes)
write_tsv(avana_Nsk_scores, "processed_data/avana_Nsk_scores.tsv")

uvm_scores <- uvm_beta_scores %>%
  filter(genes %in% common_genes)
write_tsv(uvm_scores, "processed_data/uvm_scores.tsv")
