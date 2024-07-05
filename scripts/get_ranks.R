library(dplyr)

# Load depmap data
anno <- read.table("raw_data/annotation.tsv", sep = " ", header = TRUE) %>%
  select("DepMap_ID", "primary_disease")

avana <- read.table("raw_data/avana.tsv", sep = " ",
                    header = TRUE)
colnames(avana) <- gsub("\\.", "-", colnames(avana))
avana <- na.omit(avana)

# Load UVM data
uvm_beta_scores <- read.delim("raw_data/MAGeCK_gene_corrected_beta.tsv")

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

avana_Nsk_scores <- avana_Nsk %>% 
  filter(genes %in% common_genes)

uvm_scores <- uvm_beta_scores %>% 
  filter(genes %in% common_genes)

# Get ranks
rank_scores <- function(df) {
  ranked_df <- cbind(df[, "genes"], apply(df[,-1], 2, rank)) %>%
  as.data.frame() %>% 
  rename(genes = V1) %>% 
  mutate(across(-genes, as.numeric))
  
  return(ranked_df)
}

uvm_ranks <- rank_scores(uvm_scores)
write.table(uvm_ranks, "processed_data/uvm_ranks.tsv",
            row.names = FALSE, col.names = TRUE,
            sep = '\t', quote = FALSE)

avana_sk_ranks <- rank_scores(avana_sk_scores)
write.table(avana_sk_ranks, "processed_data/avana_sk_ranks.tsv",
            row.names = FALSE, col.names = TRUE,
            sep = '\t', quote = FALSE)

avana_Nsk_ranks <- rank_scores(avana_Nsk_scores)
write.table(avana_Nsk_ranks, "processed_data/avana_Nsk_ranks",
            row.names = FALSE, col.names = TRUE,
            sep = '\t', quote = FALSE)