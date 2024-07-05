library(dplyr)
library(readr)

# Load depmap data
avana_sk_scores <- read_tsv("processed_data/avana_sk_scores.tsv")
avana_Nsk_scores <- read_tsv("processed_data/avana_Nsk_scores.tsv")

# Load UVM data
uvm_scores <- read_tsv("processed_data/uvm_scores.tsv")

# Get ranks
rank_scores <- function(df) {
  ranked_df <- cbind(df[, "genes"], apply(df[,-1], 2, rank)) %>%
  as.data.frame() %>%
  mutate(across(-genes, as.numeric))
  
  return(ranked_df)
}

uvm_ranks <- rank_scores(uvm_scores)
write_tsv(uvm_ranks, "processed_data/uvm_ranks.tsv")

avana_sk_ranks <- rank_scores(avana_sk_scores)
write_tsv(avana_sk_ranks, "processed_data/avana_sk_ranks.tsv")

avana_Nsk_ranks <- rank_scores(avana_Nsk_scores)
write_tsv(avana_Nsk_ranks, "processed_data/avana_Nsk_ranks.tsv")