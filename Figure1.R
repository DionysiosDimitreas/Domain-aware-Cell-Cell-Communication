# Load data
results_list <- readRDS("~/Downloads/Analysis/ThesisProject/results/dengue_1-026.rds")

myIsoformCellChatDB_old <- readRDS("~/Downloads/Analysis/ThesisProject/data/myIsoformCellChatDB.rds")

# Create hybrid dataset to handle cofactors
sce_hybrid <- augment_cofactors(results_list$sce_pc, results_list$sce_deeploc, myIsoformCellChatDB)

# Run CellChat on Isoform Level
cellchat_results_domain <- run_cellchat_isoform_level(sce_hybrid, "cell_type_refined")

results_list[["cc_domain"]] <- cellchat_results_domain

# Plot Figure 1A
# 1. Join to get Ligand Gene Names
# We perform a left join to keep all interactions, mapping 'ligand' (transcript) to 'gene_name'
transcript_to_gene <- data.frame(
  transcript_id = rownames(results_list$sce_baseline),
  gene_name = rowData(results_list$sce_baseline)$gene_symbol
)

interactions_mapped <- results_list$cc_domain$interactions %>%
  left_join(transcript_to_gene, by = c("ligand" = "transcript_id")) %>%
  select(-ligand) %>%
  rename(ligand = gene_name) %>%
  relocate(ligand, .after = 2)

# 2. Join to get Receptor Gene Names
# Map 'receptor' (transcript) to 'gene_name'
interactions_mapped <- interactions_mapped %>%
  left_join(transcript_to_gene, by = c("receptor" = "transcript_id")) %>%
  select(-receptor) %>%
  rename(receptor = gene_name) %>%
  relocate(receptor, .after = 3)

# 3. Group by the gene-level keys and sum the probabilities
# We group by source, target, and the new gene names.
# We also include 'pathway_name' to keep that annotation intact.
interactions_domain_genelvl <- interactions_mapped %>%
  group_by(source, target, ligand, receptor, pathway_name) %>%
  summarise(
    prob = 1 - prod(1 - prob, na.rm = TRUE),
    #prob = max(prob, na.rm = TRUE), #max instead of sum
    # For p-values, taking the minimum is a common approach when merging multiple significant hits
    pval = min(pval, na.rm = TRUE),
    .groups = "drop"
  )

interaction_counts_domain_genelvl <- interactions_domain_genelvl %>%
  # Count the number of rows for each Source-Target pair
  count(source, target, name = "Count") %>%
  # Rename columns to match your desired output format
  rename(Sender = source, Receiver = target)

counts_list <- list(
  BL = results_list$cc_baseline$interaction_counts, 
  PC = results_list$cc_pc$interaction_counts, 
  LOC = results_list$cc_deeploc$interaction_counts, 
  DOM = interaction_counts_domain_genelvl
)

counts_list <- lapply(counts_list, function(df) {
  # Identify column names as they might differ (source/target vs Sender/Receiver)
  s_col <- intersect(c("source", "Sender"), colnames(df))[1]
  r_col <- intersect(c("target", "Receiver"), colnames(df))[1]
  
  df %>% filter(
    !grepl("platelet", !!sym(s_col), ignore.case = TRUE),
    !grepl("platelet", !!sym(r_col), ignore.case = TRUE)
  )
})

# 1. Compare CondB (after) vs. CondA (before)
comp1 <- heatmap_cellchat_comparison(counts_list$BL, counts_list$PC)

# 2. Compare CondC (after) vs. CondB (before)
comp2 <- heatmap_cellchat_comparison(counts_list$BL, counts_list$LOC)

# 3. Compare CondD (after) vs. CondC (before)
comp3 <- heatmap_cellchat_comparison(counts_list$BL, counts_list$DOM)

p1 <- comp1$plot +
  labs(title = "LogFC (PC vs. BL)") +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
  )

p2 <- comp2$plot + 
  labs(title = "LogFC (LOC vs. BL)") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
  )

p3 <- comp3$plot + 
  labs(title = "LogFC (DOM vs. BL)") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
  )


p <- (p1 + p2 + p3) + 
  plot_layout(
    ncol = 3, 
    guides = "collect"
  )

# Figure 1B
combined_df <- bind_rows(
  "BL" = results_list$cc_baseline$interactions,
  "PC" = results_list$cc_pc$interactions,
  "LOC" = results_list$cc_deeploc$interactions,
  "DOM" = interactions_domain_genelvl,
  .id = "stage"
)

combined_df$stage <- factor(
  combined_df$stage, 
  levels = c("BL", "PC", "LOC", "DOM") # Enforce this order
)


results_pathways <- analyze_cell_filtering_impact(combined_df, sender = "Monocyte", receiver = "T cell", top_n = 10, focus_pathway = "all", max_pairs = 5)


print(results_pathways$plot_heatmap)