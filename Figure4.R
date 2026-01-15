library(DropletUtils)
library(Matrix)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(CellChat)

# ==============================================================================
# 1. HELPER: Fast Sparse Aggregation
# ==============================================================================
aggregate_gene_helper <- function(sce, assay_name, gene_id_col) {
  gene_ids <- rowData(sce)[[gene_id_col]]
  
  # 1. Get the matrix
  mat <- assay(sce, assay_name)
  
  # 2. Handle Sparse Matrix (The source of your error)
  if (inherits(mat, "dgCMatrix") || inherits(mat, "sparseMatrix")) {
    mat <- as.matrix(mat)
  }
  
  # 3. Efficient Aggregation (Replaces slow 'aggregate' function)
  gene_expr_mat <- rowsum(mat, group = gene_ids)
  
  return(gene_expr_mat)
}

# ==============================================================================
# 2. HELPER: Robust CellChat Wrapper (Bypasses Fragile Steps)
# ==============================================================================
run_interaction_permutation <- function(sce_deeploc,
                                        sce_baseline,
                                        cc_deeploc_results, 
                                        cc_baseline_results,
                                        n_perm = 10,
                                        seed = 42,
                                        assay_name = "counts",
                                        save_plot_path = NULL) {
  
  message("--- Setup ---")
  # We want to downsample the baseline pool to match the depth of deeploc
  target_reads <- sum(assay(sce_deeploc, assay_name))
  raw_counts_pool <- assay(sce_baseline, assay_name)
  total_available_reads <- sum(raw_counts_pool)
  prob <- min(target_reads / total_available_reads, 1)
  
  message(paste("Downsampling probability:", round(prob, 4)))
  
  # Extract pre-calculated interaction counts
  ref_counts <- cc_deeploc_results$interaction_counts # 'Real Reduced'
  base_counts <- cc_baseline_results$interaction_counts # 'Baseline'
  
  set.seed(seed)
  random_seeds <- sample(1:10000, size = n_perm)
  all_logFC <- list()
  
  counts_matrix_integer <- as(round(raw_counts_pool), "CsparseMatrix")
  
  # --- RUN PERMUTATIONS ---
  for (i in seq_len(n_perm)) {
    message(paste("Running permutation", i, "of", n_perm))
    set.seed(random_seeds[i])
    
    # Downsample the pool
    downsampled_counts <- DropletUtils::downsampleMatrix(counts_matrix_integer, prop = prob)
    
    # Create temp SCE and run CellChat
    sce_temp <- SingleCellExperiment(list(counts = downsampled_counts))
    rowData(sce_temp) <- rowData(sce_baseline) 
    colData(sce_temp) <- colData(sce_baseline)
    
    gene_expr <- aggregate_gene_helper(sce_temp, assay_name = "counts", gene_id_col = "gene_symbol")
    meta <- as.data.frame(colData(sce_temp))
    
    cc_random <- run_standard_cc(gene_expr, meta, "cell_type_refined")
    random_counts <- cc_random$interaction_counts
    
    # Calculate LogFC: (Random + 1) / (Baseline + 1)
    comparison <- full_join(base_counts, random_counts, by = c("Sender", "Receiver"), suffix = c("_base", "_random")) %>%
      mutate(
        Count_base = tidyr::replace_na(Count_base, 0),
        Count_random = tidyr::replace_na(Count_random, 0),
        logFC = log2((Count_random + 1) / (Count_base + 1)),
        Permutation = i,
        Type = "Random"
      )
    
    all_logFC[[i]] <- comparison %>% dplyr::select(Sender, Receiver, logFC, Type)
  }
  
  real_comparison <- full_join(base_counts, ref_counts, by = c("Sender", "Receiver"), suffix = c("_base", "_real")) %>%
    mutate(
      Count_base = tidyr::replace_na(Count_base, 0),
      Count_real = tidyr::replace_na(Count_real, 0),
      logFC = log2((Count_real + 1) / (Count_base + 1)),
      Type = "Real"
    ) %>%
    dplyr::select(Sender, Receiver, logFC, Type)
  
  combined_df <- bind_rows(bind_rows(all_logFC), real_comparison) %>%
    mutate(Interaction = paste(Sender, "→", Receiver)) %>%
    filter(Sender != "platelet", Receiver != "platelet")
  
  # Order by median logFC of the random permutations
  interaction_order <- combined_df %>%
    filter(Type == "Random") %>%
    group_by(Interaction) %>%
    summarise(m = median(logFC)) %>%
    arrange(desc(m)) %>%
    pull(Interaction)
  
  combined_df$Interaction <- factor(combined_df$Interaction, levels = interaction_order)
  
  # Separate for plotting
  df_random <- combined_df %>% filter(Type == "Random")
  df_real <- combined_df %>% filter(Type == "Real")
  
  p <- ggplot(df_random, aes(x = Interaction, y = logFC)) +
    geom_boxplot(fill = "gray95", color = "gray40", outlier.shape = NA) +
    geom_point(data = df_real, aes(x = Interaction, y = logFC), 
               color = "firebrick", shape = 8, size = 3, stroke = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          panel.grid.minor = element_blank()) +
    labs(title = "Interaction Significance Test",
         subtitle = "Boxplots: Random Downsampling | Red Star: Real Data (DeepLoc)",
         x = "Cell-Cell Interaction", 
         y = "log2 Fold Change (vs Baseline)")
  
  if (!is.null(save_plot_path)) ggsave(save_plot_path, plot = p, width = 14, height = 7)
  
  return(list(plot = p, data = combined_df))
}

# Load data
results_list <- readRDS("results/dengue_1-026.rds")

results_permutation <- run_interaction_permutation(
  sce_deeploc = results_list$sce_deeploc,
  sce_baseline = results_list$sce_baseline,
  cc_deeploc_results = results_list$cc_deeploc, 
  cc_baseline_results = results_list$cc_baseline, 
  n_perm = 10
)

results_permutation$plot