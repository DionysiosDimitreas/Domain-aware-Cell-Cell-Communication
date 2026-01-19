run_domain_analysis <- function(results_list, db, species = "Human", assay_name = "counts", cell_type_column = "cell_type_refined", patient_id){
  # Augment cofactors
  sce_hybrid <- augment_cofactors(results_list$sce_pc, results_list$sce_deeploc, db)
  
  # Run CellChat
  cellchat_results_domain <- run_cellchat_isoform_level(sce_hybrid, "cell_type_refined")
  
  results_list[["cc_domain"]] <- cellchat_results_domain
  
  saveRDS(results_list, 
          file=paste0("~/Downloads/Analysis/ThesisProject/results/dengue_",patient_id,".rds"))
  
  return(results_list)
}

# Run analysis for the rest 7 samples: 3 Severe Dengue and 4 Controls
sample_list = c("1-013", "1-010" , "1-036", "3-013", "3-027", "3-018", "3-006")


for (patient_id in sample_list){
  print(patient_id)
  results_list_i <- readRDS(file = paste0("~/Downloads/Analysis/ThesisProject/results/dengue_",patient_id,".rds"))
  
  results_list_i <- run_domain_analysis(results_list = results_list_i,
                                        db = myIsoformCellChatDB,
                                        patient_id = patient_id)}

# -------- Differential Communication Analysis --------
data_dir <- "~/Downloads/Analysis/ThesisProject/results"
files <- list.files(path = data_dir, pattern = "\\.rds$", full.names = TRUE)

# Lists to store data
list_bl <- list()
list_dom <- list()

for (file in files) {
  sample_name <- tools::file_path_sans_ext(basename(file))
  message(paste("Processing Sample:", sample_name))
  full_data <- readRDS(file)
  
  # Determine Condition
  condition <- ifelse(grepl("1-0", sample_name, ignore.case = TRUE), "Dengue", 
                      ifelse(grepl("3-0", sample_name, ignore.case = TRUE), "Control", "Unknown"))
  
  # --- 1. Process Baseline (Gene-level) ---
  if (!is.null(full_data$cc_baseline$interactions)) {
    df_bl <- full_data$cc_baseline$interactions %>%
      mutate(group = condition, sample = sample_name, type = "Baseline")
    list_bl[[sample_name]] <- df_bl
  }
  
  # --- 2. Process Domain (Transcript-to-Gene) ---
  if (!is.null(full_data$sce_baseline) && !is.null(full_data$cc_domain$interactions)) {
    # Create map
    t2g <- data.frame(
      transcript_id = rownames(full_data$sce_baseline),
      gene_name = rowData(full_data$sce_baseline)$gene_symbol,
      stringsAsFactors = FALSE
    )
    
    df_dom <- full_data$cc_domain$interactions %>%
      left_join(t2g, by = c("ligand" = "transcript_id")) %>%
      rename(ligand_gene = gene_name) %>%
      left_join(t2g, by = c("receptor" = "transcript_id")) %>%
      rename(receptor_gene = gene_name) %>%
      group_by(source, target, ligand = ligand_gene, receptor = receptor_gene, pathway_name) %>%
      summarise(prob = sum(prob, na.rm = TRUE), pval = min(pval, na.rm = TRUE), .groups = "drop") %>%
      mutate(group = condition, sample = sample_name, type = "Domain")
    
    list_dom[[sample_name]] <- df_dom
  }
}

# --- 3. Combine and Filter ---
keep_pattern <- "B cell|Monocyte|NK|T cell"

master_bl <- bind_rows(list_bl) %>%
  filter(grepl(keep_pattern, source, ignore.case = TRUE) & grepl(keep_pattern, target, ignore.case = TRUE))

master_dom <- bind_rows(list_dom) %>%
  filter(grepl(keep_pattern, source, ignore.case = TRUE) & grepl(keep_pattern, target, ignore.case = TRUE))

# Statistical Analysis

# Run for both
stats_bl <- calc_diff_stats(master_bl) %>% mutate(method = "Baseline")
stats_dom <- calc_diff_stats(master_dom) %>% mutate(method = "Domain")

# Final Merged Result Table
final_comparison <- bind_rows(stats_bl, stats_dom)

# -------- PLOT 2a --------
# We sum the probabilities between cell types to get total weights
dengue_weights <- master_bl %>%
  filter(group == "Dengue") %>%
  group_by(source, target) %>%
  summarise(weight = sum(prob), .groups = "drop") %>%
  pivot_wider(names_from = target, values_from = weight, values_fill = 0) %>%
  column_to_rownames("source") %>%
  as.matrix()

control_weights <- master_bl %>%
  filter(group == "Control") %>%
  group_by(source, target) %>%
  summarise(weight = sum(prob), .groups = "drop") %>%
  pivot_wider(names_from = target, values_from = weight, values_fill = 0) %>%
  column_to_rownames("source") %>%
  as.matrix()

dengue_weights_dom <- master_dom %>%
  filter(group == "Dengue") %>%
  group_by(source, target) %>%
  summarise(weight = sum(prob), .groups = "drop") %>%
  pivot_wider(names_from = target, values_from = weight, values_fill = 0) %>%
  column_to_rownames("source") %>%
  as.matrix()

control_weights_dom <- master_dom %>%
  filter(group == "Control") %>%
  group_by(source, target) %>%
  summarise(weight = sum(prob), .groups = "drop") %>%
  pivot_wider(names_from = target, values_from = weight, values_fill = 0) %>%
  column_to_rownames("source") %>%
  as.matrix()

png("results/circle_plots.png", width = 1200, height = 1200, res = 150)

# 1. Set the layout to 2x2
# mar = c(bottom, left, top, right) - we reduce these to bring plots closer
par(mfrow = c(1,1), xpd = TRUE)

# 2. Plot with increased vertex size (rel)
# The 'vertex.weight' or 'vertex.label.cex' can sometimes push the circle size.
# More importantly, netVisual_circle uses an internal 'canvas' size.

# BL - Control
netVisual_circle(control_weights, weight.scale = T, 
                 label.edge = F, title.name = "BL - Control",
                 vertex.label.cex = 0.8)


# BL - Dengue
netVisual_circle(dengue_weights, weight.scale = T, 
                 label.edge = F, title.name = "BL - Dengue",
                 vertex.label.cex = 0.8)



# DOM - Control
netVisual_circle(control_weights_dom, weight.scale = T, 
                 label.edge = F, title.name = "DOM - Control",
                 vertex.label.cex = 0.8)


# DOM - Dengue
netVisual_circle(dengue_weights_dom, weight.scale = T, 
                 label.edge = F, title.name = "DOM - Dengue",
                 vertex.label.cex = 0.8)

dev.off()

# -------- PLOT 2b --------
# Identify significant interaction IDs from BOTH methods
sig_ids <- final_comparison %>%
  filter(pval < 0.05) %>%
  pull(interaction_id) %>%
  unique()

# Filter the master stats table for these IDs
plot_data_unified <- final_comparison %>%
  filter(interaction_id %in% sig_ids)

plot_data_unified <- plot_data_unified %>%
  mutate(
    complex_label = paste0("[", pathway_name, "] ", interaction_id)
  )

# Pivot to Wide Format to compare Method A vs Method B row-by-row
pathway_logic <- plot_data_unified %>%
  select(interaction_id, pathway_name, method, log2FC) %>%
  
  pivot_wider(
    names_from = method, 
    values_from = log2FC,
    names_prefix = "logFC_" 
  ) %>%
  mutate(
    # --- CONDITION 1: EXISTENCE ---
    # Check if logFC is NA in one column but exists in the other.
    exists_in_baseline = !is.na(logFC_Baseline),
    exists_in_domain   = !is.na(logFC_Domain),
    is_missing_in_one  = exists_in_baseline != exists_in_domain,
    
    # --- CONDITION 2: DIRECTION FLIP ---
    # Only check this if it exists in BOTH 
    # Check if one is positive (+) and the other is negative (-)
    direction_flip = (exists_in_baseline & exists_in_domain) & 
      (sign(logFC_Baseline) != sign(logFC_Domain))
  )

# Aggregate by Pathway
# If a pathway has AT LEAST ONE interaction that matches either condition, PLOT IT.
target_pathways <- pathway_logic %>%
  group_by(pathway_name) %>%
  summarise(
    keep_pathway = any(is_missing_in_one | direction_flip),
    .groups = "drop"
  ) %>%
  filter(keep_pathway == TRUE) %>%
  pull(pathway_name)

# Create Final Plotting Data
# Filter the original dataset to include ONLY these pathways
plot_data_final <- df_check %>%
  filter(pathway_name %in% target_pathways) %>%
  mutate(neg_log_p = -log10(pval + 1e-10)) %>%
  mutate(method_label = if_else(method == "Baseline", "Baseline Stage (BL)", "Domain-aware Stage (DOM)"))


# PREPARE LABELS (Same as before)
plot_data_labels <- plot_data_final %>%
  filter(method == "Baseline") %>%
  separate(interaction_id, into = c("cell_pair", "lr_pair"), sep = ": ", remove = FALSE)

# CREATE THE "TABLE" PLOT (The Labels)
p_labels <- ggplot(plot_data_labels, aes(y = interaction_id)) +
  # Column 1: Cell Pair (at x = 0)
  geom_text(aes(x = 0, label = cell_pair), hjust = 0, size = 3.5, fontface = "bold") +
  # Column 2: LR Pair (at x = 1)
  geom_text(aes(x = 1, label = lr_pair), hjust = 0, size = 3.5) +
  
  scale_x_continuous(limits = c(0, 2)) +
  
  facet_grid(pathway_name ~ ., scales = "free_y", space = "free_y") +
  
  theme_void() + 
  theme(
    strip.text = element_blank(),
    panel.spacing.y = unit(5.5, "pt")
  )

# MODIFY YOUR MAIN PLOT
p_main <- ggplot(plot_data_final, aes(x = log2FC, y = interaction_id)) +
  geom_point(aes(color = method, size = neg_log_p), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  facet_grid(pathway_name ~ method_label, scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("Baseline" = "#E41A1C", "Domain" = "#377EB8")) +
  guides(color = "none") +
  scale_size_continuous(range = c(0.5, 3.5)) +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8, hjust = 0),
    strip.text.x = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95"),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "Differential Cell-Cell Communication: Severe Dengue vs Control",
    x = "Log2 Fold Change (Interaction Strength)",
    size = "-log10(p-value)"
  )

final_plot <- p_labels + p_main + plot_layout(widths = c(1, 3))
