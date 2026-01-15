library(zellkonverter)
library(SingleCellExperiment)
library(scater)
library(scran)
library(biomaRt)
library(CellChat)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(umap)
library(EnsDb.Hsapiens.v86)
library(GenomicFeatures)
library(ggplot2)
library(clusterProfiler)
library(stringr)
library(purrr)
library(glue)
library(RColorBrewer)
library(patchwork)
library(forcats)
library(ggh4x)

augment_cofactors <- function(sce_pc, sce_pc_loc, db, species = "Human", assay_name = "counts") {
  rowData(sce_pc) <- NULL
  rowData(sce_pc_loc) <- NULL
  
  # Validate inputs
  if (!assay_name %in% assayNames(sce_pc)) {
    stop(paste("Assay", assay_name, "not found in sce_pc object"))
  }
  
  # Get the counts matrix (aggregation should be on raw counts)
  counts_matrix <- assay(sce_pc, assay_name)
  
  # Connect to the Ensembl BioMart
  message("Connecting to Ensembl BioMart...")
  ensembl_mart <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "https://useast.ensembl.org"
  )
  
  # Extract unique cofactor genes
  cofactor_genes <- db$cofactor %>% 
    unlist() %>% 
    unique()
  
  genes_to_pool <- cofactor_genes[!is.na(cofactor_genes) & cofactor_genes != ""]
  
  if (length(genes_to_pool) == 0) {
    warning("No valid cofactor genes found in db$cofactor")
    return(sce_pc_loc)
  }
  
  message(paste("Querying BioMart for", length(genes_to_pool), "cofactor genes..."))
  
  # Get transcript-to-gene mappings in one batch query
  tx_map <- getBM(
    attributes = c('ensembl_transcript_id', 'hgnc_symbol'),
    filters = 'hgnc_symbol',
    values = genes_to_pool,
    mart = ensembl_mart
  )
  
  message(paste("Found", nrow(tx_map), "transcript-to-gene mappings"))
  
  # Aggregate transcripts for each gene
  pooled_rows <- list()
  genes_found <- 0
  genes_skipped <- 0
  
  message("Aggregating transcript counts to gene level...")
  
  for (gene in genes_to_pool) {
    # Find transcripts for this gene from BioMart mapping
    tx_for_gene <- tx_map$ensembl_transcript_id[tx_map$hgnc_symbol == gene]
    
    # Find which transcripts exist in the SCE object
    tx_in_sce <- intersect(tx_for_gene, rownames(counts_matrix))
    
    if (length(tx_in_sce) > 0) {
      # Extract counts for these transcripts
      tx_counts <- counts_matrix[tx_in_sce, , drop = FALSE]
      
      # Sum across transcripts - handle both single and multiple transcripts
      if (length(tx_in_sce) == 1) {
        # Single transcript: extract as vector
        pooled_expression <- as.vector(tx_counts)
      } else {
        # Multiple transcripts: sum them
        pooled_expression <- colSums(tx_counts)
      }
      
      pooled_rows[[gene]] <- pooled_expression
      genes_found <- genes_found + 1
      
      message(paste("  ✓", gene, ":", length(tx_in_sce), "transcript(s) aggregated"))
    } else {
      genes_skipped <- genes_skipped + 1
      warning(paste("  ✗", gene, ": No transcripts found in SCE object"))
    }
  }
  
  message(paste("\nSummary: Found", genes_found, "genes, skipped", genes_skipped))
  
  # Return unchanged object if no genes were pooled
  if (length(pooled_rows) == 0) {
    warning("No cofactor genes were successfully pooled. Returning original object.")
    return(sce_pc_loc)
  }
  
  # Create matrix of pooled gene-level counts
  new_rows_matrix <- do.call(rbind, pooled_rows)
  
  # --- FIX STARTS HERE ---
  
  # 1. Check if the original matrix is a DelayedArray
  # We use methods::is() to check inheritance safely
  is_delayed <- methods::is(counts_matrix, "DelayedArray")
  
  # 2. Match matrix structure
  if (is_delayed) {
    # If original is DelayedArray, the new one MUST be too
    requireNamespace("DelayedArray")
    new_rows_matrix <- DelayedArray::DelayedArray(new_rows_matrix)
  } else if (methods::is(counts_matrix, "sparseMatrix")) {
    # If original is sparse (and not Delayed), make new one sparse
    new_rows_matrix <- as(new_rows_matrix, "sparseMatrix")
  }
  
  # Match matrix type (sparse vs dense)
  if (is(counts_matrix, "sparseMatrix")) {
    new_rows_matrix <- as(new_rows_matrix, "sparseMatrix")
  }
  
  # Create temporary SCE for new rows
  sce_new_rows <- SingleCellExperiment(
    assays = stats::setNames(list(new_rows_matrix), assay_name),
    colData = colData(sce_pc_loc)
  )
  
  # Combine original and new rows
  message(paste("Combining:", nrow(sce_pc_loc), "original +", 
                nrow(sce_new_rows), "new features"))
  print(sce_pc_loc)
  print(sce_new_rows)
  sce_combined <- rbind(sce_pc_loc, sce_new_rows)
  
  message(paste("✓ Final SCE object:", nrow(sce_combined), "features"))
  message("⚠ WARNING: Re-normalize the object before downstream analysis!")
  
  return(sce_combined)
}

run_cellchat_isoform_level <- function(sce, cell_type_column, assay_type = "counts") { 
  
  # --- STEP 1: PREPARE DATA --- 
  message("Extracting data...")
  data.input <- assay(sce, assay_type)
  
  # Explicitly convert to sparse matrix (dgCMatrix)
  if (!inherits(data.input, "dgCMatrix")) {
    message("Converting data to dgCMatrix format...")
    data.input <- as(as.matrix(data.input), "dgCMatrix")
  }
  
  # --- STEP 2: LOG-NORMALIZE IF NEEDED --- 
  if (max(data.input) > 50) { 
    message("Data appears to be raw/TPM (Max value > 50). Applying Log-Normalization...") 
    data.input <- log1p(data.input)  
  } 
  
  # --- STEP 3: CREATE CELLCHAT OBJECT --- 
  meta.data <- as.data.frame(colData(sce)) 
  
  if (!cell_type_column %in% colnames(meta.data)) {
    stop(paste("Column", cell_type_column, "not found in metadata!"))
  }
  
  message("Creating CellChat object...")
  cellchat <- createCellChat(object = data.input,  
                             meta = meta.data,  
                             group.by = cell_type_column) 
  
  # --- STEP 4: SET DATABASE & RUN --- 
  if (!exists("myIsoformCellChatDB")) {
    stop("myIsoformCellChatDB not found in environment.")
  }
  cellchat@DB <- myIsoformCellChatDB 
  
  message("Running CellChat workflow...")
  cellchat <- subsetData(cellchat) 
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)  
  
  # --- FIX FOR "droplevels" ERROR ---
  # Convert the column to a factor first. This automatically sets the levels
  # based on the data present, so explicit 'droplevels' isn't strictly needed,
  # but we do it safely just in case.
  current_labels <- as.character(cellchat@meta[[cell_type_column]])
  cellchat@meta[[cell_type_column]] <- factor(current_labels)
  
  # Set identities
  cellchat@idents <- cellchat@meta[[cell_type_column]]
  
  cellchat <- identifyOverExpressedInteractions(cellchat) 
  
  # Calculate probabilities
  # Note: You can adjust trim=0.05 or 0.1 depending on sparsity
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05) 
  cellchat <- computeCommunProbPathway(cellchat) 
  cellchat <- aggregateNet(cellchat) 
  
  print("CellChat analysis complete!") 
  
  # --- STEP 5: EXTRACT RESULTS ---
  interactions <- subsetCommunication(cellchat) 
  
  if (is.null(cellchat@net$count)) {
    interaction_counts <- data.frame()
    warning("No interactions detected! Check your DB or expression levels.")
  } else {
    interaction_counts <- as.data.frame(as.table(cellchat@net$count)) 
    colnames(interaction_counts) <- c("Sender", "Receiver", "Count") 
  }
  
  return(list(cellchat = cellchat, interactions = interactions, interaction_counts = interaction_counts)) 
}

heatmap_cellchat_comparison <- function(counts_before, counts_after) {
  interaction_comparison <- full_join(counts_before, 
                                      counts_after, 
                                      by = c("Sender", "Receiver"),
                                      suffix = c("_before", "_after"))
  
  interaction_comparison <- interaction_comparison |>
    mutate(
      Count_before_pseudo = Count_before + 1,
      Count_after_pseudo = Count_after + 1
    ) |> 
    mutate(logFC = log2(Count_after_pseudo / Count_before_pseudo))
  
  p <- ggplot(interaction_comparison, aes(x = Sender, y = Receiver, fill = logFC)) +
    geom_tile(color = "white") +
    # This line ensures the blocks are square
    coord_fixed() + 
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-2.50, 2.50),
      name = "log2 Fold Change"
    ) +
    theme_minimal() +
    labs(
      title = "Log Fold Change in Cell-Cell Interactions",
      x = "Sender Cell Type",
      y = "Receiver Cell Type"
    ) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
  
  return(list(plot = p, data = interaction_comparison))
}

analyze_cell_filtering_impact <- function(data, sender = NULL, receiver = NULL, top_n = 5, focus_pathway = NULL, max_pairs = 8) {
  
  # --- 1. Setup & Logic Handling ---
  # Ensure stages are ordered
  data$stage <- factor(data$stage, levels = c("BL", "PC", "LOC", "DOM"))
  
  filter_logic <- function(column, pattern) {
    str_detect(column, regex(paste0("\\b", pattern, "\\b"), ignore_case = TRUE))
  }
  
  if (!is.null(sender) && !is.null(receiver)) {
    filtered_data <- data %>% 
      filter(filter_logic(source, sender) & filter_logic(target, receiver))
    title_prefix <- paste0(sender, " -> ", receiver)
  } else if (!is.null(sender)) {
    filtered_data <- data %>% filter(filter_logic(source, sender))
    title_prefix <- paste0("Outgoing from ", sender)
  } else if (!is.null(receiver)) {
    filtered_data <- data %>% filter(filter_logic(target, receiver))
    title_prefix <- paste0("Incoming to ", receiver)
  } else {
    stop("Please provide at least a sender or a receiver cell type.")
  }
  
  if (nrow(filtered_data) == 0) stop("No interactions found for exact cell type match.")
  
  # --- 2. Pathway Analysis ---
  pathway_stats <- filtered_data %>%
    group_by(stage, pathway_name) %>%
    summarise(total_prob = sum(prob), .groups = "drop") %>%
    pivot_wider(names_from = stage, values_from = total_prob, values_fill = 0)
  
  # Ensure BL/DOM columns exist
  for(col in c("BL", "DOM")) if(!col %in% names(pathway_stats)) pathway_stats[[col]] <- 0
  
  pathway_impact <- pathway_stats %>%
    mutate(delta = DOM - BL, abs_change = abs(DOM - BL)) %>%
    arrange(desc(abs_change)) 
  
  pathway_impact_top <- pathway_impact %>% slice_head(n = top_n)
  top_pathways <- pathway_impact_top$pathway_name
  
  # --- 3. Selection Logic ---
  if (is.null(focus_pathway)) {
    target_pathways <- top_pathways[1]
  } else if (length(focus_pathway) == 1 && focus_pathway == "all") {
    target_pathways <- top_pathways
  } else {
    target_pathways <- focus_pathway
  }
  
  # --- 4. Heatmap Data Prep ---
  heatmap_data <- filtered_data %>%
    filter(pathway_name %in% target_pathways) %>%
    mutate(lr_pair = paste0(ligand, "-", receptor)) %>%
    group_by(stage, pathway_name, lr_pair) %>%
    summarise(prob = sum(prob), .groups = "drop")
  
  top_labels <- heatmap_data %>%
    group_by(pathway_name, lr_pair) %>%
    summarise(max_p = max(prob), .groups = "drop") %>%
    group_by(pathway_name) %>%
    slice_max(max_p, n = max_pairs, with_ties = FALSE) %>%
    pull(lr_pair)
  
  heatmap_data <- heatmap_data %>% 
    filter(lr_pair %in% top_labels) %>%
    # This line ensures empty/missing combinations appear as NA
    complete(stage, nesting(pathway_name, lr_pair))
  
  # --- 5. Visualization ---
  p_heatmap <- ggplot(heatmap_data, aes(x = stage, y = lr_pair, fill = prob)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_grid(pathway_name ~ ., scales = "free_y", space = "free_y", switch = "y") +
    
    # Custom Gradient: Black (Low) -> Red (High), Grey (NA/NS)
    scale_fill_gradientn(
      # 1. Define colors: Light Pink -> Vibrant Red -> Dark Red/Black
      colors = c("#FFD9D9", "#FF0000", "#4A0000"), 
      
      # 2. Define where these colors happen (normalized 0 to 1 based on limits)
      # 0.0 = Start (0.0 interaction) -> Light Pink
      # 0.33 = One-third up (approx 0.1 interaction) -> Vibrant Red
      # 1.0 = End (0.3 interaction) -> Dark Red
      values = c(0, 0.33, 1), 
      
      # 3. Set the absolute limits of your data
      limits = c(0, 0.3), 
      
      # 4. Handle values above 0.3 safely
      oob = scales::squish, 
      
      # 5. Legend label
      name = "Interaction Score",
      
      # 6. Ensure NA (insignificant) remains Grey
      na.value = "grey96"
    ) +
    
    theme_minimal() +
    theme(
      plot.title.position = "plot", 
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, face = "bold"),
      panel.spacing.y = unit(0.3, "lines"),
      strip.background = element_rect(fill = "grey96", color = "grey90"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(hjust = 0)
    ) +
    force_panelsizes(cols = unit(6, "cm"), respect = TRUE) +
    labs(
      title = "Interaction Intensity Heatmap",
      subtitle = paste0("Focus: ", title_prefix),
      x = "Filtering Stage",
      y = "Pathway | L-R Pair",
      caption = "Grey tiles indicate no significant interaction detected."
    )
  
  # --- 6. Return List ---
  return(list(report = pathway_impact, plot_heatmap = p_heatmap))
}

calc_diff_stats <- function(df) {
  # 1. Create a mapping of interaction_id to pathway_name before pivoting
  pathway_map <- df %>%
    mutate(interaction_id = paste0(source, " → ", target, ": ", ligand, "-", receptor)) %>%
    select(interaction_id, pathway_name) %>%
    distinct()
  
  # 2. Create Matrix and Pivot
  mat <- df %>%
    mutate(interaction_id = paste0(source, " → ", target, ": ", ligand, "-", receptor)) %>%
    select(interaction_id, sample, prob) %>%
    pivot_wider(names_from = sample, values_from = prob, values_fill = 0) %>%
    column_to_rownames("interaction_id")
  
  d_samps <- colnames(mat)[grepl("1-0", colnames(mat))]
  c_samps <- colnames(mat)[grepl("3-0", colnames(mat))]
  
  # 3. Stats Calculation
  res <- data.frame(
    interaction_id = rownames(mat),
    mean_dengue = rowMeans(mat[, d_samps, drop=FALSE]),
    mean_control = rowMeans(mat[, c_samps, drop=FALSE]),
    log2FC = log2((rowMeans(mat[, d_samps]) + 1e-6) / (rowMeans(mat[, c_samps]) + 1e-6)),
    pval = apply(mat, 1, function(x) {
      if(sum(x) == 0) return(1)
      tryCatch(wilcox.test(x[d_samps], x[c_samps])$p.value, error = function(e) 1)
    }),
    stringsAsFactors = FALSE
  ) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    # Re-attach the pathway name
    left_join(pathway_map, by = "interaction_id")
  
  return(res)
}