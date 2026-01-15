library(ggpubr)
library(ggplot2)
library(scales)

# 1. Define a helper function to extract data
get_counts_data <- function(sce_obj, step_name, cell_type_col = "cell_type") {
  total_counts <- colSums(assay(sce_obj, "counts"))
  
  # Extract cell type annotations
  cell_types <- colData(sce_obj)[[cell_type_col]]
  
  data.frame(
    Transcript_Counts = total_counts,
    Cell_Type = cell_types,
    Stage = step_name
  )
}

# 2. Extract data from your 3 objects
df1 <- get_counts_data(results_list$sce_baseline, 
                       "BL", 
                       cell_type_col = "cell_type_refined")
df2 <- get_counts_data(results_list$sce_pc, 
                       "PC", 
                       cell_type_col = "cell_type_refined")
df3 <- get_counts_data(results_list$sce_deeploc, 
                       "LOC", 
                       cell_type_col = "cell_type_refined")

# 3. Combine into one large data frame
plot_data <- rbind(df1, df2, df3)

plot_data$Stage <- factor(plot_data$Stage, levels = c("BL", "PC", "LOC"))

plot_data <- plot_data %>%
  filter(Cell_Type != "platelet")

my_comparisons <- list( 
  c("BL", "PC"),
  c("BL", "LOC")
)

p <- ggplot(plot_data, aes(x = Stage, y = Transcript_Counts, fill = Stage)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.6) +
  
  geom_jitter(width = 0.2, alpha = 0.25, size = 0.8) +
  
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
  
  facet_wrap(~Cell_Type, scales = "free_y") +
  
  scale_y_log10(
    labels = scales::comma, 
    expand = expansion(mult = c(0.05, 0.15)) 
  ) +
  
  scale_fill_viridis_d(option = "viridis", begin = 0.2, end = 0.8) +
  
  # --- STATISTICAL TEST LAYER ---
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.05, 1), 
      symbols = c("*", "ns")
    )
  ) +
  # ------------------------------

theme_bw() +
  labs(
    title = "Smart-seq2 Library Size per Cell Type",
    y = "Total Mapped Reads (Log Scale)",
    x = "Filtering Stage",
    caption = "ns: p > 0.05; *: p <= 0.05"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"), # Left-align and style it
    
    # Title settings
    plot.title = element_text(size = 12, face = "bold", hjust = 0))

# Summarize total library size per stage
total_loss_summary <- plot_data %>%
  group_by(Stage) %>%
  summarise(
    Total_Reads = sum(Transcript_Counts),
    Mean_Reads_Per_Cell = mean(Transcript_Counts),
    Cell_Count = n()
  ) %>%
  mutate(Percent_Remaining = (Total_Reads / first(Total_Reads)) * 100)


# Create the global loss plot
p_loss <- ggplot(total_loss_summary, aes(x = Stage, y = Total_Reads, fill = Stage)) +
  geom_col(color = "black", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = comma(Total_Reads)), vjust = -0.5, fontface = "bold") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_viridis_d(option = "viridis", begin = 0.2, end = 0.8) +
  theme_classic() +
  labs(
    title = "Total Read Loss Across Filtering Stages",
    subtitle = "Global sample view (all cell types combined)",
    y = "Aggregate Sum of Mapped Reads",
    x = "Filtering Stage"
  )

ggsave("results/global-read-loss.png", plot = p_loss, width = 8, height = 6, dpi = 400)