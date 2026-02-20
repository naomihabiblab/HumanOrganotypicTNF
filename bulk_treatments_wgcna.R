library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(viridis)
library(MASS)

# =============================================================================
# 1. LDA Analysis (Panel A)
# =============================================================================

lda_scores<-read.csv("LDA_treatments.csv")
treatments <- c("Control", "BLM", "H2O2", "LPS", "TNF")

lda_df <- data.frame(
  LD1 = lda_scores[,1],
  LD2 = lda_scores[,2],
  Treatment = meta_data$Treatment
)

p_lda <- ggplot(lda_df, aes(x = LD1, y = LD2, color = Treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c(
    "Control" = "#0072B2",     
    "BLM" = "#D55E00",       
    "H2O2" = "#009E73",       
    "LPS" = "#E69F00",        
    "TNF" = "#CC79A7"        
  )) +
  stat_ellipse(level = 0.68, alpha = 0.3) +
  theme_classic(base_size = 14) +
  labs(x = "LD1", y = "LD2") +
  theme(legend.position = "right", legend.title = element_blank())

# =============================================================================
# 2. WGCNA Program-Trait Heatmap (Panel B) 
# =============================================================================

moduleTraitCor <- read.csv("WGCNA_program_trait_correlations.csv")
moduleTraitPvalue <- read.csv("WGCNA_program_trait_p_values.csv")

rownames(moduleTraitCor) <- moduleTraitCor$X
moduleTraitCor <- moduleTraitCor[,-1]
rownames(moduleTraitPvalue) <- moduleTraitPvalue$X
moduleTraitPvalue <- moduleTraitPvalue[,-1]

selected_cols <- c("Control", "BLM", "H2O2", "LPS", "TNF")
selected_rows <- c("Prog.1","Prog.2","Prog.3","Prog.4","Prog.5") 
moduleTraitPvalue <- moduleTraitPvalue[selected_rows,selected_cols]
moduleTraitCor <- moduleTraitCor[selected_rows,selected_cols]

stars_matrix <- matrix("", nrow = nrow(moduleTraitPvalue), ncol = ncol(moduleTraitPvalue))
stars_matrix[moduleTraitPvalue < 0.001] <- "***"
stars_matrix[moduleTraitPvalue < 0.01 & moduleTraitPvalue >= 0.001] <- "**"
stars_matrix[moduleTraitPvalue < 0.05 & moduleTraitPvalue >= 0.01] <- "*"
rownames(stars_matrix) <- rownames(moduleTraitPvalue)
colnames(stars_matrix) <- colnames(moduleTraitPvalue)

col_fun <- colorRamp2(c(-0.75, 0, 0.75), c("navy", "white", "darkred"))

ht <- Heatmap(
  as.matrix(moduleTraitCor),
  name = "Correlation",
  col = col_fun,
  heatmap_legend_param = list(at = c(-0.75, 0, 0.75), labels = c("-0.75", "0", "0.75")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (stars_matrix[i, j] != "") {
      grid.text(stars_matrix[i, j], x, y, gp = gpar(fontsize = 14, fontface = "bold"))
    }
  },
  show_row_names = TRUE, show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  show_heatmap_legend = FALSE, border = TRUE, border_gp = gpar(col = "black", lwd = 1)
)

lgd_sig <- Legend(labels = c("* <0.05", "** <0.01", "*** <0.001"), type = "points", pch = "*", 
                  title = "Significance", title_gp = gpar(fontsize = 12, fontface = "bold"))
lgd_col <- Legend(col_fun = col_fun, title = "Correlation", at = c(-0.75, 0, 0.75), 
                  labels = c("-0.75", "0", "0.75"), title_gp = gpar(fontsize = 12, fontface = "bold"))

draw(ht, annotation_legend_list = list(lgd_col, lgd_sig))

# =============================================================================
# 3. PATHWAY DOTPLOT (PANEL C)
# =============================================================================
pathway_dir <- "."  # UPDATE TO YOUR PATHWAY FOLDER
module_files <- list.files(pathway_dir, pattern = "*.csv", full.names = TRUE)
module_names <- gsub(".csv$", "", basename(module_files))

combined_df <- bind_rows(
  lapply(seq_along(module_files), function(i) {
    df <- read_csv(module_files[i])
    df$Program <- module_names[i]
    df
  })
)

# selected pathways from enrichment analysis per program
selected_pathways <- list(
  Prog.1 = c("positive regulation of apoptotic process","p53 signaling pathway", 
             "Apoptosis", "Cell cycle","Interleukin-10 signaling",
             "NF-kappa B signaling pathway","regulation of JNK cascade",
             "cytokine-mediated signaling pathway"),
  `Prog.2` = c("collagen metabolic process","Fibronectin matrix formation", 
               "positive regulation of cell migration","response to wounding",
               "regulation of cellular response to growth factor stimulus"),
  `Prog.3` = c("response to lipopolysaccharide", "TNF signaling pathway", 
               "chemotaxis","cell activation","NF-kappa B signaling pathway",
               "response to steroid hormone","cell-cell adhesion"),
  `Prog.4` = c("Interferon Signaling", "cellular response to cytokine stimulus", 
               "cytolysis in another organism"),
  `Prog.5` = c("inflammatory response", "positive regulation of cell activation", 
               "response to copper ion")
)

module_order <- c("Prog.1", "Prog.2", "Prog.3", "Prog.4", "Prog.5")

filtered_df <- combined_df %>%
  rowwise() %>%
  filter(Description %in% selected_pathways[[Program]])


description_order_df <- filtered_df %>%
  group_by(Program, Description) %>%
  summarise(max_LogP = max(-LogP, na.rm = TRUE), .groups = "drop") %>%
  arrange(factor(Program, levels = module_order), desc(max_LogP))
ordered_descriptions <- description_order_df %>% pull(Description) %>% unique()

filtered_df_ordered <- filtered_df %>%
  ungroup() %>%
  mutate(
    Program = factor(Program, levels = module_order),
    Description = factor(Description, levels = ordered_descriptions)
  )

# plotting
p_dot <- ggplot(filtered_df_ordered, aes(x = Program, y = Description)) +
  geom_point(aes(color = -LogP, size = Count), alpha = 0.8) +
  scale_color_gradient(low = "white", high = "darkred", limits = c(0, 15), 
                       name = expression(-log[10](p~value))) +
  scale_size_continuous(range = c(3, 10), name = "Gene Count") +
  scale_y_discrete(limits = rev(levels(filtered_df_ordered$Description))) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.title = element_blank(),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(p_lda)
print(p_dot)


