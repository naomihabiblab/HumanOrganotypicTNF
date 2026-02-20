
# =============================================================================
# Single-Nucleus RNA-seq Visualization
# =============================================================================

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(pheatmap)
library(ComplexHeatmap)
library(gridExtra)
library(viridis)
library(Nebulosa)

# ---- CONFIGURATION ----
# Update these paths to your data
f_tnf <- "..."                    # Main TNF object folder
f_tnf_de <- "data/de_results/"          # DE results folder
out_tnf_f <- "figures/"                 # Output figures folder

# Create output directory
dir.create(out_tnf_f, recursive = TRUE, showWarnings = FALSE)

# Color palettes (from your script)
new_colors_celltype <- c("#1F78B4","#E69F00","#B2DF8A","#33A02C","#A6CEE3",
                         "#CC6677","#D694B9","#FDBF6F","#AA4499")
colors_tnf_nb_umap <- c("Steel Blue", "Sandy Brown")

# Celltype markers 
Molecular_markers <- list(
  Astrocytes = c("SLC1A2","SLC1A3","ALDH1L1","SLC4A4","AQP4","CD44","GFAP"),
  Endothelial = c("FLT1", "CLDN5","ABCB1","ATP10A"),
  Ex.neurons = c("SNAP25","CAMK2A","SLC17A7", "CUX2","RBFOX3"),
  In.neurons = c("GAD1","GAD2","VIP","SST","ADARB2","LAMP5","LHX6","PVALB","CHRNA7"),
  Microglia = c("TMEM119","P2RY12","CX3CR1","PTPRC","TREM2","ITGAX"),
  Oligos = c("MBP", "MOG", "MAG"),
  OPCs = c("PCDH15", "MEGF11", "VCAN","PDGFRA"),
  Perivascular_cells = c("PDGFRB","DCN","FN1","COL3A1","COL5A2","COL6A1","CEMIP","FBLN1")
)
gene_list <-c("NAMPT","SOD2","WTAP","IRAK1","NFKB1","PARP14","TNFAIP2","TNIP1")

# =============================================================================
# FIGURE 1: Cell Types Analysis
# =============================================================================

#' Figure 1A: UMAP by Cell Types
figure1a_umap_celltype <- function(obj) {
  p <- DimPlot(obj, group.by = "CellType", pt.size = 0.5, cols = new_colors_celltype) +
    ggtitle("Cell Types UMAP") +
    theme_classic(base_size = 14) +
    theme(legend.position = "right")
  return(p)
}

#' Figure 1C: Cell Type Markers DotPlot
figure1b_dotplot_markers <- function(obj) {
  p <- DotPlot(obj, features = unlist(Molecular_markers), 
               group.by = "CellType", cols = c("yellow", "blue")) +
    RotatedAxis() + ggtitle("Cell Types") +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
  return(p + theme_classic())
}

#' Figure 1D: Cell Type Proportions by Sample
figure1c_proportions <- function(obj) {
  prop_data <- prop.table(table(obj$orig.ident, obj$CellType), margin = 1)
  prop_df <- as.data.frame(prop_data)
  colnames(prop_df) <- c("Sample", "CellType", "Proportion")
  
  p <- ggplot(prop_df, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = new_colors_celltype) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(y = "Proportion", x = "Sample", title = "Cell Type Proportions")
  return(p)
}

# =============================================================================
# FIGURE 2: Treatment Analysis (TNF vs Control)
# =============================================================================

#' Figure 2A: UMAP by Treatment (TNF vs Control)
figure2a_umap_treatment <- function(tnf) {
  p <- DimPlot(tnf, group.by = "treatment", pt.size = 0.5, 
               cols = colors_tnf_nb_umap) +
    ggtitle("TNF vs Control") +
    theme_classic(base_size = 14)
  return(p)
}

#' Figure 2B: NFKB1 Density Plot
figure2b_nfkb1_density <- function(tnf) {
  plot_density(tnf, features = "NFKB1", joint = FALSE) +
    ggtitle("NFKB1") +
    theme_classic(base_size = 14) +
    scale_color_viridis_c(option = "magma")
}

#' Figure 2C: Shared Main Genes Heatmap 
figure2c_shared_heatmap <- function(tnf, des.tnf.vs.ctr) {
  # shared genes
  shared_genes <- intersect(gene_list, rownames(tnf))

  # Average expression matrix
  avg_exp <- AverageExpression(tnf, features = shared_genes, 
                               group.by = "CellType")$RNA
  heatmap_data <- t(scale(t(avg_exp)))
  
  # Order cell types
  celltype_order <- c("Ex.neurons", "In.neurons", "Astrocytes", "Microglia", 
                      "OPCs", "Oligos", "Endothelial", "Perivascular_cells")
  celltype_order <- celltype_order[celltype_order %in% colnames(heatmap_data)]
  heatmap_data <- heatmap_data[, celltype_order]
  
  # Create annotation for cell types
  annotation_col <- data.frame(
    CellType = factor(colnames(heatmap_data), levels = celltype_order),
    row.names = colnames(heatmap_data)
  )
  annotation_colors <- list(
    CellType = setNames(rainbow(length(unique(annotation_col$CellType))), 
                        unique(annotation_col$CellType))
  )
  
# plot
  p <- pheatmap(heatmap_data,
                scale = "none", 
                cluster_rows = TRUE,
                cluster_cols = FALSE,  
                show_rownames = TRUE,
                show_colnames = TRUE,
                fontsize_row = 10,
                fontsize_col = 11,
                angle_col = 0,
                main = "C. Shared TNF Response Genes",
                annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                color = colorRampPalette(c("white", "navy"))(100),
                border_color = "white",
                cellwidth = 25,
                cellheight = 12,
                treeheight_col = 0,
                legend = TRUE)
  
  return(p)
}

#' Figure 2D-F: Volcano Plots (Astrocytes, Microglia, OPCs)
figure2d_volcano <- function(des.tnf.vs.ctr) {
  glia_types <- c("microglia", "astrocytes", "opcs")
  plots <- list()
  
  for (i in seq_along(glia_types)) {
    ctype <- glia_types[i]
    data <- des.tnf.vs.ctr[des.tnf.vs.ctr$sample == ctype, ]
    data$sig <- ifelse(data$p_val_adj < 0.05, 
                       ifelse(data$avg_log2FC > 0.2, "Upregulated", "Downregulated"), 
                       "Non-significant")
    
    p <- ggplot(data, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                          color = sig, size = abs(avg_log2FC))) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("Upregulated" = "Sandy Brown", "Non-significant" = "black","Downregulated"="Steel Blue")) +
      scale_size_continuous(range = c(0.5, 2)) +
      theme_classic(base_size = 11) +
      labs(x = "Log2FC", y = "-Log10(p-adj)") +
      theme(legend.position = "none")
    
    if (i == 1) p <- p + ggtitle("D") + theme(plot.title = element_text(hjust = 0.5))
    plots[[i]] <- p
  }
  return(plots)
}

#' Figure 2G-I: Pathway Barplots (Astrocytes, Microglia, OPCs)
figure2e_pathway_barplots <- function(fs_tnf_de_pathway, paths_per_type) {
  glia_types <- c("astrocytes", "microglia", "opcs")
  plots <- list()
  
  for (obj_n in names(fs_tnf_de_pathway)) {
    if (any(grepl(glia_types, obj_n))) {
      path <- read.csv(file.path(f_tnf_de, fs_tnf_de_pathway[[obj_n]]))
      
            # Selected pathways
      path_list <- intersect(path$Description, paths_per_type[[obj_n]])
      path_subset <- path[path$Description %in% path_list, ]
      
      p <- ggplot(path_subset, aes(x = reorder(Description, -log10(p.adjust)), 
                                   y = GeneRatioFrac, fill = p.adjust)) +
        geom_bar(stat = "identity") + coord_flip() +
        theme_classic() + ggtitle(obj_n) +
        labs(y = "-log10(p.adjusted)", x = "", fill = "p.adjusted")
      plots[[obj_n]] <- p
    }
  }
  return(plots)
}

# =============================================================================
# MAIN PIPELINE - Generate All Figures
# =============================================================================

main_snuc_pipeline <- function() {
  cat("Loading SNUC objects...\n")
  
  # Load main objects
  obj <- readRDS(file.path(f_tnf, 'TNF.rds'))

  cat("Generating Figure 1...\n")
  # Figure 1
  p1a <- figure1a_umap_celltype(obj)
  p1b <- figure1b_dotplot_markers(obj)
  p1c <- figure1c_proportions(obj)
  
  # Save Figure 1
  combined_fig1 <- plot_grid(p1a, p1b, p1c, ncol = 3, labels = c("A", "B", "C"))
  ggsave(file.path(out_tnf_f, "Figure1_SNUC_CellTypes.pdf"), combined_fig1, 
         width = 18, height = 6, dpi = 300)
  
  cat("Generating Figure 2...\n")
  # Figure 2  
  p2a <- figure2a_umap_treatment(tnf)
  p2b <- figure2b_nfkb1_density(tnf)
  p2c <- figure2c_shared_heatmap(tnf, des.tnf.vs.ctr)
  
  # Volcano and pathway plots
  volcano_plots <- figure2d_volcano_plots(des.tnf.vs.ctr)
  pathway_plots <- figure2e_pathway_barplots(fs_tnf_de_pathway, paths_per_type)
  
  cat("all figures saved to", out_tnf_f)

}

# =============================================================================
# RUN PIPELINE
# =============================================================================
if (!interactive()) {
  main_snuc_pipeline()
} else {
  cat("SNUC Analysis Pipeline Loaded!\n")
  cat("Run: main_snuc_pipeline()\n")
}
