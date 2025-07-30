# Load necessary libraries
library(circlize)
library(dplyr)
library(RColorBrewer)

# Output directory
base_dir <- "C:/Users/lenov/Downloads/GSEDATASETS"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
setwd(base_dir)

gene_pathway <- tribble(
  ~gene,      ~pathway,
  "CXCL17",   "hsa04060 Cytokine-cytokine receptor interaction",
  "CXCL1",    "hsa04060 Cytokine-cytokine receptor interaction",
  "CXCL2",    "hsa04060 Cytokine-cytokine receptor interaction",
  "CXCL3",    "hsa04060 Cytokine-cytokine receptor interaction",
  "CXCL8",    "hsa04060 Cytokine-cytokine receptor interaction",
  "CXCL14",   "hsa04060 Cytokine-cytokine receptor interaction",
  
  "CXCL1",    "hsa04062 Chemokine signaling pathway",
  "CXCL2",    "hsa04062 Chemokine signaling pathway",
  "CXCL3",    "hsa04062 Chemokine signaling pathway",
  "CXCL8",    "hsa04062 Chemokine signaling pathway",
  "CXCL14",   "hsa04062 Chemokine signaling pathway",
  
  "CXCL1",    "hsa04061 Viral protein interaction with cytokine and cytokine receptor",
  "CXCL2",    "hsa04061 Viral protein interaction with cytokine and cytokine receptor",
  "CXCL3",    "hsa04061 Viral protein interaction with cytokine and cytokine receptor",
  "CXCL8",    "hsa04061 Viral protein interaction with cytokine and cytokine receptor",
  "CXCL14",   "hsa04061 Viral protein interaction with cytokine and cytokine receptor",
  
  "CXCL1",    "hsa04936 Alcoholic liver disease",
  "CXCL2",    "hsa04936 Alcoholic liver disease",
  "CXCL3",    "hsa04936 Alcoholic liver disease",
  "CXCL8",    "hsa04936 Alcoholic liver disease",
  
  "CXCL1",    "hsa05417 Lipid and atherosclerosis",
  "CXCL2",    "hsa05417 Lipid and atherosclerosis",
  "CXCL3",    "hsa05417 Lipid and atherosclerosis",
  "CXCL8",    "hsa05417 Lipid and atherosclerosis",
  
  "CXCL1",    "hsa04621 NOD-like receptor signaling pathway",
  "CXCL2",    "hsa04621 NOD-like receptor signaling pathway",
  "CXCL3",    "hsa04621 NOD-like receptor signaling pathway",
  "CXCL8",    "hsa04621 NOD-like receptor signaling pathway",
  
  "CXCL1",    "hsa05167 Kaposi sarcoma-associated herpesvirus infection",
  "CXCL2",    "hsa05167 Kaposi sarcoma-associated herpesvirus infection",
  "CXCL3",    "hsa05167 Kaposi sarcoma-associated herpesvirus infection",
  "CXCL8",    "hsa05167 Kaposi sarcoma-associated herpesvirus infection",
  
  "CXCL1",    "hsa05323 Rheumatoid arthritis",
  "CXCL2",    "hsa05323 Rheumatoid arthritis",
  "CXCL3",    "hsa05323 Rheumatoid arthritis",
  "CXCL8",    "hsa05323 Rheumatoid arthritis",
  
  "CXCL1",    "hsa05134 Legionellosis",
  "CXCL2",    "hsa05134 Legionellosis",
  "CXCL3",    "hsa05134 Legionellosis",
  "CXCL8",    "hsa05134 Legionellosis",
  
  "CXCL1",    "hsa05120 Epithelial cell signaling in Helicobacter pylori infection",
  "CXCL2",    "hsa05120 Epithelial cell signaling in Helicobacter pylori infection",
  "CXCL3",    "hsa05120 Epithelial cell signaling in Helicobacter pylori infection",
  "CXCL8",    "hsa05120 Epithelial cell signaling in Helicobacter pylori infection",
  
  "CXCL1",    "hsa04064 NF-kappa B signaling pathway",
  "CXCL2",    "hsa04064 NF-kappa B signaling pathway",
  "CXCL3",    "hsa04064 NF-kappa B signaling pathway",
  "CXCL8",    "hsa04064 NF-kappa B signaling pathway",
  
  "CXCL1",    "hsa05146 Amoebiasis",
  "CXCL2",    "hsa05146 Amoebiasis",
  "CXCL3",    "hsa05146 Amoebiasis",
  "CXCL8",    "hsa05146 Amoebiasis",
  
  "CXCL1",    "hsa04657 IL-17 signaling pathway",
  "CXCL2",    "hsa04657 IL-17 signaling pathway",
  "CXCL3",    "hsa04657 IL-17 signaling pathway",
  "CXCL8",    "hsa04657 IL-17 signaling pathway",
  
  "CXCL1",    "hsa04668 TNF signaling pathway",
  "CXCL2",    "hsa04668 TNF signaling pathway",
  "CXCL3",    "hsa04668 TNF signaling pathway"
)


# 2. Custom colors
gene_colors <- c("CXCL17"="#DAA520", "CXCL1"="#E41A1C", "CXCL2"="#377EB8", "CXCL3"="#4DAF4A", "CXCL8"="#984EA3","CXCL14"="#FF1493")
all_pathways <- unique(gene_pathway$pathway)
pathway_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(all_pathways))
pathway_colors <- setNames(pathway_palette, all_pathways)
all_colors <- c(gene_colors, pathway_colors)

# 3. Split KEGG ID and wrap pathway name helper
split_and_wrap_hsaid_label <- function(label, wrap_width = 24) {
  m <- regexpr("^hsa\\d+\\s+", label)
  if (m[1] > 0) {
    id <- substr(label, 1, attr(m, "match.length") - 1)
    txt <- substr(label, attr(m, "match.length") + 1, nchar(label))
    wrapped_txt <- paste(strwrap(txt, width = wrap_width), collapse = "\n")
    return(paste0(id, "\n", wrapped_txt))
  } else {
    return(label)
  }
}

# 4. Plot chord diagram
png(filename = file.path(base_dir, "CXCL_ChordDiagram_split_and_wrap.png"), width = 2400, height = 2000, res = 220)
par(mar = c(2, 2, 2, 2))  # Extra large margin for text

chordDiagram(
  gene_pathway,
  grid.col = all_colors,
  transparency = 0.35,
  annotationTrack = c("grid"),
  preAllocateTracks = 1,
  directional = 0,
  link.sort = TRUE,
  link.decreasing = FALSE
)

circos.trackPlotRegion(
  track.index = 1, 
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    if (sector.name %in% names(gene_colors)) {
      circos.text(mean(xlim), ylim[1], sector.name, facing="clockwise", niceFacing=TRUE,
                  adj=c(0,0.5), font=2, cex=1.0)
    } else {
      split_label <- split_and_wrap_hsaid_label(sector.name, wrap_width = 22)
      circos.text(mean(xlim), ylim[1], split_label, facing="clockwise", niceFacing=TRUE,
                  adj=c(0,0.5), font=1, cex=0.85)
    }
  }, bg.border = NA
)

dev.off()