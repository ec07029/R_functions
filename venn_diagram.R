library("VennDiagram")

vd <- venn.diagram(
  x = list(
    "LT: \nWT vs KO" = lt.sig,
    "LT vs others &\n WT vs KO" = lt.vs.others.sig
  ),
  filename = NULL, 
  col = "transparent",
  fill = c("cornflowerblue", "green"),
  alpha = 0.50,
  scaled = TRUE,
  ext.text = TRUE,
  ext.line.lwd = 2,
  ext.dist = -0.15,
  ext.length = 0.9,
  ext.pos = -4,
  inverted = FALSE,
  cex = 2,
  fontfamily = "sans",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen"),
  cat.cex = 2,
  cat.fontfamily = "sans",
  rotation.degree = 0
)

library(grDevices)

png(file=paste0( data.name, "_", name.type, "_DESeq2_significant_genes_comparison.png" ), 800, 800 )
par(mar=c(5,5,5,5))
grid.draw(vd)
dev.off()


# Old version:
# 10/25/19: Started to save all venn diagrams as 36MB image files, which takes up too much memory and too long to load and view. Fixed with version above using library(grDevices)
#   
# library("VennDiagram")
# 
# vd <- venn.diagram(
#   x = list(
#     "LSK" = filtered.edit.lsk$site.id,
#     "LSC" = filtered.edit.lsc$site.id
#   ),
#   filename = "./editSiteComparison.tiff",
#   col = "transparent",
#   fill = c("cornflowerblue", "green"),
#   alpha = 0.50,
#   label.col = c("blue", "blue", "blue"),
#   scaled = TRUE,
#   ext.text = TRUE,
#   ext.line.lwd = 2,
#   ext.dist = -0.15,
#   ext.length = 0.9,
#   ext.pos = -4,
#   inverted = FALSE,
#   cex = 1.5,
#   fontfamily = "sans",
#   fontface = "bold",
#   cat.col = c("darkblue", "darkgreen"),
#   cat.cex = 1.5,
#   cat.fontfamily = "sans",
#   rotation.degree = 0
# )