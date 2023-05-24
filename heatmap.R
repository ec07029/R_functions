# Heatmap with color scale control

library(pheatmap)
#library(dichromat)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

#########################
# pheatmap
# With color scale kept the same for any plot and transitions from blue to red at 0.1
pheatmap( log2(lsc.heatmap.input+1), 
          main = paste0("LSC-unique gene targets\np.adj < 0.1; ", nrow(lsc.filter), " / ", nrow(lsc), " significant gene targets\n",
                        "( log2(diff.freq+1) )"),
          color = colorRampPalette(c("blue", "white", "red"))(10), # number of colors must equal number of breaks
          breaks = c( seq(-0.3, 0.1, length.out=5), seq(0.1+1/15, 1, length.out=5) ),
          cex=1.5,
          angle_col=0,
          labels_col=c("\nLSC\ndifferential edit frequency", "\nLSK\ndifferential edit frequency"),
          show_rownames = FALSE, cluster_cols = FALSE )

# Adding annotation to rows
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

# define the annotation
annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6))),
  AdditionalAnnotation = c(rep("random1", 10), rep("random2", 10))
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")

pheatmap(test, annotation_row = annotation_row)



#########################
# ComplexHeatmap package
# ht_opt(heatmap_row_names_gp = gpar(fontface = "italic")) # italicize genes
Heatmap(heatmap.diff.freq.input,
        col=colorRamp2(c(min(heatmap.diff.freq.input),
                         0.1,
                         max(heatmap.diff.freq.input)),
                       c("blue", "white", "red")),
        #column_title = "Differential Editing Frequency",
        name="Differential Editing Frequency",
        show_row_names = FALSE,
        cluster_columns = FALSE) +
  # Heatmap(heatmap.ann,
  #         col=c("red", "blue", "green", "gray"),
  #         name="Cell Type",
  #         show_row_names = FALSE) +
  Heatmap(heatmap.gene.expression.mig.input,
          col=colorRamp2(c(gene.x.min, gene.x.middle, gene.x.max),
                         c("blue", "white", "red")),
          column_title = "Gene Expression in MIG",
          name="MIG: VST read counts",
          show_row_names = FALSE,
          cluster_columns = FALSE) +
  Heatmap(heatmap.gene.expression.adar.input,
          col=colorRamp2(c(gene.x.min, gene.x.middle, gene.x.max),
                         c("blue", "white", "red")),
          column_title = "Gene Expression in MSI2-ADA",
          name="MSI2-ADA: VST read counts",
          show_row_names = TRUE,
          cluster_columns = FALSE)

# Add borders
decorate_heatmap_body(c("Differential Editing Frequency"), {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  
})

decorate_heatmap_body(c("MIG: VST read counts"), {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  
})

decorate_heatmap_body(c("MSI2-ADA: VST read counts"), {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  
})