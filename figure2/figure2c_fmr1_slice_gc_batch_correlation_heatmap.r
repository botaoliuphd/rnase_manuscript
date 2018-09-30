#module load R/3.4.0
#R
#Heatmap of GC contents of all hippocampal slice samples #

gc <- read.table("slice_gc_summary.txt", row.names=1, header=T)

library("pheatmap")
library("RColorBrewer")
library("grid")


draw_colnames_0 <- function (coln, gaps, ...) {
    coord <- pheatmap:::find_coordinates(length(coln), gaps)
    x     <- coord$coord - 0.5 * coord$size
    res   <- grid::textGrob(
      coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
      vjust = 1, hjust = 0.5, rot = 0, gp = grid::gpar(...)
    )
    return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_0",
  ns = asNamespace("pheatmap")
)

output_name="slice_gc_summary_clustering_heatmap.pdf"
pdf(output_name,width=5,height=15)
par(mar=c(5,5,5,5),cex=5)
colors <- colorRampPalette(brewer.pal(9, "Blues") )(255)
pheatmap(gc, cutree_rows = 2, col=colors, fontsize=25,cellwidth = 50, cellheight = 50,cluster_col=FALSE)
dev.off()