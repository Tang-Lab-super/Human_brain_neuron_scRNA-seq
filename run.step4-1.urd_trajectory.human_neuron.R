###################################################################################################
# human neuron scRNA-seq data URD trajectory analysis                                             #
###################################################################################################

library(glue)
library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(RColorBrewer)
library(cowplot)
library(stringr)
require("GenomicRanges")
require("Signac")
options(future.globals.maxSize = 30 * 1024^3)

#################################################################################################################################################
cell_colors <- c("Astro"="#081D58","Oligo"="#253494","RG"="#F46D43","EOMES_IPC"="#FEE090","ASCL1_IPC"="#FFFFBF","EN"="#41AB5D","IN"="#54278F","IN_Precursor"="#6A51A3")

obj <- readRDS("data_merged.rds")
obj$cellType <- NA
obj$cellType[obj$celltype == "EOMES+ IPC"] <- "EOMES_IPC"
obj$cellType[obj$celltype == "ASCL1+ IPC"] <- "ASCL1_IPC"
obj$cellType[obj$celltype == "Astro"] <- "Astro"
obj$cellType[obj$celltype == "Oligo"] <- "Oligo"
obj$cellType[obj$celltype == "RG"] <- "RG"
obj$cellType[obj$celltype == "EN"] <- "EN"
obj$cellType[obj$celltype == "IN"] <- "IN"
obj$cellType[obj$celltype == "IN_Precursor"] <- "IN_Precursor"

pdf(file="DimPlot.combine_by_cellType.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap", pt.size=1.2, label=TRUE, label.size = 6, group.by='cellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="DimPlot.combine_by_sample.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

developmentalTrajectories_obj <- subset(obj, subset=cellType %in% c("Astro","Oligo"), invert=TRUE)
pdf(file="DimPlot.combine_by_cellType.pdf", height=10,width=10)
DimPlot(developmentalTrajectories_obj, reduction = "umap", pt.size=1.2, label=TRUE, label.size = 6, group.by='cellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()

URDobj <- developmentalTrajectories_obj
URDobj <- RunTSNE(URDobj, reduction = "harmony", assay = "RNA", dims = 1:50, reduction.name="tsne")
pdf(file="DimPlot_tsne.combine_by_cellType.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "tsne", pt.size=1.2, label=TRUE, label.size = 6, group.by='cellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="DimPlot.URDobj.combine_by_sample.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "tsne", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()
saveRDS(URDobj, file="URDobj.rds")


#################################################################################################################################################
#######################################################################################
# conda create -n R4.2
# source activate R4.2
# conda install -c conda-forge r-base
# conda install udunits2
# conda install -c conda-forge libgit2
# source("https://raw.githubusercontent.com/farrellja/URD/master/URD-Install.R")

#######################################################################################
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))
knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()

cell_colors <- c("Astro"="#081D58","Oligo"="#253494","RG"="#F46D43","EOMES_IPC"="#FEE090","ASCL1_IPC"="#FFFFBF","EN"="#41AB5D","IN"="#54278F","IN_Precursor"="#6A51A3")

#######################################################################################
## Create URD object
URDobj <- readRDS("URDobj.rds")
seurat.object <- URDobj
# Copy over data
ds <- new("URD")
ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
ds@count.data <- as(as.matrix(seurat.object@assays$RNA@counts[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
ds@logupx.data <- ds@logupx.data[-grep("XIST",rownames(ds@logupx.data)),]
ds@count.data <- ds@count.data[-grep("XIST",rownames(ds@count.data)),]
get.data <- as.data.frame(seurat.object@meta.data) 
get.data <- get.data[,c("nCount_RNA","nFeature_RNA","sample","seurat_clusters","cellType")]
ds@meta <- get.data
ds@group.ids <- get.data[,c("sample","seurat_clusters","cellType")]
# Move over tSNE projection
ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
# Move over PCA results
ds@pca.load <- as.data.frame(seurat.object@reductions$harmony@feature.loadings)
ds@pca.scores <- as.data.frame(seurat.object@reductions$harmony@cell.embeddings)
ds@pca.sdev <- seurat.object@reductions$harmony@stdev
ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
axial <- ds
#######################################################################################
## check tSNE (CellType and sample)
pdf(file="tSNE_CellType.pdf", width=10, height=10)
plotDim(axial, "cellType", plot.title = "tSNE: Stage") + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="tSNE_sample.pdf", width=10, height=10)
plotDim(axial, "sample", plot.title = "tSNE: Stage")
dev.off()
#######################################################################################
## Calculate variable genes
# (Normally would do this for each stage, but there are not very many cells in this subset of the data)
# diffCV.cutoff can be varied to include more or fewer genes.
# Copy stage from @meta to @group.ids 
axial@group.ids$stage <- as.character(axial@group.ids$cellType)
# Find a list of cells from each stage.
stages <- c("IN_Precursor","IN","EN","ASCL1_IPC","RG","EOMES_IPC")
cells.each.stage <- lapply(stages, function(stage) rownames(axial@meta)[which(axial@meta$cellType == stage)])
# Compute variable genes for each stage.
var.by.stage <- lapply(1:length(stages), function(n) findVariableGenes(axial, cells.fit = cells.each.stage[[n]], set.object.var.genes = F, diffCV.cutoff = 0.6, mean.min = 0.005, mean.max = 100, main.use = stages[[n]], do.plot = T))
# Combine the results from each group of stages into a single list of variable genes and load into the URD object
var.genes <- sort(unique(unlist(var.by.stage)))
axial@var.genes <- var.genes
#######################################################################################
## Calculate Diffusion Map
# In this case, knn=100 (larger than sqrt(n.cells)) works well because there are not many cell types.
# Sigma 16 is slightly smaller than the sigma auto-determined by using NULL parameter.
axial <- calcDM(axial, knn = 200, sigma=10)
saveRDS(axial,file="axial_calcDM.rds")

# axial <- readRDS("axial_calcDM.rds")
pdf(file="DiffusionMap.pdf", width=10, height=10)
plotDimArray(axial, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map (Sigma 16, 100 NNs): Stage", label="stage", plot.title="", legend=F)
dev.off()
pdf(file="tSNE_transitions.pdf", width=10, height=10)
plotDim(axial, "cellType", transitions.plot = 10000, plot.title="Developmental stage (with transitions)") + scale_color_manual(values=c(cell_colors))
dev.off()


axial <- readRDS("axial_calcDM.rds")
pdf(file="tSNE_transitions_cellType.pdf", width=10, height=10)
plotDim(axial, "cellType", transitions.plot = 10000, plot.title="Developmental stage (with transitions)") + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="tSNE_transitions_sample.pdf", width=10, height=10)
plotDim(axial, "sample", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")
dev.off()

df <- axial@group.ids
axialGW11 <- urdSubset(axial, rownames(df[which(df$sample=="GW11_d20"),]))
axialGW20 <- urdSubset(axial, rownames(df[which(df$sample=="GW20_d20"),]))

pdf(file="tSNE_transitions_GW11.pdf", width=10, height=10)
plotDim(axialGW11, "cellType", transitions.plot = 10000, plot.title="Developmental stage (with transitions)") + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="tSNE_transitions_GW20.pdf", width=10, height=10)
plotDim(axialGW20, "cellType", transitions.plot = 10000, plot.title="Developmental stage (with transitions)") + scale_color_manual(values=c(cell_colors))
dev.off()


####################################################################################################
# Select transition according to the direction of differentiation between cell subtypes

# plot function
plot_selected_transitions <- function(obj_urd, transition_file, output_file) {
    data.plot = obj_urd@dm@eigenvectors
    dim.x = 1
    dim.y = 2
    dim.x <- paste0("DC", dim.x)
    dim.y <- paste0("DC", dim.y)
    transitions.alpha = 0.5

    transitions.df <- read.table(transition_file, row.names=1, header=T, sep='\t')
    transitions.df$x1 <- data.plot[transitions.df$from, dim.x]
    transitions.df$x2 <- data.plot[transitions.df$to, dim.x]
    transitions.df$y1 <- data.plot[transitions.df$from, dim.y]
    transitions.df$y2 <- data.plot[transitions.df$to, dim.y]
    transitions.df$alpha <- transitions.df$weight/max(transitions.df$weight) *
        transitions.alpha
    transitions.df$from_ct = NULL
    transitions.df$to_ct = NULL


    pdf(file=output_file, width=10, height=10)
    p <- plotDim(obj_urd, "cellType", transitions.plot = 10000, plot.title="Developmental stage (with transitions)", point.size=2, transitions.df=transitions.df) + scale_color_manual(values=c(cell_colors))
    print(p)
    dev.off()
}

# write all transition information 
gw11_transion = edgesFromDM(axialGW11)
write.table(gw11_transion, 'axialGW11.allTransition.txt', sep='\t', quote=F)

gw20_transion = edgesFromDM(axialGW20)
write.table(gw20_transion, 'axialGW20.allTransition.txt', sep='\t', quote=F)

# > head(gw11_transion, n=3)
#                          from                          to       weight
# 1 GW11_d20_CGATTGACAACTTGAC-1 GW11_d20_GTGTGCGGTAAGGGAA-1 0.0004752729
# 2 GW11_d20_CGATTGACAACTTGAC-1 GW11_d20_TCAGCAATCCGCATAA-1 0.0004050893
# 3 GW11_d20_CGATTGACAACTTGAC-1 GW11_d20_TTCTTAGGTGACTCAT-1 0.0009065816

# write cell subtypes information
write.table(axialGW11@meta[, 'cellType', drop=F], 'axialGW11.meta.txt', sep='\t', quote=F)
write.table(axialGW20@meta[, 'cellType', drop=F], 'axialGW20.meta.txt', sep='\t', quote=F)


# use a custom python script to select the transitions
# this step will generate 3 files contains selected transitions of gw11, gw20 and merge, respectively
system('python run.step4-2.select_transition.human_neuron.py')


plot_selected_transitions(axialGW11, 'axialGW11.UsedTransition.txt', 'tSNE_transitions_GW11.v2.pdf')
plot_selected_transitions(axialGW20, 'axialGW20.UsedTransition.txt', 'tSNE_transitions_GW20.v2.pdf')
plot_selected_transitions(axial, 'axialMerge.UsedTransition.txt', 'tSNE_transitions_Merge.v2.pdf')


######################################################## feature plot 
for (gene in c('ASCL1', 'DLX2', 'MKI67', 'NPY', 'SCGN')) {
    pdf(file=glue("Expr.{gene}_merge.pdf"), width=10, height=10)
    p <- plotDim(axial, gene, plot.title=gene, colors=rev(paletteer_c("grDevices::Spectral", 50)), point.size=2)
    print(p)
    dev.off()

    pdf(file=glue("Expr.{gene}_GW11.pdf"), width=10, height=10)
    p <- plotDim(axialGW11, gene, plot.title=gene, colors=rev(paletteer_c("grDevices::Spectral", 50)), point.size=2)
    print(p)
    dev.off()

    pdf(file=glue("Expr.{gene}_GW20.pdf"), width=10, height=10)
    p <- plotDim(axialGW20, gene, plot.title=gene, colors=rev(paletteer_c("grDevices::Spectral", 50)), point.size=2)
    print(p)
    dev.off()
}