###################################################################################################
# Pre-processing human neuron scRNA-seq data                                                      #
###################################################################################################
rm(list=ls())

# load required packages
library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
library(harmony)
library(RColorBrewer)
library(future)
# it is tricky to set options
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 30 * 1024^3)

########################################################################################################################
### run SoupX estimation ###
library(SoupX)

GW11B3_d20_sc <- load10X("/data/ProcessedData/scRNA-seq/9_GW11B3_d20merge_5scRNA_20220317_B1S2/outs/")
#Loading raw count data
#Loading cell-only count data
#Loading extra analysis data where available
GA20W_d24_sc <- load10X("/data/ProcessedData/scRNA-seq/GA20Wmerge_301_d24_5scRNA_20220117_B1S2/outs/")
#Loading raw count data
#Loading cell-only count data
#Loading extra analysis data where available

GW11B3_d20_sc = setContaminationFraction(GW11B3_d20_sc, 0.02)
GA20W_d24_sc = setContaminationFraction(GA20W_d24_sc, 0.02)

GW11B3_d20_out = adjustCounts(GW11B3_d20_sc)
GA20W_d24_out = adjustCounts(GA20W_d24_sc)

# define function for read data
create_obj <- function(count_data, sample_name, condition, tissue_source, week, scrna_kit, project) {
  data_obj <- CreateSeuratObject(counts = count_data, min.cells = 3, min.features = 200)

  data_obj <- RenameCells(object = data_obj, add.cell.id = sample_name)
  data_obj[["sample"]] <- sample_name
  data_obj[["condition"]] <- condition
  data_obj[["tissue_source"]] <- tissue_source
  data_obj[["week"]] <- week
  data_obj[["scRNA_Kit"]] <- scrna_kit
  data_obj[["project"]] <- project

  data_obj <- PercentageFeatureSet(data_obj, pattern = "^MT-", col.name = "percent.mt")
  pdf(file=paste0("QC_metrics_raw.", sample_name, ".pdf"), height=10, width=30)
  print(VlnPlot(data_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()

# filter cell with high percentage of mt reads
  data_obj <- subset(data_obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# filter mt and ribo genes
  genes_use <- rownames(data_obj[['RNA']]@counts)
  genes_use <- genes_use[which(!grepl('^MT-', genes_use))]
  genes_use <- genes_use[which(!grepl('^RPS', genes_use))]
  genes_use <- genes_use[which(!grepl('^RPL', genes_use))]
  data_obj <- subset(data_obj, features=genes_use)

# normalize data
  data_obj <- NormalizeData(data_obj)

  return(data_obj)
}

# read data
GW11B3_d20 <- create_obj(count_data = GW11B3_d20_out, sample_name = "GW11_d20", condition="D20", tissue_source="Cortex", week="Week_11", scrna_kit="scRNA_5end_R2_only", project='Human_Neuron')
GA20W_d24 <- create_obj(count_data = GA20W_d24_out, sample_name = "GW20_d20", condition="D20", tissue_source="Cortex", week="Week_20", scrna_kit="scRNA_5end_R2_only", project='Human_Neuron')

###################################################################################################
##### for GW11B3_d20 sample data pre-processing                                               #####
###################################################################################################
GW11B3_d20 <- FindVariableFeatures(GW11B3_d20, selection.method = "vst", nfeatures = 2000)
GW11B3_d20 <- ScaleData(GW11B3_d20,features = row.names(GW11B3_d20))

# run PCA and UMAP
GW11B3_d20 <- RunPCA(object = GW11B3_d20, assay = "RNA", npcs = 50, verbose = TRUE)
GW11B3_d20 <- RunUMAP(GW11B3_d20, dims = 1:50)
GW11B3_d20 <- FindNeighbors(GW11B3_d20, dims = 1:50)
GW11B3_d20 <- FindClusters(GW11B3_d20,resolution=1.7)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(GW11B3_d20)))+1
jBrewColors <- SpatialColors(num)

pdf(file="Human_neuron_GW11B3_d20_DimPlot.pdf", height=12,width=13)
DimPlot(GW11B3_d20, reduction = "umap", pt.size=1.2, label=T, cols=jBrewColors)
dev.off()
#saveRDS(GW11B3_d20, file="Human_neuron_GW11B3_d20.rds")

### subset Right data
df <- GW11B3_d20@meta.data
df$seurat_clusters_new <- "middle"
df[which(df$seurat_clusters %in% c(2,3,5,7,8,13,27,18,16,23,24,22,15,19,29)),]$seurat_clusters_new <- "Left"
df[which(df$seurat_clusters %in% c(10,14,20,1,9,0,4,26,11,17,12,25)),]$seurat_clusters_new <- "Right"
GW11B3_d20@meta.data <- df

GW11B3_d20_Right <- subset(GW11B3_d20, subset=seurat_clusters_new=='Right')
GW11B3_d20_Right <- ScaleData(GW11B3_d20_Right,features = row.names(GW11B3_d20_Right))
GW11B3_d20_Right <- FindVariableFeatures(GW11B3_d20_Right, selection.method = "vst", nfeatures = 2000)
GW11B3_d20_Right <- RunPCA(object = GW11B3_d20_Right, assay = "RNA", npcs = 50, verbose = TRUE)
GW11B3_d20_Right <- RunUMAP(GW11B3_d20_Right, dims = 1:50)
GW11B3_d20_Right <- FindNeighbors(GW11B3_d20_Right, dims = 1:50)
GW11B3_d20_Right <- FindClusters(GW11B3_d20_Right,resolution=1.7)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(GW11B3_d20_Right)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))

pdf(file="GW11B3_d20_Right.umap.plot.pdf", height=12,width=13.5)
DimPlot(GW11B3_d20_Right, reduction = "umap", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 8)
dev.off()

### remove predicted doublets
library(DoubletFinder)

# define the expected number of doublet cells.
# loading 20,000 cells, expect 10% doublets
nExp <- round(ncol(GW11B3_d20_Right) * 0.10)
GW11B3_d20_Right <- doubletFinder_v3(GW11B3_d20_Right, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:50)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(GW11B3_d20_Right@meta.data)[grepl("DF.classification", colnames(GW11B3_d20_Right@meta.data))]
DF.name
#[1] "DF.classifications_0.25_0.09_1564"

pdf("GW11B3_d20_Right.doublet.umap.plot.pdf",height=5,width=11)
cowplot::plot_grid(ncol = 2, DimPlot(GW11B3_d20_Right, group.by = "sample") + NoAxes(),
                             DimPlot(GW11B3_d20_Right, group.by = DF.name) + NoAxes())
dev.off()

# 过滤双细胞
GW11B3_d20_Right <- GW11B3_d20_Right[, GW11B3_d20_Right@meta.data[, DF.name] == "Singlet"]

saveRDS(GW11B3_d20_Right, file="Human_neuron_GW11B3_d20_Right.rds")

###################################################################################################
##### for GA20W_d24 sample data pre-processing                                                #####
###################################################################################################
GA20W_d24 <- FindVariableFeatures(GA20W_d24, selection.method = "vst", nfeatures = 2000)
GA20W_d24 <- ScaleData(GA20W_d24,features = row.names(GA20W_d24))

# run PCA and UMAP
GA20W_d24 <- RunPCA(object = GA20W_d24, assay = "RNA", npcs = 50, verbose = TRUE)
GA20W_d24 <- RunUMAP(GA20W_d24, dims = 1:50)
GA20W_d24 <- FindNeighbors(GA20W_d24, dims = 1:50)
GA20W_d24 <- FindClusters(GA20W_d24)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(GA20W_d24)))+1
jBrewColors <- SpatialColors(num)

pdf(file="Human_neuron_GA20W_d24_DimPlot.pdf", height=12,width=13)
DimPlot(GA20W_d24, reduction = "umap", pt.size=1.2, label=T, cols=jBrewColors)
dev.off()

### remove predicted doublets
library(DoubletFinder)

# define the expected number of doublet cells.
# loading 5000 cells, expect 3% doublets
nExp <- round(ncol(GA20W_d24) * 0.03)
GA20W_d24 <- doubletFinder_v3(GA20W_d24, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:30)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(GA20W_d24@meta.data)[grepl("DF.classification", colnames(GA20W_d24@meta.data))]
DF.name
#[1] "DF.classifications_0.25_0.09_148"

pdf("GA20W_d24.doublet.umap.plot.pdf",height=5,width=11)
cowplot::plot_grid(ncol = 2, DimPlot(GA20W_d24, group.by = "sample") + NoAxes(),
                             DimPlot(GA20W_d24, group.by = DF.name) + NoAxes())
dev.off()

# 过滤双细胞
GA20W_d24 <- GA20W_d24[, GA20W_d24@meta.data[, DF.name] == "Singlet"]

saveRDS(GA20W_d24, file="Human_neuron_GA20W_d24.rds")
