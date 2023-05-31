###################################################################################################
# human neuron scRNA-seq data integration                                                         #
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
# load pre-processed data
GW11B3_d20_Right <- readRDS("Human_neuron_GW11B3_d20_Right.rds")
GA20W_d24 <- readRDS("Human_neuron_GA20W_d24.rds")

##### data combined and harmony integration
# downsampling GW11B3_d20_Right to 5000 cells
cells <- Cells(GW11B3_d20_Right)
set.seed(2023)
cells_ds <- sample(cells,size=5000,replace = FALSE)

GW11B3_d20_Right_ds <- subset(GW11B3_d20_Right, cells = cells_ds)
GW11B3_d20_Right_ds
#An object of class Seurat
#24423 features across 5000 samples within 1 assay
#Active assay: RNA (24423 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

#saveRDS(GW11B3_d20_Right_ds, file="Human_neuron_GW11B3_d20_Right_ds.rds")

GA20W_d24
#An object of class Seurat
#21424 features across 4780 samples within 1 assay
#Active assay: RNA (21424 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

# data combined
data_merged <- merge(GA20W_d24,y=GW11B3_d20_Right_ds)

# get variable genes
genes_var <- unique(c(VariableFeatures(GA20W_d24),VariableFeatures(GW11B3_d20_Right_ds)))
length(genes_var)
#[1] 3179

# only scaled data can be applied for PCA
data_merged <- ScaleData(data_merged)

# run PCA
data_merged <- RunPCA(
        object = data_merged,
        features = genes_var,
        assay = "RNA",
        npcs = 30,
        verbose = TRUE)

pdf(file='RunHarmony.plot.final.pdf',height=12,width=13)
data_merged <- RunHarmony(data_merged,
       group.by.vars=c("sample","condition"),
       assay.use="RNA",
       reduction="pca",
       dims.use=1:30,
       theta=c(2,1),
       max.iter.harmony = 40, max.iter.cluster=200,
       kmeans_init_nstart=2, kmeans_init_iter_max=1000,
       return_object = TRUE,
       plot_convergence = TRUE)
dev.off()

data_merged <- RunUMAP(data_merged, dims = 1:30, reduction = "harmony")
data_merged <- FindNeighbors(data_merged, reduction = "harmony", dims = 1:30)
data_merged <- FindClusters(data_merged,resolution=1.5)

pdf("data_combined.harmony.pdf",width=16,height=7.2)
DimPlot(data_merged, group.by = c("sample", "ident"), ncol = 2, label=T)
dev.off()

data_merged <- subset(data_merged,idents=c(8,18),invert=T)

### celltype annotation ###
##### data reference mapping #####
ref <- readRDS("InVitro_Full_STICR.rds")
ref
#An object of class Seurat
#61023 features across 121290 samples within 3 assays
#Active assay: RNA (33538 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 2 dimensional reductions calculated: pca, umap
table(Idents(ref))
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
#6573 5888 5491 5408 5025 4913 4602 4589 4577 4537 4476 4337 4275 4029 3849 3801
#  17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32
#3777 3623 3334 3264 3255 3051 3034 2947 2911 2805 2415 2278 2241 2011 1917  898
#  33   34
# 616  543
pdf(file="reference.umap.cluster.plot.pdf", height=7,width=8.2)
DimPlot(ref, reduction = "umap", label=TRUE, pt.size=0.6,label.size = 4)
dev.off()

ref <- RenameIdents(object = ref,
'3'='DLX2+ IPCs', '4'='DLX2+ IPCs', '6'='DLX2+ IPCs', '23'='DLX2+ IPCs', '29'='DLX2+ IPCs', '1'='Glial', '7'='Glial', '8'='Glial', '10'='Glial', '22'='Glial', '24'='Glial', '25'='Glial', '30'='Glial', '31'='EOMES+ IPCs', '16'='EN', '21'='EN', '26'='EN', '2'='IN', '5'='IN','27'='IN', '28'='IN','32'='IN', '14'='IN','18'='IN', '11'='IN','12'='IN', '13'='IN','19'='IN', '9'='IN','20'='IN', '15'='IN','17'='IN', '33'='IN','34'='Other')
table(Idents(ref))
# DLX2+ IPCs       Glial EOMES+ IPCs          EN          IN       Other
#      21087       31221        1917        9861       56661         543
ref$celltype <- Idents(ref)
pdf(file="reference.umap.celltype.plot.pdf", height=7,width=8)
DimPlot(ref, reduction = "umap", label=TRUE, pt.size=0.6,label.size = 4)
dev.off()

DefaultAssay(ref) <- "integrated"
ref.anchors <- FindTransferAnchors(reference = ref, query = data_merged, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = ref.anchors, refdata = ref$celltype, dims = 1:30)
head(predictions)
#                            predicted.id prediction.score.Glial
#GW20_d20_AAACCTGAGCAGCCTC-1   DLX2+ IPCs             0.00000000
#GW20_d20_AAACCTGAGCTATGCT-1   DLX2+ IPCs             0.00000000
#GW20_d20_AAACCTGAGGGCTTGA-1           EN             0.00000000

data_merged <- AddMetaData(data_merged, metadata = predictions)
table(data_merged$predicted.id)
# DLX2+ IPCs          EN EOMES+ IPCs       Glial          IN
#       2254        4685         105         524        1683

pdf(file="data_combined.harmony.umap.celltype.reference.plot.pdf", height=7,width=8.2)
DimPlot(data_merged, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3,
    repel = TRUE) + ggtitle("Reference annotations")
dev.off()

### rename cluster idents
data_merged <- RenameIdents(object = data_merged,
'0'='EN', '1'='EN', '3'='EN', '14'='EN', '15'='EN', '21'='EN', '2'='EOMES+ IPC', '16'='Astro', '13'='RG', '20'='RG', '17'='Oligo', '9'='ASCL1+ IPC', '6'='IN_Precursor', '7'='IN_Precursor', '12'='IN', '19'='IN', '5'='IN', '11'='IN', '4'='IN', '10'='IN', '22'='IN', '23'='IN')
table(Idents(data_merged))
#          EN   EOMES+ IPC        Astro           RG        Oligo   ASCL1+ IPC
#        3770          871          202          353          193          382
#IN_Precursor           IN
#         943         2537

data_merged$celltype <- Idents(data_merged)
table(data_merged$celltype)

cluster_colors <- rev(RColorBrewer::brewer.pal(8,"Paired"))[c(1,6,3,8,4,5,7,2)]
pdf(file="data_combined.harmony.umap.celltype.plot.pdf", height=7,width=8)
DimPlot(data_merged, reduction = "umap", label=TRUE, cols=cluster_colors, pt.size=0.5, label.size = 4)
dev.off()

sample_colors <- c('#d73027','#4575b4')
pdf(file="data_combined.harmony.umap.sample.plot.pdf", height=7,width=8)
DimPlot(data_merged, reduction = "umap", group.by="sample", label=TRUE, cols=sample_colors, pt.size=0.5, label.size = 4)
dev.off()

# featureplot visualization
#RG：HOPX
#NPC:EOMES
#Pro_NPC:MKI67
#OLIG:OLIG1
#EN:NEUROD2,SLA
#IN:GAD1,GAD2
#ASTR:SLC1A3
pdf("data_combined.harmony.celtype_FeaturePlot.pdf",height=8.5,width=10.5)
FeaturePlot(data_merged, features = c("HOPX", "EOMES","MKI67","OLIG1","NEUROD2","SLA","ASCL1","GAD2","SLC1A3"),ncol=3,min.cutoff="q5",max.cutoff="q95",pt.size=0.1,order=T)
dev.off()

# dotplot visualization
pdf("data_combined.harmony.celtype_DotPlot.pdf",height=5,width=6)
DotPlot(data_merged, features = c("NEUROD2", "NEUROD6","EOMES","SLC1A3","ATP1A2","HOPX","VIM","MKI67","TOP2A","SOX2","OLIG1","OLIG2","GAD2","DLX1","DLX2","ASCL1"), col.min=0, col.max=1) + coord_flip() + RotatedAxis()
dev.off()

### sample-cluster stats
sample_stats <- table(data_merged$sample, Idents(data_merged))
library(reshape2)
sample_stats <- melt(sample_stats)
colnames(sample_stats) <- c("Sample","Cluster","Value")
write.table(sample_stats,"All_sample_cluster_stats_by_sample_number.txt",sep="\t",quote=F)

sample_colors <- c('#d73027','#4575b4')
pdf("All_sample_cluster_stats_by_cluster.pdf",height=6,width=6.5)
ggplot(sample_stats,aes(Cluster,Value,fill=Sample)) + geom_bar(stat="identity",position="fill") + geom_hline(yintercept=0.5,color="black",linetype=2) + scale_fill_manual(values=sample_colors) + scale_y_continuous(expand = c(0,0)) +
theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
#    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
dev.off()

pdf("All_sample_cluster_stats_by_sample.pdf",height=5.5,width=5)
ggplot(sample_stats,aes(Sample,Value,fill=Cluster)) + geom_bar(stat="identity",position="fill") + scale_fill_manual(values=cluster_colors) + scale_y_continuous(expand = c(0,0)) +
theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
#    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
dev.off()

# save metadata information
meta_data <- data_merged@meta.data
head(meta_data)
write.table(meta_data,"All_sample_meta_data.txt",sep="\t",quote=F)

saveRDS(data_merged,"data_merged.rds")
