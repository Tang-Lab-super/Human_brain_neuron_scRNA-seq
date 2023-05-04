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
GW11B3_d20_Right
#An object of class Seurat
#24423 features across 13962 samples within 1 assay
#Active assay: RNA (24423 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

GA20W_d24 <- readRDS("Human_neuron_GA20W_d24.rds")
GA20W_d24
#An object of class Seurat
#21424 features across 4785 samples within 1 assay
#Active assay: RNA (21424 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

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
#21424 features across 4785 samples within 1 assay
#Active assay: RNA (21424 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

# data combined
data_merged <- merge(GA20W_d24,y=GW11B3_d20_Right_ds)

# get variable genes
genes_var <- unique(c(VariableFeatures(GA20W_d24),VariableFeatures(GW11B3_d20_Right_ds)))
#length(genes_var)
#[1] 3162

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

data_merged <- RunUMAP(data_merged, dims = 1:20, reduction = "harmony")
data_merged <- FindNeighbors(data_merged, reduction = "harmony", dims = 1:20)
data_merged <- FindClusters(data_merged,resolution=0.7)

pdf("data_combined.harmony.pdf",width=15,height=7.5)
DimPlot(data_merged, group.by = c("sample", "ident"), ncol = 2, label=T)
dev.off()

### rename cluster idents
data_merged <- RenameIdents(object = data_merged,
'0'='EN', '1'='EN', '11'='EN', '12'='EN', '18'='EN', '3'='EOMES+ IPC', '14'='Astro', '15'='Astro', '9'='RG', '16'='RG', '13'='Oligo', '7'='ASCL1+ IPC', '5'='IN_Precursor', '8'='IN_Precursor', '2'='IN', '4'='IN', '10'='IN', '6'='IN', '17'='IN')
table(Idents(data_merged))
#          EN   EOMES+ IPC        Astro           RG        Oligo   ASCL1+ IPC
#        4200          769          273          424          182          416
#IN_Precursor           IN
#         963         2558

cluster_colors <- rev(RColorBrewer::brewer.pal(8,"Paired"))[c(1,6,3,8,4,5,7,2)]
pdf(file="data_combined.harmony.umap.celltype.noLegend.plot.pdf", height=7,width=7.3)
DimPlot(data_merged, reduction = "umap", label=TRUE, cols=cluster_colors, pt.size=0.6, label.size = 4) + NoLegend()
dev.off()

sample_colors <- c('#d73027','#4575b4')
pdf(file="data_combined.harmony.umap.sample.noLegend.plot.pdf", height=7,width=7.3)
DimPlot(data_merged, reduction = "umap", group.by="sample", label=TRUE, cols=sample_colors, pt.size=0.6, label.size = 4) + NoLegend()
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
data_merged$cellType <- Idents(data_merged)
meta_data <- data_merged@meta.data
head(meta_data)
#                                orig.ident nCount_RNA nFeature_RNA    sample
#GA20W_d24_AAACCTGAGCAGCCTC-1 SeuratProject      14249         5047 GA20W_d24
#GA20W_d24_AAACCTGAGCTATGCT-1 SeuratProject      10937         4251 GA20W_d24
#GA20W_d24_AAACCTGAGGGCTTGA-1 SeuratProject       3399         2231 GA20W_d24
write.table(meta_data,"All_sample_meta_data.txt",sep="\t",quote=F)

saveRDS(data_merged,"data_merged.rds")
