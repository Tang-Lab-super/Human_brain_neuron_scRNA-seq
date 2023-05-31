###################################################################################################
# human neuron scRNA-seq data differential analysis                                               #
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
# load data
data_merged <- readRDS("data_merged.rds")

library(MAST)

DefaultAssay(data_merged)
#[1] "RNA"

# All celltype differential analysis
data_merged.markers <- FindAllMarkers(data_merged, assay="RNA", slot="data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
write.table(data_merged.markers,"data_merged_celltype_DEs.txt",quote=F,sep="\t")

library(dplyr)
top20 <- data_merged.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(top20,"data_merged_celltype_DEs_top20.txt",quote=F,sep="\t")

pdf("data_combined.harmony.celtype_Heatmap_top20.pdf",height=5,width=9)
DoHeatmap(data_merged, features = top20$gene, group.colors=cluster_colors) + viridis::scale_fill_viridis()
dev.off()

# GW11B3_d20 vs. GA20W_d24 differential analysis for each celltype
All_DEs <- data.frame()

for(cluster in c("EN","EOMES+ IPC","Astro","RG","ASCL1+ IPC","IN_Precursor","IN")){
	data_cluster <- subset(data_merged,idents=cluster)
	Idents(data_cluster) <- "sample"
	cluster_des <- FindMarkers(data_cluster, assay="RNA", slot="data", min.pct =0, logfc.threshold =0, test.use = "MAST", ident.1="GW11B3_d20", ident.2="GA20W_d24")
	cluster_des$cluster <- cluster
	cluster_des$gene <- rownames(cluster_des)
	# remove XIST gene
	cluster_des <- cluster_des[grep("XIST",cluster_des$gene,invert=T),]
	write.table(cluster_des,paste0("DEs_GW11B3_d20vsGA20W_d24_",cluster,"_all.txt"),sep="\t",quote=F)
	# get significant genes
	cluster_des_sig <- subset(cluster_des, p_val_adj < 0.05 & abs(avg_log2FC) >= 0.5)
	write.table(cluster_des_sig,paste0("DEs_GW11B3_d20vsGA20W_d24_",cluster,"_sig.txt"),sep="\t",quote=F)
	All_DEs <- rbind(All_DEs,cluster_des)
}

write.table(All_DEs,"DEs_GW11B3_d20vsGA20W_d24_all_celltype.txt",sep="\t",quote=F)
