library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(nloptr)
library(hdf5r)

IR<- readRDS("~/IR.rds")
IR$log10GenesPerUMI <- log10(IR$nFeature_RNA) / log10(IR$nCount_RNA)
metadata <- IR@meta.data
metadata$orig.ident <- factor(metadata$orig.ident, 
                              levels = c("ShamA","ShamB","IRH4","IRH12","IRD1","IRD5A","IRD5B","IRD5C","IRD14","IRD28"))

metadata <- metadata %>%
  dplyr::rename(
    nUMI = nCount_RNA,
    nGene = nFeature_RNA)

metadata %>% 
  ggplot(aes(color=orig.ident, x=nUMI, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  ylab("Cell density") +
  geom_vline(xintercept = 45000)+
  geom_vline(xintercept = 500)


metadata %>% 
  ggplot(aes(color=orig.ident, x=nGene, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 8000)


metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

DimPlot(IR, reduction = "tsne",label =T,group.by = "orig.ident")
DimPlot(IR, reduction = "tsne",label =T,group.by = "labels2")
