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

celltype_marker<-c(
  "Slc27a2",#PT
  "Havcr1",#PT-Injured
  "Slc12a1",#LOH
  "Slc12a3",#DCT
  "Slc4a1",#CD-IC
  "Aqp2",#CD-PC
  "Flt1",#Endo
  "Chil3",#Fib
  "Col1a2",#Myo
  "Myh11",#Peri
  "Igkc",#B lymph
  "Trbc2",#T lymph
  "Cd209a",#cDC
  "Ccl12",#Macro
  "Ifitm1",#Neutrophil
  "Siglech")#pDC

IR$labels2<-factor(IR$labels2,levels = c("PT","PT-Injured","LOH","DCT","CD-IC","CD-PC","Endo","Fib","Myo","Peri","B lymph","T lymph","cDC","Macro","Neutrophil","pDC"))

p<-VlnPlot(IR,features = celltype_marker,split.by = "orig.ident",group.by = "labels2")
p+ theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5))
