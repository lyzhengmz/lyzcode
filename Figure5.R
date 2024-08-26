library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(mindr)
library(patchwork)
library(NMF)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(future)

IR<- readRDS("~/IR.rds")

dim(IR@active.ident)
str(IR@active.ident)
table(IR@active.ident)
plan(sequential)
DefaultAssay(IR) <-"RNA"
IR <- NormalizeData(IR)

data.input  <- IR@assays$RNA@data
data.input[1:4,1:4]
Idents(IR) <- "labels2" 
identity = data.frame(group =IR@active.ident, row.names = names(IR@active.ident))

cellchat <- createCellChat(object = data.input, meta = identity, group.by = "group")
cellchat <- addMeta(cellchat, meta = identity)
cellchat <- setIdent(cellchat, ident.use = "group") 
levels(cellchat@idents)


CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM receptor")

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse) 
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


df.net <- subsetCommunication(cellchat)
df.net

df.net <- subsetCommunication(cellchat, signaling = c("PDGF")) 

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",targets.use = c("Myo","Peri"))


mat <- cellchat@net$count
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


cellchat@netP$pathways
head(cellchat@LR$LRsig)

cellchat@netP$pathways
levels(cellchat@idents)

pathways.show <- c("PDGF")

#circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat,signaling = "PDGF",layout = "circle" ,edge.width.max =10)
netVisual_aggregate(cellchat,signaling = cellchat@netP$pathways,layout = "circle",targets.use = c("Myo") )


bb1<-netVisual_bubble(cellchat.D5, sources.use = c("Myo"),grid.on = F,signaling =c("PDGF"),
                      angle.x = 45,font.size = 15,remove.isolate = FALSE,)
bb2<-netVisual_bubble(cellchat.D5, targets.use = c("Myo"), grid.on = F,signaling =c("PDGF"),
                      angle.x = 45,font.size = 15,remove.isolate = FALSE)
bb1+bb2



