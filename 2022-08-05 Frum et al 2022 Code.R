library(Seurat)
library(SeuratWrappers)
library(CellChat)
library(slingshot)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
#Figure 1
#Integrate Fetal Data, Extract Epithelium, Identify Populations for Cell Chat Analysis
f59ddistal.data <- Read10X(data.dir = "./")
f59dairway.data <- Read10X(data.dir = "./")
f132ddistal.data <- Read10X(data.dir = "./")
f80ddistal.data <- Read10X(data.dir = "./")
f80dairway.data <- Read10X(data.dir = "./")
f103ddistal.data <- Read10X(data.dir = "./")
f103dairway.data <- Read10X(data.dir = "./")


f103dairway <- CreateSeuratObject(counts = f103dairway.data, project = "f103dairway", min.cells = 3, min.features = 200)
f103ddistal <- CreateSeuratObject(counts = f103ddistal.data, project = "f103ddistal", min.cells = 3, min.features = 200)
f132ddistal <- CreateSeuratObject(counts = f132ddistal.data, project = "f132ddistal", min.cells = 3, min.features = 200)
f59dairway <- CreateSeuratObject(counts = f59dairway.data, project = "f59dairway", min.cells = 3, min.features = 200)
f59ddistal <- CreateSeuratObject(counts = f59ddistal.data, project = "f59ddistal", min.cells = 3, min.features = 200)
f80dairway  <- CreateSeuratObject(counts = f80dairway.data, project = "f80dairway", min.cells = 3, min.features = 200)
f80ddistal <- CreateSeuratObject(counts = f80ddistal.data, project = "f80ddistal", min.cells = 3, min.features = 200)


f103dairway <- RenameCells(object = f103dairway, add.cell.id = "f103dairway")
f103ddistal <- RenameCells(object = f103ddistal, add.cell.id = "f103ddistal")
f132ddistal <- RenameCells(object = f132ddistal, add.cell.id = "f132ddistal")
f59dairway <- RenameCells(object = f59dairway, add.cell.id = "f59dairway")
f59ddistal <- RenameCells(object = f59ddistal, add.cell.id = "f59ddistal")
f80dairway  <- RenameCells(object = f80dairway, add.cell.id = "f80dairway")
f80ddistal <- RenameCells(object = f80ddistal, add.cell.id = "f80ddistal")


f103dairway[["percent.mt"]] <- PercentageFeatureSet(f103dairway, pattern = "^MT-") 
f103ddistal[["percent.mt"]] <- PercentageFeatureSet(f103ddistal, pattern = "^MT-") 
f132ddistal[["percent.mt"]] <- PercentageFeatureSet(f132ddistal, pattern = "^MT-") 
f59dairway[["percent.mt"]] <- PercentageFeatureSet(f59dairway, pattern = "^MT-") 
f59ddistal[["percent.mt"]] <- PercentageFeatureSet(f59ddistal, pattern = "^MT-") 
f80dairway[["percent.mt"]]  <- PercentageFeatureSet(f80dairway, pattern = "^MT-") 
f80ddistal[["percent.mt"]] <- PercentageFeatureSet(f80ddistal, pattern = "^MT-")

fetal.combined <- merge(f103dairway, y = c(f103ddistal, f132ddistal, f59dairway,f59ddistal, f80dairway, f80ddistal), project = "fetallung")

fetal.combined <- subset(fetal.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

fetal.combined <- SplitObject(fetal.combined, split.by = "orig.ident")
for (i in seq_along(fetal.combined)) {
  fetal.combined[[i]] <- NormalizeData(fetal.combined[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(fetal.combined)
for (i in seq_along(along.with = fetal.combined)) {
  fetal.combined[[i]] <- ScaleData(fetal.combined[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(fetal.combined, anchor.features = features, reduction = "rpca", dims = 1:30)
fetal.combined.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(fetal.combined.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fetal.combined.integrated <- CellCycleScoring(fetal.combined.integrated, s.features = s.genes, g2m.features = g2m.genes)
fetal.combined.integrated <- ScaleData(fetal.combined.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
fetal.combined.integrated <- RunPCA(fetal.combined.integrated)
ElbowPlot(fetal.combined.integrated, ndims = 50)
fetal.combined.integrated <- RunUMAP(fetal.combined.integrated, dims = 1:18, reduction.name = "umap", return.model = TRUE)
fetal.combined.integrated <- FindNeighbors(fetal.combined.integrated, dims = 1:18)
fetal.combined.integrated <- FindClusters(fetal.combined.integrated, resolution = 0.5)
DimPlot(fetal.combined.integrated, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)

DefaultAssay(fetal.combined.integrated) <- "RNA"
fetal.combined.integrated <- NormalizeData(fetal.combined.integrated)

fetal.combined.simplified <- subset(fetal.combined.integrated, idents = c(0, 1, 5, 21, 19, 9, 7)) #subset populations of interest: Bud-Tip Progenitors (9), RSPO2+ Mesenchyme (0, 1), Differentiating Epithelium (7,21,19), Myofibroblasts (5)

DefaultAssay(fetal.combined.simplified) <- "integrated"
fetal.combined.simplified <- RunPCA(fetal.combined.simplified, verbose = FALSE)
ElbowPlot(fetal.combined.simplified, ndims = 50)
fetal.combined.simplified <- RunUMAP(fetal.combined.simplified, reduction = "pca", dims = 1:12)
fetal.combined.simplified <- FindNeighbors(fetal.combined.simplified, dims = 1:12)
fetal.combined.simplified <- FindClusters(fetal.combined.simplified, resolution = 0.15) 
DimPlot(fetal.combined.simplified, label = TRUE, pt.size = 2)

fetal.combined.simplified <- subset(fetal.combined.simplified, idents = c(0, 1, 2, 3)) #clean out weird small clusters
DefaultAssay(fetal.combined.simplified) <- "integrated"
fetal.combined.simplified <- RunPCA(fetal.combined.simplified, verbose = FALSE)
ElbowPlot(fetal.combined.simplified, ndims = 50)
fetal.combined.simplified <- RunUMAP(fetal.combined.simplified, reduction = "pca", dims = 1:10)
fetal.combined.simplified <- FindNeighbors(fetal.combined.simplified, dims = 1:10)
fetal.combined.simplified <- FindClusters(fetal.combined.simplified, resolution = 0.15) 
DimPlot(fetal.combined.simplified, label = TRUE, pt.size = 2)                                    

DefaultAssay(fetal.combined.simplified) <- "RNA"
fetal.combined.simplified <- NormalizeData(fetal.combined.simplified)   

fetal.combined.simplified.renamed <- RenameIdents(object = fetal.combined.simplified, '0' = "RSPO2+", '1' = "Myofibroblasts", '2' = "Differentiating Epithelium", '3' = "Bud Tip Progenitors")
fetal.combined.simplified.renamed.cellchat <- createCellChat(fetal.combined.simplified.renamed, group.by = "ident", assay = "RNA")
#CellChat Analysis
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling


fetal.combined.simplified.renamed.cellchat@DB <- CellChatDB.use
fetal.combined.simplified.renamed.cellchat <- subsetData(fetal.combined.simplified.renamed.cellchat)
fetal.combined.simplified.renamed.cellchat <- identifyOverExpressedGenes(fetal.combined.simplified.renamed.cellchat)
fetal.combined.simplified.renamed.cellchat <- identifyOverExpressedInteractions(fetal.combined.simplified.renamed.cellchat)
fetal.combined.simplified.renamed.cellchat <- projectData(fetal.combined.simplified.renamed.cellchat, PPI.human)
fetal.combined.simplified.renamed.cellchat <- computeCommunProb(fetal.combined.simplified.renamed.cellchat, raw.use = FALSE)


fetal.combined.simplified.renamed.cellchat <- filterCommunication(fetal.combined.simplified.renamed.cellchat, min.cells = 10)
fetal.combined.simplified.renamed.cellchat <- computeCommunProbPathway(fetal.combined.simplified.renamed.cellchat, thresh = 0.05)
fetal.combined.simplified.renamed.cellchat <- aggregateNet(fetal.combined.simplified.renamed.cellchat)
fetal.combined.simplified.renamed.cellchat <- netAnalysis_computeCentrality(fetal.combined.simplified.renamed.cellchat, slot.name = "netP")

groupSize <- as.numeric(table(fetal.combined.simplified.renamed.cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(fetal.combined.simplified.renamed.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(fetal.combined.simplified.renamed.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#visualize communication between populations
mat <- fetal.combined.simplified.renamed.cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
fetal.combined.simplified.renamed.cellchat@netP$pathways


fetal.combined.simplified.renamed.cellchat <- netAnalysis_computeCentrality(fetal.combined.simplified.renamed.cellchat, slot.name = "netP")

#Figure 1A
pdf(file.path("./", paste0("Fetal TGFb", ".pdf")), w=6, h=6)
pathways.show <- c("TGFb")
netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.001)
dev.off()
pdf(file.path("./", paste0("Fetal BMP", ".pdf")), w=6, h=6)
pathways.show <- c("BMP")
netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.001)
dev.off()
pdf(file.path("./", paste0("Fetal Legend", ".pdf")), w=6, h=6)

netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.001, edge.weight.max = 100)
dev.off()
pdf(file.path("./", paste0("Fetal WNT", ".pdf")), w=6, h=6)
pathways.show <- c("WNT")
netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,vertex.label.cex = 0.001 )
dev.off()
pdf(file.path("./", paste0("Fetal FGF", ".pdf")), w=6, h=6)
pathways.show <- c("FGF")
netVisual_aggregate(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.label.cex = 0.001)
dev.off()

##Figure S1A
pdf(file.path("./", paste0("Fetal TGFb LR", ".pdf")), w=6, h=6)
pathways.show <- c("TGFb")
netAnalysis_contribution(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show, font.size = 14, font.size.title = 16)
dev.off()

pdf(file.path("./", paste0("Fetal BMP LR", ".pdf")), w=6, h=6)
pathways.show <- c("BMP")
netAnalysis_contribution(fetal.combined.simplified.renamed.cellchat, signaling = pathways.show,font.size = 14, font.size.title = 16)
dev.off()

#Figure 1B 
#process each distal fetal dataset separately to identify and extract budtips

VlnPlot(f103ddistal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f103ddistal <- subset(f103ddistal, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
f103ddistal <- CellCycleScoring(f103ddistal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
f103ddistal <- SCTransform(f103ddistal, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
f103ddistal <- RunPCA(f103ddistal, features = VariableFeatures(object = f103ddistal))
ElbowPlot(f103ddistal, ndims = 50)
f103ddistal <- RunUMAP(f103ddistal, dims = 1:16, reduction.name = "umap", return.model = TRUE)
f103ddistal <- FindNeighbors(f103ddistal, dims = 1:16)
f103ddistal <- FindClusters(f103ddistal, resolution = 0.5)
DimPlot(f103ddistal, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
pdf(file.path("./", paste0("FullData Compartment", ".pdf")), w=11, h=8.5)
FeaturePlot(f103ddistal, features = compartment, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f103ddistal, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f103ddistal, label = TRUE, pt.size = 0.5)
dev.off()
FeaturePlot(f103ddistal, features = compartment, pt.size = 0.2, order = TRUE)

f103ddistal.epithelial <- subset(f103ddistal, idents = c(6, 7))
DefaultAssay(f103ddistal.epithelial) <- "SCT"
f103ddistal.epithelial <- RunPCA(f103ddistal.epithelial, verbose = FALSE)
ElbowPlot(f103ddistal.epithelial, ndims = 50)
f103ddistal.epithelial <- RunUMAP(f103ddistal.epithelial, reduction = "pca", dims = 1:10)
f103ddistal.epithelial <- FindNeighbors(f103ddistal.epithelial, dims = 1:10)
f103ddistal.epithelial <- FindClusters(f103ddistal.epithelial, resolution = 0.5)
DimPlot(f103ddistal.epithelial, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)



pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f103ddistal.epithelial, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()


DefaultAssay(f103ddistal.epithelial) <- "RNA"
f103ddistal.epithelial <- NormalizeData(f103ddistal.epithelial)
pdf(file.path("./", paste0("Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f103ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()

zero.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 3, min.pct = 0.25)
four.markers <- FindMarkers(f103ddistal.epithelial, ident.1 = 4, min.pct = 0.25)

write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")
write.csv(four.markers, "Cluster4.csv")

f103ddistal.budtip.cellids <- WhichCells(f103ddistal.epithelial, idents = c(1, 4, 2))
f103ddistal.bta.cellids <- WhichCells(f103ddistal.epithelial, idents = c(0))
f103ddistal.budtips <- subset(f103ddistal.epithelial, idents = c(1, 4, 2))
#f132distal
VlnPlot(f132ddistal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f132ddistal <- subset(f132ddistal, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
f132ddistal <- CellCycleScoring(f132ddistal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
f132ddistal <- SCTransform(f132ddistal, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
f132ddistal <- RunPCA(f132ddistal, features = VariableFeatures(object = f132ddistal))
ElbowPlot(f132ddistal, ndims = 50)
f132ddistal <- RunUMAP(f132ddistal, dims = 1:16, reduction.name = "umap", return.model = TRUE)
f132ddistal <- FindNeighbors(f132ddistal, dims = 1:16)
f132ddistal <- FindClusters(f132ddistal, resolution = 0.5)
DimPlot(f132ddistal, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
pdf(file.path("./", paste0("FullData Compartment", ".pdf")), w=11, h=8.5)
FeaturePlot(f132ddistal, features = compartment, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f132ddistal, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f132ddistal, label = TRUE, pt.size = 0.5)
dev.off()
FeaturePlot(f132ddistal, features = compartment, pt.size = 0.2, order = TRUE)

f132ddistal.epithelial <- subset(f132ddistal, idents = c(7, 8))
DefaultAssay(f132ddistal.epithelial) <- "SCT"
f132ddistal.epithelial <- RunPCA(f132ddistal.epithelial, verbose = FALSE)
ElbowPlot(f132ddistal.epithelial, ndims = 50)
f132ddistal.epithelial <- RunUMAP(f132ddistal.epithelial, reduction = "pca", dims = 1:10)
f132ddistal.epithelial <- FindNeighbors(f132ddistal.epithelial, dims = 1:10)
f132ddistal.epithelial <- FindClusters(f132ddistal.epithelial, resolution = 0.5)

pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f132ddistal.epithelial, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()
pdf(file.path("./", paste0("UMAPbySample", ".pdf")), w=11, h=8.5)
DimPlot(f132ddistal.epithelial, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()    

DefaultAssay(f132ddistal.epithelial) <- "RNA"
f132ddistal.epithelial <- NormalizeData(f132ddistal.epithelial)
pdf(file.path("./", paste0("Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f132ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()


zero.markers <- FindMarkers(f132ddistal.epithelial, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(f132ddistal.epithelial, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(f132ddistal.epithelial, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(f132ddistal.epithelial, ident.1 = 3, min.pct = 0.25)


write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")



f132ddistal.budtip.cellids <- WhichCells(f132ddistal.epithelial, idents = c(1,3))
f132ddistal.bta.cellids <- WhichCells(f132ddistal.epithelial, idents = c(0))
f132ddistal.mixed.cellids <- WhichCells(f132ddistal.epithelial, idents = c(2))
f132ddistal.budtip <- subset(f132ddistal.epithelial, idents = c(1, 3))

#80 day distal fetal lung
VlnPlot(f80ddistal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f80ddistal <- subset(f80ddistal, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
f80ddistal <- CellCycleScoring(f80ddistal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
f80ddistal <- SCTransform(f80ddistal, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
f80ddistal <- RunPCA(f80ddistal, features = VariableFeatures(object = f80ddistal))
ElbowPlot(f80ddistal, ndims = 50)
f80ddistal <- RunUMAP(f80ddistal, dims = 1:22, reduction.name = "umap", return.model = TRUE)
f80ddistal <- FindNeighbors(f80ddistal, dims = 1:22)
f80ddistal <- FindClusters(f80ddistal, resolution = 0.5)
DimPlot(f80ddistal, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
pdf(file.path("./", paste0("FullData Compartment", ".pdf")), w=11, h=8.5)
FeaturePlot(f80ddistal, features = compartment, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f80ddistal, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f80ddistal, label = TRUE, pt.size = 0.5)
dev.off()

f80ddistal.epithelial <- subset(f80ddistal, idents = c(7, 13, 10))
DefaultAssay(f80ddistal.epithelial) <- "SCT"
f80ddistal.epithelial <- RunPCA(f80ddistal.epithelial, verbose = FALSE)
ElbowPlot(f80ddistal.epithelial, ndims = 50)
f80ddistal.epithelial <- RunUMAP(f80ddistal.epithelial, reduction = "pca", dims = 1:10)
f80ddistal.epithelial <- FindNeighbors(f80ddistal.epithelial, dims = 1:10)
f80ddistal.epithelial <- FindClusters(f80ddistal.epithelial, resolution = 0.5)
DimPlot(f80ddistal.epithelial, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
FeaturePlot(f80ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)

pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f80ddistal.epithelial, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()
pdf(file.path("./", paste0("UMAPbySample", ".pdf")), w=11, h=8.5)
DimPlot(f80ddistal.epithelial, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()    

DefaultAssay(f80ddistal.epithelial) <- "RNA"
f80ddistal.epithelial <- NormalizeData(f80ddistal.epithelial)
pdf(file.path("./", paste0("Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f80ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()




zero.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 3, min.pct = 0.25)
four.markers <- FindMarkers(f80ddistal.epithelial, ident.1 = 4, min.pct = 0.25)


write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")
write.csv(four.markers, "Cluster4.csv")

f80ddistal.budtip.cellids <- WhichCells(f80ddistal.epithelial, idents = c(0))
f80ddistal.bta.cellids <- WhichCells(f80ddistal.epithelial, idents = c(1))
f80ddistal.budtip <- subset(f80ddistal.epithelial, idents = c(0,3))

##59 day distal fetal lung
VlnPlot(f59ddistal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f59ddistal <- subset(f59ddistal, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
f59ddistal <- CellCycleScoring(f59ddistal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
f59ddistal <- SCTransform(f59ddistal, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
f59ddistal <- RunPCA(f59ddistal, features = VariableFeatures(object = f59ddistal))
ElbowPlot(f59ddistal, ndims = 50)
f59ddistal <- RunUMAP(f59ddistal, dims = 1:22, reduction.name = "umap", return.model = TRUE)
f59ddistal <- FindNeighbors(f59ddistal, dims = 1:22)
f59ddistal <- FindClusters(f59ddistal, resolution = 0.5)
DimPlot(f59ddistal, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
pdf(file.path("./", paste0("FullData Compartment", ".pdf")), w=11, h=8.5)
FeaturePlot(f59ddistal, features = compartment, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f59ddistal, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()
pdf(file.path("./", paste0("FullData Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f59ddistal, label = TRUE, pt.size = 0.5)
dev.off()

f59ddistal.epithelial <- subset(f59ddistal, idents = c(10, 13))
DefaultAssay(f59ddistal.epithelial) <- "SCT"
f59ddistal.epithelial <- RunPCA(f59ddistal.epithelial, verbose = FALSE)
ElbowPlot(f59ddistal.epithelial, ndims = 50)
f59ddistal.epithelial <- RunUMAP(f59ddistal.epithelial, reduction = "pca", dims = 1:10)
f59ddistal.epithelial <- FindNeighbors(f59ddistal.epithelial, dims = 1:10)
f59ddistal.epithelial <- FindClusters(f59ddistal.epithelial, resolution = 1.45)
DimPlot(f59ddistal.epithelial, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)

pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(f59ddistal.epithelial, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()
pdf(file.path("./", paste0("UMAPbySample", ".pdf")), w=11, h=8.5)
DimPlot(f59ddistal.epithelial, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()    

DefaultAssay(f59ddistal.epithelial) <- "RNA"
f59ddistal.epithelial <- NormalizeData(f59ddistal.epithelial)
pdf(file.path("./", paste0("Epithelial Markers", ".pdf")), w=11, h=8.5)
FeaturePlot(f59ddistal.epithelial, features = epithelial.markers, pt.size = 0.2, order = TRUE)
dev.off()


zero.markers <- FindMarkers(f59ddistal.epithelial, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(f59ddistal.epithelial, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(f59ddistal.epithelial, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(f59ddistal.epithelial, ident.1 = 3, min.pct = 0.25)


write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")


f59ddistal.budtip.cellids <- WhichCells(f59ddistal.epithelial, idents = c(1))
f59ddistal.bta.cellids <- WhichCells(f59ddistal.epithelial, idents = c(3))
f59ddistal.budtip <- subset(f59ddistal.epithelial, idents = c(1))

allbudtips <- merge(f59ddistal.budtip, c(f80ddistal.budtip, f103ddistal.budtips, f132ddistal.budtip))

allbudtips$orig.ident <- factor(allbudtips$orig.ident, levels = c("f132ddistal", "f103ddistal", "f80ddistal", "f59ddistal"))
DefaultAssay(allbudtips) <- "RNA"
allbudtips <- NormalizeData(allbudtips)


Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE, group.by = "orig.ident") + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}

##Figure 1B
pdf(file.path("./", paste0("Fetal BudTip Only BMP4 ID2 SOX2 SOX9 dotplot and others", ".pdf")), w=8.0, h=7.5)
Dotplot_Zhiwei_Version(allbudtips, c("BMP2", "BMP4", "BMP5", "ID2" ,"TGFB2", "TGFB3", "SOX2", "SOX9", "SFTPC"))
dev.off()

pdf(file.path("./", paste0("Fetal BudTip Only BMP4 ID2 SOX2 SOX9 dotplot and others Aug 2022", ".pdf")), w=8.0, h=7.5)
Dotplot_Zhiwei_Version(allbudtips, c("BMP2", "BMP4", "BMP5","GDF5", "ID2" ,"TGFB1","TGFB2", "TGFB3", "SOX2", "SOX9", "SFTPC"))
dev.off()

##Figure 1C
rspo2 <- subset(fetal.combined.simplified.renamed, idents = "RSPO2+")
rspo2.59ddistal <- subset(rspo2, subset = orig.ident == "f59ddistal")
rspo2.80ddistal <- subset(rspo2, subset = orig.ident == "f80ddistal")
rspo2.103ddistal <- subset(rspo2, subset = orig.ident == "f103ddistal")
rspo2.132ddistal <- subset(rspo2, subset = orig.ident == "f132ddistal")

rspo2.merge <- merge(rspo2.59ddistal, c(rspo2.80ddistal, rspo2.103ddistal, rspo2.132ddistal))
DefaultAssay(rspo2.merge) <- "RNA"
rspo2.merge <- NormalizeData(rspo2.merge)

rspo2.merge$orig.ident <- factor(rspo2.merge$orig.ident, levels = c("f132ddistal", "f103ddistal", "f80ddistal", "f59ddistal"))

pdf(file.path("./", paste0("Fetal RSPO2 only BMP4 ID2 SOX2 SOX9 dotplot and others Aug 22", ".pdf")), w=8.0, h=7.5)
Dotplot_Zhiwei_Version(rspo2.merge, c("BMP2", "BMP4", "BMP5","GDF5", "ID2" ,"TGFB1", "TGFB2", "TGFB3", "RSPO2", "LIFR"))
dev.off()

#Explant Processing, Figure S1D-G, Figure 1I-K. 

daythree.data <- Read10X(data.dir = "./")
daysix.data <- Read10X(data.dir = "./")
daynine.data <- Read10X(data.dir = "./")
daytwelve.data <- Read10X(data.dir = "./")


daythree <- CreateSeuratObject(counts = daythree.data, project = "daythree", min.cells = 3, min.features = 200)
daysix <- CreateSeuratObject(counts = daysix.data, project = "daysix", min.cells = 3, min.features = 200)
daynine <- CreateSeuratObject(counts = daynine.data, project = "daynine", min.cells = 3, min.features = 200)
daytwelve <- CreateSeuratObject(counts = daytwelve.data, project = "daytwelve", min.cells = 3, min.features = 200)

daythree <- RenameCells(object = daythree, add.cell.id = "daythree")
daysix <- RenameCells(object = daysix, add.cell.id = "daysix")
daynine <- RenameCells(object = daynine, add.cell.id = "daynine")
daytwelve <- RenameCells(object = daytwelve, add.cell.id = "daytwelve")


daythree[["percent.mt"]] <- PercentageFeatureSet(daythree, pattern = "^MT-") 
daysix[["percent.mt"]] <- PercentageFeatureSet(daysix, pattern = "^MT-") 
daynine[["percent.mt"]] <- PercentageFeatureSet(daynine, pattern = "^MT-") 
daytwelve[["percent.mt"]] <- PercentageFeatureSet(daytwelve, pattern = "^MT-")

explant.combined = merge(daythree, y = c(daysix, daynine, daytwelve))
explant.combined <- AddMetaData(explant.combined, "explant", col.name = "group")
VlnPlot(explant.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

explant.combined <- subset(explant.combined, subset = nFeature_RNA > 1000 & nFeature_RNA < 12000 & percent.mt < 10)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
explant.combined <- CellCycleScoring(explant.combined, s.features = s.genes, g2m.features = g2m.genes)
explant.combined.list <- SplitObject(explant.combined, split.by = "orig.ident")
explant.combined.list <- lapply(X = explant.combined.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = explant.combined.list)
explant.combined.list <- lapply(X = explant.combined.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE, vars.to.regress = c("S.Score", "G2M.Score"))
  x <- RunPCA(x, features = features, verbose = TRUE)
})
explant.anchors <- FindIntegrationAnchors(object.list = explant.combined.list, anchor.features = features, reduction = "rpca",  k.anchor = 5)
explant.combined.rpca <- IntegrateData(anchorset = explant.anchors)
DefaultAssay(explant.combined.rpca) <- "integrated"
explant.combined.rpca <- ScaleData(explant.combined.rpca, verbose = FALSE, vars.to.regress = c("S.Score", "G2M.Score"))
explant.combined.rpca <- RunPCA(explant.combined.rpca, npcs = 50, verbose = FALSE)
ElbowPlot(explant.combined.rpca, ndims = 50)
DefaultAssay(explant.combined.rpca) <- "integrated"

explant.combined.rpca <- RunUMAP(explant.combined.rpca, reduction = "pca", dims = 1:18)
explant.combined.rpca <- FindNeighbors(explant.combined.rpca, reduction = "pca", dims = 1:18)
explant.combined.rpca <- FindClusters(explant.combined.rpca, resolution = 0.5)

DimPlot(explant.combined.rpca, label = TRUE, pt.size = 2)

DefaultAssay(explant.combined.rpca) <- "RNA"
explant.combined.rpca <- NormalizeData(explant.combined.rpca)

zero.markers <- FindMarkers(explant.combined.rpca, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(explant.combined.rpca, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(explant.combined.rpca, ident.1 = 2, min.pct = 0.25)
three.markers <- FindMarkers(explant.combined.rpca, ident.1 = 3, min.pct = 0.25)
four.markers <- FindMarkers(explant.combined.rpca, ident.1 = 4, min.pct = 0.25)
five.markers <- FindMarkers(explant.combined.rpca, ident.1 = 5, min.pct = 0.25)
six.markers <- FindMarkers(explant.combined.rpca, ident.1 = 6, min.pct = 0.25)
seven.markers <- FindMarkers(explant.combined.rpca, ident.1 = 7, min.pct = 0.25)
eight.markers <- FindMarkers(explant.combined.rpca, ident.1 = 8, min.pct = 0.25)
nine.markers <- FindMarkers(explant.combined.rpca, ident.1 = 9, min.pct = 0.25)
ten.markers <- FindMarkers(explant.combined.rpca, ident.1 = 10, min.pct = 0.25)
eleven.markers <- FindMarkers(explant.combined.rpca, ident.1 = 11, min.pct = 0.25)
twelve.markers <- FindMarkers(explant.combined.rpca, ident.1 = 12, min.pct = 0.25)
thirteen.markers <- FindMarkers(explant.combined.rpca, ident.1 = 13, min.pct = 0.25)
fourteen.markers <-FindMarkers(explant.combined.rpca, ident.1 = 14, min.pct = 0.25)
write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")
write.csv(three.markers, "Cluster3.csv")
write.csv(four.markers, "Cluster4.csv")
write.csv(five.markers, "Cluster5.csv")
write.csv(six.markers, "Cluster6.csv")
write.csv(seven.markers, "Cluster7.csv")
write.csv(eight.markers, "Cluster8.csv")
write.csv(nine.markers, "Cluster9.csv")
write.csv(ten.markers, "Cluster10.csv")
write.csv(eleven.markers, "Cluster11.csv")
write.csv(twelve.markers, "Cluster12.csv")
write.csv(thirteen.markers, "Cluster13.csv")
write.csv(fourteen.markers, "Cluster14.csv")
explant.combined.rpca.renamed <- RenameIdents(explant.combined.rpca, '0' = "RSPO2+ Mesenchyme", '1' = "AT1", '2' = "Proliferative Epithelium", '3' = "Myofibroblast", '4' = "Pericyte", '5' = "Matrix Fibroblast", '6' = "Proliferative Fibroblast 1", '7' = "AT2", '8' = "Proliferative Fibroblast 2", '9' = "Endothelial", '10' = "Macrophages", '11' = "Airway", '12' = "Lymphatic", '13' = "Mesothelial", '14' = "T-Cell")

#Figure S1D
pdf(file.path("./", paste0("Louvain no label", ".pdf")), w=11, h=8.5)
DimPlot(explant.combined.rpca, reduction = "umap", label = FALSE, pt.size = 2)
dev.off()

#Figure S1E
levels(explant.combined.rpca) <- c("13","12", "9", "10", "14","8", "6", "5", "4","3", "0" ,"7", "2", "1", "11")
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}

pdf(file.path("./", paste0("Full Data Dotplot ", ".pdf")), w=18, h=11)
Dotplot_Zhiwei_Version(explant.combined.rpca, c("EPCAM", "CDH1", "SOX2", "AQP4",  "AGER", "SFTPC" , "LIFR","RSPO2", "PDGFRA", "ACTA2", "THY1", "PDGFRB","MGP", "TOP2A", "PCNA", "PTPRC", "CD3G","CD86", "PECAM1", "CA4","PROX1", "WT1", "UPK3B"))
dev.off()

#Data for Figure S1F (graph made in prism)
write.csv(table(Idents(explant.combined.rpca), explant.combined.rpca$orig.ident))


#Figure 1I
explant.epithelium.rpca <- subset(explant.combined.rpca, idents = c(1,7,2,11))

DefaultAssay(explant.epithelium.rpca) <- "integrated"
explant.epithelium.rpca <- RunPCA(explant.epithelium.rpca, verbose = FALSE)
ElbowPlot(explant.epithelium.rpca, ndims = 50)
explant.epithelium.rpca <- RunUMAP(explant.epithelium.rpca, reduction = "pca", dims = 1:7)
explant.epithelium.rpca <- FindNeighbors(explant.epithelium.rpca, dims = 1:7)
explant.epithelium.rpca <- FindClusters(explant.epithelium.rpca, resolution = 0.5)
DimPlot(explant.epithelium.rpca, label = TRUE, pt.size = 2)
epithelium.cluster5.cellids  <- WhichCells(epithelium.cluster5.cellids, idents = c(5)) ##remove cluster 5 which has low NKX2-1 expression

explant.epithelium.rpca <- subset(explant.epithelium.rpca, idents = c(0,1,2,3,4,6))
explant.epithelium.rpca <- RunPCA(explant.epithelium.rpca, verbose = FALSE)
ElbowPlot(explant.epithelium.rpca, ndims = 50)
explant.epithelium.rpca <- RunUMAP(explant.epithelium.rpca, reduction = "pca", dims = 1:7)
explant.epithelium.rpca <- FindNeighbors(explant.epithelium.rpca, dims = 1:7)
explant.epithelium.rpca <- FindClusters(explant.epithelium.rpca, resolution = 0.4)
pdf(file.path("./", paste0("Louvain no Label", ".pdf")), w=11, h=8.5)
DimPlot(explant.epithelium.rpca, reduction = "umap",  pt.size = 2)
dev.off()

#Figure 1J
DefaultAssay(explant.epithelium.rpca) <- "RNA"
explant.epithelium.rpca <- NormalizeData(explant.epithelium.rpca)

ExplantEpithelialCluster.markers <- c("SFTPC","SFTPB", "NAPSA", "ABCA3", "SLC34A2" ,"DMBT1","SFTPA1", "HOPX", "PDPN", "AGER", "RTKN2", "SPOCK2" ,"SCGB3A2", "SOX2","FOXJ1","TP63", "SCGB1A1", "SPDEF", "SOX9", "TESC", "SOX11", "CPM", "TOP2A", "MKI67")

levels(explant.epithelium.rpca) <- c("0", "5", "3", "2", "4", "1")
pdf(file.path("./", paste0("Epithelial Marker DotPlot v2", ".pdf")), w=17, h=5)
Dotplot_Zhiwei_Version(explant.epithelium.rpca, ExplantEpithelialCluster.markers)
dev.off()


#Figure 1K
pdf(file.path("./", paste0("Ligands of interest dotplot.pdf")), w=6, h=5)
Dotplot_Zhiwei_Version(explant.epithelium.rpca, c("BMP2", "BMP4", "BMP5", "ID2" ,"TGFB2", "TGFB3")) 
dev.off()



#Figure S1G
#pull AT2 cells from Travaglini et al., 2020 data
krasnow.data <- read.table(file=paste0("krasnow_hlca_10x_UMIs.csv"), sep = ",", row.names = 1, header = TRUE) #import data from Travaglini et al., 2020
krasnow <- CreateSeuratObject(krasnow.data, project = "krasnow", min.cells = 3, min.features = 200)
krasnow <- subset(krasnow, subset = nFeature_RNA > 500 & nFeature_RNA < 10000)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
krasnow <- CellCycleScoring(krasnow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
krasnow <- SCTransform(krasnow, vars.to.regress = c("S.Score", "G2M.Score"))
krasnow <- RunPCA(krasnow, features = VariableFeatures(object = krasnow))
ElbowPlot(krasnow, ndims = 50)
krasnow <- RunUMAP(krasnow, dims = 1:25, reduction.name = "umap", return.model = TRUE)
krasnow <- FindNeighbors(krasnow, dims = 1:25)
krasnow <- FindClusters(krasnow, resolution = 0.5)

#extract epithelium from full data
krasnow.epithelial <- subset(krasnow, idents = c(16, 22, 14, 21, 4))
DefaultAssay(krasnow.epithelial) <- "SCT"
krasnow.epithelial <- RunPCA(krasnow.epithelial, verbose = FALSE)
ElbowPlot(krasnow.epithelial, ndims = 50)
krasnow.epithelial <- RunUMAP(krasnow.epithelial, reduction = "pca", dims = 1:20)
krasnow.epithelial <- FindNeighbors(krasnow.epithelial, dims = 1:20)
krasnow.epithelial <- FindClusters(krasnow.epithelial, resolution = 0.5)
DimPlot(krasnow.epithelial, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)

#extract at2 from epithelium
krasnow.at2 <- subset(krasnow.epithelial, idents = c(0, 1, 3))
krasnow.at2 <- AddMetaData(krasnow.at2, "Primary", col.name = "Sample")
#extract at2s from explant epithelium
explant.at2 <- subset(explant.epithelium.rpca, idents = c(1, 4))
explant.at2 <- AddMetaData(explant.at2, "Explant", col.name = "Sample")

#extract btps from fetal epithelium
fetal.btps <- merge(f59ddistal.budtip, f80ddistal.budtip, f103ddistal.budtips, f132distal.budtip)
fetal.btps <- AddMetaData(explant.at2, "Fetal", col.name = "Sample")
#merge adult AT2s with BTPs from fetal data and AT2s from explants
at2.combined <- merge(krasnow.at2, explant.at2, fetal.btps)

Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE, group.by = "Sample") + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}


pdf(file.path("./", paste0("Fetal Explant Adult ATII marker DotPlot", ".pdf")), w=8, h=5)
Dotplot_Zhiwei_Version(at2.combined, c("SFTPC", "SFTPB",  "SLC34A2","MUC1", "NAPSA","LAMP3", "ABCA3", "SFTPA1"))
dev.off()

##Figure 3
frum.aec2.curated <- c("SFTPC", "SFTPB", "SFTPA1", "LAMP3", "NAPSA", "SLC34A2")
frum.btp.curated <- frum.btp.curated <- c("SOX9", "SOX2", "TESC", "SOX11", "CPM", "NKX2-1")




twoefabday21.data <- Read10X(data.dir = "./")

threef.data <- Read10X(data.dir = "./") #BTP-Organoids under maintenance condition, 'day 0'

twoefabday1.data <- Read10X(data.dir = "./")
twoefabday6.data <- Read10X(data.dir = "./")


threef <- CreateSeuratObject(counts = threef.data, project = "threef", min.cells = 3, min.features = 200)
twoefabday21 <- CreateSeuratObject(counts = twoefabday21.data, project = "twoefabday21", min.cells = 3, min.features = 200)
twoefabday1 <- CreateSeuratObject(counts = twoefabday1.data, project = "2fabday1", min.cells = 3, min.features = 200)
twoefabday6 <- CreateSeuratObject(counts = twoefabday6.data, project = "2fabday6", min.cells = 3, min.features = 200)
threef <- RenameCells(object = threef, add.cell.id = "3F")
threef[["percent.mt"]] <- PercentageFeatureSet(threef, pattern = "^MT-")
twoefabday21 <- RenameCells(object = twoefabday21, add.cell.id = "2FABDAY21")

twoefabday1 <- RenameCells(object = twoefabday1, add.cell.id = "2FABDAY1")
twoefabday6 <- RenameCells(object = twoefabday6, add.cell.id = "2FABDAY6")
twoefabday1[["percent.mt"]] <- PercentageFeatureSet(twoefabday1, pattern = "^MT-")
twoefabday6[["percent.mt"]] <- PercentageFeatureSet(twoefabday6, pattern = "^MT-")
twoefabday21[["percent.mt"]] <- PercentageFeatureSet(twoefabday21, pattern = "^MT-")


threef <- AddMetaData(threef, "Day 0", col.name = "group") 
twoefabday6 <- AddMetaData(twoefabday6, "Day 6", col.name = "group")
twoefabday1 <- AddMetaData(twoefabday1, "Day 1", col.name = "group")
twoefabday21 <- AddMetaData(twoefabday21, "Day 21", col.name = "group")

threef <- AddMetaData(threef, 0, col.name = "timepoint") 
twoefabday6 <- AddMetaData(twoefabday6, 1, col.name = "timepoint")
twoefabday1 <- AddMetaData(twoefabday1, 6, col.name = "timepoint")
twoefabday21 <- AddMetaData(twoefabday21, 21, col.name = "timepoint")


twoefab <- merge(threef, y = c(twoefabday1, twoefabday6, twoefabday21), project = "2FAB")

VlnPlot(twoefab, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
twoefab <- subset(twoefab, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)


DefaultAssay(twoefab) <- "RNA"
twoefab <- NormalizeData(twoefab, verbose = TRUE)

#Figure 3C
twoefab$group <- factor(twoefab$group, levels = c("Day 0", "Day 1", "Day 6", "Day 21"))
p <- VlnPlot(twoefab, features = frum.aec2.curated, group.by = "group",  ncol = 6, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN AT2 marker curated set", ".pdf")), w=5.7, h=1.8)
p
dev.off()

p <- VlnPlot(twoefab, features = frum.btp.curated, group.by = "group",  ncol = 6, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN BTP marker curated set + NKX21", ".pdf")), w=5.7, h=1.8)
p
dev.off()


##Figure S3
###Process each 2FAB timepoint individually and identify the best AT2s on the basis of AT2 marker expression
####Day 1
twoefabday1 <- CreateSeuratObject(counts = twoefabday1.data, project = "twoefabday1", min.cells = 3, min.features = 200)
twoefabday1 <- AddMetaData(twoefabday1, "2FAB Day 1", col.name = "group")
twoefabday1[["percent.mt"]] <- PercentageFeatureSet(twoefabday1, pattern = "^MT-")
twoefabday1 <- RenameCells(object = twoefabday1, add.cell.id = "2FABDAY1")
VlnPlot(twoefabday1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
twoefabday1 <- subset(twoefabday1, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
twoefabday1 <- CellCycleScoring(twoefabday1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
twoefabday1 <- SCTransform(twoefabday1, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
twoefabday1 <- RunPCA(twoefabday1, features = VariableFeatures(object = twoefabday1))

ElbowPlot(twoefabday1, ndims = 50)
DefaultAssay(twoefabday1) <- "SCT"
twoefabday1 <- RunUMAP(twoefabday1, dims = 1:18)
twoefabday1 <- FindNeighbors(twoefabday1, dims = 1:18)
twoefabday1 <- FindClusters(twoefabday1, resolution = 0.5)
DimPlot(twoefabday1, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)


DefaultAssay(twoefabday1) <- "RNA"
twoefabday1 <- NormalizeData(twoefabday1)

####Day 6
twoefabday6 <- CreateSeuratObject(counts = twoefabday6.data, project = "twoefabday6", min.cells = 3, min.features = 200)
twoefabday6 <- AddMetaData(twoefabday6, "2FAB Day 6", col.name = "group")
twoefabday6[["percent.mt"]] <- PercentageFeatureSet(twoefabday6, pattern = "^MT-")
twoefabday6 <- RenameCells(object = twoefabday6, add.cell.id = "2FABDAY6")
VlnPlot(twoefabday6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
twoefabday6 <- subset(twoefabday6, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
twoefabday6 <- CellCycleScoring(twoefabday6, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
twoefabday6 <- SCTransform(twoefabday6, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
twoefabday6 <- RunPCA(twoefabday6, features = VariableFeatures(object = twoefabday6))

ElbowPlot(twoefabday6, ndims = 50)
DefaultAssay(twoefabday6) <- "SCT"
twoefabday6 <- RunUMAP(twoefabday6, dims = 1:18)
twoefabday6 <- FindNeighbors(twoefabday6, dims = 1:18)
twoefabday6 <- FindClusters(twoefabday6, resolution = 0.5)
DimPlot(twoefabday6, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)


DefaultAssay(twoefabday6) <- "RNA"
twoefabday6 <- NormalizeData(twoefabday6)
####Day 21
twoefabday21 <- CreateSeuratObject(counts = twoefabday21.data, project = "twoefabday21", min.cells = 3, min.features = 200)
twoefabday21 <- AddMetaData(twoefabday21, "2FAB Day 21", col.name = "group")
twoefabday21[["percent.mt"]] <- PercentageFeatureSet(twoefabday21, pattern = "^MT-")
twoefabday21 <- RenameCells(object = twoefabday21, add.cell.id = "2FABDAY21")
VlnPlot(twoefabday21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
twoefabday21 <- subset(twoefabday21, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
twoefabday21 <- CellCycleScoring(twoefabday21, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
twoefabday21 <- SCTransform(twoefabday21, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
twoefabday21 <- RunPCA(twoefabday21, features = VariableFeatures(object = twoefabday21))

ElbowPlot(twoefabday21, ndims = 50)
DefaultAssay(twoefabday21) <- "SCT"
twoefabday21 <- RunUMAP(twoefabday21, dims = 1:18)
twoefabday21 <- FindNeighbors(twoefabday21, dims = 1:18)
twoefabday21 <- FindClusters(twoefabday21, resolution = 0.5)
DimPlot(twoefabday21, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

DefaultAssay(twoefabday21) <- "RNA"
twoefabday21 <- NormalizeData(twoefabday21)

#Figure S3B
pdf(file.path("./", paste0("Day 1 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday1, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 1 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday1, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 1 SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday1, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 1 TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday1, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("Day 1 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(twoefabday1, pt.size = 2)
dev.off()
pdf(file.path("./", paste0("Day 1 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(twoefabday1, cells.highlight = twoefabday1.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()


#Figure S3C
pdf(file.path("./", paste0("Day 6 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday6, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 6 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday6, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 6 SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday6, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 6 TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday6, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("Day 6 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(twoefabday6, pt.size = 2)
dev.off()
pdf(file.path("./", paste0("Day 6 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(twoefabday6, cells.highlight = twoefabday6.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()

#Figure S3D

pdf(file.path("./", paste0("Day 21 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday21, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday21, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday21, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabday21, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(twoefabday21, pt.size = 2)
dev.off()
pdf(file.path("./", paste0("Day 21Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(twoefabday21, cells.highlight = twoefabday21.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()

##Figure 3 Continued
###Integrate BTP-Organoids (Day 0) with Day 1, 6, and 21 
twoefab.fastmnn <- twoefab
set.seed(888)
DefaultAssay(twoefab.fastmnn) <- "RNA"
twoefab.fastmnn <- NormalizeData(twoefab.fastmnn,
                                 normalization.method = "LogNormalize", scale.factor =10000)
twoefab.fastmnn <-  FindVariableFeatures(twoefab.fastmnn)
twoefab.fastmnn <- RunFastMNN(object.list = SplitObject(twoefab.fastmnn, split.by = "orig.ident"))
ElbowPlot(twoefab.fastmnn, ndims = 50)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
twoefab.fastmnn <- CellCycleScoring(twoefab.fastmnn, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
twoefab.fastmnn <- RunUMAP(twoefab.fastmnn, dims = 1:30,reduction = "mnn")
twoefab.fastmnn <- FindNeighbors(twoefab.fastmnn, dims = 1:30,reduction = "mnn")
twoefab.fastmnn <- FindClusters(twoefab.fastmnn, resolution = 0.3, algorithm = 1)
twoefab.fastmnn$group <- factor(twoefab.fastmnn$group, levels = c("Day 0", "Day 1", "Day 6", "Day 21"))

DimPlot(twoefab.fastmnn, group.by = c("group", "ident"), label = TRUE, pt.size = 2)

##PrctCellExpringGene and calc_helper function from : https://github.com/satijalab/seurat/issues/371#issuecomment-486384854
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}
##Data for Figure 3D (graph made in prism)
PrctCellExpringGene(twoefab.fastmnn ,genes =c("SFTPC", "SFTPA1", "SFTPB"), group.by = "group")

##Figure 3E
pdf(file.path("./", paste0("2FAB FastMNN by Phase", ".pdf")), w=11, h=8.5)
DimPlot(twoefab.fastmnn, reduction = "umap",  pt.size = 2, group.by = "Phase")
dev.off()
##Data for Figure 3F (graph made in prism)
write.csv(table(twoefab.fastmnn$group, twoefab.fastmnn$Phase), "RunbyPhase.csv")
##Figure 3G
pdf(file.path("./", paste0("2FAB FastMNN by Sample", ".pdf")), w=11, h=8.5)
DimPlot(twoefab.fastmnn, reduction = "umap",  pt.size = 2, group.by = "group")
dev.off()
##Figure 3H
DefaultAssay(twoefab.fastmnn) <- "RNA"
twoefab.fastmnn <- NormalizeData(twoefab.fastmnn)
pdf(file.path("./", paste0("2FAB FASTMNN SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefab.fastmnn, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("2FAB FASTMNN SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefab.fastmnn, features = "SFTPB", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("2FAB FASTMNN SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefab.fastmnn, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("2FAB FASTMNN LAMP3", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefab.fastmnn, features = "LAMP3", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("2FAB FASTMNN KRT5", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefab.fastmnn, features = "KRT5", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
##Figure 3I
pdf(file.path("./", paste0("2FAB FastMNN", ".pdf")), w=11, h=8.5)
DimPlot(twoefab.fastmnn, reduction = "umap",  pt.size = 2, label = FALSE)
dev.off()
##Figure 3J
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}



levels(twoefab.fastmnn) <- c("6", "1", "0", "5", "3" ,"2", "4")

twoefab.fastmnn.markers <- c("SOX9", "SOX2", "TESC", "SFTPC","NAPSA", "SLC34A2","HOPX", "AGER", "PDPN", "TP63", "FOXJ1" , "SCGB1A1","SPDEF", "ASCL1", "CHGA", "TOP2A", "PCNA")
pdf(file.path("./", paste0("2FAB FASTMNN Dotplot", ".pdf")), w=12.5, h=5.75)
Dotplot_Zhiwei_Version(twoefab.fastmnn, twoefab.fastmnn.markers)
dev.off()

##Figure 3K
twoefab.fastmnn <- FindClusters(twoefab.fastmnn, resolution = 1, algorithm = 1)

sce.twoefab.nosct <- as.SingleCellExperiment(twoefab.fastmnn)
sce.twoefab.nosct <- slingshot(sce.twoefab.nosct, clusterLabels = sce.twoefab.nosct@colData@listData[["ident"]]
                               , reducedDim = "UMAP",
                               start.clus = "6", allow.breaks = FALSE, stretch = 0, omega = TRUE, extend = "n")
lnes.twoefab <- getLineages(reducedDim(sce.twoefab.nosct,"UMAP"), sce.twoefab.nosct@colData@listData[["ident"]], start.clus = "6")
pt.twoefab.nosct <- slingPseudotime(sce.twoefab.nosct, na=TRUE)
pt.twoefab.nosct[is.na(pt.twoefab.nosct)] <- 0

colors <- pal[cut(pt.twoefab.nosct[,2], breaks = 100)]

pdf(file.path("./", paste0("2FAB Psuedotime with Lineage", ".pdf")), w=11, h=8.5)

plot(reducedDim(sce.twoefab.nosct, "UMAP"), col = colors, pch = 16, cex = 1)

plot(SlingshotDataSet(sce.twoefab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 1, add = TRUE)
plot(SlingshotDataSet(sce.twoefab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 3, add = TRUE)
plot(SlingshotDataSet(sce.twoefab.nosct), lwd =6, type = 'lineage', col = c("gray"), linInd = 5, add = TRUE)
plot(SlingshotDataSet(sce.twoefab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 4, add = TRUE)
plot(SlingshotDataSet(sce.twoefab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 6, add = TRUE)
plot(SlingshotDataSet(sce.twoefab.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 7, add = TRUE)
plot(SlingshotDataSet(sce.twoefab.nosct), lwd=8, type = 'lineage', col = c("red"), linInd = 2, add = TRUE)

dev.off()
##Figure 3L
pdf(file.path("./", paste0("2FAB FastMNN All Highlight", ".pdf")), w=11, h=8.5)
DimPlot(twoefab.fastmnn, cells.highlight = list(twoefabday21.bestat2.cellids, twoefabday6.bestat2.cellids, twoefabday1.bestat2.cellids, threef.cellids), cols.highlight = c("purple", "blue", "green", "yellow"), pt.size = 2, sizes.highlight = 2, order = TRUE)
dev.off()
##Figure 3M
##Compare Primary AT2s (see Figure 5) to BTP organoids to identify genes expected to onset during alveolar differentiation

threef <- CreateSeuratObject(counts = threef.data, project = "threef", min.cells = 3, min.features = 200)
threef <- RenameCells(object = threef, add.cell.id = "3F")
threef[["percent.mt"]] <- PercentageFeatureSet(threef, pattern = "^MT-")
adultnof.data  <- Read10X(data.dir = "./") ##see figure 5
adultnof <- CreateSeuratObject(counts = threef.data, project = "adult", min.cells = 3, min.features = 200)
adultnof <- RenameCells(object = adultnof, add.cell.id = "adult")
adultnof[["percent.mt"]] <- PercentageFeatureSet(adultnof, pattern = "^MT-")
adultnof.bestat2 <- subset(adultnof, cells = adultnof.bestat2.cellids)

btptoat2.invitro <-merge(threef, y = adultnof.bestat2)
btptoat2.invitro <- subset(btptoat2.invitro, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 20)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
btptoat2.invitro <- CellCycleScoring(btptoat2.invitro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
btptoat2.invitro <- SCTransform(btptoat2.invitro, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
setwd("~/scRNAseq/Output/2021-10-13 HT313 Long Term Alveolar Differentiation/2021-10-18 BTPtoAT2InVitro")
btptoat2.invitro <- RunPCA(btptoat2.invitro, features = VariableFeatures(object = btptoat2.invitro))
ElbowPlot(btptoat2.invitro, ndims = 50)
DefaultAssay(btptoat2.invitro) <- "SCT"

btptoat2.invitro <- RunUMAP(btptoat2.invitro, dims = 1:15)
btptoat2.invitro <- FindNeighbors(btptoat2.invitro, dims = 1:15)
btptoat2.invitro <- FindClusters(btptoat2.invitro, resolution = 0.1)
DimPlot(btptoat2.invitro, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

DefaultAssay(btptoat2.invitro) <- "RNA"
btptoat2.invitro <- NormalizeData(btptoat2.invitro, verbose = FALSE)

Adultat2vBTP <- FindMarkers(btptoat2.invitro, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
write.csv(Adultat2vBTP, "In vitro Diff.csv")
#top 200 primary AT2 organoids enirched from primary AT2 organoid v BTP organoid comparison
invitro.at2.genes.gained.nolist <- list(c("SFTPC",
                                          "SFTPA1",
                                          "SFTPA2",
                                          "SFTPB",
                                          "NAPSA",
                                          "SPINK5",
                                          "SLC34A2",
                                          "SERPIND1",
                                          "HPGD",
                                          "SCGB3A1",
                                          "PIGR",
                                          "AQP1",
                                          "LRRK2",
                                          "SFTA2",
                                          "HLA-DRA",
                                          "HHIP",
                                          "CEACAM6",
                                          "CD74",
                                          "AQP5",
                                          "VEPH1",
                                          "NUPR1",
                                          "FTL",
                                          "SFTPD",
                                          "SERPINA1",
                                          "MFSD2A",
                                          "HLA-B",
                                          "LAMP3",
                                          "MT-ND4L",
                                          "ADGRF5",
                                          "HLA-DPA1",
                                          "HLA-DPB1",
                                          "SLPI",
                                          "S100A9",
                                          "CD36",
                                          "SUSD2",
                                          "DBI",
                                          "XIST",
                                          "TMEM213",
                                          "FASN",
                                          "ALOX15B",
                                          "CTSH",
                                          "CCND2",
                                          "C3",
                                          "MALL",
                                          "ALPL",
                                          "RASGRF1",
                                          "SLC22A31",
                                          "S100A14",
                                          "HOPX",
                                          "SFTA1P",
                                          "PHLDA2",
                                          "CYB5A",
                                          "CTSD",
                                          "LMO3",
                                          "CD59",
                                          "C19orf33",
                                          "SLC39A8",
                                          "C2",
                                          "MT-CO2",
                                          "CACNA2D2",
                                          "LGI3",
                                          "SDR16C5",
                                          "HLA-DRB1",
                                          "MICAL2",
                                          "TFPI",
                                          "TNC",
                                          "S100A6",
                                          "POLR2L",
                                          "KRT7",
                                          "SNX30",
                                          "HIP1",
                                          "SFRP5",
                                          "ABCA3",
                                          "CAT",
                                          "MYO1B",
                                          "IL18",
                                          "MT-ND5",
                                          "MPZL2",
                                          "CHI3L1",
                                          "DUOX1",
                                          "CDC25B",
                                          "SDC4",
                                          "AP000357.2",
                                          "SFN",
                                          "RPS26",
                                          "SELENOW",
                                          "B2M",
                                          "SELENOP",
                                          "ARPC1B",
                                          "C1orf116",
                                          "FTH1",
                                          "CSTB",
                                          "MT-ND3",
                                          "LPCAT1",
                                          "DRAM1",
                                          "ETV1",
                                          "NFIC",
                                          "NNMT",
                                          "NQO1",
                                          "NPC2",
                                          "ITGB6",
                                          "RAB27A",
                                          "HSPH1",
                                          "EPHX1",
                                          "LGALS3",
                                          "MSMO1",
                                          "FBP1",
                                          "GPX4",
                                          "BCL2L1",
                                          "NFIX",
                                          "HLA-C",
                                          "BRI3",
                                          "ICAM1",
                                          "SCD",
                                          "HLA-DMA",
                                          "C16orf89",
                                          "FGGY",
                                          "AK1",
                                          "FOLR1",
                                          "MSN",
                                          "ISCU",
                                          "HPCAL1",
                                          "LIPH",
                                          "MGST1",
                                          "CAPN8",
                                          "IL1R1",
                                          "PLXND1",
                                          "MEGF6",
                                          "MBNL1",
                                          "GGT5",
                                          "PDXK",
                                          "CXCL17",
                                          "CA2",
                                          "MT-ND4",
                                          "HIF1A",
                                          "TMEM238",
                                          "UQCR10",
                                          "MT-ATP6",
                                          "SCNN1A",
                                          "HHIP-AS1",
                                          "STC1",
                                          "CBR1",
                                          "MTUS1",
                                          "ACSL4",
                                          "RAB27B",
                                          "SNHG7",
                                          "LY6E",
                                          "MEGF9",
                                          "CEBPA",
                                          "SLC6A20",
                                          "SNX25",
                                          "ANXA1",
                                          "CYP1B1",
                                          "SDC1",
                                          "TGFBR2",
                                          "ADIPOR1",
                                          "DCXR",
                                          "CAPN2",
                                          "MMP28",
                                          "TSTD1",
                                          "ACSS2",
                                          "MT-CO1",
                                          "GGTLC1",
                                          "DCBLD2",
                                          "TGM2",
                                          "HSPB8",
                                          "CPM",
                                          "PID1",
                                          "RHOBTB2",
                                          "EPDR1",
                                          "QPRT",
                                          "SULT1A1",
                                          "ADGRF1",
                                          "TAOK3",
                                          "MID1IP1",
                                          "AQP4",
                                          "ALDH2",
                                          "DGKD",
                                          "ZMAT3",
                                          "CD9",
                                          "GALNT10",
                                          "MEG3",
                                          "IFI16",
                                          "CREG1",
                                          "HLA-A",
                                          "SELENBP1",
                                          "DHCR24",
                                          "ARHGDIB",
                                          "TST",
                                          "STOM",
                                          "ACSL1",
                                          "SLC66A1L",
                                          "CTSS",
                                          "FAH",
                                          "PARM1",
                                          "COX17",
                                          "TMEM125",
                                          "HLA-DOA",
                                          "MCUR1"))

#gene set module scoring of BTP-organoids and all TGFbi/BMPa treatment timepoints
twoefab <- merge(threef, y = c(twoefabday1, twoefabday6, twoefabday21), project = "2FAB")
VlnPlot(twoefab, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
twoefab <- subset(twoefab, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
DefaultAssay(twoefab) <- "RNA"
twoefab <- NormalizeData(twoefab, verbose = TRUE)
y <- twoefab
y$group <- factor(y$group, levels = c("Day 0", "Day 1", "Day 6", "Day 21"))
y <- AddModuleScore(y, features = invitro.at2.genes.gained.nolist, name = "At2invitro", search = TRUE, assay = "RNA")
pdf(file.path("./", paste0("In Vitro Alveolar Differentation Scoring")), w=11, h=8.5)
print(VlnPlot(y, features = "At2invitro1", pt.size = 0, group.by = "group"))
dev.off()

##Figure 5
twoefabnof.data <- Read10X(data.dir = "./") #TGFbi/BMPa induced organoids for 21 days, expanded in primary AT2 organoid media without FGF10 for 120 days

twoefabnof <- CreateSeuratObject(counts = twoefabnof.data, project = "twoefabnof", min.cells = 3, min.features = 200)
twoefabnof <- RenameCells(object = twoefabnof, add.cell.id = "twoefabnof")
twoefabnof[["percent.mt"]] <- PercentageFeatureSet(twoefabnof, pattern = "^MT-")

VlnPlot(twoefabnof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
twoefabnof <- subset(twoefabnof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)


twoefabnof <- SCTransform(twoefabnof)
twoefabnof <- RunPCA(twoefabnof)
ElbowPlot(twoefabnof, ndims = 50)
twoefabnof <- RunUMAP(twoefabnof, dims = 1:16, reduction.name = "umap")
twoefabnof <- FindNeighbors(twoefabnof, dims = 1:16)
twoefabnof <- FindClusters(twoefabnof, resolution = 0.5)
DimPlot(twoefabnof, label = TRUE)

##Figure 5A
pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(twoefabnof, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(twoefabnof) <- "RNA"
twoefabnof <- NormalizeData(twoefabnof)
twoefabnof  <- CellCycleScoring(twoefabnof, s.features = s.genes, g2m.features = g2m.genes)

PrctCellExpringGene(twoefabnof, "MUC5AC")
FeaturePlot(twoefabnof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))

pdf(file.path("./", paste0("Day 140 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabnof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 140 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabnof, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 140 SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabnof, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 140 TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabnof, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 140 MUC5AC", ".pdf")), w=11, h=8.5)
FeaturePlot(twoefabnof, features = "MUC5AC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

#Figure 5B
pdf(file.path("./", paste0("AT2 Marker Marker Dotplot", ".pdf")), w=11, h=8.5)
Dotplot_Zhiwei_Version(twoefabnof, c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "PGC", "NAPSA", "SFTA2", "SFTA3", "MUC5AC"))
dev.off()

#Figure S5
adultnof.data  <- Read10X(data.dir = "./") #primary AT2 organoids cultured in primary AT2 organoid media for 90 days and then switched to primary AT2 organoid with no FGF10 for 30 days before sequencing

adultnof <- CreateSeuratObject(counts = adultnof.data, project = "adultnof", min.cells = 3, min.features = 200)
adultnof <- RenameCells(object = adultnof, add.cell.id = "adultnof")
adultnof[["percent.mt"]] <- PercentageFeatureSet(adultnof, pattern = "^MT-")

VlnPlot(adultnof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
adultnof <- subset(adultnof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)


adultnof <- SCTransform(adultnof)
adultnof <- RunPCA(adultnof)
ElbowPlot(adultnof, ndims = 50)
adultnof <- RunUMAP(adultnof, dims = 1:16, reduction.name = "umap")
adultnof <- FindNeighbors(adultnof, dims = 1:16)
adultnof <- FindClusters(adultnof, resolution = 0.5)
DimPlot(adultnof, label = TRUE)

##Figure S5A
pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(adultnof, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(adultnof) <- "RNA"
adultnof <- NormalizeData(adultnof)



pdf(file.path("./", paste0("Adult NOF SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult NOF SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult NOF SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult NOF TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult NOF MUC5AC", ".pdf")), w=11, h=8.5)
FeaturePlot(adultnof, features = "MUC5AC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

#Figure S5B
pdf(file.path("./", paste0("AT2 Marker Marker Dotplot", ".pdf")), w=11, h=8.5)
Dotplot_Zhiwei_Version(adultnof, c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "PGC", "NAPSA", "SFTA2", "SFTA3", "MUC5AC"))
dev.off()

##Figure S5C
fadmanof.data  <- Read10X(data.dir = "./") #CKDCI induced organoids for 21 days, expanded in primary AT2 organoid media without FGF10 for 120 days
fadmanof <- CreateSeuratObject(counts = fadmanof.data, project = "fadmanof", min.cells = 3, min.features = 200)
fadmanof <- RenameCells(object = fadmanof, add.cell.id = "fadmanof")
fadmanof[["percent.mt"]] <- PercentageFeatureSet(fadmanof, pattern = "^MT-")

VlnPlot(fadmanof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fadmanof <- subset(fadmanof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)
setwd("~/scRNAseq/Output/2022-03-16 Frum 2022 Figure 5/Fadma Nof")


fadmanof <- SCTransform(fadmanof)
fadmanof <- RunPCA(fadmanof)
ElbowPlot(fadmanof, ndims = 50)
fadmanof <- RunUMAP(fadmanof, dims = 1:16, reduction.name = "umap")
fadmanof <- FindNeighbors(fadmanof, dims = 1:16)
fadmanof <- FindClusters(fadmanof, resolution = 0.5)
DimPlot(fadmanof, label = TRUE)
write.csv(table(Idents(fadmanof), fadmanof$orig.ident), "ClusterbyRun.csv") #number of cells in each cluster by run

pdf(file.path("./", paste0("Louvain", ".pdf")), w=11, h=8.5)
DimPlot(fadmanof, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(fadmanof) <- "RNA"
fadmanof <- NormalizeData(fadmanof)

FeaturePlot(fadmanof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))

pdf(file.path("./", paste0("CKDCI NOF SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmanof, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("CKDCI NOF SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmanof, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("CKDCI NOF SFTPB", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmanof, features =  "SFTPB" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("CKDCI NOF TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmanof, features = "TOP2A", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("CKDCI NOF MUC5AC", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmanof, features = "MUC5AC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

##Figure S5D
pdf(file.path("./", paste0("AT2 Marker Marker Dotplot", ".pdf")), w=11, h=8.5)
Dotplot_Zhiwei_Version(fadmanof, c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "PGC", "NAPSA", "SFTA2", "SFTA3", "MUC5AC"))
dev.off()

#Figure S5F
day140.merge <- merge(twoefabnof, c(fadmanof, adultnof))
DefaultAssay(day140.merge) <- "RNA"
day140.merge <- NormalizeData(day140.merge)

otherlineages.markers <- c("AGER", "TP63", "FOXJ1" , "SCGB1A1", "ASCL1", "CHGA", "SPDEF", "MUC5AC")

VlnPlot(day140.merge, otherlineages.markers, group.by = "orig.ident")
p <- VlnPlot(day140.merge, features = otherlineages.markers, group.by = "orig.ident",  ncol = 4, pt.size = 0.5)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[7]] <- p[[7]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[8]] <- p[[8]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())

pdf(file.path("./", paste0("VLN Bulk Day 140 Alternative CellTypes Pt", ".pdf")), w=5.7, h=3.5)
p
dev.off()

p <- VlnPlot(day140.merge, features = otherlineages.markers, group.by = "orig.ident",  ncol = 4, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[7]] <- p[[7]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[8]] <- p[[8]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN Bulk Day 140 Alternative CellTypes", ".pdf")), w=5.7, h=3.5)
p
dev.off()

##Data for Figure 5C (graph made in prism)
PrctCellExpringGene(day140.merge, "SFTPC", group.by = "orig.ident")
PrctCellExpringGene(day140.merge, "SFTPA1", group.by = "orig.ident")
PrctCellExpringGene(day140.merge, "SFTPB", group.by = "orig.ident")

##Figure 5D
day140.merge$orig.ident <- factor(day140.merge$orig.ident, levels = c("adultnof", "twoefabnof", "fadmanof"))

VlnPlot(day140.merge, c("SFTPC", "SFTPA1", "SFTPB"), pt.size = 0, group.by = "orig.ident")

p <- VlnPlot(day140.merge, features = c("SFTPC", "SFTPA1", "SFTPB"), group.by = "orig.ident",  ncol = 4, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())

pdf(file.path("./", paste0("VLN Bulk Compare SFTPC A1 B", ".pdf")), w=5.7, h=1.5)
p
dev.off()

##Figure SE
PrctCellExpringGene(day140.merge, "MUC5AC", group.by = "orig.ident")

##Figure 5F - H
##Reference Based Mapping of organoid samples to single cell sequencing of microdissected distal and proximal lung (Murthy et al., 2022)
#load Murthy et al., 2022
epi.sub_SAE10x3.integrated_v2 <- readRDS("~/filepath/epi.sub_SAE10x3.integrated_v2.RDS")
epi.sub_SAE10x3.integrated_v2 <- UpdateSeuratObject(epi.sub_SAE10x3.integrated_v2)
epi.sub_SAE10x3.integrated_v2 <- RunUMAP(object = epi.sub_SAE10x3.integrated_v2, dims = 1:20,   min.dist = 1, n.neighbors = 20, spread = 1.5, local.connectivity = 10, return.model = TRUE)
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA" #data already normalized

#Map primary AT2 organoids cultured in primary AT2 organoid media for 90 days and then switched to primary AT2 organoid with no FGF10 for 30 days before sequencing
adultnof.data  <- Read10X(data.dir = "./")
adultnof <- CreateSeuratObject(counts = adultnof.data, project = "adultnof", min.cells = 3, min.features = 200)
adultnof <- RenameCells(object = adultnof, add.cell.id = "adultnof")
adultnof[["percent.mt"]] <- PercentageFeatureSet(adultnof, pattern = "^MT-")
VlnPlot(adultnof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
adultnof <- subset(adultnof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)
adultnof <- NormalizeData(adultnof)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
adultnof <- CellCycleScoring(adultnof, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
adultnof <- FindVariableFeatures(adultnof, selection.method = "vst", nfeatures = 2000)
adultnof <- ScaleData(adultnof)
adultnof <- RunPCA(adultnof, features = VariableFeatures(adultnof))
ElbowPlot(adultnof, ndims = 50)
adultnof <- RunUMAP(adultnof, dims = 1:20, reduction.name = "umap", return.model = TRUE)
adultnof <- FindNeighbors(adultnof, dims = 1:20)
adultnof <- FindClusters(adultnof, resolution = 0.3)
DimPlot(adultnof, label = TRUE, pt.size = 1)

common.features <- intersect(rownames(epi.sub_SAE10x3.integrated_v2), rownames(adultnof))
length(x = common.features)
epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = adultnof,
  k.filter = NA,
  reference.reduction = "pca",
  features = common.features,
  dims = 1:20
)


adultnof.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = adultnof,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)
adultnof.label.at2 <- subset(adultnof.mapped.tata, predicted.id == "AT2")
adultnof.label.at2.cellids <- WhichCells(adultnof.label.at2)
adultnof.label.pat2 <- subset(adultnof.mapped.tata, predicted.id == "Proliferating AT2")
adultnof.label.pat2.cellids <- WhichCells(adultnof.label.pat2)
adultnof.label.ne <- subset(adultnof.mapped.tata, predicted.id == "Neuroendocrine")
adultnof.label.ne.cellids <- WhichCells(adultnof.label.ne)
adultnof.label.1a1neg <- subset(adultnof.mapped.tata, predicted.id == "SFTPB+ SCGB3A2+ SCGB1A1-")
adultnof.label.1a1neg.cellids <- WhichCells(adultnof.label.1a1neg)
adultnof.label.trb <- subset(adultnof.mapped.tata, predicted.id == "SFTPC+ SCGB3A2+")
adultnof.label.trb.cellids <- WhichCells(adultnof.label.trb)

adultnof.mapped.tata@reductions[["umap"]] <- adultnof.mapped.tata@reductions[["ref.umap"]]
adultnof.visual <- merge(epi.sub_SAE10x3.integrated_v2, adultnof.mapped.tata, merge.dr = "umap")

predictions.adultnof.tata <- table(adultnof.mapped.tata$predicted.id)

#Map TGFbi/BMPa induced organoids for 21 days, expanded in primary AT2 organoid media without FGF10 for 120 days
twoefabnof.data <- Read10X(data.dir = "./")
twoefabnof <- CreateSeuratObject(counts = twoefabnof.data, project = "twoefabnof", min.cells = 3, min.features = 200)
twoefabnof <- RenameCells(object = twoefabnof, add.cell.id = "twoefabnof")
twoefabnof[["percent.mt"]] <- PercentageFeatureSet(twoefabnof, pattern = "^MT-")
VlnPlot(twoefabnof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
twoefabnof <- subset(twoefabnof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)

setwd("~/scRNAseq/Output/2022-03-16 Frum 2022 Figure 5/Label Transfer/twoefabnof")

twoefabnof <- NormalizeData(twoefabnof)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
twoefabnof <- CellCycleScoring(twoefabnof, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
twoefabnof <- FindVariableFeatures(twoefabnof, selection.method = "vst", nfeatures = 2000)
twoefabnof <- ScaleData(twoefabnof)
twoefabnof <- RunPCA(twoefabnof, features = VariableFeatures(twoefabnof))
ElbowPlot(twoefabnof, ndims = 50)
twoefabnof <- RunUMAP(twoefabnof, dims = 1:20, reduction.name = "umap", return.model = TRUE)
twoefabnof <- FindNeighbors(twoefabnof, dims = 1:20)
twoefabnof <- FindClusters(twoefabnof, resolution = 0.3)
DimPlot(twoefabnof, label = TRUE, pt.size = 1)

common.features <- intersect(rownames(epi.sub_SAE10x3.integrated_v2), rownames(twoefabnof))
length(x = common.features)
epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = twoefabnof,
  k.filter = NA,
  reference.reduction = "pca",
  features = common.features,
  dims = 1:20
)


twoefabnof.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = twoefabnof,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)

twoefabnof.label.at2 <- subset(twoefabnof.mapped.tata, predicted.id == "AT2")
twoefabnof.label.at2.cellids <- WhichCells(twoefabnof.label.at2)
twoefabnof.label.at1 <- subset(twoefabnof.mapped.tata, predicted.id == "AT1")
twoefabnof.label.at1.cellids <- WhichCells(twoefabnof.label.at1)
twoefabnof.label.pat2 <- subset(twoefabnof.mapped.tata, predicted.id == "Proliferating AT2")
twoefabnof.label.pat2.cellids <- WhichCells(twoefabnof.label.pat2)
twoefabnof.label.diffbas <- subset(twoefabnof.mapped.tata, predicted.id == "Differentiating Basal")
twoefabnof.label.diffbas.cellids <- WhichCells(twoefabnof.label.diffbas)
twoefabnof.label.ne <- subset(twoefabnof.mapped.tata, predicted.id == "Neuroendocrine")
twoefabnof.label.ne.cellids <- WhichCells(twoefabnof.label.ne)
twoefabnof.label.1a1neg <- subset(twoefabnof.mapped.tata, predicted.id == "SFTPB+ SCGB3A2+ SCGB1A1-")
twoefabnof.label.1a1neg.cellids <- WhichCells(twoefabnof.label.1a1neg)
twoefabnof.label.1a1pos <- subset(twoefabnof.mapped.tata, predicted.id == "SFTPB+ SCGB3A2+ SCGB1A1+")
twoefabnof.label.1a1pos.cellids <- WhichCells(twoefabnof.label.1a1pos)
twoefabnof.label.goblet <- subset(twoefabnof.mapped.tata, predicted.id == "MUC5AC+ MUC5B+")
twoefabnof.label.goblet.cellids <- WhichCells(twoefabnof.label.goblet)
twoefabnof.mapped.tata@reductions[["umap"]] <- twoefabnof.mapped.tata@reductions[["ref.umap"]]
twoefabnof.visual <- merge(epi.sub_SAE10x3.integrated_v2, twoefabnof.mapped.tata, merge.dr = "umap")



predictions.twoefabnof.tata <- table(twoefabnof.mapped.tata$predicted.id)


#Map CKDCI induced organoids for 21 days, expanded in primary AT2 organoid media without FGF10 for 120 days
fadmanof.data  <- Read10X(data.dir = "./")
fadmanof <- CreateSeuratObject(counts = fadmanof.data, project = "fadmanof", min.cells = 3, min.features = 200)
fadmanof <- RenameCells(object = fadmanof, add.cell.id = "fadmanof")
fadmanof[["percent.mt"]] <- PercentageFeatureSet(fadmanof, pattern = "^MT-")
VlnPlot(fadmanof, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fadmanof <- subset(fadmanof, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 20)



fadmanof <- NormalizeData(fadmanof)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fadmanof <- CellCycleScoring(fadmanof, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
fadmanof <- FindVariableFeatures(fadmanof, selection.method = "vst", nfeatures = 2000)
fadmanof <- ScaleData(fadmanof)
fadmanof <- RunPCA(fadmanof, features = VariableFeatures(fadmanof))
ElbowPlot(fadmanof, ndims = 50)
fadmanof <- RunUMAP(fadmanof, dims = 1:20, reduction.name = "umap", return.model = TRUE)
fadmanof <- FindNeighbors(fadmanof, dims = 1:20)
fadmanof <- FindClusters(fadmanof, resolution = 0.3)
DimPlot(fadmanof, label = TRUE, pt.size = 1)

common.features <- intersect(rownames(epi.sub_SAE10x3.integrated_v2), rownames(fadmanof))
length(x = common.features)
epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = fadmanof,
  k.filter = NA,
  reference.reduction = "pca",
  features = common.features,
  dims = 1:20
)


fadmanof.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = fadmanof,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)




fadmanof.mapped.tata@reductions[["umap"]] <- fadmanof.mapped.tata@reductions[["ref.umap"]]
fadmanof.visual <- merge(epi.sub_SAE10x3.integrated_v2, fadmanof.mapped.tata, merge.dr = "umap")

predictions.fadmanof.tata <- table(fadmanof.mapped.tata$predicted.id)

##Figure 5F
cols <- c("AT1" = "#f5fc86", "AT2" = "#7CAE00", "Proliferating" = "#ED68ED", "Proliferating AT2" ="#FF61CC", "Differentiating Basal" = "#8494FF", "Neuroendocrine" = "#ab467a" , "SFTPB+ SCGB3A2+ SCGB1A1-" = "#ABA300", "SFTPB+ SCGB3A2+ SCGB1A1+" = "#E68613", "SFTPB- KRT5+ Basal" = "#00C19A",  "Proliferating" = "#ED68ED", "MUC5AC+ MUC5B+" = "#00BE67", "Immature AT1" ="#0CB702", "Ciliated" ="#F8766D", "FOXJ1+ Secretory" = "#CD9600", "MUC5B+" = "#8494FF", "SFTPB+ KRT5_low Basal" = "#0CB702", "SFTPB+ KRT5- Basal" = "#00A9FF", "Deuterosomal" = "#FF61CC", "SFTPC+ SCGB3A2+" = "#00B8E7") 
pdf(file.path("./", paste0("twoefabnof Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(twoefabnof.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()
pdf(file.path("./", paste0("fadmanof Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(fadmanof.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()
pdf(file.path("./", paste0("adultnof Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(adultnof.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()
pdf(file.path("./", paste0("Color Mapped Reference", ".pdf")), w=8, h=8)
DimPlot(epi.sub_SAE10x3.integrated_v2, group.by = "cell_type", cols = cols, pt.size = 1, order = FALSE, label = TRUE, repel = TRUE, label.size = 6) + NoAxes() + NoLegend()
dev.off()
##Data for Figure 5 G and H (graphs made in prism)
write.csv(predictions.fadmanof.tata, "fadmanof to tata Predictions.csv")
write.csv(predictions.twoefabnof.tata, "twoefabnof to tata Predictions.csv")
write.csv(predictions.adultnof.tata, "twoefabnof to tata Predictions.csv")





##Merge Cells from each organoid mapping to AT2 cluster for comparison
adult.label.at2 <- subset(adultnof.mapped.tata, predicted.id == "AT2")
twoefab.label.at2 <- subset(twoefabnof.mapped.tata, predicted.id == "AT2")
fadma.label.at2 <- subset(fadmanof.mapped.tata, predicted.id == "AT2")

label.at2.merge <- merge(adult.label.at2, c(twoefab.label.at2, fadma.label.at2))
label.at2.merge.fastmnn <- label.at2.merge
set.seed(888)
DefaultAssay(label.at2.merge.fastmnn) <- "RNA"
label.at2.merge.fastmnn <- NormalizeData(label.at2.merge.fastmnn,
                                         normalization.method = "LogNormalize", scale.factor =10000)
label.at2.merge.fastmnn <-  FindVariableFeatures(label.at2.merge.fastmnn)
label.at2.merge.fastmnn <- RunFastMNN(object.list = SplitObject(label.at2.merge.fastmnn, split.by = "orig.ident"))
ElbowPlot(label.at2.merge.fastmnn, ndims = 50)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
label.at2.merge.fastmnn <- CellCycleScoring(label.at2.merge.fastmnn, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
label.at2.merge.fastmnn <- RunUMAP(label.at2.merge.fastmnn, dims = 1:30,reduction = "mnn")
label.at2.merge.fastmnn <- FindNeighbors(label.at2.merge.fastmnn, dims = 1:30,reduction = "mnn")
label.at2.merge.fastmnn <- FindClusters(label.at2.merge.fastmnn, resolution = 0.2, algorithm = 1)
label.at2.merge.fastmnn$orig.ident <- factor(label.at2.merge.fastmnn$orig.ident, levels = c("adultnof", "twoefabnof", "fadmanof"))

DimPlot(label.at2.merge.fastmnn, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 2)

##Figure 5I
pdf(file.path("./", paste0("Integrated Labeled AT2s by Sample", ".pdf")), w=11, h=8.5)
DimPlot(label.at2.merge.fastmnn, reduction = "umap",  pt.size = 2, group.by = "orig.ident")
dev.off()
#Figure 5J
pdf(file.path("./", paste0("Integrated Labeled AT2s FastMNN", ".pdf")), w=11, h=8.5)
DimPlot(label.at2.merge.fastmnn, reduction = "umap",  pt.size = 2)
dev.off()

#Figure 5K
DefaultAssay(label.at2.merge.fastmnn) <- "RNA"
label.at2.merge.fastmnn <- NormalizeData(label.at2.merge.fastmnn)
kotton.aec2.differentiation.set <- c("SFTPC", "CLDN18", "LAMP3", "SFTPB", "SFTPD", "NAPSA", "SLC34A2", "CXCL8")
kotton.aec2.maturation.set <- c("SFTPA2", "SFTPA1", "PGC", "SLPI", "CXCL5")

p <- VlnPlot(label.at2.merge.fastmnn, features = kotton.aec2.differentiation, group.by = "orig.ident",  ncol = 4, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[6]] <- p[[6]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[7]] <- p[[7]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[8]] <- p[[8]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())

pdf(file.path("./", paste0("AT2 Mapped Differentiation Kotton set", ".pdf")), w=5.7, h=1.8)
p
dev.off()

p <- VlnPlot(label.at2.merge.fastmnn, features = kotton.aec2.maturation, group.by = "orig.ident",  ncol = 4, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())
p[[4]] <- p[[4]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[5]] <- p[[5]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())

pdf(file.path("./", paste0("AT2 Mapped Maturation Kotton set", ".pdf")), w=5.7, h=1.8)
p
dev.off()





label.at2.merge.fastmnn.markers <- FindAllMarkers(label.at2.merge.fastmnn, min.pct = 0.25, logfc.threshold = 0.25)
label.at2.merge.fastmnn <- ScaleData(label.at2.merge.fastmnn)
label.at2.merge.fastmnn.markers %>% group_by(cluster)  %>% top_n(n = 20, wt = avg_log2FC) -> top10
pdf(file.path("./", paste0("Top10 Feature Heatmap v2", ".pdf")), w=22, h=17)
DoHeatmap(label.at2.merge.fastmnn, features = top10$gene) + NoLegend()
dev.off()

#Supplemental Table
zero.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = 0, min.pct = 0.25)
one.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = 1, min.pct = 0.25)
two.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = 2, min.pct = 0.25)

write.csv(zero.markers, "Cluster0.csv")
write.csv(one.markers, "Cluster1.csv")
write.csv(two.markers, "Cluster2.csv")

Idents(label.at2.merge.fastmnn) <- "orig.ident"
label.at2.merge.fastmnn.markers <- FindAllMarkers(label.at2.merge.fastmnn, min.pct = 0.25, logfc.threshold = 0.25)

#Supplemental Table

write.csv(label.at2.merge.fastmnn.markers, "AT2 label transfer enrichment markers.csv")


adultnof.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "adultnof", min.pct = 0.25)
twoefabnof.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "twoefabnof", min.pct = 0.25)
fadmanof.markers <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "fadmanof", min.pct = 0.25)

write.csv(adultnof.markers, "AdultEnriched.csv")
write.csv(twoefabnof.markers, "TwoefabEnriched.csv")
write.csv(fadmanof.markers, "FADMAEnriched.csv")

#Figure 5L
krasnow.at2.markers.2.list <- list(c("SFTPC",
                                     "NPC2",
                                     "SFTPA1",
                                     "NAPSA",
                                     "SFTPA2",
                                     "CTSH",
                                     "PGC",
                                     "SFTPD",
                                     "LAMP3",
                                     "ABCA3",
                                     "CHI3L2",
                                     "CA2",
                                     "DBI",
                                     "SERPINA1",
                                     "WIF1",
                                     "LRRK2",
                                     "C11orf96",
                                     "NRGN",
                                     "SFTA2",
                                     "PLA2G1B",
                                     "HHIP",
                                     "PEBP4",
                                     "CPB2",
                                     "NNMT",
                                     "MFSD2A",
                                     "SFTPB",
                                     "HLA-DPB1",
                                     "CRTAC1",
                                     "C2",
                                     "MALL",
                                     "MID1IP1",
                                     "SLC22A31",
                                     "FTL",
                                     "P3H2",
                                     "HLA-DPA1",
                                     "AK1",
                                     "LHFPL3-AS2",
                                     "C4BPA",
                                     "C16orf89",
                                     "CACNA2D2",
                                     "HLA-DRA",
                                     "TMEM243",
                                     "DRAM1",
                                     "LPCAT1",
                                     "TMEM163",
                                     "CXCL2",
                                     "SFTA3",
                                     "HHIP-AS1",
                                     "CD74",
                                     "ETV5",
                                     "SLC34A2",
                                     "DMBT1",
                                     "FGGY",
                                     "HLA-DRB1",
                                     "RASGRF1",
                                     "NECAB1",
                                     "SELENBP1",
                                     "SDR16C5",
                                     "TFPI",
                                     "HOPX",
                                     "RND1",
                                     "CD36",
                                     "FABP5",
                                     "MUC1",
                                     "RGS16",
                                     "ALPL",
                                     "ALOX15B",
                                     "LRRC36",
                                     "KCNJ15",
                                     "CSF3R",
                                     "SCD",
                                     "LGALSL",
                                     "LANCL1-AS1",
                                     "PPP1R1B",
                                     "SLC46A2",
                                     "DCXR",
                                     "C3",
                                     "NFKBIA",
                                     "BMP2",
                                     "DUSP6",
                                     "SLC6A14",
                                     "HLA-DRB5",
                                     "AREG",
                                     "GKN2",
                                     "CAT",
                                     "EDNRB",
                                     "CEBPD",
                                     "KCNJ8",
                                     "MSMO1",
                                     "KIAA1324L",
                                     "SNX30",
                                     "ETV1",
                                     "PARM1",
                                     "ZNF385B",
                                     "FASN",
                                     "FBP1",
                                     "HMOX1",
                                     "CITED2",
                                     "PLD3",
                                     "PMM1",
                                     "CDC42EP1",
                                     "ODC1",
                                     "ORM1",
                                     "HLA-DMA",
                                     "SPRY4",
                                     "SMAGP",
                                     "ACADL",
                                     "B3GNT8",
                                     "AGPAT2",
                                     "ESAM",
                                     "ASRGL1",
                                     "EPHX1",
                                     "LPL",
                                     "QDPR",
                                     "CISH",
                                     "MTRR",
                                     "CHI3L1",
                                     "LGMN",
                                     "CD44",
                                     "HLA-DQB1",
                                     "S100A14",
                                     "MSN",
                                     "MLPH",
                                     "GADD45B",
                                     "MBIP",
                                     "SOCS2",
                                     "GSTA4",
                                     "EP300-AS1",
                                     "TTN",
                                     "ACSL4",
                                     "ZDHHC3",
                                     "HP",
                                     "PID1",
                                     "AQP1",
                                     "HSD17B4",
                                     "SNX25",
                                     "PEBP1",
                                     "FGG",
                                     "C1orf21",
                                     "SPTSSA",
                                     "IDI1",
                                     "CXCL17",
                                     "AKAP13",
                                     "BTG1",
                                     "FMO5",
                                     "FDPS",
                                     "RAB27A",
                                     "TMSB4X",
                                     "ASAH1",
                                     "BLVRB",
                                     "FLRT3",
                                     "SECISBP2L",
                                     "CDK2AP2",
                                     "RBPMS-AS1",
                                     "TSC22D1",
                                     "MRPL14",
                                     "CHCHD7",
                                     "STC1",
                                     "ADI1",
                                     "ATP6V0E1",
                                     "NTN4",
                                     "LDHA",
                                     "TIFA",
                                     "GEM",
                                     "SAT2",
                                     "STEAP4",
                                     "SLC25A5",
                                     "TMEM41A",
                                     "POLR2C",
                                     "IFITM2",
                                     "LTA4H",
                                     "TPD52L1",
                                     "ZFP36",
                                     "ENO1",
                                     "SCP2",
                                     "CKS2",
                                     "HLA-DMB",
                                     "CYP51A1",
                                     "CHP1",
                                     "CREB3L1",
                                     "GSPT1",
                                     "CD83",
                                     "BRI3",
                                     "HMGCS1",
                                     "PTP4A3",
                                     "SOCS3",
                                     "SERPINB1",
                                     "AZGP1",
                                     "PNRC1",
                                     "SOD2",
                                     "XBP1",
                                     "GADD45G",
                                     "CSF3",
                                     "NFKBIZ",
                                     "TXNIP",
                                     "CXCL3",
                                     "PLIN2",
                                     "SEC61G",
                                     "MED24"))


epi.sub_SAE10x3.integrated_v2 <- readRDS("~/scRNAseq/Data/Cell/TATA Lab TRB-SC/epi.sub_SAE10x3.integrated_v2.RDS")
epi.sub_SAE10x3.integrated_v2 <- UpdateSeuratObject(epi.sub_SAE10x3.integrated_v2)
epi.sub_SAE10x3.integrated_v2 <- RunUMAP(object = epi.sub_SAE10x3.integrated_v2, dims = 1:20,   min.dist = 1, n.neighbors = 20, spread = 1.5, local.connectivity = 10, return.model = TRUE)
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA"




tata.at2s <- subset(epi.sub_SAE10x3.integrated_v2, idents = "AT2")
tata.mcs <- subset(epi.sub_SAE10x3.integrated_v2, idents = "Ciliated")

tata.mcs <- AddMetaData(tata.mcs, "ciliated", col.name = "orig.ident")
label.at2.merge.tata.threef.cil <- merge(tata.at2s, c(adult.label.at2, twoefab.label.at2, fadma.label.at2, threef, tata.mcs))

label.at2.merge.tata.threef.cil$orig.ident <- factor(label.at2.merge.tata.threef.cil$orig.ident, levels = c("SeuratProject","adultnof", "twoefabnof", "fadmanof", "threef", "ciliated"))


DefaultAssay(label.at2.merge.tata.threef.cil) <- "RNA"
label.at2.merge.tata.threef.cil <- NormalizeData(label.at2.merge.tata.threef.cil)

y <- label.at2.merge.tata.threef.cil
y <- AddModuleScore(y, features = krasnow.at2.markers.2.list, name = "At2KrasnowMarkers", search = TRUE, assay = "RNA")
p <- VlnPlot(y, features = "At2KrasnowMarkers1", pt.size = 0, group.by = "orig.ident")
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())

pdf(file.path("./", paste0("Alveolar Differentation Scoring 199 Krasnow Published TATA Threef Cil Included")), w=8, h=2.5)
p
dev.off()   

##Data for 5M
adultnofvtwoefab.markers.kat2.features <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "adultnof", ident.2 = "twoefabnof", features = krasnow.at2.markers.2, min.pct = 0.25)
adultnofvfadma.markers.kat2.features <- FindMarkers(label.at2.merge.fastmnn, ident.1 = "adultnof", ident.2 = "fadmanof", features = krasnow.at2.markers.2, min.pct = 0.25)


##Figure 6
fadmaday1.data <- Read10X(data.dir = "./") #BTP organoids treated with CKDCI for 1 day
fadmaday6.data <- Read10X(data.dir = "./") #BTP organoids treated with CKDCI for 6 days
fadmaday21.data <- Read10X(data.dir = "./")  #BTP organoids treated with CKDCI for 21 days

fadmaday1 <- CreateSeuratObject(counts = fadmaday1.data, project = "fadmaday1", min.cells = 3, min.features = 200)
fadmaday1 <- AddMetaData(fadmaday1, "FADMA Day 1", col.name = "group")
fadmaday1[["percent.mt"]] <- PercentageFeatureSet(fadmaday1, pattern = "^MT-")
fadmaday6 <- CreateSeuratObject(counts = fadmaday6.data, project = "fadmaday6", min.cells = 3, min.features = 200)
fadmaday6 <- AddMetaData(fadmaday6, "FADMA Day 1", col.name = "group")
fadmaday6[["percent.mt"]] <- PercentageFeatureSet(fadmaday6, pattern = "^MT-")
fadmaday21 <- CreateSeuratObject(counts = fadmaday21.data, project = "fadmaday21", min.cells = 3, min.features = 200)
fadmaday21 <- AddMetaData(fadmaday21, "FADMA Day 1", col.name = "group")
fadmaday21[["percent.mt"]] <- PercentageFeatureSet(fadmaday21, pattern = "^MT-")
fadma <- merge(threef, c(fadmaday1, fadmaday6, fadmaday21))
fadma.fastmnn <- subset(fadma, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
set.seed(888)
DefaultAssay(fadma.fastmnn) <- "RNA"
fadma.fastmnn <- NormalizeData(fadma.fastmnn,
                               normalization.method = "LogNormalize", scale.factor =10000)
fadma.fastmnn <-  FindVariableFeatures(fadma.fastmnn)
fadma.fastmnn <- RunFastMNN(object.list = SplitObject(fadma.fastmnn, split.by = "orig.ident"))
ElbowPlot(fadma.fastmnn, ndims = 50)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fadma.fastmnn <- CellCycleScoring(fadma.fastmnn, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
fadma.fastmnn <- RunUMAP(fadma.fastmnn, dims = 1:30,reduction = "mnn")
fadma.fastmnn <- FindNeighbors(fadma.fastmnn, dims = 1:30,reduction = "mnn")
fadma.fastmnn <- FindClusters(fadma.fastmnn, resolution = 0.3, algorithm = 1)
DimPlot(fadma.fastmnn, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 2)



##Data for Figure 6A (graph made in prism)
PrctCellExpringGene(fadma.fastmnn, genes = c("SFTPC", "SFTPA1", "SFTPB"), group.by = "orig.ident")


##Figure 6B
pdf(file.path("./", paste0("fadma FastMNN by Sample", ".pdf")), w=11, h=8.5)
DimPlot(fadma.fastmnn, reduction = "umap",  pt.size = 2, group.by = "group")
dev.off()
##Figure 6C
pdf(file.path("./", paste0("fadma FastMNN", ".pdf")), w=11, h=8.5)
DimPlot(fadma.fastmnn, reduction = "umap",  pt.size = 2, label = FALSE)
dev.off()
##


DefaultAssay(fadma.fastmnn) <- "RNA"
fadma.fastmnn <- NormalizeData(fadma.fastmnn)

##Figure 6D
pdf(file.path("./", paste0("FADMA Lineage Markers from Fig 3J Dotplot", ".pdf")), w=12.5, h=5.75)
Dotplot_Zhiwei_Version(fadma.fastmnn, twoefab.fastmnn.markers)
dev.off()
##Figure S6A
write.csv(table(fadma.fastmnn$group, fadma.fastmnn$Phase), "RunbyPhase.csv")
 
##Figure S6B
twoefab.fastmnn.cluster.6 <- subset(twoefab.fastmnn, ident = 6)
fadma.fastmnn.cluster.5 <- subset(fadma.fastmnn, ident = 5)
diff.neuroendocrine.compare <- merge(twoefab.fastmnn.cluster.6, fadma.fastmnn.cluster.5)

DefaultAssay(diff.neuroendocrine.compare) <- "RNA"
diff.neuroendocrine.compare <- NormalizeData(diff.neuroendocrine.compare)

p <- VlnPlot(diff.neuroendocrine.compare, features = c("SFTPC", "ASCL1", "CHGA"),  ncol = 3, pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
p[[2]] <- p[[2]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
p[[3]] <- p[[3]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.y = element_blank(), axis.title.x = element_blank(),  axis.text.x = element_blank())


pdf(file.path("./", paste0("Neuroendocrine Comparison", ".pdf")), w=5.7, h=1.8)
p
dev.off()

##Figure S6C
pdf(file.path("./", paste0("FADMA FastMNN Dotplot Cluster 3 and 6", ".pdf")), w=12.5, h=5.75)
Dotplot_Zhiwei_Version(fadma.fastmnn, c("SFTPC", "SFTPA1", "SFTPB", "NDRG1", "SLC5A3", "SLC2A3", "SLC6A6", "CXCL8", "TF", "LCN2"))
dev.off()

pdf(file.path("./", paste0("fadma FastMNN All Highlight", ".pdf")), w=11, h=8.5)
DimPlot(fadma.fastmnn, cells.highlight = list(fadmaday21.bestat2.cellids, fadmaday6.bestat2.cellids, fadmaday1.bestat2.cellids, threef.cellids), cols.highlight = c("purple", "blue", "green", "yellow"), pt.size = 2, sizes.highlight = 2, order = TRUE)
dev.off()

##Figure 6E
fadma.fastmnn <- FindClusters(fadma.fastmnn, resolution = 1, algorithm = 2)
DimPlot(fadma.fastmnn, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 2)

DefaultAssay(fadma.fastmnn) <- "RNA"
fadma.fastmnn <- NormalizeData(fadma.fastmnn)

sce.fadma.nosct <- as.SingleCellExperiment(fadma.fastmnn)
sce.fadma.nosct <- slingshot(sce.fadma.nosct, clusterLabels = sce.fadma.nosct@colData@listData[["ident"]]
                             , reducedDim = "UMAP",
                             start.clus = "6", allow.breaks = FALSE, stretch = 0, omega = TRUE, extend = "n")
lnes.fadma <- getLineages(reducedDim(sce.fadma.nosct,"UMAP"), sce.fadma.nosct@colData@listData[["ident"]], start.clus = "6")
pt.fadma.nosct <- slingPseudotime(sce.fadma.nosct, na = TRUE)
pt.fadma.nosct[is.na(pt.fadma.nosct)] <- 0
colors <- pal[cut(pt.fadma.nosct[,1], breaks = 100)]

pdf(file.path("./", paste0("FADMA Psuedotime with Lineage", ".pdf")), w=11, h=8.5)
plot(reducedDim(sce.fadma.nosct, "UMAP"), col = colors, pch = 16, cex = 1)
plot(SlingshotDataSet(sce.fadma.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 2, add = TRUE)
plot(SlingshotDataSet(sce.fadma.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 3, add = TRUE)
plot(SlingshotDataSet(sce.fadma.nosct), lwd =6, type = 'lineage', col = c("gray"), linInd = 5, add = TRUE)
plot(SlingshotDataSet(sce.fadma.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 4, add = TRUE)
plot(SlingshotDataSet(sce.fadma.nosct), lwd=6, type = 'lineage', col = c("gray"), linInd = 6, add = TRUE)
plot(SlingshotDataSet(sce.fadma.nosct), lwd=8, type = 'lineage', col = c("red"), linInd = 1, add = TRUE)
dev.off()

##Figure 6F
pdf(file.path("./", paste0("FADMA FastMNN SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(fadma.fastmnn, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()


pdf(file.path("./", paste0("FADMA FastMNN SCGB3A2", ".pdf")), w=11, h=8.5)
FeaturePlot(fadma.fastmnn, features = "SCGB3A2", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("FADMA FastMNN RNASE1", ".pdf")), w=11, h=8.5)
FeaturePlot(fadma.fastmnn, features = "RNASE1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

##Figure 6G
alvdiff.fastmnn$group <- factor(alvdiff.fastmnn$group, c("2FAB Day 21", "2FAB Day 6", "2FAB Day 1", "3F", "FADMA Day 1", "FADMA Day 6", "FADMA Day 21"))
p <- VlnPlot(alvdiff.fastmnn, features = "SCGB3A2", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN SCGB3A2 set", ".pdf")), w=8, h=4)
p
dev.off()
p <- VlnPlot(alvdiff.fastmnn, features = "SFTPC", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN SFTPC set", ".pdf")), w=6, h=4)
p
dev.off()
p <- VlnPlot(alvdiff.fastmnn, features = "RNASE1", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN RNASE1 set", ".pdf")), w=6, h=4)
p
dev.off()

##Figure 6H
alvdiff.combined <- merge(threef, c(twoefabday1, twoefabday6, twoefabday21, fadmaday1, fadmaday6, fadmaday21))
alvdiff.combined <- subset(alvdiff.combined, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
DefaultAssay(alvdiff.combined) <- "RNA"
alvdiff.combined <- NormalizeData(alvdiff.combined)
p <- VlnPlot(alvdiff.combined, features = "LAMP3", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN LAMP3 set", ".pdf")), w=6, h=4)
p
dev.off()
p <- VlnPlot(alvdiff.combined, features = "SFTPB", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN LAMP3 set", ".pdf")), w=6, h=4)
p
dev.off()
p <- VlnPlot(alvdiff.combined, features = "SFTPA1", group.by = "group",  pt.size = 0)
p[[1]] <- p[[1]] + theme(text = element_text(size=8), axis.text=element_text(size = 8),axis.title.x = element_blank(), axis.text.x = element_blank())
pdf(file.path("./", paste0("VLN SFTPA1 set", ".pdf")), w=6, h=4)
p
dev.off()


##Figure 6K - L
###Label transfer cells from all timepoints of TGFbi/BMPa and CK+DCI differentiations to same reference dataset as used in figure 5 (Murthy et al., 2022)
epi.sub_SAE10x3.integrated_v2 <- readRDS("~/scRNAseq/Data/Cell/TATA Lab TRB-SC/epi.sub_SAE10x3.integrated_v2.RDS")
epi.sub_SAE10x3.integrated_v2 <- UpdateSeuratObject(epi.sub_SAE10x3.integrated_v2)
epi.sub_SAE10x3.integrated_v2 <- RunUMAP(object = epi.sub_SAE10x3.integrated_v2, dims = 1:20,   min.dist = 1, n.neighbors = 20, spread = 1.5, local.connectivity = 10, return.model = TRUE)
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA"


fadma <- merge(fadmaday1, c(fadmaday6, fadmaday21))

fadma <- subset(fadma, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
fadma <- UpdateSeuratObject(fadma)

fadma <- NormalizeData(fadma)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fadma <- CellCycleScoring(fadma, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
fadma <- FindVariableFeatures(fadma, selection.method = "vst", nfeatures = 2000)
fadma <- ScaleData(fadma)
fadma <- RunPCA(fadma, features = VariableFeatures(fadma))
ElbowPlot(fadma, ndims = 50)
fadma <- RunUMAP(fadma, dims = 1:20, reduction.name = "umap", return.model = TRUE)
fadma <- FindNeighbors(fadma, dims = 1:20)
fadma <- FindClusters(fadma, resolution = 0.3)
DimPlot(fadma, label = TRUE, pt.size = 1)

common.features <- intersect(rownames(epi.sub_SAE10x3.integrated_v2), rownames(fadma))
common.features.2 <- intersect(rownames(common.features, rownames(twoefab)))
length(x = common.features.2)
epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = fadma,
  reference.reduction = "pca",
  features = common.features.2,
  dims = 1:20
)


fadma.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = fadma,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)



fadma.mapped.tata@reductions[["umap"]] <- fadma.mapped.tata@reductions[["ref.umap"]]
fadma.visual <- merge(epi.sub_SAE10x3.integrated_v2, fadma.mapped.tata, merge.dr = "umap")




epi.sub_SAE10x3.integrated_v2 <- readRDS("~/scRNAseq/Data/Cell/TATA Lab TRB-SC/epi.sub_SAE10x3.integrated_v2.RDS")
epi.sub_SAE10x3.integrated_v2 <- UpdateSeuratObject(epi.sub_SAE10x3.integrated_v2)
epi.sub_SAE10x3.integrated_v2 <- RunUMAP(object = epi.sub_SAE10x3.integrated_v2, dims = 1:20,   min.dist = 1, n.neighbors = 20, spread = 1.5, local.connectivity = 10, return.model = TRUE)
DefaultAssay(epi.sub_SAE10x3.integrated_v2) <- "RNA"


twoefab <- merge(twoefabday1, c(twoefabday6, twoefabday21))

twoefab <- subset(twoefab, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
twoefab <- UpdateSeuratObject(twoefab)

twoefab <- NormalizeData(twoefab)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
twoefab <- CellCycleScoring(twoefab, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
twoefab <- FindVariableFeatures(twoefab, selection.method = "vst", nfeatures = 2000)
twoefab <- ScaleData(twoefab)
twoefab <- RunPCA(twoefab, features = VariableFeatures(twoefab))
ElbowPlot(twoefab, ndims = 50)
twoefab <- RunUMAP(twoefab, dims = 1:20, reduction.name = "umap", return.model = TRUE)
twoefab <- FindNeighbors(twoefab, dims = 1:20)
twoefab <- FindClusters(twoefab, resolution = 0.3)
DimPlot(twoefab, label = TRUE, pt.size = 1)


epi.sub_SAE10x3.integrated_v2.anchors <- FindTransferAnchors(
  reference = epi.sub_SAE10x3.integrated_v2,
  normalization.method = "LogNormalize",
  query = twoefab,
  reference.reduction = "pca",
  features = common.features.2,
  dims = 1:20
)


twoefab.mapped.tata <- MapQuery(
  anchorset = epi.sub_SAE10x3.integrated_v2.anchors,
  query = twoefab,
  reference = epi.sub_SAE10x3.integrated_v2,
  refdata = "cell_type",
  reference.reduction = "pca", 
  transferdata.args = list(k.weight = 10),
  integrateembeddings.args = list(k.weight = 10),
  reduction.model = "umap"
)

twoefab.mapped.tata@reductions[["umap"]] <- twoefab.mapped.tata@reductions[["ref.umap"]]
twoefab.visual <- merge(epi.sub_SAE10x3.integrated_v2, twoefab.mapped.tata, merge.dr = "umap")

#Figure 6J
cols <- c("AT1" = "#f5fc86", "AT2" = "#7CAE00", "Proliferating" = "#ED68ED", "Proliferating AT2" ="#FF61CC", "Differentiating Basal" = "#8494FF", "Neuroendocrine" = "#ab467a" , "SFTPB+ SCGB3A2+ SCGB1A1-" = "#ABA300", "SFTPB+ SCGB3A2+ SCGB1A1+" = "#E68613", "SFTPB- KRT5+ Basal" = "#00C19A",  "Proliferating" = "#ED68ED", "MUC5AC+ MUC5B+" = "#00BE67", "Immature AT1" ="#0CB702", "Ciliated" ="#F8766D", "FOXJ1+ Secretory" = "#CD9600", "MUC5B+" = "#8494FF", "SFTPB+ KRT5_low Basal" = "#0CB702", "SFTPB+ KRT5- Basal" = "#00A9FF", "Deuterosomal" = "#FF61CC", "SFTPC+ SCGB3A2+" = "#00B8E7") 


pdf(file.path("./", paste0("fadma Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(fadma.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()
pdf(file.path("./", paste0("Color Mapped Reference", ".pdf")), w=8, h=8)
DimPlot(epi.sub_SAE10x3.integrated_v2, group.by = "cell_type", cols = cols, pt.size =1, order = FALSE, label = TRUE, label.size = 6, repel = TRUE) + NoAxes() + NoLegend()
dev.off()

pdf(file.path("./", paste0("twoefab Label Transfer to TATA colored", ".pdf")), w=8, h=8)
DimPlot(twoefab.visual, group.by = "predicted.id", pt.size = 1, order = FALSE, cols = cols, na.value = "grey68") + NoAxes() + NoLegend()
dev.off()


##Data for Figure 6J and K
predictions.twoefab.tata <- table(twoefab.mapped.tata$predicted.id)
write.csv(predictions.twoefab.tata, "Twoefab to tata Predictions.csv")
predictions.fadma.tata <- table(fadma.mapped.tata$predicted.id)
write.csv(predictions.fadma.tata, "Fadma to tata Predictions.csv")

##Figure 6L


twoefab.mapped.tata$group<- factor(twoefab.mapped.tata$group, levels = c("Day 1", "Day 6", "Day 21")) #for some reason using orig.ident would not work here.
fadma.mapped.tata$orig.ident<- factor(fadma.mapped.tata$orig.ident, levels = c("fadmaday1", "fadmaday6", "fadmaday21"))
pdf(file.path("./", paste0("FADMA LT Assembly", ".pdf")), w=8, h=8)
DimPlot(fadma.mapped.tata, reduction = "ref.umap", group.by = "orig.ident", label = FALSE, na.value = "grey68",
              label.size = 3, repel = TRUE, pt.size = 1) + NoLegend() + NoAxes()

dev.off()

pdf(file.path("./", paste0("TWOEFAB LT Assembly", ".pdf")), w=8, h=8)
DimPlot(twoefab.visual,  group.by = "group", cols = c("Day 1" = "#F8766D", "Day 6" = "#00BA38", "Day 21" = "#619CFF"), label = FALSE,
        na.value = "grey68", label.size = 3, repel = TRUE, pt.size = 1)  + NoLegend()+ NoAxes()
dev.off()
pdf(file.path("./", paste0("FADMA LT Assembly", ".pdf")), w=8, h=8)
DimPlot(fadma.visual,  group.by = "orig.ident", cols = c("fadmaday1" = "#F8766D", "fadmaday6" = "#00BA38", "fadmaday21" = "#619CFF"), label = FALSE,
na.value = "grey68", label.size = 3, repel = TRUE, pt.size = 1)  + NoLegend()+ NoAxes()
dev.off()








##Figure S6 D-F
##Analyze each timepoint of CKDCI differentiation seperately
###Day 1
fadmaday1 <- CreateSeuratObject(counts = fadmaday1.data, project = "fadmaday1", min.cells = 3, min.features = 200)
fadmaday1 <- AddMetaData(fadmaday1, "FADMA Day 1", col.name = "group")
fadmaday1[["percent.mt"]] <- PercentageFeatureSet(fadmaday1, pattern = "^MT-")
fadmaday1 <- RenameCells(object = fadmaday1, add.cell.id = "FADMADAY1")
VlnPlot(fadmaday1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fadmaday1 <- subset(fadmaday1, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fadmaday1 <- CellCycleScoring(fadmaday1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
fadmaday1 <- SCTransform(fadmaday1, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
fadmaday1 <- RunPCA(fadmaday1, features = VariableFeatures(object = fadmaday1))
ElbowPlot(fadmaday1, ndims = 50)
DefaultAssay(fadmaday1) <- "SCT"
fadmaday1 <- RunUMAP(fadmaday1, dims = 1:18)
fadmaday1 <- FindNeighbors(fadmaday1, dims = 1:18)
fadmaday1 <- FindClusters(fadmaday1, resolution = 0.5)
DimPlot(fadmaday1, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

##Figure S6D
pdf(file.path("./", paste0("Day 1 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(fadmaday1, pt.size = 2, label = FALSE)
dev.off()
DefaultAssay(fadmaday1) <- "RNA"
fadmaday1 <- NormalizeData(fadmaday1)

pdf(file.path("./", paste0("Day 1 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmaday1, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 1 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmaday1, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
fadmaday1.bestat2 <- subset(fadmaday1, idents = 5)

pdf(file.path("./", paste0("Day 1 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(fadmaday1, cells.highlight = fadmaday1.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()

###Day 6
fadmaday6 <- CreateSeuratObject(counts = fadmaday6.data, project = "fadmaday6", min.cells = 3, min.features = 200)
fadmaday6 <- AddMetaData(fadmaday6, "FADMA Day 6", col.name = "group")
fadmaday6[["percent.mt"]] <- PercentageFeatureSet(fadmaday6, pattern = "^MT-")
fadmaday6 <- RenameCells(object = fadmaday6, add.cell.id = "FADMADAY6")
VlnPlot(fadmaday6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fadmaday6 <- subset(fadmaday6, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fadmaday6 <- CellCycleScoring(fadmaday6, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
fadmaday6 <- SCTransform(fadmaday6, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
fadmaday6 <- RunPCA(fadmaday6, features = VariableFeatures(object = fadmaday6))
ElbowPlot(fadmaday6, ndims = 50)
DefaultAssay(fadmaday6) <- "SCT"
fadmaday6 <- RunUMAP(fadmaday6, dims = 1:18)
fadmaday6 <- FindNeighbors(fadmaday6, dims = 1:18)
fadmaday6 <- FindClusters(fadmaday6, resolution = 0.5)
DimPlot(fadmaday6, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

##Figure S6E
pdf(file.path("./", paste0("Day 6 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(fadmaday6, pt.size = 2, label = FALSE)
dev.off()
DefaultAssay(fadmaday6) <- "RNA"
fadmaday6 <- NormalizeData(fadmaday6)

pdf(file.path("./", paste0("Day 6 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmaday6, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 6 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmaday6, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
fadmaday6.bestat2 <- subset(fadmaday6, idents = 1)

pdf(file.path("./", paste0("Day 6 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(fadmaday6, cells.highlight = fadmaday6.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()

##Day 21
fadmaday21 <- CreateSeuratObject(counts = fadmaday21.data, project = "fadmaday21", min.cells = 3, min.features = 200)
fadmaday21 <- AddMetaData(fadmaday21, "FADMA Day 21", col.name = "group")
fadmaday21[["percent.mt"]] <- PercentageFeatureSet(fadmaday21, pattern = "^MT-")
fadmaday21 <- RenameCells(object = fadmaday21, add.cell.id = "FADMADAY21")
VlnPlot(fadmaday21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fadmaday21 <- subset(fadmaday21, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 20)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fadmaday21 <- CellCycleScoring(fadmaday21, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
fadmaday21 <- SCTransform(fadmaday21, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
fadmaday21 <- RunPCA(fadmaday21, features = VariableFeatures(object = fadmaday21))
ElbowPlot(fadmaday21, ndims = 50)
DefaultAssay(fadmaday21) <- "SCT"
fadmaday21 <- RunUMAP(fadmaday21, dims = 1:18)
fadmaday21 <- FindNeighbors(fadmaday21, dims = 1:18)
fadmaday21 <- FindClusters(fadmaday21, resolution = 0.5)
DimPlot(fadmaday21, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 2)

##Figure S6F

pdf(file.path("./", paste0("Day 21 Louvain", ".pdf")), w=11, h=8.5)
DimPlot(fadmaday21, pt.size = 2, label = FALSE)
dev.off()
DefaultAssay(fadmaday21) <- "RNA"
fadmaday21 <- NormalizeData(fadmaday21)

pdf(file.path("./", paste0("Day 21 SFTPC", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmaday21, features = "SFTPC", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Day 21 SFTPA1", ".pdf")), w=11, h=8.5)
FeaturePlot(fadmaday21, features = "SFTPA1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
fadmaday21.bestat2 <- subset(fadmaday21, idents = 0)

pdf(file.path("./", paste0("Day 21 Louvain Best AT2s", ".pdf")), w=11, h=8.5)
DimPlot(fadmaday21, cells.highlight = fadmaday21.bestat2.cellids, pt.size = 2, sizes.highlight = 2, cols.highlight = "lightblue")
dev.off()




##Figure S6G

fadmaday1.bestat2.cellids <- WhichCells(fadmaday1, idents = 5)
fadmaday6.bestat2.cellids <- WhichCells(fadmaday6, idents = 1)
fadmaday21.bestat2.cellids <- WhichCells(fadmaday21, idents = 0)


pdf(file.path("./", paste0("fadma FastMNN All Highlight", ".pdf")), w=11, h=8.5)
DimPlot(fadma.fastmnn, cells.highlight = list(fadmaday21.bestat2.cellids, fadmaday6.bestat2.cellids, fadmaday1.bestat2.cellids, threef.cellids), cols.highlight = c("purple", "blue", "green", "yellow"), pt.size = 2, sizes.highlight = 2, order = TRUE)
dev.off()

