library(Seurat)
library(patchwork)



#read data and create two seurat objects
a.data <- Read10X(data.dir = "D:/210626CKD_5vs3/c/")
b.data <- Read10X(data.dir = "D:/210619CKD-3vs1/c/")
a <- CreateSeuratObject(counts = a.data, project = "silke", min.cells = 5)
b <- CreateSeuratObject(counts = b.data, project = "silke", min.cells = 5)


#read data and create two further seurat objects
c.data <- Read10X(data.dir = "D:/210626CKD_1vs6/ckd/")
d.data <- Read10X(data.dir = "D:/210626CKD_1vs2/ckd/")
c <- CreateSeuratObject(counts = c.data, project = "silke", min.cells = 5)
d <- CreateSeuratObject(counts = d.data, project = "silke", min.cells = 5)


#read data and create two further seurat objects
e.data <- Read10X(data.dir = "D:/210626CKD_5vs3/ckd/")
f.data <- Read10X(data.dir = "D:/210626CKD_1vs4/ckd/")
e <- CreateSeuratObject(counts = e.data, project = "silke", min.cells = 5)
f <- CreateSeuratObject(counts = f.data, project = "silke", min.cells = 5)





#normalize, find variable features, subset, mito!!!!
a[["percent.mt"]] <- PercentageFeatureSet(a, pattern = "^mt-")
a <- subset(a, subset = nFeature_RNA > 600 & percent.mt < 5)
a <- NormalizeData(a, verbose = FALSE)
a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000)


b[["percent.mt"]] <- PercentageFeatureSet(b, pattern = "^mt-")
b <- subset(b, subset = nFeature_RNA > 600 & percent.mt < 5)
b <- NormalizeData(b, verbose = FALSE)
b <- FindVariableFeatures(b, selection.method = "vst", nfeatures = 2000)

c[["percent.mt"]] <- PercentageFeatureSet(c, pattern = "^mt-")
c <- subset(c, subset = nFeature_RNA > 600 & percent.mt < 5)
c <- NormalizeData(c, verbose = FALSE)
c <- FindVariableFeatures(c, selection.method = "vst", nfeatures = 2000)



d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^mt-")
d <- subset(d, subset = nFeature_RNA > 600 & percent.mt < 5)
d <- NormalizeData(d, verbose = FALSE)
d <- FindVariableFeatures(d, selection.method = "vst", nfeatures = 2000)

e[["percent.mt"]] <- PercentageFeatureSet(e, pattern = "^mt-")
e <- subset(e, subset = nFeature_RNA > 600 & percent.mt < 5)
e <- NormalizeData(e, verbose = FALSE)
e <- FindVariableFeatures(e, selection.method = "vst", nfeatures = 2000)


f[["percent.mt"]] <- PercentageFeatureSet(f, pattern = "^mt-")
f <- subset(f, subset = nFeature_RNA > 600 & percent.mt < 5)
f <- NormalizeData(f, verbose = FALSE)
f <- FindVariableFeatures(f, selection.method = "vst", nfeatures = 2000)





#create a merged object of three seurat objects (a,c,e)
ac.combined <- merge(a, y = c(b, c), add.cell.ids = c("A", "B", "C"), project = "3merge")
ac.combined

#create a merged object of three seurat objects (b,d,f)
ebdf.combined <- merge(d, y = c(e, f), add.cell.ids = c("D", "E", "F"), project = "3merge")
ebdf.combined



#this
ac.combined$ebdf.combined <- "ac.combined"
ebdf.combined$ebdf.combined <- "ebdf.combined"


#create a merged object of two seurat objects, this
ckd <- merge(ac.combined, y = ebdf.combined, add.cell.ids = c("ac.combined", "ebdf.combined"), project = "mergeall")
ckd

ckd[["percent.mt"]] <- PercentageFeatureSet(ckd, pattern = "^mt-")
ckd <- subset(ckd, subset = nFeature_RNA > 600 & percent.mt < 5)



#this
head(colnames(ckd))
table(ckd$orig.ident)

#then this
ckd.list <- SplitObject(ckd, split.by = "ebdf.combined")
ckd.list <- lapply(X = ckd.list, FUN = function(x) {x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})

#this
features <- SelectIntegrationFeatures(object.list = ckd.list)

END.anchors <- FindIntegrationAnchors(object.list = ckd.list, anchor.features = features)
END.combined <- IntegrateData(anchorset = END.anchors)
DefaultAssay(END.combined) <- "integrated"


#scale, PCA, clusters, UMAP of the merged object
END.combined <- ScaleData(END.combined, verbose = FALSE)
END.combined <- RunPCA(END.combined, npcs = 30, verbose = FALSE)

#determine number of PCAs
ElbowPlot(END.combined)
END.combined <- FindNeighbors(END.combined, reduction = "pca", dims = 1:12)
END.combined <- FindClusters(END.combined, resolution = 0.4)
END.combined <- RunUMAP(END.combined, reduction = "pca", dims = 1:12)
p1 <- DimPlot(END.combined, reduction = "umap", group.by = "ebdf.combined")
p2 <- DimPlot(END.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2









#find conserved markers of each cluster
DefaultAssay(END.combined) <- "RNA"

a.markers <- FindConservedMarkers(END.combined, ident.1 = 0, grouping.var = "ebdf.combined", verbose = FALSE)
head(a.markers)

b.markers <- FindConservedMarkers(END.combined, ident.1 = 1, grouping.var = "ebdf.combined", verbose = FALSE)
head(b.markers)

c.markers <- FindConservedMarkers(END.combined, ident.1 = 2, grouping.var = "ebdf.combined", verbose = FALSE)
head(c.markers)

d.markers <- FindConservedMarkers(END.combined, ident.1 = 3, grouping.var = "ebdf.combined", verbose = FALSE)
head(d.markers)

e.markers <- FindConservedMarkers(END.combined, ident.1 = 4, grouping.var = "ebdf.combined", verbose = FALSE)
head(e.markers)

f.markers <- FindConservedMarkers(END.combined, ident.1 = 5, grouping.var = "ebdf.combined", verbose = FALSE)
head(f.markers)


g.markers <- FindConservedMarkers(END.combined, ident.1 = 6, grouping.var = "ebdf.combined", verbose = FALSE)
head(g.markers)

h.markers <- FindConservedMarkers(END.combined, ident.1 = 7, grouping.var = "ebdf.combined", verbose = FALSE)
head(h.markers)

i.markers <- FindConservedMarkers(END.combined, ident.1 = 8, grouping.var = "ebdf.combined", verbose = FALSE)
head(i.markers)

j.markers <- FindConservedMarkers(END.combined, ident.1 = 9, grouping.var = "ebdf.combined", verbose = FALSE)
head(j.markers)

k.markers <- FindConservedMarkers(END.combined, ident.1 = 10, grouping.var = "ebdf.combined", verbose = FALSE)
head(k.markers)

l.markers <- FindConservedMarkers(END.combined, ident.1 = 11, grouping.var = "ebdf.combined", verbose = FALSE)
head(l.markers)

m.markers <- FindConservedMarkers(END.combined, ident.1 = 12, grouping.var = "ebdf.combined", verbose = FALSE)
head(m.markers)

n.markers <- FindConservedMarkers(END.combined, ident.1 = 13, grouping.var = "ebdf.combined", verbose = FALSE)
head(n.markers)

o.markers <- FindConservedMarkers(END.combined, ident.1 = 14, grouping.var = "ebdf.combined", verbose = FALSE)
head(o.markers)

p.markers <- FindConservedMarkers(END.combined, ident.1 = 15, grouping.var = "ebdf.combined", verbose = FALSE)
head(p.markers)

q.markers <- FindConservedMarkers(END.combined, ident.1 = 16, grouping.var = "ebdf.combined", verbose = FALSE)
head(q.markers)

r.markers <- FindConservedMarkers(END.combined, ident.1 = 17, grouping.var = "ebdf.combined", verbose = FALSE)
head(r.markers)




DimPlot(END.combined, reduction = "umap", split.by = "ebdf.combined")

p1 <- DimPlot(END.combined, reduction = "umap", group.by = "ebdf.combined")
p1
p2 <- DimPlot(END.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


#heatmap, did not work
DoHeatmap(END.combined, features = NULL, cells = NULL, group.by = "bdf.combined", group.bar = TRUE, group.colors = NULL, disp.min = -2.5,disp.max = NULL, slot = "scale.data", assay = NULL, label = TRUE, size = 5.5, hjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02, combine = TRUE)

#try...
KD.combined %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(KD.combined, features = top10$gene) + NoLegend()




#detect and identify clusters, celltypes....
VlnPlot(KD.combined, features = c("Tmem119"))
VlnPlot(KD.combined, features = c("Sox11"))
VlnPlot(KD.combined, features = c("Gja1"))
VlnPlot(KD.combined, features = c("Cldn11"))
VlnPlot(KD.combined, features = c("Tmeff2"))
VlnPlot(KD.combined, features = c("Kcnj8"))
VlnPlot(pro.combined, features = c("Syt1"))


VlnPlot(pro.combined, features = c("Kcnj13"))
VlnPlot(pro.combined, features = c("Ccdc153"))
VlnPlot(pro.combined, features = c("Acta2"))

VlnPlot(pro.combined, features = c("Slc47a1"))
VlnPlot(pro.combined, features = c("Clec2d"))
VlnPlot(pro.combined, features = c("Map2"))





#dorsal root ganglia

#ac
VlnPlot(END.combined, features = c("F3"))
VlnPlot(END.combined, features = c("Gfap"))
VlnPlot(END.combined, features = c("Aldh1l1"))
VlnPlot(END.combined, features = c("Aqp4"))
VlnPlot(END.combined, features = c("Atp1a2"))
VlnPlot(END.combined, features = c("Gja1"))
VlnPlot(END.combined, features = c("Slc1a2"))


VlnPlot(pro.combined, features = c("S100b"))


VlnPlot(pro.combined, features = c("Foxj1"))
VlnPlot(END.combined, features = c("Sox11"))




VlnPlot(END.combined, features = c("Sox11"))
VlnPlot(END.combined, features = c("Kcnj8"))


#ab hier; Neurons
VlnPlot(END.combined, features = c("Tubb3"))
VlnPlot(END.combined, features = c("Mef2c"))
VlnPlot(END.combined, features = c("Snap25"))
VlnPlot(END.combined, features = c("Syp"))
VlnPlot(END.combined, features = c("Rbfox3"))
VlnPlot(END.combined, features = c("Snhg11"))
VlnPlot(END.combined, features = c("Reln"))

VlnPlot(END.combined, features = c("Stmn2"))


#tcell like
VlnPlot(END.combined, features = c("Ms4a4b"))
VlnPlot(END.combined, features = c("Cd3d"))

#endoth cells
VlnPlot(END.combined, features = c("Apold1"))

#vlmc cell
VlnPlot(END.combined, features = c("Ogn"))
VlnPlot(END.combined, features = c("Lum"))
VlnPlot(END.combined, features = c("Dcn"))


#neuronal precursors
VlnPlot(END.combined, features = c("Gad1"))
VlnPlot(END.combined, features = c("Neurod2"))

#IL cells
VlnPlot(END.combined, features = c("Pomc"))
VlnPlot(END.combined, features = c("Pcsk2"))

#pituicyte
VlnPlot(END.combined, features = c("Rax"))
VlnPlot(END.combined, features = c("Adm"))


#epithelial  like cells
VlnPlot(END.combined, features = c("Krt18"))
VlnPlot(END.combined, features = c("Krt8"))
VlnPlot(END.combined, features = c("Clu"))

# mitotic cells
VlnPlot(END.combined, features = c("Ube2c"))
VlnPlot(END.combined, features = c("Top2a"))

VlnPlot(END.combined, features = c("Casp3"))
VlnPlot(END.combined, features = c("Il1b"))
VlnPlot(END.combined, features = c("Nlrp3"))
VlnPlot(END.combined, features = c("Casp8"))
VlnPlot(END.combined, features = c("Ctsc"))





#OC
VlnPlot(END.combined, features = c("Olig1"))
VlnPlot(END.combined, features = c("Plp1"))
VlnPlot(END.combined, features = c("Mbp"))
VlnPlot(END.combined, features = c("Mobp"))
VlnPlot(END.combined, features = c("Mog"))
VlnPlot(END.combined, features = c("Ugt8"))
VlnPlot(END.combined, features = c("Ermn"))


#OPCs
VlnPlot(END.combined, features = c("Pdgfra"))
VlnPlot(END.combined, features = c("C1ql1"))
VlnPlot(END.combined, features = c("Cspg4"))
VlnPlot(END.combined, features = c("Gpr17"))
VlnPlot(END.combined, features = c("Pdgfra"))

#vasc.cells
VlnPlot(END.combined, features = c("Rgs5"))
VlnPlot(END.combined, features = c("Flt1"))
VlnPlot(END.combined, features = c("Pecam1"))
VlnPlot(END.combined, features = c("Tek"))
VlnPlot(END.combined, features = c("Myl9"))
VlnPlot(END.combined, features = c("Pdgfrb"))


# epith..cells
VlnPlot(END.combined, features = c("Ttr"))

#emedp.cells
VlnPlot(END.combined, features = c("Tmem212"))
VlnPlot(END.combined, features = c("Kcnj13"))
VlnPlot(END.combined, features = c("Ccdc153"))

#macrophages/microglia
VlnPlot(END.combined, features = c("Tmem119"))
VlnPlot(END.combined, features = c("Ccl4"))
VlnPlot(END.combined, features = c("C1qa"), split.by = "ebdf.combined", split.plot = TRUE)
VlnPlot(END.combined, features = c("Ctss"))
VlnPlot(END.combined, features = c("Itgam"))
VlnPlot(END.combined, features = c("Ptprc"))

#pericytes
VlnPlot(END.combined, features = c("Kcnj8"))
VlnPlot(END.combined, features = c("Sox17"))

#Schwann
VlnPlot(END.combined, features = c("Mpz"))
VlnPlot(END.combined, features = c("Pmp22"))
VlnPlot(END.combined, features = c("Prx"))

#Menigeal cells
VlnPlot(END.combined, features = c("Dcn"))
VlnPlot(END.combined, features = c("Col3a1"))
VlnPlot(END.combined, features = c("Igf2"))



#does not work..
VlnPlot(m.markers, features = c("Casp1"), split.by = "ebdf.combined", split.plot = TRUE)



#rename clusters
END.combined <- RenameIdents(END.combined, `0` = "Neu1", `1` = "Neu2", `2` = "Neu3", `3` = "Neu4", `4` = "Neu5", `5` = "OC1", `6` = "AC1", `7` = "AC2", `8` = "OC2", `9` = "OPCs", `10` = "OC",  `11` = "VascCells/Peri", `12`= "MG", `13`= "Mtt", `14`= "x", `15`= "zz", `16`= "jj", `17`= "xi", `18` = "ee")

DimPlot(END.combined, label = TRUE)


#
p1 <- DimPlot(END.combined, reduction = "umap", group.by = "ebdf.combined")
p2 <- DimPlot(END.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


#determine cell number of each cluster, worked!!
n_cells <- FetchData(END.combined, vars = c("ident", "orig.ident")) %>%
dplyr::count(ident, orig.ident) %>%
tidyr::spread(ident, n)

n_cells <- FetchData(END.combined, vars = c("ident", "orig.ident", split.by = "ebdf.combined")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv2(n_cells, file="n-counts1.csv")


#did not work, i want to get average expression of each gene within each cluster..to create heatmaps...
avg.mNeu1$gene <- rownames(avg.mNeu1)
cells.0 <- subset(END.combined, idents = "a")
Idents(cells.0) <- ""
avg.cells.0 <- as.data.frame(log1p(AverageExpression(cells.0, verbose = FALSE)$RNA))
avg.cells.0$gene <- rownames(avg.cells.0)


#this worked
END.combined$celltype.ebdf.combined <- paste(Idents(END.combined), END.combined$ebdf.combined, sep = "_")
END.combined$celltype <- Idents(END.combined)
Idents(END.combined) <- "celltype.ebdf.combined"


View(END.combined)


#not yet tried, DotPlot
Idents(END.combined) <- factor(Idents(END.combined), levels = c("PrecNeurons1", "Neurons1", "PrecNeurons2", "Astrocytes", "Neurons2", "Oligodendrocytes1", "Astrocytes", "Ilcells", "Oligodendrocytes2", "OPCs", "VascCells/peri", "Microglia"))
markers.to.plot <- c("Tubb3", "Gad1", "Aldh1l1", "Olig1", "Plp1", "Pdgfra", "C1qa", "Ctss", "Dcn", "Prx")
DotPlot(END.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "bd.combined") +
  RotatedAxis()



#do this!!!
Neu1.response <- FindMarkers(END.combined, ident.1 = "Neu1_ebdf.combined", ident.2 = "Neu1_ac.combined", verbose = FALSE)
head(Neu1.response, n = 15)
rownames(Neu1.response[Neu1.response$avg_log2FC < -0.4 | Neu1.response$avg_log2FC > 0.4  ,]) # deutliche Expression


#perform this
write.csv2(Neu1.response, file="VARNeu1.csv")

Neu2.response <- FindMarkers(END.combined, ident.1 = "Neu2_ebdf.combined", ident.2 = "Neu2_ac.combined", verbose = FALSE)
head(Neu2.response, n = 15)
rownames(Neu2.response[Neu2.response$avg_log2FC < -0.4 | Neu2.response$avg_log2FC > 0.4  ,]) # deutliche Expression


#perform this
write.csv2(Neu2.response, file="VARNeu2.csv")


Neu3.response <- FindMarkers(END.combined, ident.1 = "Neu3_ebdf.combined", ident.2 = "Neu3_ac.combined", verbose = FALSE)
head(Neu3.response, n = 15)
rownames(Neu3.response[Neu3.response$avg_log2FC < -0.4 | Neu3.response$avg_log2FC > 0.4  ,]) # deutliche Expression


#perform this
write.csv2(Neu3.response, file="VARNeu3.csv")

Neu4.response <- FindMarkers(END.combined, ident.1 = "Neu4_ebdf.combined", ident.2 = "Neu4_ac.combined", verbose = FALSE)
head(Neu4.response, n = 15)
rownames(Neu4.response[Neu4.response$avg_log2FC < -0.4 | Neu4.response$avg_log2FC > 0.4  ,]) # deutliche Expression


#perform this
write.csv2(Neu4.response, file="VARNeu4.csv")



#MG
MG.response <- FindMarkers(END.combined, ident.1 = "MG_ebdf.combined", ident.2 = "MG_ac.combined", verbose = FALSE)
head(MG.response, n = 15)
rownames(MG.response[MG.response$avg_log2FC < -0.4 | MG.response$avg_log2FC > 0.4  ,]) # deutliche Expression

#then this, dafe it!!
write.csv2(MG.response, file="VAR.MG.csv")

Neu5.response <- FindMarkers(END.combined, ident.1 = "Neu5_ebdf.combined", ident.2 = "Neu5_ac.combined", verbose = FALSE)
head(Neu5.response, n = 15)
rownames(Neu5.response[Neu5.response$avg_log2FC < -0.4 | Neu5.response$avg_log2FC > 0.4  ,]) # deutliche Expression


#perform this
write.csv2(Neu5.response, file="VARNeu5.csv")



cluster.averages3 <- AverageExpression(object = END.combined)
#order the columns
cluster.averages3[["RNA"]] <- cluster.averages3[["RNA"]][,order(colnames(cluster.averages3[["RNA"]]))]
head(x = cluster.averages3[["RNA"]][, 1:5])  # all genes
write.table(cluster.averages3[["RNA"]], paste0("cluster_avg_expr.txt"))


#split by conditions
cluster.averages4 <- AverageExpression(object = END.combined, add.ident = "ebdf.combined")
#order the columns
cluster.averages4[["RNA"]] <- cluster.averages4[["RNA"]][,order(colnames(cluster.averages4[["RNA"]]))]
head(x = cluster.averages4[["RNA"]][, 1:5])  # all genes
write.table(cluster.averages4[["RNA"]], paste0("cluster_avg_expr_diet.txt"))





DefaultAssay(END.combined) <- "RNA"



#this, do it!
MGaverage <- AverageExpression(END.combined, ident.1 = 0, grouping.var = "ebdf.combined", verbose = FALSE)
head(MGaverage)

install.packages("writexl")
lapply(MGaverage, function(x) write.table( data.frame(x), 'testMGaverage1.csv'  , append= T, sep=',' ))

write.csv2(aver.expr., file="aver.expr.csv")
write.table(aver.expr., file="aver.expr.txt", dec =",", sep = "\t")


MGaverage <- AverageExpression(KD.combined, ident.1 = 0, grouping.var = "nx", verbose = FALSE)
head(MGaverage)

install.packages("writexl")
lapply(MGaverage, function(x) write.table( data.frame(x), 'testMGaverage.csv'  , append= T, sep=',' ))

write.csv2(aver.expr., file="aver.expr.csv")
write.table(aver.expr., file="aver.expr.txt", dec =",", sep = "\t")








GetAssayData(object = END.combined, slot = "counts")
END.combined
write.csv2(GetAssayData, file="aver.expr.csv")


Neu1 <- FetchData(object = END.combined, vars = 'Neu1')
head(x = Neu1)
head(x = FetchData(object = END.combined, vars = c('groups', 'ident')))




GENES<-rownames(END.combined)
CELLs<-colnames(END.combined)
names(END.combined)
dim(x = END.combined)
FindMyGene<-which(GENES=="myGene")
FindMyCell<-which(GENES=="myCell")

CountMyGeneMyCell<-GetAssayData(object = END.combined, slot = 'counts')[FindMyGene, FindMyCell]
write.csv2(CountMyGeneMyCell, file="ckd")

#Raw counts can be found in object@raw.data, with gene names as rownames and cell names as colnames
END.combined@raw.data



write.table(as.matrix(GetAssayData(object = END.combined, slot = "counts")), 
            'E:/ckd/counts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)




save.image("211222-4er.1,5,6vs2,3,4.RData")


