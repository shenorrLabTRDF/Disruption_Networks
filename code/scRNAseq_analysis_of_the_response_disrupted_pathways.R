#---
#title: "scRNAseq analysis of the response disrupted pathways"
#author: "Shiran Vainberg"
#Proj : "A personalized network framework reveals predictive axis of anti-TNF response across diseases"
#---

# Loading packages and data

library(Seurat)
library(SingleR)
library(MTGOsc)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(caret)
library(glmnet)
library(ROCR)
library(pROC)
library(nsROC)
library(RColorBrewer)
library(huex10sttranscriptcluster.db)
library(GEOquery)
library(igraph)
library(psych)
library(Hmisc)
library(data.table)
library(CellMix)
library(dplyr)    
library(plyr)

options(stringsAsFactors = F)

code_dir = "./"
source(paste0(code_dir,'Disruption_Networks_functions.R'))

#Read data

#define the data directory path:
dataDir = '...'
setwd(dataDir)


V1.module.disrupted.nodes.wilcoxon <- readRDS(paste(dataDir, "/filtered.esetALL.rds", sep=""))

#######Analysis of 10x scRNAseq#######

pbmcR39.data <- Read10X(data.dir = paste(dataDir, "/singleCell_10X_IFX/Resp39_counts/filtered_feature_bc_matrix", sep=""))
pbmcNR32.data <- Read10X(data.dir = paste(dataDir, "/singleCell_10X_IFX/Resp39_counts/filtered_feature_bc_matrix", sep=""))

pbmcR39 <- CreateSeuratObject(counts = pbmcR39.data, project = "IFX.R", min.cells = 3)
pbmcNR32 <- CreateSeuratObject(counts = pbmcNR32.data, project = "IFX.NR", min.cells = 3)

pbmc.IFX <- merge(pbmcR39, y = pbmcNR32, add.cell.ids = c("R", "NR"), project = "IFX.R.NR")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc.IFX@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc.IFX@assays$RNA@data[mito.genes, ])/Matrix::colSums(pbmc.IFX@assays$RNA@data)
pbmc.IFX <- AddMetaData(object = pbmc.IFX, metadata = percent.mito, col.name = "percent.mito")
pbmc.IFX <- subset(pbmc.IFX, subset = nFeature_RNA > 200 & percent.mito < 0.25)
pbmc.IFX <- NormalizeData(object = pbmc.IFX)
pbmc.IFX <- FindVariableFeatures(object = pbmc.IFX)
pbmc.IFX <- ScaleData(object = pbmc.IFX, vars.to.regress = c("nUMI", "percent.mito"))
pbmc.IFX <- RunPCA(object = pbmc.IFX)


#Determine statistically significant principal components
pbmc.IFX <- JackStraw(object = pbmc.IFX, num.replicate = 100, dims = 50, reduction = "pca", verbose = TRUE)
pbmc.IFX<- ScoreJackStraw(object = pbmc.IFX, reduction = "pca", dims = 1:40,  do.plot = FALSE)
JackStrawPlot(pbmc.IFX, dims = 1:40, xmax=1, ymax = 1)
ElbowPlot(object = pbmc.IFX, ndims=50)

pbmc.IFX <- FindClusters(object = pbmc.IFX)
pbmc.IFX <- RunTSNE(object = pbmc.IFX, dims.use = 1:30)
DimPlot(object = pbmc.IFX, reduction = "tsne")

#SingleR cell annotation

hpca.se <- HumanPrimaryCellAtlasData()
ImmCells <- DatabaseImmuneCellExpressionData()
MonacoImm <- MonacoImmuneData()	
BlueprintEncodeData <- BlueprintEncodeData()

singleR.pbmc <- SingleR(test = as.SingleCellExperiment(pbmc.IFX), ref = ImmCells, labels = ImmCells$label.main, 
                        method = "single",
                        clusters = NULL, genes = "de", quantile = 0.8, fine.tune = TRUE,
                        tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                        assay.type.test = "logcounts", assay.type.ref = "logcounts",
                        check.missing = TRUE)


singleR.pbmc.hpca <- SingleR(test = as.SingleCellExperiment(pbmc.IFX), ref = hpca.se, labels = hpca.se$label.fine, 
                             method = "single",
                             clusters = NULL, genes = "de", quantile = 0.8, fine.tune = TRUE,
                             tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                             assay.type.test = "logcounts", assay.type.ref = "logcounts",
                             check.missing = TRUE)

singleR.pbmc.MonacoImm <- SingleR(test = as.SingleCellExperiment(pbmc.IFX), ref = MonacoImm, labels = MonacoImm$label.fine, 
                                  method = "single",
                                  clusters = NULL, genes = "de", quantile = 0.8, fine.tune = TRUE,
                                  tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                                  assay.type.test = "logcounts", assay.type.ref = "logcounts",
                                  check.missing = TRUE)

singleR.pbmc.BlueprintEncodeData <- SingleR(test = as.SingleCellExperiment(pbmc.IFX), ref = BlueprintEncodeData, labels = BlueprintEncodeData$label.fine, 
                                            method = "single",
                                            clusters = NULL, genes = "de", quantile = 0.8, fine.tune = TRUE,
                                            tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                                            assay.type.test = "logcounts", assay.type.ref = "logcounts",
                                            check.missing = TRUE)

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
pbmc.IFX$pruned.clusters <- singleR.pbmc$pruned.labels
pbmc.IFX$clusters <- singleR.pbmc$labels

pbmc.IFX$pruned.clusters.hpca <- singleR.pbmc.hpca$pruned.labels
pbmc.IFX$clusters.hpca <- singleR.pbmc.hpca$labels

pbmc.IFX$pruned.clusters.MonacoImm <- singleR.pbmc.MonacoImm$pruned.labels
pbmc.IFX$clusters.MonacoImm <- singleR.pbmc.MonacoImm$labels

pbmc.IFX$pruned.clusters.Blueprint <- singleR.pbmc.BlueprintEncodeData$pruned.labels
pbmc.IFX$clusters.Blueprint <- singleR.pbmc.BlueprintEncodeData$labels

pbmc.IFX$clusters.Blueprint.main <- ifelse(grepl("CD4+|Tregs", pbmc.IFX$clusters.Blueprint), "CD4+", 
                                           ifelse(grepl("CD8+", pbmc.IFX$clusters.Blueprint), "CD8+",
                                                  ifelse(grepl("B-cells|Plasma cells", pbmc.IFX$clusters.Blueprint), "B-cells", 
                                                         pbmc.IFX$clusters.Blueprint)))

Idents(pbmc.IFX) = "clusters.Blueprint.main" 
pbmc.IFX.sub <- subset(pbmc.IFX, idents = c("Monocytes", "CD8+", "CD4+", "B-cells", "NK cells"))


#Plot Fig. 4c, left
#tSNE plot representing cell types identities

DimPlot(pbmc.IFX, group.by="clusters.MonacoImm", reduction="tsne")
DimPlot(pbmc.IFX, group.by="clusters", reduction="tsne")
DimPlot(pbmc.IFX, group.by="clusters.Blueprint", reduction="tsne", pt.size=0.03, label = TRUE)

DimPlot(pbmc.IFX, group.by="clusters.MonacoImm", reduction="tsne", size=0.03,
        cells.highlight=findCells(pbmc.IFX, 'clusters.MonacoImm', 'Classical monocytes'))

DimPlot(pbmc.IFX, group.by="clusters.MonacoImm", reduction="tsne",size=0.03,
        cells.highlight=findCells(pbmc.IFX, 'clusters.MonacoImm', 'Non classical monocytes'))

DimPlot(pbmc.IFX, group.by="clusters.MonacoImm", reduction="tsne",size=0.03,
        cells.highlight=findCells(pbmc.IFX, 'clusters.MonacoImm', 'Intermediate monocytes'))


#Plot Fig. 4c, right
#tSNE plot colored by the expended fiber organization scaled expression

Baseline.diff.genes <- V1.module.disrupted.nodes.wilcoxon[V1.module.disrupted.nodes.wilcoxon$FDR<0.2,'SYMBOL']
fiber.organization.nodes <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% "GO_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION",'SYMBOL']
Baseline.diff.fiber.organization.nodes <- intersect(Baseline.diff.genes, fiber.organization.nodes)


fiber.organization.signature.score.res <- do.call('rbind', lapply(unique(Idents(pbmc.IFX)), function(cur.clust) {
  print(paste("choosing clust", cur.clust, sep="."))
  pbmc.IFX.sub.clust <- subset(pbmc.IFX, idents = cur.clust)
  features.in.cur.path <-Baseline.diff.fiber.organization.nodes 
  pbmc.IFX.sub.clust.path <- subset(pbmc.IFX.sub.clust, features=features.in.cur.path)
  print("start calculating scaled score")
  submat <- pbmc.IFX.sub.clust.path@assays$RNA@scale.data
  submat.scaled <- data.frame(t(submat), stringsAsFactors = F)
  submat.scaled$module.score <-apply(submat.scaled, 1, mean)
  submat.scaled$group <- pbmc.IFX.sub.clust.path@meta.data[match(row.names(submat.scaled), row.names(pbmc.IFX.sub.clust.path@meta.data)), 'orig.ident']
  submat.scaled$clust <- cur.clust
  submat.scaled$module.score <- as.numeric(as.character(submat.scaled$module.score))
  return(data.frame(Ident=row.names(pbmc.IFX.sub.clust.path@meta.data), clust=submat.scaled$clust, Module.score=submat.scaled$module.score, group=submat.scaled$group, stringsAsFactors = F))
}))

pbmc.IFX@meta.data$Ident.inf <- row.names(pbmc.IFX@meta.data)

pbmc.IFX@meta.data$fiber.V1.diff.score <-fiber.organization.signature.score.res[match(row.names(pbmc.IFX@meta.data), fiber.organization.signature.score.res$Ident), 'Module.score']

for.tsne.plot <- data.frame(pbmc.IFX@reductions$tsne@cell.embeddings, stringsAsFactors = F)
for.tsne.plot$fiber.V1.diff.score <-fiber.organization.signature.score.res[match(row.names(for.tsne.plot), fiber.organization.signature.score.res$Ident), 'Module.score']


ggplot(for.tsne.plot,
       aes(x = tSNE_1, y = tSNE_2,  color=fiber.V1.diff.score)) +
  geom_point(size=0.1)+
  scale_colour_gradientn(colours = c("#2B83BA", "#ABDDA4", "#FDAE61", "#D7191C"),
                         breaks = c(0,0.5,1,1.5), labels = c(0,0.5,1,1.5), limits=c(-0.6, 2))+
  theme_bw()+theme(strip.text = element_text(size=8), panel.grid.major = element_blank(),
                   axis.text.x=element_text(size=8,hjust=1),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))


saveRDS(pbmc.IFX, paste(dataDir, "/pbmc.IFX.rds", sep=""))

#Plot Supp. Fig 6
#Relative fiber-organization pathwaט score in each cell subset

fiber.organization.signature.score.res.sub <- fiber.organization.signature.score.res[fiber.organization.signature.score.res$clust %in% c("Monocytes", "CD8+", "CD4+", "B-cells", "NK cells"),]
fiber.organization.signature.score.res.sub$group <- factor(fiber.organization.signature.score.res.sub$group, levels=c("IFX.R", "IFX.NR"))
fiber.organization.signature.score.res.sub$clust <- factor(fiber.organization.signature.score.res.sub$clust, levels=c("Monocytes", "CD8+", "CD4+", "B-cells", "NK cells"))
ggplot(fiber.organization.signature.score.res.sub, aes(x = clust, y = Module.score, fill=group, color="black"))+
  geom_boxplot(outlier.color = NA)+
  scale_color_manual(values=c("black"="black"))+
  scale_fill_manual(values=c("IFX.R"="cyan3", "IFX.NR"="brown2"))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("Relative pathway score")+
  stat_compare_means(method = "wilcox",size=3, group="clust", comparisons=list(c(1,2), c(1,3), c(1,4), c(1,5)), label='p.format')+
  stat_compare_means(data=fiber.organization.signature.score.res[fiber.organization.signature.score.res$clust %in% "Monocytes",], method = "wilcox",size=3, group="group", label='p.format', label.y=3.5)


#Expand fiber organization differential signal in monocytes 
#using stringDB (knowledge based) and scRNA-seq (data-driven) nets

fiber.organization.core.genes <- V1.module.disrupted.nodes.wilcoxon[V1.module.disrupted.nodes.wilcoxon$FDR<0.1,'gene']

#generate scRNAseq networks within monocytes and in each monocyte subset seperately
geneSets.with.Go <- do.call('append', list(readList("h.all.v7.0.symbols.gmt"), readList("c2.cp.v7.0.symbols.gmt")))
geneSets.with.Go <- do.call('append', list(geneSets.with.Go, readList("c5.bp.v7.0.symbols.gmt")))

geneSets.with.Go <- lapply(geneSets.with.Go, function(l) data.frame(l, stringsAsFactors = F))

pathways.gene.df <-data.frame(rbindlist(geneSets.with.Go, use.names=T, idcol="pathway"), stringsAsFactors = F)
colnames(pathways.gene.df) <- c("pathway", "gene")


fiber.enrich.func <- function(cur.clust ,cur.cond, clust.ref, percentile, selected.genes) {
  Idents(pbmc.IFX) <- 'orig.ident'
  pbmc.R <- subset(pbmc.IFX, idents = cur.cond)
  Idents(pbmc.R) <- clust.ref
  cells.selected = Idents(pbmc.R) == cur.clust
  print(paste("selected.cells:", cur.clust, sep=" "))
  my.pbmc.R = pbmc.R
  my.pbmc.R = my.pbmc.R[, cells.selected]
  
  #Define MTGOsc folder to save temporary files and final results.
  root = paste(dataDir, "/singleCell_10X_IFX/MTGOres", sep="")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  
  #building a genes-pathways dictionary
  #use KEGG, Reactom, and Hallmark
  dict = write.dictionary(genes=pathways.genes.df$genes,
                          terms = pathways.genes.df$pathway, outfolder = root, overwrite=T)
  
  coexp = write.coexpressionMatrix(geneExpression = t(my.pbmc.R@assays$RNA@data[VariableFeatures(my.pbmc.R),]), fun= cor, method="spearman", outfolder = root, overwrite=T)

  edges = write.edges(coexpression = coexp, outfolder = root,
                      keep.weights = FALSE, fun = thinning_percentile, top_percentile = percentile, overwrite=T)
  
  #writing a parameter file, useful for MGTO
  write.paramFile(outfolder = root, overwrite=T)
  
  #Extract relevant sub-network:
  edges$gene1 <- as.character(edges$gene1)
  edges$gene2 <- as.character(edges$gene2)
  
  edges.sub <- edges[edges$gene1 %in% selected.genes| edges$gene2 %in% selected.genes,]
  coexp.sub <- coexp[coexp$gene1 %in% selected.genes| coexp$gene2 %in% selected.genes,]
  
  #the list of all genes involved in the cluster
  genes = unique(c(as.character(edges.sub$gene1), as.character(edges.sub$gene2)))
  
  #translating gene names to ENTREZID via org.Mm.eg.db database
  genes = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  #the actual enrichment
  enriched <-  enrichPathway(gene=genes$ENTREZID, pvalueCutoff=0.05, readable=T, organism = "human", pAdjustMethod = "BH")
  enriched@result <- enriched@result[enriched@result$p.adjust<0.05,]
  enriched@result <- enriched@result[enriched@result$Count>2,]
  
  g <- ReactomePA::dotplot(enriched, showCategory=20) + ggtitle(paste("Enrichment analysis", cur.cond, sep="."))
  if(is.null(enriched)) {
    return(list(net=edges.sub))} else {
      return(list(plot=g, data=enriched, net=edges.sub, coexp=coexp.sub))}
  return(list(net=edges.sub))}

res.R.NR.Monocytes <- fiber.enrich.func("Monocytes", c('IFX.R', 'IFX.NR'), "clusters.Blueprint.main", 0.05, fiber.organization.cor.genes)
res.R.NR.Intermediate.monocytes <- fiber.enrich.func("Intermediate monocytes", c('IFX.R', 'IFX.NR'), "clusters.MonacoImm", 0.05, fiber.organization.core.genes)
res.R.NR.Classical.monocytes <- fiber.enrich.func("Classical monocytes", c('IFX.R', 'IFX.NR'), "clusters.MonacoImm", 0.05, fiber.organization.core.genes)
res.R.NR.non.Classical.monocytes <- fiber.enrich.func("Non classical monocytes", c('IFX.R', 'IFX.NR'), "clusters.MonacoImm", 0.05, fiber.organization.core.genes)

res.R.NR.combined.monocytes.subsets <- unique(rbind(res.R.NR.Intermediate.monocytes$net[,colnames(res.R.NR.Intermediate.monocytes$net) %in% c("gene1", "gene2")], 
                                                    res.R.NR.Classical.monocytes$net[,colnames(res.R.NR.Classical.monocytes$net) %in% c("gene1", "gene2")],
                                                    res.R.NR.non.Classical.monocytes$net[,colnames(res.R.NR.non.Classical.monocytes$net) %in% c("gene1", "gene2")])) 


#Test knowledge based features based on stringdb

converter <- read.table(gzfile(paste(dataDir, "/human.name_2_string.tsv.gz", sep=""), open="r"), header=F)

converter$V2 <- as.character(converter$V2)
converter$V3 <- as.character(converter$V3)

string.edges <- read.table(gzfile(paste(dataDir, "/9606.protein.links.full.v11.0.txt.gz", sep=""), open="r"), header=T)

string.edges$protein1 <- as.character(string.edges$protein1)
string.edges$protein2 <- as.character(string.edges$protein2)

string.edges$SYMBOL1 <- converter[match(string.edges$protein1, converter$V3),'V2']
string.edges$SYMBOL2 <- converter[match(string.edges$protein2, converter$V3),'V2']

string.edges.sub <- string.edges[string.edges$SYMBOL1 %in% c(fiber.organization.core.genes, "TREM1", "CCR2", "CCL7")|string.edges$SYMBOL2 %in% c(fiber.organization.core.genes, "TREM1", "CCR2", "CCL7"),]
string.edges.sub <- string.edges.sub[string.edges.sub$combined_score>900,]

single.cell.net.monocytes <- res.R.NR.Monocytes$net

#put the original network genes in gene1
i <- as.character(single.cell.net.monocytes$gene2) %in% c(fiber.organization.core.genes, "TREM1", "CCR2", "CCL7")
tmp <- single.cell.net.monocytes$gene2[i]
single.cell.net.monocytes$gene2[i] <- single.cell.net.monocytes$gene1[i]
single.cell.net.monocytes$gene1[i] <- tmp

#use expanded list as the intersection between data driven and knowledge based networks
fiber.expansion.by.string.monocytes.genes <-  merge(single.cell.net.monocytes, string.edges.sub, by.x=c("gene1", "gene2"), by.y=c("SYMBOL1", "SYMBOL2"))
fiber.expansion.by.string.monocytes.genes <- unique(c(fiber.expansion.by.string.monocytes.genes$gene1, fiber.expansion.by.string.monocytes.genes$gene2))

saveRDS(fiber.expansion.by.string.monocytes.genes, paste(dataDir, "/fiber.expansion.by.string.monocytes.genes.rds", sep=""))

Idents(pbmc.IFX.monocytes) <- "orig.ident"

#Test differential expression between response groups for the expanded list
DE.summary.pbmc.IFX.monocytes.clusters <- do.call('rbind', lapply(unique(pbmc.IFX.monocytes@meta.data$clusters.MonacoImm), function(cur.clust) {
  Idents(pbmc.IFX.monocytes) <- "clusters.MonacoImm"
  subpbmc <- subset(pbmc.IFX.monocytes, idents = cur.clust)
  Idents(subpbmc) <- "orig.ident"
  DE.cur.clust <- FindMarkers(subpbmc, ident.1 =  "IFX.R" , ident.2 = "IFX.NR",  min.pct=0, logfc.threshold=0, test.use="wilcox")
  DE.cur.clust$Clust <- cur.clust
  DE.cur.clust$gene <- row.names(DE.cur.clust)
  return(DE.cur.clust)
}))

DE.summary.pbmc.IFX.monocytes.clusters.sub <- DE.summary.pbmc.IFX.monocytes.clusters[DE.summary.pbmc.IFX.monocytes.clusters$gene %in% fiber.expansion.by.string.monocytes.genes,]

DE.summary.pbmc.IFX.monocytes.clusters.sub <- do.call('rbind', lapply(unique(DE.summary.pbmc.IFX.monocytes.clusters.sub$Clust), function(cur.clust) {
  cur.clust.DE.summary.pbmc.IFX.monocytes.clusters.sub <- DE.summary.pbmc.IFX.monocytes.clusters.sub[DE.summary.pbmc.IFX.monocytes.clusters.sub$Clust %in% cur.clust,]
  cur.clust.DE.summary.pbmc.IFX.monocytes.clusters.sub$FDR <- p.adjust(cur.clust.DE.summary.pbmc.IFX.monocytes.clusters.sub$p_val, method="BH")
  return(cur.clust.DE.summary.pbmc.IFX.monocytes.clusters.sub)
}))

DE.summary.pbmc.IFX.monocytes.clusters.sub <- DE.summary.pbmc.IFX.monocytes.clusters.sub[DE.summary.pbmc.IFX.monocytes.clusters.sub$gene %in% unique(DE.summary.pbmc.IFX.monocytes.clusters.sub[DE.summary.pbmc.IFX.monocytes.clusters.sub$FDR<0.005,'gene']),]


#Plot Fig 4d, top
#Average expression of the selected expanded fiber-organization signature
#in monocyte-subsets

pbmc.IFX@meta.data$class <- paste(pbmc.IFX@meta.data$clusters.MonacoImm, pbmc.IFX@meta.data$orig.ident, sep=".")
Idents(pbmc.IFX) <- "class"

cluster.averages <- AverageExpression(pbmc.IFX, use.scale = T)

cluster.averages.df <- cluster.averages[["RNA"]]
cluster.averages.df.monocytes <- cluster.averages.df[,grepl("monocytes", colnames(cluster.averages.df))]
cluster.averages.df.monocytes.sub <- cluster.averages.df.monocytes[row.names(cluster.averages.df.monocytes) %in% fiber.expansion.by.string.monocytes.genes,]
cluster.averages.df.monocytes.sub <- data.frame(t(cluster.averages.df.monocytes.sub), stringsAsFactors = F)
cluster.averages.df.monocytes.sub$class <- row.names(cluster.averages.df.monocytes.sub)

cluster.averages.df.monocytes.sub.long <- reshape2::melt(cluster.averages.df.monocytes.sub, id.vars='class')
cluster.averages.df.monocytes.sub.long$subset <- sub(".IFX.R|.IFX.NR", "",cluster.averages.df.monocytes.sub.long$class )

cluster.averages.df.monocytes.sub.long$variable <- as.character(cluster.averages.df.monocytes.sub.long$variable)

cluster.averages.df.monocytes.sub.long.merged <- merge(cluster.averages.df.monocytes.sub.long, DE.summary.pbmc.IFX.monocytes.clusters.sub, by.x=c("variable", "subset"), by.y=c("gene", "Clust"), all.x=T)

cluster.averages.df.monocytes.sub.long.merged$mea.exp <- ifelse(cluster.averages.df.monocytes.sub.long.merged$FDR<0.05,1,0)*cluster.averages.df.monocytes.sub.long.merged$value

cluster.averages.df.monocytes.sub.long.merged.wide <- reshape2::dcast(class ~ variable, data= cluster.averages.df.monocytes.sub.long.merged, value.var='mea.exp')
row.names(cluster.averages.df.monocytes.sub.long.merged.wide) <- cluster.averages.df.monocytes.sub.long.merged.wide$class
cluster.averages.df.monocytes.sub.long.merged.wide$class <- NULL

anno.df = data.frame(clust = cluster.averages.df.monocytes.sub.long.merged.wide$subset <- sub(".IFX.R|.IFX.NR", "",row.names(cluster.averages.df.monocytes.sub.long.merged.wide) ),
                     stringsAsFactors = F)

anno.df$clust <- factor(anno.df$clust, levels=c("Intermediate monocytes", "Classical monocytes", "Non classical monocytes"))

ha1 <-   rowAnnotation(df =anno.df ,
                       col = list(clust =c("Classical monocytes" = "coral", "Intermediate monocytes" = "orange", "Non classical monocytes"="Gold")))

cluster.averages.df.monocytes.sub.long.merged.wide.final <- cluster.averages.df.monocytes.sub.long.merged.wide[,!colnames(cluster.averages.df.monocytes.sub.long.merged.wide) %in% "subset"]

i <- (colSums(cluster.averages.df.monocytes.sub.long.merged.wide.final, na.rm=T) != 0)
cluster.averages.df.monocytes.sub.long.merged.wide.final <- cluster.averages.df.monocytes.sub.long.merged.wide.final[,i]

ht <- Heatmap(cluster.averages.df.monocytes.sub.long.merged.wide.final, 
              col=colorRamp2(c( -1.2,-0.5, 0, 0.5,1.2), c("#084594", "#9ECAE1", "white", "#f0e3ab", "#b23030")),
              right_annotation = ha1,rect_gp = gpar(col = "black", lwd = 0.7), row_title = NULL,
              split=factor(anno.df$clust),cluster_row_slices = FALSE,
              name="scaled mean log2 TP10k+1", row_names_gp = gpar(fontsize = 10), 
              column_names_gp = gpar(fontsize = 10), show_row_names = TRUE, show_column_names = T,
              column_names_rot = 90, column_names_max_height = unit(10,"cm"))

draw(ht, padding = unit(c(2, 10, 2, 30), "mm"), heatmap_legend_side = "bottom") #bottom, left, top, right paddings


#Plot fig. 4d , bottom
#Module score for monocyte subsets

Idents(pbmc.IFX) <- "clusters.Blueprint"
pbmc.IFX.monocytes <- subset(pbmc.IFX, idents = c("Monocytes"))
Idents(pbmc.IFX.monocytes) <- "clusters.MonacoImm"
pbmc.IFX.monocytes <- subset(pbmc.IFX.monocytes, idents = c("Classical monocytes", "Intermediate monocytes", "Non classical monocytes"))

fiber.organization.signature.score.res.monocytes.subsets <-   do.call('rbind', lapply(unique(pbmc.IFX.monocytes@meta.data$clusters.MonacoImm), function(cur.clust) {
  print(paste("choosing clust", cur.clust, sep="."))
  pbmc.IFX.sub <- subset(pbmc.IFX.monocytes, idents = cur.clust)
  features.in.cur.path <- fiber.expansion.by.string.monocytes.genes
  pbmc.IFX.sub.path <- subset(pbmc.IFX.sub, features=features.in.cur.path)
  print("start calculating scaled score")
  submat <- pbmc.IFX.sub.path@assays$RNA@scale.data
  submat.scaled <- data.frame(t(submat), stringsAsFactors = F)
  submat.scaled$module.score <-apply(submat.scaled, 1, mean)
  submat.scaled$group <- pbmc.IFX.sub.path@meta.data[match(row.names(submat.scaled), row.names(pbmc.IFX.sub.path@meta.data)), 'orig.ident']
  submat.scaled$clust <- cur.clust
  return(data.frame(clust=submat.scaled$clust, Module.score=submat.scaled$module.score, group=submat.scaled$group, stringsAsFactors = F))
}))


fiber.organization.signature.score.res.monocytes.subsets$group <- factor(fiber.organization.signature.score.res.monocytes.subsets$group, levels=c("IFX.R", "IFX.NR"))

ggplot(fiber.organization.signature.score.res.monocytes.subsets, 
       aes(x = group, y = as.numeric(Module.score), fill=group, colour="black"))+
  geom_boxplot(outlier.colour = NA)+
  facet_wrap(~clust, scales = "free_y", ncol=7)+
  scale_color_manual(values=c("black"="black"))+
  scale_fill_manual(values=c("IFX.R"="cyan3", "IFX.NR"="brown2"))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("Relative pathway score")+
  stat_compare_means(method = "wilcox",size=3, comparisons=list(c(1,2)))


