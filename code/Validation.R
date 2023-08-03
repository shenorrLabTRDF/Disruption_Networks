#---
#title: "Predictive signature and validation"
#author: "Shiran Vainberg"
#Proj : "A personalized network framework reveals predictive axis of anti-TNF response across diseases"
#---

# Loading packages and functions

library(GEOquery)
library(affycoretools)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(reshape2)
library(limma)
library(glmnet)
library(glmnetcr)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ComplexHeatmap)
library(dplyr)    
library(plyr)
library(shenR)
library(CellMix)
library(gtools)
library(rstatix)
library(scales)
library(VennDiagram)
library(grid)
library(CyTNF)
library(Hmisc)
library(data.table)
library(fgsea)
library(ggbeeswarm)
library(gage)
library(xbioc)
library(omics)
library(AnnotationDbi)
library(networkD3)
library(Biobase)
library(data.table)
library(igraph)
library(NetPathMiner)
library(psych)
library(foreach)
library(doParallel)
library(ggforce)
library(visNetwork)
library(rlist)
library(circlize)
library(gridExtra)
library(grid)
library(caret)
library(glmnet)
library(ROCR)
library(pROC)
library(nsROC)
library(Seurat)
library(MTGOsc)
library(ReactomePA)
library(clusterProfiler)
library(RColorBrewer)
library(huex10sttranscriptcluster.db)
library(dplyr)    
library(plyr)


options(stringsAsFactors = F)

code_dir = "./"
source(paste0(code_dir,'Disruption_Networks_functions.R'))

#define the data directory path:
dataDir = '...'
setwd(dataDir)

# Read data

filtered.esetALL <- readRDS(paste(dataDir, "/filtered.esetALL.rds", sep=""))
fiber.expansion.by.string.monocytes.genes <- readRDS(paste(dataDir, "/fiber.expansion.by.string.monocytes.genes.rds", sep=""))
pbmc.IFX <- readRDS(paste(dataDir, "/pbmc.IFX.rds", sep=""))

#Predictive signature identification
#perform variable selection for extended differential fiber-organization related genes

#lasso for variable selection

Baseline.DEGs <- fiber.expansion.by.string.monocytes.genes
esetMat <- filtered.esetALL[fData(filtered.esetALL)$SYMBOL %in% Baseline.DEGs & fData(filtered.esetALL)$dataType %in% "AGX", pData(filtered.esetALL)$Visit %in% "V1"]
esetMat.exprs <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
esetMat.exprs <- esetMat.exprs[!row.names(esetMat.exprs) %in% 'featureID',] 
esetMat.exprs$group <- pData(filtered.esetALL)[match(row.names(esetMat.exprs), pData(filtered.esetALL)$sampleID),'group']
esetMat.exprs$group <- ifelse( esetMat.exprs$group %in% "IM", "NR",  esetMat.exprs$group)
esetMat.exprs$sampleID <- row.names(esetMat.exprs)
colnames(esetMat.exprs)[!colnames(esetMat.exprs) %in% c("group", "sampleID")] <- as.character(fData(filtered.esetALL)[match(colnames(esetMat.exprs[,!colnames(esetMat.exprs) %in% c("group", "sampleID")]), fData(filtered.esetALL)$featureID),'SYMBOL'])

x <- as.matrix(esetMat.exprs[,!colnames(esetMat.exprs) %in% c("group", "sampleID")])
y <- factor(esetMat.exprs$group, levels=c("R", "NR"))


#Classifier using the Caret package with glmnet

my_trainControl <- trainControl(
  method='repeatedcv', 
  number=2, 
  repeats=100,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  verboseIter = FALSE,
  savePredictions = TRUE)

# Build the model

lambda = seq(0.0001, 1, length = 20)
alpha <- seq(0.5, 1, length = 6)

lasso <- train(
  x = x, y = y, method = "glmnet",
  trControl =my_trainControl,
  tuneGrid = expand.grid(alpha = alpha, lambda = lambda)
)

plot(lasso)


# Model coefficients
coeff_df <- coef(lasso$finalModel, lasso$bestTune$lambda)
coeff_df = data.frame(as.matrix(coeff_df), stringsAsFactors = F)
coeff_df$SYMBOL <- fData(filtered.esetALL)[match(row.names(coeff_df), fData(filtered.esetALL)$featureID),'displayName.final']

#####################################################
#calculate final classifier AUC based on primary cohort
#####################################################

coeff_df$gene <- row.names(coeff_df)
coeff_df <- coeff_df[!coeff_df$gene %in% "(Intercept)",]
Baseline.DEGs <- coeff_df[!coeff_df$X1==0,'gene']

esetMat <- filtered.esetALL[fData(filtered.esetALL)$SYMBOL %in% Baseline.DEGs & fData(filtered.esetALL)$dataType %in% "AGX", pData(filtered.esetALL)$Visit %in% "V1"]
V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
V1.module.disrupted.nodes.exp.mat$group <- pData(filtered.esetALL)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(filtered.esetALL)$sampleID),'group']
V1.module.disrupted.nodes.exp.mat$group <- ifelse(V1.module.disrupted.nodes.exp.mat$group %in% "IM", "NR", V1.module.disrupted.nodes.exp.mat$group)
V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))
V1.module.disrupted.nodes.exp.mat$sampleID <- row.names(V1.module.disrupted.nodes.exp.mat)

#Test performance on primary cohort using bootstraping
res.samp1000 <- sapply(1:1000, simplify=T, function(k) {
  samp <- sample(seq(nrow(V1.module.disrupted.nodes.exp.mat)), nrow(V1.module.disrupted.nodes.exp.mat), replace = T)
})

AUC.boot <- do.call('rbind', lapply(1:ncol(res.samp1000), function(cur.samp) {
  cur.V1.module.disrupted.nodes.exp.mat <- V1.module.disrupted.nodes.exp.mat[res.samp1000[,cur.samp],]
  if(ncol(cur.V1.module.disrupted.nodes.exp.mat)>3) {
    cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")] <- apply(cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group",  "sampleID")], 2, scale)
  } else {
    cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")] <- scale( cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")])
  }
  cur.V1.module.disrupted.nodes.exp.mat$scaled.module.score <- apply(cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")],1,mean)
  Gmean.difference.df <- data.frame(cur.V1.module.disrupted.nodes.exp.mat[,colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("scaled.module.score", "group", "sampleID")], stringsAsFactors = F)
  Gmean.difference.df$group <- as.character(pData(filtered.esetALL)[match(Gmean.difference.df$sampleID, as.character(pData(filtered.esetALL)$sampleID)),'group'])
  Gmean.difference.df$group <- ifelse(Gmean.difference.df$group %in% "IM", "NR", Gmean.difference.df$group)
  Gmean.difference.df$group <- factor(Gmean.difference.df$group, levels=c("R", "NR"))
  res <- try(roc(Gmean.difference.df$group, as.numeric(as.character(Gmean.difference.df$scaled.module.score)), quiet = TRUE)$auc)
  if(inherits(res, "try-error"))
  {
    print("error")
  } else {
    AUC <- roc(Gmean.difference.df$group, as.numeric(as.character(Gmean.difference.df$scaled.module.score)), quiet = TRUE)$auc
    return(AUC)
  }
}))

#test CI
AUC <- mean(AUC.boot)
CI2.5 <- quantile(AUC.boot, 0.025)
CI97.5 <- quantile(AUC.boot, 0.975)

#ROC of the relative predictive score
cur.V1.module.disrupted.nodes.exp.mat <- V1.module.disrupted.nodes.exp.mat
cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")] <- apply(cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group",  "sampleID")], 2, scale)
View(cur.V1.module.disrupted.nodes.exp.mat)
cur.V1.module.disrupted.nodes.exp.mat$scaled.module.score <- apply(cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")],1,mean)
Gmean.difference.df <- data.frame(cur.V1.module.disrupted.nodes.exp.mat[,colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("scaled.module.score", "group", "sampleID")], stringsAsFactors = F)
Gmean.difference.df$group <- as.character(pData(filtered.esetALL)[match(Gmean.difference.df$sampleID, as.character(pData(filtered.esetALL)$sampleID)),'group'])
Gmean.difference.df$group <- ifelse(Gmean.difference.df$group %in% "IM", "NR", Gmean.difference.df$group)
Gmean.difference.df$group <- factor(Gmean.difference.df$group, levels=c("R", "NR"))
View(Gmean.difference.df)
rocobj=plot.roc(Gmean.difference.df$group, Gmean.difference.df$scaled.module.score, main="ROC curve with 95% CIs",
                percent=TRUE, ci=TRUE, print.auc=TRUE, col="red")

roc(Gmean.difference.df$group, Gmean.difference.df$scaled.module.score)

## Calculate CI of sensitivity at select set of
## specificities and form a 'band' (might take a bit):
ciobj=ci.se(rocobj,specificities=seq(0, 100, 2))
plot(ciobj, type="shape", col="lightgrey") # blue band

# permute labels and calculate AUC for the optimal model
perm.label <- do.call('cbind', lapply(1:10000, function(i) {
  test <-  sample(seq(nrow(V1.module.disrupted.nodes.exp.mat)), nrow(V1.module.disrupted.nodes.exp.mat), replace = FALSE)
  return(test)
}))

perm.AUC <- do.call('rbind', lapply(1:ncol(perm.label), function(cur.perm) {
  V1.module.disrupted.nodes.exp.mat.scaled <- V1.module.disrupted.nodes.exp.mat
  V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")] <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group",  "sampleID")], 2, scale)
  V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")],1,mean)
  V1.module.disrupted.nodes.exp.mat.scaled$perm.group <- V1.module.disrupted.nodes.exp.mat.scaled[perm.label[,cur.perm],'group']
  AUC <-  roc(V1.module.disrupted.nodes.exp.mat.scaled$perm.group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score)), direction=">", quiet=T)$auc
  return(AUC)
}))

p.value <- sum(perm.AUC>AUC)/ncol(perm.label)


#Plot ROC
#Supp. Fig 5c, left

res.samp200 <- sapply(1:200, simplify=T, function(k) {
  samp <- sample(seq(nrow(V1.module.disrupted.nodes.exp.mat)), nrow(V1.module.disrupted.nodes.exp.mat), replace = T)
})

#Plot 
newcol <- colorRampPalette(c("grey80", "grey84"))
ncols <- 200
bluecols2 <- newcol(ncols)
cur.model <- Baseline.DEGs

cur.V1.module.disrupted.nodes.exp.mat <- V1.module.disrupted.nodes.exp.mat[res.samp100[,cur.samp],]
cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")] <- apply(cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group",  "sampleID")], 2, scale)
cur.V1.module.disrupted.nodes.exp.mat$scaled.module.score <- apply(cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")],1,mean)
plot(roc(cur.V1.module.disrupted.nodes.exp.mat$group, as.numeric(as.character(cur.V1.module.disrupted.nodes.exp.mat$scaled.module.score)),  smooth=T), col = "grey",lwd=1)

AUC.boot <- do.call('rbind', lapply(c(1:ncol(res.samp200)), function(cur.samp) {
  print(cur.samp)
  cur.V1.module.disrupted.nodes.exp.mat <- V1.module.disrupted.nodes.exp.mat[res.samp200[,cur.samp],]
  cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")] <- apply(cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group",  "sampleID")], 2, scale)
  cur.V1.module.disrupted.nodes.exp.mat$scaled.module.score <- apply(cur.V1.module.disrupted.nodes.exp.mat[,!colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("group", "sampleID")],1,mean)
  
  Gmean.difference.df <- data.frame(cur.V1.module.disrupted.nodes.exp.mat[,colnames(cur.V1.module.disrupted.nodes.exp.mat) %in% c("scaled.module.score", "group", "sampleID")], stringsAsFactors = F)
  Gmean.difference.df$group <- as.character(pData(filtered.esetALL)[match(Gmean.difference.df$sampleID, as.character(pData(filtered.esetALL)$sampleID)),'group'])
  Gmean.difference.df$group <- ifelse(Gmean.difference.df$group %in% "IM", "NR", Gmean.difference.df$group)
  Gmean.difference.df$group <- factor(Gmean.difference.df$group, levels=c("R", "NR"))
  res <- try(plot(roc(Gmean.difference.df$group, as.numeric(as.character(Gmean.difference.df$scaled.module.score)), smooth=T), col = bluecols2[cur.samp], add = TRUE,lwd=1))
  if(inherits(res, "try-error"))
  {
    print("error")
  } else {
    plot(roc(Gmean.difference.df$group, as.numeric(as.character(Gmean.difference.df$scaled.module.score)), smooth=T), col = bluecols2[cur.samp], add = TRUE,lwd=1)
  }
}))

plot(roc(Gmean.difference.df$group, Gmean.difference.df$scaled.module.score,  smooth=F), add=T, col = "red",lwd=2)


#Plot relative axis score
#Supp. Fig 5c, left

esetMat <- filtered.esetALL[fData(filtered.esetALL)$SYMBOL %in% Baseline.DEGs & fData(filtered.esetALL)$dataType %in% "AGX", pData(filtered.esetALL)$Visit %in% "V1" & !pData(filtered.esetALL)$sampleID %in% "HR-37-V1"]
V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
V1.module.disrupted.nodes.exp.mat$group <- pData(filtered.esetALL)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(filtered.esetALL)$sampleID),'group']
V1.module.disrupted.nodes.exp.mat$group <- ifelse(V1.module.disrupted.nodes.exp.mat$group %in% "IM", "NR", V1.module.disrupted.nodes.exp.mat$group)
V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))
V1.module.disrupted.nodes.exp.mat$sampleID <- row.names(V1.module.disrupted.nodes.exp.mat)
V1.module.disrupted.nodes.exp.mat.scaled <- V1.module.disrupted.nodes.exp.mat
V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")] <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group",  "sampleID")], 2, scale)
V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")],1,mean)


ggplot(V1.module.disrupted.nodes.exp.mat.scaled, aes(x = group, y = as.numeric(scaled.module.score), fill=group, colour="black"))+
  geom_boxplot()+
  scale_color_manual(values=c("black"="black"))+
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme_bw()+theme(strip.text = element_text(size=12), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("Relative axis score")+
  stat_compare_means(method = "wilcox",size=4, comparisons=list(c(1,2)), method.args = list(c("greater")))


#Test pathway enrichment in intermediate monocytes network 
#associated with Baseline.DEGs
#Plot Supp. Fig 7

cell.specific.net.enrich.func <- function(cur.clust ,cur.cond, clust.ref, percentile, selected.genes, enrich.dataset) {
  Idents(pbmc.IFX) <- 'orig.ident'
  pbmc.R <- subset(pbmc.IFX, idents = cur.cond)
  Idents(pbmc.R) <- clust.ref
  cells.selected = Idents(pbmc.R) == cur.clust
  print(paste("selected.cells:", cur.clust, sep=" "))
  my.pbmc.R = pbmc.R
  my.pbmc.R = my.pbmc.R[, cells.selected]
  
  #MTGOsc needs a folder where to save temporary files and final results.
  
  root = "~/Shiran/Remicade.analysis.July2018/all.features.together/IFX_paper/singleCell/MTGOres" #change this to your preferred local path
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
  
  #the list of all genes involved in the cluster
  genes = unique(c(as.character(edges.sub$gene1), as.character(edges.sub$gene2)))
  
  #translating gene names to ENTREZID via org.Mm.eg.db database
  genes = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  #the actual enrichment
  if(enrich.dataset %in% "Reactome") {
    enriched <-  enrichPathway(gene=genes$ENTREZID, pvalueCutoff=0.05, readable=T, organism = "human", pAdjustMethod = "BH")
    enriched@result <- enriched@result[enriched@result$p.adjust<0.05,]
    enriched@result <- enriched@result[enriched@result$Count>2,]
    g <- ReactomePA::dotplot(enriched, showCategory=20) + ggtitle(paste("Enrichment analysis", cur.cond, sep="."))
  } else {
    if(enrich.dataset %in% "KEGG") {
      enriched <- enrichKEGG(gene=genes$ENTREZID, pvalueCutoff=0.05, organism = "hsa", pAdjustMethod = "BH")
      enriched@result <- enriched@result[enriched@result$p.adjust<0.05,]
      enriched@result <- enriched@result[enriched@result$Count>2,]
      g <-clusterProfiler::dotplot(enriched, showCategory=20) + ggtitle(paste("Enrichment analysis", cur.cond, sep="."))
    } else {
      if(enrich.dataset %in% "GO") {
        enriched <-   enrichGO(gene=genes$ENTREZID, pvalueCutoff=0.05, OrgDb= org.Hs.eg.db, ont= "BP", pAdjustMethod = "BH", readable= TRUE)
        enriched@result <- enriched@result[enriched@result$p.adjust<0.05,]
        enriched@result <- enriched@result[enriched@result$Count>2,]
        g <-clusterProfiler::dotplot(enriched, showCategory=20) + ggtitle(paste("Enrichment analysis", cur.cond, sep="."))
      } else {
        print("choose enrichment database KEGG or Reactome")
      }}}
  if(is.null(enriched)) {
    return(list(net=edges.sub))} else {
      return(list(plot=g, data=enriched, net=edges.sub))}
  return(list(net=edges.sub))}


Intermediate.monocytes.network.KEGG <- cell.specific.net.enrich.func("Intermediate monocytes", c('IFX.R', 'IFX.NR'), "clusters.MonacoImm", 0.1, Baseline.DEGs, "KEGG")
Intermediate.monocytes.network.Reactome <- cell.specific.net.enrich.func("Intermediate monocytes", c('IFX.R', 'IFX.NR'), "clusters.MonacoImm", 0.1, Baseline.DEGs, "Reactome")
Intermediate.monocytes.network.GO <- cell.specific.net.enrich.func("Intermediate monocytes", c('IFX.R', 'IFX.NR'), "clusters.MonacoImm", 0.1, Baseline.DEGs, "GO")

#plot Heatmap for the differences between R-NR in the different subsets

Idents(pbmc.IFX) <- "clusters.MonacoImm"
pbmc.IFX.monocytes <- subset(pbmc.IFX, idents = c("Classical monocytes", "Intermediate monocytes", "Non classical monocytes"))

pathways <- Intermediate.monocytes.network.Reactome$data@result$Description
pathways <- pathways[order(pathways$qvalue), ]
pathways <- pathways[1:20,]
pathways <- unique(c(Intermediate.monocytes.network.Reactome$data@result$Description))
module.score.res <-   do.call('rbind', lapply(unique(Idents(pbmc.IFX.monocytes)), function(cur.clust) {
  tmp <- do.call('rbind', lapply(pathways, function(cur.path) {
    print(paste("choosing clust", cur.clust, sep="."))
    pbmc.IFX.monocytes.sub <- subset(pbmc.IFX.monocytes, idents = cur.clust)
    print(paste("choosing path related features", cur.path, sep=" "))
    features.in.cur.path.mono <- res.R.NR.Monocytes.sub[res.R.NR.Monocytes.sub$Description.x %in% cur.path, ]
    features.in.cur.path.momo <- unique(c(unlist(strsplit(features.in.cur.path.mono$geneID.x, "/")) ,unlist(strsplit(features.in.cur.path.mono$geneID.y, "/")) ))
    
    features.in.cur.path.class <- res.R.NR.classical.monocytes.sub[res.R.NR.classical.monocytes.sub$Description.x %in% cur.path, ]
    features.in.cur.path.class <- unique(c(unlist(strsplit(features.in.cur.path.class$geneID.x, "/")) ,unlist(strsplit(features.in.cur.path.class$geneID.y, "/")) ))
    
    features.in.cur.path.interm <- res.R.NR.Intermediate.monocytes.sub[res.R.NR.Intermediate.monocytes.sub$Description.x %in% cur.path, ]
    features.in.cur.path.interm<- unique(c(unlist(strsplit(features.in.cur.path.interm$geneID.x, "/")) ,unlist(strsplit(features.in.cur.path.interm$geneID.y, "/")) ))
    
    features.in.cur.path.Nonclass <- res.R.NR.Nonclassical.monocytes.sub[res.R.NR.Nonclassical.monocytes.sub$Description.x %in% cur.path, ]
    features.in.cur.path.Nonclass<- unique(c(unlist(strsplit(features.in.cur.path.Nonclass$geneID.x, "/")) ,unlist(strsplit(features.in.cur.path.Nonclass$geneID.y, "/")) ))
    
    features.in.cur.path <- unique(c(features.in.cur.path.momo,features.in.cur.path.class, features.in.cur.path.interm, features.in.cur.path.Nonclass))
    
    pbmc.IFX.monocytes.sub.path <- subset(pbmc.IFX.monocytes.sub, features=features.in.cur.path)
    print("start calculating scaled score")
    submat <- pbmc.IFX.monocytes.sub.path@assays$RNA@scale.data
    submat.scaled <- data.frame(t(submat), stringsAsFactors = F)
    module.score <- data.frame(apply(submat.scaled, 1, function(x) mean(x, na.rm=T)), stringsAsFactors = F)
    colnames(module.score) <- "mean.module.score"
    module.score$group <- pbmc.IFX.monocytes.sub.path@meta.data[match(row.names(module.score), row.names(pbmc.IFX.monocytes.sub.path@meta.data)), 'orig.ident']
    module.score$module.name <-cur.path
    module.score$clust <- cur.clust
    return(module.score)
  }))
  return(tmp)
}))


module.score.res$class <- paste(module.score.res$module.name, module.score.res$clust, sep='.')

module.score.wilcoxon <-  do.call('rbind', lapply(unique(module.score.res$class), function(cur.class) {
  submat <- module.score.res[module.score.res$class %in% cur.class,]
  mean.R <- mean(submat[submat$group %in% "IFX.R",'mean.module.score'], na.rm=T)
  mean.NR <- mean(submat[submat$group %in% "IFX.NR",'mean.module.score'], na.rm=T)
  p.value=wilcox.test(mean.module.score ~ group, paired=F, data=submat)$p.value 
  df <- data.frame(clust=cur.class, p.value=p.value, mean.R=mean.R, mean.NR=mean.NR, stringsAsFactors = F)
  return(df)
}))

module.score.wilcoxon$mono.clust <- ifelse(grepl("Classical monocytes", module.score.wilcoxon$clust), "Classical monocytes", 
                                           ifelse(grepl("Intermediate monocytes", module.score.wilcoxon$clust), "Intermediate monocytes",  
                                                  ifelse(grepl("Non classical monocytes", module.score.wilcoxon$clust), "Non classical monocytes", NA )))


module.score.wilcoxon.final <- do.call('rbind', lapply(unique(module.score.wilcoxon$mono.clust), function(cur.clust) {
  submat <- module.score.wilcoxon[module.score.wilcoxon$mono.clust %in% cur.clust,]
  submat$FDR <- p.adjust(submat$p.value, method="BH")
  return(submat)
}))

module.score.wilcoxon.final$clust.name <- sub("\\..*", "", module.score.wilcoxon.final$clust)
module.score.wilcoxon.final <- module.score.wilcoxon.final[order(module.score.wilcoxon.final$FDR),]

top.pathways <- unique(module.score.wilcoxon.final[1:20, 'clust.name'])

module.score.wilcoxon.final.for.heatmap <- module.score.wilcoxon.final[module.score.wilcoxon.final$clust.name %in% top.pathways,]
module.score.wilcoxon.final.for.heatmap$mean.R <- ifelse(module.score.wilcoxon.final.for.heatmap$FDR<0.05, module.score.wilcoxon.final.for.heatmap$mean.R, 0)
module.score.wilcoxon.final.for.heatmap$mean.NR <- ifelse(module.score.wilcoxon.final.for.heatmap$FDR<0.05, module.score.wilcoxon.final.for.heatmap$mean.NR, 0)
module.score.wilcoxon.final.for.heatmap$path <- gsub(".Classical monocytes|.Intermediate monocytes|.Non classical monocytes", "", module.score.wilcoxon.final.for.heatmap$clust)


module.score.wilcoxon.final.for.heatmap <- module.score.wilcoxon.final.for.heatmap[!module.score.wilcoxon.final.for.heatmap$path %in% c("PD-1 signaling", "Translocation of ZAP-70 to Immunological synapse"),]
module.score.wilcoxon.final.for.heatmap.wide.R <- reshape2::dcast(module.score.wilcoxon.final.for.heatmap, mono.clust ~ path, value.var="mean.R")
module.score.wilcoxon.final.for.heatmap.wide.R$group <- "R"

module.score.wilcoxon.final.for.heatmap.wide.NR <- reshape2::dcast(module.score.wilcoxon.final.for.heatmap, mono.clust ~ path, value.var="mean.NR")
module.score.wilcoxon.final.for.heatmap.wide.NR$group <- "NR"

module.score.wilcoxon.final.for.heatmap.wide.R.NR <- rbind.fill(module.score.wilcoxon.final.for.heatmap.wide.R, module.score.wilcoxon.final.for.heatmap.wide.NR)
row.names(module.score.wilcoxon.final.for.heatmap.wide.R.NR) <- paste(module.score.wilcoxon.final.for.heatmap.wide.R.NR$mono.clust, module.score.wilcoxon.final.for.heatmap.wide.R.NR$group, sep=".")




module.score.wilcoxon.final.for.heatmap.wide.R.NR_t <- data.frame(t(module.score.wilcoxon.final.for.heatmap.wide.R.NR), stringsAsFactors = F)

pathway.Intermediate.net.for.heatmap <- module.score.wilcoxon.final.for.heatmap.wide.R.NR_t[!row.names(module.score.wilcoxon.final.for.heatmap.wide.R.NR_t) %in% c("mono.clust", "group"),]
pathway.Intermediate.net.for.heatmap.final <- apply(pathway.Intermediate.net.for.heatmap, 2, as.numeric)
row.names(pathway.Intermediate.net.for.heatmap.final) <- row.names(pathway.Intermediate.net.for.heatmap)
colnames(pathway.Intermediate.net.for.heatmap.final) <- c("R", "R", "R", "NR", "NR", "NR")

pathway.Intermediate.net.for.heatmap.final <- pathway.Intermediate.net.for.heatmap.final[,c(1,4,2,5,3,6)]

anno.df <- anno.df[c(1,4,2,5,3,6),]
anno.df <- data.frame(mono.clust=anno.df, stringsAsFactors = F)

ha1 <-  HeatmapAnnotation(df =anno.df ,
                          col = list(mono.clust =c("Classical monocytes" = "coral", "Intermediate monocytes" = "orange", "Non classical monocytes"="Gold")))




ht <- Heatmap(pathway.Intermediate.net.for.heatmap.final, 
              top_annotation = ha1, rect_gp = gpar(col = "black", lwd = 2), row_title = NULL,
              cluster_columns =F, 
              col=colorRamp2(c(-1, -0.5, 0,0.5,1), c( "#084594", "#9ECAE1", "white", "#f0e3ab", "#b23030")),
              name="module score scaled log2 TP10k+1", row_names_gp = gpar(fontsize = 12), 
              column_names_gp = gpar(fontsize = 10), show_row_names = TRUE, show_column_names = T,
              column_names_rot = 0, column_names_max_height = unit(10,"cm"), row_names_max_width = unit(30,"cm"))

draw(ht, padding = unit(c(2, 10, 2, 20), "mm"), heatmap_legend_side = "left", annotation_legend_side = "left") #bottom, left, top, right paddings


#Test mTNF on intermediate monocytes as measured by CyTOF
#Figure 4D

Monocytes.homing.molecules <- read.csv(paste(dataDir, "/Monocytes.subsets.homing.molecules.csv", sep=""),stringsAsFactors = F)
Monocytes.homing.molecules.long <- reshape2::melt(Monocytes.homing.molecules,  id.vars=c("marker", "cell.subset", "sampleID"))
Monocytes.homing.molecules.long$variable <- gsub("\\.", "-", Monocytes.homing.molecules.long$variable)
Monocytes.homing.molecules.long$group <- pData(filtered.esetALL)[match(Monocytes.homing.molecules.long$variable, pData(filtered.esetALL)$sampleID),'group']
Monocytes.homing.molecules.long$group <- ifelse(Monocytes.homing.molecules.long$group %in% "IM", "NR", Monocytes.homing.molecules.long$group)

Monocytes.homing.molecules.long <- Monocytes.homing.molecules.long[!is.na(Monocytes.homing.molecules.long$group),]

Monocytes.homing.molecules.long$Visit <- pData(filtered.esetALL)[match(Monocytes.homing.molecules.long$variable, pData(filtered.esetALL)$sampleID),'Visit']

Monocytes.homing.molecules.long$class <- paste(Monocytes.homing.molecules.long$group, Monocytes.homing.molecules.long$Visit, sep=".")
Monocytes.homing.molecules.long$marker <- ifelse(Monocytes.homing.molecules.long$marker %in% "Humira", "mTNF", Monocytes.homing.molecules.long$marker)

Monocytes.homing.molecules.long$marker <- factor(Monocytes.homing.molecules.long$marker, levels=c("mTNF", "IntegrinB7", "Remicade", "CCR9"))

Monocytes.homing.molecules.long$class <- paste(Monocytes.homing.molecules.long$cell.subset, Monocytes.homing.molecules.long$group, sep=".")
Monocytes.homing.molecules.long$class <- factor(Monocytes.homing.molecules.long$class, levels=c("classical.monocytes.R", "classical.monocytes.NR", "intermediate.monocytes.R", "intermediate.monocytes.NR", "non-classical.monocytes.R", "non-classical.monocytes.NR"))
Monocytes.homing.molecules.long$group <- factor(Monocytes.homing.molecules.long$group, levels=c("R", "NR"))
Monocytes.homing.molecules.long.sub <- Monocytes.homing.molecules.long[Monocytes.homing.molecules.long$marker %in% "mTNF",]

ggplot(Monocytes.homing.molecules.long.sub, aes(x = cell.subset, y = as.numeric(value), fill=group)) +
  geom_boxplot() +
  ylab("CyTOF intensity")+
  facet_wrap(~marker, scales = "free", ncol=6)+
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0), 
        strip.text = element_text(size=10), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))+ 
  stat_compare_means(data=Monocytes.homing.molecules.long.sub ,method = "wilcox.test", size=4, comparisons=list(c(2,1), c(2,3)), aes(group = group), method.args = list(c("greater")), label.y = c(19,21))+
  stat_compare_means(data=Monocytes.homing.molecules.long.sub , aes(group = group, label=..p.format..),
                     method = "wilcox.test", method.args = list(alternative = "less"))


#Plot disrupted pathways dynamics - fiber organization pathway
#Fig 4b

V1.module.disrupted.nodes.wilcoxon.sig <- V1.module.disrupted.nodes.wilcoxon[V1.module.disrupted.nodes.wilcoxon$FDR<0.2,]
cur.path <- "GO_POSITIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION"
node.list <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
node.list.sig <- node.list[node.list %in% V1.module.disrupted.nodes.wilcoxon.sig$gene]
esetMat.eset <- filtered.esetALL[fData(filtered.esetALL)$featureID %in% node.list.sig,!pData(filtered.esetALL)$Patient.code %in% "37"]
expMat <- data.frame(exprs(esetMat.eset), stringsAsFactors = F)
expMat <- data.frame(t(expMat), stringsAsFactors = F)
expMat.scaled <- apply(expMat, 2, scale)
row.names(expMat.scaled) <- row.names(expMat)
row.names(expMat.scaled) <- gsub("\\.", "-", row.names(expMat.scaled))
expMat.scaled <- data.frame(expMat.scaled, stringsAsFactors = F)
expMat.scaled$group <- pData(esetMat.eset)[match(row.names(expMat.scaled), pData(esetMat.eset)$sampleID),'group']
expMat.scaled$Visit <- pData(esetMat.eset)[match(row.names(expMat.scaled), pData(esetMat.eset)$sampleID),'Visit']
expMat.scaled$group <- ifelse(expMat.scaled$group %in% "IM", "NR", expMat.scaled$group)
expMat.scaled$module.score <- apply(expMat.scaled[,!colnames(expMat.scaled) %in% c("group", "Visit")], 1,  function(x) mean(x, na.rm=T))

fiber.organization.pathway.score <- data.frame(path=cur.path, module.score=expMat.scaled$module.score, group=expMat.scaled$group, Visit=expMat.scaled$Visit, sampleID=row.names(expMat.scaled),stringsAsFactors = F)

ggplot(fiber.organization.pathway.score, aes(x = Visit, y = module.score, colour = group, group=group, width=0.2)) + 
  stat_summary(fun = "mean", geom = "line", size =1.5)+ # adding connecting lines
  stat_summary(fun = "mean", geom = "point") +        # adding data points
  stat_summary(fun.data = "mean_se", geom = "errorbar", size=1)+
  theme_bw()+theme(strip.text = element_text(size=12), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=12),
                   axis.text.x=element_text(size=12,hjust=1),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("Scaled axis score")+
  stat_compare_means(data=fiber.organization.pathway.score[fiber.organization.pathway.score$group %in% "R",], method = "wilcox.test", method.args = list(alternative = "greater"), comparisons=list(c(1,2), c(1,3)), label.y=c(1.7, 1.9), color="cyan3")+
  stat_compare_means(data=fiber.organization.pathway.score[fiber.organization.pathway.score$group %in% "NR",], method = "wilcox.test", method.args = list(alternative = "less"), comparisons=list(c(1,3), c(1,2)), label.y=c(1, 1.2), color="Indianred")+
  stat_compare_means(data=fiber.organization.pathway.score[fiber.organization.pathway.score$Visit %in% "V1",] ,aes(group = group, label=sub("p = ", "", as.numeric(..p.format..))), method = "wilcox.test", method.args = list(alternative = "greater"),
                     label.y = c(0.8))



#Test prediction of the fiber-organization signature in RA 

processing.function <- function(gset) {
  gset <- getGEO(gset, GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
  
  if(!is_logscale(exprs(gset))){
    exprs(gset) <- log2(exprs(gset))
  }
  
  #test if data was quantile normelized
  test <- data.frame(exprs(gset), stringsAsFactors=F)
  test$featureID <- row.names(test)
  test.long <- reshape2::melt(test, id.vars="featureID")
  
  
  g <- ggplot(test.long, aes(x = variable, y = as.numeric(value))) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
    theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),                                                                  
                     strip.background = element_blank(),
                     panel.border = element_rect(colour = "black"), 
                     axis.ticks.length=unit(.15, "cm"))
  plot(g)
  
  return(gset)
}


processing.function.imputation <- function(gset) {
  gset.exprs <- exprs(gset)
  gset.exprs <- data.frame(gset.exprs, stringsAsFactors = F)
  gset.exprs$SYMBOL <- fData(gset)[match(row.names(gset.exprs), fData(gset)$ID),'Gene symbol']
  gset.exprs$mean <- as.numeric(apply(gset.exprs[,!colnames(gset.exprs) %in% 'SYMBOL'], 1, function(x) mean(x, na.rm=T)))
  gset.exprs$ID <- row.names(gset.exprs)
  
  gset.exprs <- gset.exprs[!gset.exprs$SYMBOL %in% "",]
  
  gset.exprs.filtered <- gset.exprs %>%
    group_by(SYMBOL) %>% 
    dplyr::slice(which.max(mean))
  
  gset.exprs.filtered <- as.data.frame(gset.exprs.filtered, stringsAsFactors=F)
  
  gset <- gset[fData(gset)$ID %in% gset.exprs.filtered$ID,]
  
  
  NA.perc = apply(exprs(gset), 1, function(x) sum(is.na(x))/length(x)*100)
  gset<-gset[!NA.perc>20,]
  
  #For the remaining NAs put average value of the feature
  
  Features.with.NA = which(apply(exprs(gset), 1, function(x) any(is.na(x))))
  test <- do.call('rbind',lapply(Features.with.NA, function(row.index) {
    exprs(gset)[row.index,is.na(exprs(gset)[row.index,])] <- mean(exprs(gset)[row.index,], na.rm=T)
    return(exprs(gset)[row.index,])
  }))
  
  exprs(gset)[Features.with.NA,]<- test
  return(gset)
  
}

cellMix.deconv.and.adjustment.function <- function(gset) {
  fData(gset)$hgu133.prob_id <- ENTREZ.probeid.hgu133[match(fData(gset)$`Gene ID`,ENTREZ.probeid.hgu133$gene_id),'probe_id']
  for.deconv.sig.comp <- exprs(gset)
  names <- fData(gset)[match(row.names(for.deconv.sig.comp),fData(gset)$ID),'hgu133.prob_id']
  row.names(for.deconv.sig.comp) <- names
  res <- gedBlood(for.deconv.sig.comp, rescale=T)
  
  P.cellmix <- res@fit@H
  P.cellmix <- rbind(P.cellmix, P.cellmix['neutro',])
  row.names(P.cellmix)[nrow(P.cellmix)] <- 'Gran'
  P.cellmix <- rbind(P.cellmix, P.cellmix['mono',]+P.cellmix['mono act',])
  row.names(P.cellmix)[nrow(P.cellmix)] <- 'Mono'
  
  #adjust data by major celltypes
  
  # select and compute major cell types
  major <- c("Th", "Tc", "B", "Mono", "Gran", "NK")
  P.P.cellmix.sub <- P.cellmix[major, ]
  
  # correct gene expression within each visit separately
  
  # for V1
  gset.exprs <- exprs(gset)
  P.cellmix.sub.V1 <- P.P.cellmix.sub[,colnames(P.P.cellmix.sub) %in% colnames(gset.exprs)]
  
  plot.abundance <- data.frame(t(P.cellmix.sub.V1), stringsAsFactors = F)
  plot.abundance$group <- pData(gset)[match(row.names(plot.abundance), pData(gset)$geo_accession),'group']
  
  plot.abundance$group <- factor(plot.abundance$group, levels=c("R", "NR"))
  
  g1 <- ggplot(plot.abundance, aes(x = group, y = as.numeric(Mono), fill=group)) +
    geom_boxplot( outlier.shape = NA) +
    scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
    theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),                                                                  
                     strip.background = element_blank(),
                     panel.border = element_rect(colour = "black"), 
                     axis.ticks.length=unit(.15, "cm"))+
    stat_compare_means(method = "wilcox.test", comparisons=list(c(1,2)))
  
  
  #using cellmix#
  AG.V1.cellmix <- ged(gset, as.matrix(P.cellmix.sub.V1) ,method = 'correction')
  AG.V1.cellmix <- log2(exprs(AG.V1.cellmix))
  
  ##incorporate to expressionset
  gxeset_adjust <- gset
  featureNames(gxeset_adjust) <- sub("GX.", "", featureNames(gxeset_adjust))
  row.names(AG.V1.cellmix) <- sub("^GX\\.", "", row.names(AG.V1.cellmix))
  exprs(gxeset_adjust) <- AG.V1.cellmix[featureNames(gxeset_adjust), sampleNames(gxeset_adjust)]
  fData(gxeset_adjust)$dataType <- 'AGX'
  fData(gset)$dataType <- 'GX'
  featureNames(gset) <- sub("^GX\\.", "",  featureNames(gset))
  
  ## combine all data types
  ##combine.expression.data
  exprs.data <- rbind.fill(data.frame(exprs(gset), stringsAsFactors = FALSE), 
                           data.frame(exprs(gxeset_adjust), stringsAsFactors = FALSE))
  
  
  row.names(exprs.data) <- c( paste("GX", row.names(exprs(gset)), sep="."),
                              paste("AG", row.names(exprs(gxeset_adjust)), sep="."))
  
  
  
  exprs.data <- as.matrix(exprs.data) 
  colnames(exprs.data) <- gsub("\\.", "-",  colnames(exprs.data))
  
  ##combine fData
  
  fData(gset)$featureID <- paste("GX", row.names(exprs(gset)), sep=".")
  fData(gxeset_adjust)$featureID <- paste("AG", row.names(exprs(gxeset_adjust)), sep=".")
  
  feature.data <- rbind.fill(data.frame(fData(gset), stringsAsFactors = FALSE), 
                             data.frame(fData(gxeset_adjust), stringsAsFactors = FALSE))
  
  row.names(feature.data) <- feature.data$featureID
  
  gset.cellmix <- ExpressionSet(exprs.data, phenoData = AnnotatedDataFrame(pData(gset))
                                , featureData = AnnotatedDataFrame(feature.data)
                                , annotation = 'Analysis')
  
  
  return(list(eset=gset.cellmix, plot=g1))
  
}



#load GEO20690

gset.GEO20690 <- processing.function("GSE20690")

#data should be normalized
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

test <- quantile_normalisation(2^exprs(gset.GEO20690))
test <- log2(test)

test <- data.frame(test, stringsAsFactors=F)
test$featureID <- row.names(test)
test.long <- reshape2::melt(test, id.vars="featureID")


ggplot(test.long, aes(x = variable, y = as.numeric(value))) +
  geom_boxplot(alpha=0.4, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))

test$featureID <- NULL
test <- as.matrix(test)
exprs(gset.GEO20690) <- test


pData(gset.GEO20690)$group <- ifelse(grepl("Residual Inflammation", pData(gset.GEO20690)$title), "NR", "R")

gset.GEO20690 <- processing.function.imputation(gset.GEO20690)
esetALL.GEO20690 <- cellMix.deconv.and.adjustment.function(gset.GEO20690)$eset


#plots for adjusted GX the expression of signature-related genes
#Fig 5b

nodes.in.relevant.pathways <- Baseline.DEGs
esetALL.GEO20690.nodes <- esetALL.GEO20690[fData(esetALL.GEO20690)$Gene.symbol %in% as.character(nodes.in.relevant.pathways) & fData(esetALL.GEO20690)$dataType %in% "AGX",]
esetALL.GEO20690.nodes.long <- reshape2::melt(exprs(esetALL.GEO20690.nodes))
esetALL.GEO20690.nodes.long$SYMBOL <- fData(esetALL.GEO20690)[match(esetALL.GEO20690.nodes.long$Var1, fData(esetALL.GEO20690)$featureID),'Gene.symbol']
esetALL.GEO20690.nodes.long$dataType <- fData(esetALL.GEO20690)[match(esetALL.GEO20690.nodes.long$Var1, fData(esetALL.GEO20690)$featureID),'dataType']
esetALL.GEO20690.nodes.long$Var1displayName <- paste(esetALL.GEO20690.nodes.long$SYMBOL, esetALL.GEO20690.nodes.long$dataType, sep=".")
esetALL.GEO20690.nodes.long$group <- pData(esetALL.GEO20690)[match(esetALL.GEO20690.nodes.long$Var2, pData(esetALL.GEO20690)$geo_accession),'group']
esetALL.GEO20690.nodes.long$group <- factor(esetALL.GEO20690.nodes.long$group, levels=c("R", "NR"))

esetALL.GEO20690.nodes.long$Var1displayName <- factor(esetALL.GEO20690.nodes.long$Var1displayName, levels=c("PAK1.AGX", "LYN.AGX", "ICAM1.AGX", "FCGR3A.AGX", "IL1B.AGX"))

ggplot(esetALL.GEO20690.nodes.long, aes(x = group, y = as.numeric(value), fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  #coord_fixed()+
  facet_wrap(~Var1displayName, scales = "free_y", ncol=6)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox.test", comparisons=list(c(1,2)), method.args = list(c("greater")))

#Plot relative signature score and ROC
#Fig 5b, Right

nodes.in.relevant.pathways <- Baseline.DEGs
esetMat <- esetALL.GEO20690[fData(esetALL.GEO20690)$'Gene.symbol' %in% nodes.in.relevant.pathways & fData(esetALL.GEO20690)$dataType %in% "AGX", ]
V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
V1.module.disrupted.nodes.exp.mat$group <- pData(esetALL.GEO20690)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(esetALL.GEO20690)$geo_accession),'group']
V1.module.disrupted.nodes.exp.mat$group <- ifelse(V1.module.disrupted.nodes.exp.mat$group %in% "IM", "NR", V1.module.disrupted.nodes.exp.mat$group)
V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))
V1.module.disrupted.nodes.exp.mat$sampleID <- row.names(V1.module.disrupted.nodes.exp.mat)
V1.module.disrupted.nodes.exp.mat.scaled <- V1.module.disrupted.nodes.exp.mat
V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")] <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group",  "sampleID")], 2, scale)
V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")],1,mean)

ggplot(V1.module.disrupted.nodes.exp.mat.scaled, aes(x = group, y = as.numeric(scaled.module.score), fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  ylab("Relative pathway score")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "greater"), comparisons=list(c(1,2)))

roc.data.GEO20690 <- roc(V1.module.disrupted.nodes.exp.mat.scaled$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score)),  grid=TRUE, plot=F)
roc.data.GEO20690.tp.fp <- coords(roc.data.GEO20690,  ret=c("tp", "fp", "tn", "fn"), transpose = F)
roc.data.GEO20690.tp.fp <- data.frame(roc.data.GEO20690.tp.fp, stringsAsFactors = F)
Author <- c(rep("GEO20690", nrow(roc.data.GEO20690.tp.fp)))
roc.data.GEO20690.tp.fp <- cbind(Author, roc.data.GEO20690.tp.fp)

plot(roc(V1.module.disrupted.nodes.exp.mat.scaled$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score)),  smooth=T))


#load GSE33377

gset.GSE33377 <- processing.function("GSE33377")
Annot <- data.frame(SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=","),
                    DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=","),
                    ENSEMBLID=sapply(contents(huex10sttranscriptclusterENSEMBL), paste, collapse=","), 
                    ENTREZID=sapply(contents(huex10sttranscriptclusterENTREZID), paste, collapse=","))


fData(gset.GSE33377)$'Gene symbol' <- Annot[match(fData(gset.GSE33377)$ID, row.names(Annot)) , 'SYMBOL']
fData(gset.GSE33377)$'Gene ID' <- Annot[match(fData(gset.GSE33377)$ID, row.names(Annot)) , 'ENTREZID']


gset.GSE33377 <- processing.function.imputation(gset.GSE33377)

pData(gset.GSE33377)$group <- ifelse(pData(gset.GSE33377)$'response:ch1' %in% "anti-TNF non-responder", "NR", "R")
gset.GSE33377.cellmix <- cellMix.deconv.and.adjustment.function(gset.GSE33377)

gset.GSE33377.cellmix.final <- gset.GSE33377.cellmix$eset
esetALL.GSE33377 <- gset.GSE33377.cellmix.final

gset.GSE33377.cellmix.nodes <- gset.GSE33377.cellmix.final[fData(gset.GSE33377.cellmix.final)$Gene.symbol %in% as.character(nodes.in.relevant.pathways) & fData(gset.GSE33377.cellmix.final)$dataType %in% "AGX",]
gset.GSE33377.cellmix.nodes.long <- reshape2::melt(exprs(gset.GSE33377.cellmix.nodes))
gset.GSE33377.cellmix.nodes.long$SYMBOL <- fData(gset.GSE33377.cellmix.final)[match(gset.GSE33377.cellmix.nodes.long$Var1, fData(gset.GSE33377.cellmix.final)$featureID),'Gene.symbol']
gset.GSE33377.cellmix.nodes.long$dataType <- fData(gset.GSE33377.cellmix.final)[match(gset.GSE33377.cellmix.nodes.long$Var1, fData(gset.GSE33377.cellmix.final)$featureID),'dataType']
gset.GSE33377.cellmix.nodes.long$Var1displayName <- paste(gset.GSE33377.cellmix.nodes.long$SYMBOL, gset.GSE33377.cellmix.nodes.long$dataType, sep=".")
gset.GSE33377.cellmix.nodes.long$group <- pData(gset.GSE33377.cellmix.final)[match(gset.GSE33377.cellmix.nodes.long$Var2, pData(gset.GSE33377.cellmix.final)$geo_accession),'group']
gset.GSE33377.cellmix.nodes.long$group <- factor(gset.GSE33377.cellmix.nodes.long$group, levels=c("R", "NR"))

ggplot(gset.GSE33377.cellmix.nodes.long, aes(x = group, y = as.numeric(value), fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  facet_wrap(~Var1displayName, scales = "free_y", ncol=10)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox.test", size=3, comparisons=list(c(1,2)), method.args = list(c("greater")))


esetMat <- gset.GSE33377.cellmix.final[fData(gset.GSE33377.cellmix.final)$'Gene.symbol' %in% nodes.in.relevant.pathways & fData(gset.GSE33377.cellmix.final)$dataType %in% "AGX", ]
V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
V1.module.disrupted.nodes.exp.mat$group <- pData(gset.GSE33377.cellmix.final)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(gset.GSE33377.cellmix.final)$geo_accession),'group']
V1.module.disrupted.nodes.exp.mat$group <- ifelse(V1.module.disrupted.nodes.exp.mat$group %in% "IM", "NR", V1.module.disrupted.nodes.exp.mat$group)
V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))
V1.module.disrupted.nodes.exp.mat$sampleID <- row.names(V1.module.disrupted.nodes.exp.mat)
V1.module.disrupted.nodes.exp.mat.scaled <- V1.module.disrupted.nodes.exp.mat
V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")] <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group",  "sampleID")], 2, scale)
V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")],1,mean)

ggplot(V1.module.disrupted.nodes.exp.mat.scaled, aes(x = group, y = as.numeric(scaled.module.score), fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox.test", size=3, comparisons=list(c(1,2)), method.args = list(alternative = "greater") )

plot(roc(V1.module.disrupted.nodes.exp.mat.scaled$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score)),  smooth=T))

roc.data.GSE33377 <- roc(V1.module.disrupted.nodes.exp.mat.scaled$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score)),  grid=TRUE, plot=F)
roc.data.GSE33377.tp.fp <- coords(roc.data.GSE33377,  ret=c("tp", "fp", "tn", "fn"), transpose = F)
roc.data.GSE33377.tp.fp <- data.frame(roc.data.GSE33377.tp.fp, stringsAsFactors = F)
Author <- c(rep("GSE33377", nrow(roc.data.GSE33377.tp.fp)))
roc.data.GSE33377.tp.fp <- cbind(Author, roc.data.GSE33377.tp.fp)


#GSE42296

gset.GSE42296 <- processing.function("GSE42296")

gset.GSE42296 <- processing.function.imputation(gset.GSE42296)
gset.GSE42296.sub <- gset.GSE42296[,pData(gset.GSE42296)$'disease state:ch1' %in% "Rheumatoid arthritis" & pData(gset.GSE42296)$'treatment:ch1' %in% "Before treatment"]
pData(gset.GSE42296.sub)$group <- ifelse(pData(gset.GSE42296.sub)$'response:ch1' %in% "NR - non-responder", "NR", "R")

#Deconvolution
gset.GSE42296.cellmix <- cellMix.deconv.and.adjustment.function(gset.GSE42296.sub)
esetALL.GSE42296 <- gset.GSE42296.cellmix$eset

gset.GSE42296.cellmix.nodes <- esetALL.GSE42296[fData(esetALL.GSE42296)$Gene.symbol %in% as.character(nodes.in.relevant.pathways) & fData(esetALL.GSE42296)$dataType %in% "AGX",]
gset.GSE42296.cellmix.nodes.long <- reshape2::melt(exprs(gset.GSE42296.cellmix.nodes))
gset.GSE42296.cellmix.nodes.long$SYMBOL <- fData(esetALL.GSE42296)[match(gset.GSE42296.cellmix.nodes.long$Var1, fData(esetALL.GSE42296)$featureID),'Gene.symbol']
gset.GSE42296.cellmix.nodes.long$dataType <- fData(esetALL.GSE42296)[match(gset.GSE42296.cellmix.nodes.long$Var1, fData(esetALL.GSE42296)$featureID),'dataType']
gset.GSE42296.cellmix.nodes.long$Var1displayName <- paste(gset.GSE42296.cellmix.nodes.long$SYMBOL, gset.GSE42296.cellmix.nodes.long$dataType, sep=".")
gset.GSE42296.cellmix.nodes.long$group <- pData(esetALL.GSE42296)[match(gset.GSE42296.cellmix.nodes.long$Var2, pData(esetALL.GSE42296)$geo_accession),'group']
gset.GSE42296.cellmix.nodes.long$group <- factor(gset.GSE42296.cellmix.nodes.long$group, levels=c("R", "NR"))

ggplot(gset.GSE42296.cellmix.nodes.long, aes(x = group, y = as.numeric(value), fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  facet_wrap(~Var1displayName, scales = "free_y", ncol=7)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox.test", comparisons=list(c(1,2)), method.args = list(alternative = "greater"))

esetMat <- esetALL.GSE42296[fData(esetALL.GSE42296)$'Gene.symbol' %in% nodes.in.relevant.pathways & fData(esetALL.GSE42296)$dataType %in% "AGX", ]
V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
V1.module.disrupted.nodes.exp.mat$group <- pData(esetALL.GSE42296)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(esetALL.GSE42296)$geo_accession),'group']
V1.module.disrupted.nodes.exp.mat$group <- ifelse(V1.module.disrupted.nodes.exp.mat$group %in% "IM", "NR", V1.module.disrupted.nodes.exp.mat$group)
V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))
V1.module.disrupted.nodes.exp.mat$sampleID <- row.names(V1.module.disrupted.nodes.exp.mat)
V1.module.disrupted.nodes.exp.mat.scaled <- V1.module.disrupted.nodes.exp.mat
V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")] <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group",  "sampleID")], 2, scale)
V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")],1,mean)

ggplot(V1.module.disrupted.nodes.exp.mat.scaled, aes(x = group, y = as.numeric(scaled.module.score), fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "greater"), comparisons=list(c(1,2)))

plot(roc(V1.module.disrupted.nodes.exp.mat.scaled$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score)),  smooth=T))

roc.data.GSE42296 <- roc(V1.module.disrupted.nodes.exp.mat.scaled$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score)),  grid=TRUE, plot=F)
roc.data.GSE42296.tp.fp <- coords(roc.data.GSE42296,  ret=c("tp", "fp", "tn", "fn"), transpose = F)
roc.data.GSE42296.tp.fp <- data.frame(roc.data.GSE42296.tp.fp, stringsAsFactors = F)
Author <- c(rep("GSE42296", nrow(roc.data.GSE42296.tp.fp)))
roc.data.GSE42296.tp.fp <- cbind(Author, roc.data.GSE42296.tp.fp)


#Combine all ROC into one metaROC
#Fig 5c

combined.tp.fp <- rbind(roc.data.GEO20690.tp.fp,roc.data.GSE33377.tp.fp, roc.data.GSE42296.tp.fp)
colnames(combined.tp.fp) <- c("Author", "TP", "FP", "TN", "FN")
output1 <- metaROC(combined.tp.fp, plot.bands=TRUE, plot.Author=F, col.curve='black', col.bands='lightgrey', col.border='black', cex.Author=0, lwd.Author=0.01)
summary.curve <- data.frame(fpr=output1$t, tpr=output1$sRA)
Author <- c(rep("Summary", nrow(summary.curve)))
summary.curve <- cbind(Author, summary.curve)
summary.curve$CI_upper <- output1$sRA+output1$se.RA
summary.curve$CI_lower <- output1$sRA-output1$se.RA

#add to ggplot
combined.tp.fp.final <- rbind.fill(summary.curve, output1$data)
combined.tp.fp.final$size <- factor(ifelse(combined.tp.fp.final$Author %in% "Summary", 1, 0.1))
combined.tp.fp.final$CI_upper <- ifelse(is.na(combined.tp.fp.final$CI_upper), 0, combined.tp.fp.final$CI_upper)
combined.tp.fp.final$CI_lower <- ifelse(is.na(combined.tp.fp.final$CI_lower), 0, combined.tp.fp.final$CI_lower)

ggplot(combined.tp.fp.final, aes(x = fpr, y = tpr, fill=Author, color=Author)) +
  geom_line(aes(size = size))+
  scale_size_manual(values = c(0.1,1), guide = "none")+
  xlab("False-Positive Rate")+ylab("True-Positive Rate")+
  geom_ribbon(data=combined.tp.fp.final[combined.tp.fp.final$Author %in% "Summary",], aes(ymin=CI_lower, ymax=CI_upper), fill="black", alpha=0.1, colour = NA) +
  scale_color_manual(values=c("Summary"="brown2", "GEO20690"="darkgrey", "GSE33377"="darkgrey", "GSE42296"="darkgrey"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0), 
        strip.text = element_text(size=10), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))


#Validation with additional IBD- real-life independent cohort using qPCR
#Fig 5a


#Use CD14 (monocytes) normalized expression values 

RT.PCR.HCK.RAC.PAK <- read.csv(paste(dataDir, "/IFX_CD14.normelized_signature_qPCR_average_R.csv", sep=""), stringsAsFactors = F)
RT.PCR.HCK.RAC.PAK.IFX <- RT.PCR.HCK.RAC.PAK[RT.PCR.HCK.RAC.PAK$Drug %in% "IFX",]

RT.PCR.HCK.RAC.PAK.IFX$module.score <-  apply(apply(RT.PCR.HCK.RAC.PAK.IFX[,!colnames(RT.PCR.HCK.RAC.PAK.IFX) %in% c("sampleID", "group", "Patien.code", "Drug")], 2, scale),1,mean)

RT.PCR.HCK.RAC.PAK.IFX.long <- reshape2::melt(RT.PCR.HCK.RAC.PAK.IFX, id.vars=c("group", "sampleID", "Patien.code", "Drug"))

RT.PCR.HCK.RAC.PAK.IFX.long$group <- factor(RT.PCR.HCK.RAC.PAK.IFX.long$group, levels=c("R", "NR"))

RT.PCR.HCK.RAC.PAK.IFX.long$variable <- factor(RT.PCR.HCK.RAC.PAK.IFX.long$variable, levels=c("PAK1", "LYN", "ICAM1", "FCGR3a", "IL1..", "RAC1", "CD14", "cd66b", "module.score"))
ggplot(RT.PCR.HCK.RAC.PAK.IFX.long[!RT.PCR.HCK.RAC.PAK.IFX.long$variable %in% "module.score",], aes(x = group, y = as.numeric(value), fill=group, colour="black"))+
  facet_wrap(~variable, scales = "free_y", ncol=7)+
  stat_summary(fun.data=function(...) mean_se(..., mult=1), 
               geom='errorbar', size=0.2, width=0.4, position="dodge")+
  stat_summary(fun=mean,  geom='bar', position="dodge")+
  scale_color_manual(values=c("black"="black"))+
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  scale_y_continuous(name="Relative expression", limits=c(-1,6.5), breaks=c(0, 0.5, 1.0, 1.5, 2.0))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("Relative expression")+
  stat_compare_means(method = "wilcox.test",size=3.5, method.args = list(alternative = "greater"), comparisons=list(c(1,2)))


ggplot(RT.PCR.HCK.RAC.PAK.IFX.long[RT.PCR.HCK.RAC.PAK.IFX.long$variable %in% c("module.score"),], aes(x = group, y = as.numeric(value), fill=group, colour="black"))+
  facet_wrap(~variable, scales = "free_y", ncol=7)+
  geom_boxplot(outlier.color = NA)+
  scale_color_manual(values=c("black"="black"))+
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("Relative pathway score")+
  stat_compare_means(method = "wilcox.test",size=3.5, method.args = list(alternative = "greater"), comparisons=list(c(1,2)), label.y=1.2)



##Test ROC by scaled module score

roc(controls=RT.PCR.HCK.RAC.PAK.IFX$module.score[RT.PCR.HCK.RAC.PAK.IFX$group=="R"], cases=RT.PCR.HCK.RAC.PAK.IFX$module.score[RT.PCR.HCK.RAC.PAK.IFX$group=="NR"])$auc
plot(roc(controls=RT.PCR.HCK.RAC.PAK.IFX$module.score[RT.PCR.HCK.RAC.PAK.IFX$group=="R"], cases=RT.PCR.HCK.RAC.PAK.IFX$module.score[RT.PCR.HCK.RAC.PAK.IFX$group=="NR"], smooth=T))


#Test prediction in Vedolizumab whole blood cohorts
#Supp Fig 8

setwd("~/Shiran/Vedolizumab.all.pipline030518/second batch")

Vedo.combined.eset.final <- readRDS(paste(dataDir, "/Vedo.combined.eset.final.rds", sep=""))
nodes.in.relevant.pathways <- Baseline.DEGs
esetMat <- Vedo.combined.eset.final[fData(Vedo.combined.eset.final)$SYMBOL %in% nodes.in.relevant.pathways & fData(Vedo.combined.eset.final)$dataType %in% "AGX", pData(Vedo.combined.eset.final)$Visit %in% "V1"]
V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
V1.module.disrupted.nodes.exp.mat$group <- pData(Vedo.combined.eset.final)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(Vedo.combined.eset.final)$X),'response']
V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))
V1.module.disrupted.nodes.exp.mat <- V1.module.disrupted.nodes.exp.mat[!is.na(V1.module.disrupted.nodes.exp.mat$group),]
V1.module.disrupted.nodes.exp.mat$sampleID <- row.names(V1.module.disrupted.nodes.exp.mat)
V1.module.disrupted.nodes.exp.mat$Run <- pData(Vedo.combined.eset.final)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(Vedo.combined.eset.final)$X),'Run']
V1.module.disrupted.nodes.exp.mat.scaled <- V1.module.disrupted.nodes.exp.mat
V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID", "Run")] <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group",  "sampleID", "Run")], 2, scale)
V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")],1,mean)

ggplot(V1.module.disrupted.nodes.exp.mat.scaled, aes(x = group, y = as.numeric(scaled.module.score), fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  facet_wrap(~Run)+
  ylab("Relative score")+xlab("Vedolizumab response groups")+
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox",size=4, comparisons=list(c(1,2)), method.args = list(c("greater")))
#  stat_compare_means(method = "wilcox.test", size=3, comparisons=list(c(1,2)), method.args = list(alternative = "greater") )

plot(roc(V1.module.disrupted.nodes.exp.mat.scaled$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score)),  smooth=T))


#Comparison to existing markers

esetMat <- filtered.esetALL[fData(filtered.esetALL)$SYMBOL %in% Baseline.DEGs & fData(filtered.esetALL)$dataType %in% "AGX", pData(filtered.esetALL)$Visit %in% "V1"]
V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
V1.module.disrupted.nodes.exp.mat$group <- pData(filtered.esetALL)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(filtered.esetALL)$sampleID),'group']
V1.module.disrupted.nodes.exp.mat$group <- ifelse(V1.module.disrupted.nodes.exp.mat$group %in% "IM", "NR", V1.module.disrupted.nodes.exp.mat$group)
V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))
V1.module.disrupted.nodes.exp.mat$sampleID <- row.names(V1.module.disrupted.nodes.exp.mat)

V1.module.disrupted.nodes.exp.mat.scaled <- V1.module.disrupted.nodes.exp.mat
V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID", "Run")] <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group",  "sampleID", "Run")], 2, scale)
V1.module.disrupted.nodes.exp.mat.scaled$scaled.module.score <- apply(V1.module.disrupted.nodes.exp.mat.scaled[,!colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("group", "sampleID")],1,mean)


esetMat.current.markers <- filtered.esetALL[fData(filtered.esetALL)$SYMBOL %in% c("TREM1", "OSM"), pData(filtered.esetALL)$Visit %in% "V1"]
esetMat.current.markers.exp.mat <- data.frame(t(exprs(esetMat.current.markers)), stringsAsFactors = F)
colnames(esetMat.current.markers.exp.mat) <- paste(fData(filtered.esetALL)[match(colnames(esetMat.current.markers.exp.mat), fData(filtered.esetALL)$featureID), 'SYMBOL'], fData(filtered.esetALL)[match(colnames(esetMat.current.markers.exp.mat), fData(filtered.esetALL)$featureID), 'dataType'], sep=".")

esetMat.current.markers.exp.mat.merged <- merge(esetMat.current.markers.exp.mat, V1.module.disrupted.nodes.exp.mat.scaled[,colnames(V1.module.disrupted.nodes.exp.mat.scaled) %in% c("sampleID", "scaled.module.score")], by=0)
esetMat.current.markers.exp.mat.merged$Row.names <- NULL
esetMat.current.markers.exp.mat.merged.long <- reshape2::melt(esetMat.current.markers.exp.mat.merged, id.vars=c("scaled.module.score", "sampleID"))
esetMat.current.markers.exp.mat.merged.long$group <- pData(filtered.esetALL)[match(esetMat.current.markers.exp.mat.merged.long$sampleID, pData(filtered.esetALL)$sampleID), 'group']
esetMat.current.markers.exp.mat.merged.long$group <- ifelse(esetMat.current.markers.exp.mat.merged.long$group %in% "IM", "NR", esetMat.current.markers.exp.mat.merged.long$group)
esetMat.current.markers.exp.mat.merged.long <- esetMat.current.markers.exp.mat.merged.long[!esetMat.current.markers.exp.mat.merged.long$sampleID %in% "HR-37-V1",]

ggplot(esetMat.current.markers.exp.mat.merged.long,aes(x=scaled.module.score, y=value, group=1, fill=group)) + 
  geom_point(shape=21, size=4)+
  facet_wrap(~variable, scale='free')+
  #geom_text(aes(label=label, size=2.5),hjust=0, vjust=0)+
  #geom_path(aes(group=patient), arrow=arrow(length = unit(0.3, "cm"), type="closed"))  +
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border =element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_smooth(method = "lm", aes(group="1"), fullrange = TRUE, color="black", se=F)+
  stat_cor(method = "spearman", aes(group="1"))

esetMat.current.markers.exp.mat.merged.long$group <- factor(esetMat.current.markers.exp.mat.merged.long$group, levels=c("R", "NR"))

ggplot(esetMat.current.markers.exp.mat.merged.long, aes(x = group, y = as.numeric(value), fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  facet_wrap(~variable, scale="free")+
  ylab("Relative score")+xlab("Vedolizumab response groups")+
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox",size=4, comparisons=list(c(1,2)), method.args = list(c("greater")))

#Test TREM1 biomarker in RA

TREM1.RA <- lapply(c(esetALL.GEO20690, esetALL.GSE33377, esetALL.GSE42296), function(cur.eset) {
  
  esetMat <- cur.eset[fData(cur.eset)$'Gene.symbol' %in% "TREM1" & fData(cur.eset)$dataType %in% "GX", ]
  V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
  V1.module.disrupted.nodes.exp.mat$group <- pData(esetMat)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(esetMat)$geo_accession),'group']
  V1.module.disrupted.nodes.exp.mat$group <- ifelse(V1.module.disrupted.nodes.exp.mat$group %in% "IM", "NR", V1.module.disrupted.nodes.exp.mat$group)
  V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))
  V1.module.disrupted.nodes.exp.mat$sampleID <- row.names(V1.module.disrupted.nodes.exp.mat)

  g <- ggplot(V1.module.disrupted.nodes.exp.mat, aes(x = group, y = as.numeric(V1.module.disrupted.nodes.exp.mat[,1]), fill=group)) +
    geom_boxplot( outlier.shape = NA) +
    scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
    theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),                                                                  
                     strip.background = element_blank(),
                     panel.border = element_blank(), 
                     axis.ticks.length=unit(.15, "cm"))+
    stat_compare_means(method = "wilcox.test", comparisons=list(c(1,2)))
  
  cur.ROC <- plot(roc(V1.module.disrupted.nodes.exp.mat$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat[,1])),  smooth=T))
  cur.auc <- roc(V1.module.disrupted.nodes.exp.mat$group, as.numeric(as.character(V1.module.disrupted.nodes.exp.mat[,1])))$auc
  
  return(list(MAT=V1.module.disrupted.nodes.exp.mat, boxplot=g, ROC=cur.ROC, auc=cur.auc))
})

names(TREM1.RA) <- c("GEO20690", "GSE33377", "GSE42296")

TREM1.RA.data.for.boxplots <- do.call('rbind', lapply(names(TREM1.RA), function(cur.eset) {
  submat <- TREM1.RA[[cur.eset]]$MAT
  submat$eset <- cur.eset
  colnames(submat) <- c("TREM1", "group", "sampleID", "eset")
  return(submat)
}))


ggplot(TREM1.RA.data.for.boxplots, aes(x = group, y = TREM1, fill=group)) +
  geom_boxplot( outlier.shape = NA) +
  facet_wrap(~eset, scale="free")+
  scale_fill_manual(values=c("R"="cyan3", "NR"="brown2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox.test", comparisons=list(c(1,2)))


#Show TREM1 dynamics 

esetMat.TREM.dynamics <- filtered.esetALL[fData(filtered.esetALL)$SYMBOL %in% c("TREM1") & fData(filtered.esetALL)$dataType %in% "GX", ]
esetMat.TREM.dynamics.exp.mat <- data.frame(t(exprs(esetMat.TREM.dynamics)), stringsAsFactors = F)
colnames(esetMat.TREM.dynamics.exp.mat) <- paste(fData(filtered.esetALL)[match(colnames(esetMat.TREM.dynamics.exp.mat), fData(filtered.esetALL)$featureID), 'SYMBOL'], fData(filtered.esetALL)[match(colnames(esetMat.TREM.dynamics.exp.mat), fData(filtered.esetALL)$featureID), 'dataType'], sep=".")

esetMat.TREM.dynamics.exp.mat$group <- pData(filtered.esetALL)[match(row.names(esetMat.TREM.dynamics.exp.mat), pData(filtered.esetALL)$sampleID), 'group']
esetMat.TREM.dynamics.exp.mat$group <- ifelse(esetMat.TREM.dynamics.exp.mat$group %in% "IM", "NR", esetMat.TREM.dynamics.exp.mat$group)
esetMat.TREM.dynamics.exp.mat <- esetMat.TREM.dynamics.exp.mat[!row.names(esetMat.TREM.dynamics.exp.mat) %in% c("HR-37-V1", "HR-37-V2", "HR-37-V3"),]
esetMat.TREM.dynamics.exp.mat$Visit <- pData(filtered.esetALL)[match(row.names(esetMat.TREM.dynamics.exp.mat), pData(filtered.esetALL)$sampleID), 'Visit']


ggplot(esetMat.TREM.dynamics.exp.mat, aes(x = Visit, y = TREM1.GX, colour = group, group=group, width=0.2)) + 
  stat_summary(fun = "mean", geom = "line", size =1.5)+ # adding connecting lines
  stat_summary(fun = "mean", geom = "point") +        # adding data points
  stat_summary(fun.data = "mean_se", geom = "errorbar", size=1)+
  theme_bw()+theme(strip.text = element_text(size=12), 
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=12),
                   axis.text.x=element_text(size=12,hjust=1),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("TREM1.GX")+
  stat_compare_means(data=esetMat.TREM.dynamics.exp.mat[esetMat.TREM.dynamics.exp.mat$group %in% "R",], method = "wilcox.test", method.args = list(alternative = "greater"), comparisons=list(c(1,2), c(1,3)), label.y=c(11.4, 11.5), color="cyan3")+
  stat_compare_means(data=esetMat.TREM.dynamics.exp.mat[esetMat.TREM.dynamics.exp.mat$group %in% "NR",], method = "wilcox.test", method.args = list(alternative = "less"), comparisons=list(c(1,3), c(1,2)), label.y=c(11.1, 11.3), color="Indianred")

