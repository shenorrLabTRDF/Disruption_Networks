
#---
#title: "Data pre-processing"
#author: "Shiran Vainberg"
#Proj : "A personalized network framework reveals predictive axis of anti-TNF response across diseases"
#---

# Loading packages and functions

library(oligoClasses)
library(oligo)
library(pd.clariom.s.human)
library(clariomshumantranscriptcluster.db)
library(affycoretools)
library(genefilter)
library(sva)
library(plyr)
library(CellMix)
library(hgu133plus2.db)
library(citrus)
library(ComplexHeatmap)
library(data.table)
library(XLConnect)
library(reshape2)
library(foreach)
library(doParallel)
library(permute)
library(plyr)
library(data.table)
library(gridExtra)

code_dir = "./"
source(paste0(code_dir,'Disruption_Networks_functions.R'))

options(stringsAsFactors = F)


#####Gene expression pre-processing #####
#Read data

#define the data directory path:
dataDir = '...'
setwd(dataDir)

# Read in the CEL files in the directory
celFiles <- list.celfiles(dataDir)
affyRaw <- read.celfiles(celFiles)

#pre-processing using rma
Remicade_eset <- rma(affyRaw)

#annotation
#add prob annotation including SYMBOL and description

Remicade_eset <- annotateEset(Remicade_eset, clariomshumantranscriptcluster.db)
fData(Remicade_eset)$dataType <- rep("GX", nrow(fData(Remicade_eset)))
sampleID = colnames(Remicade_eset)
sampleID = sub("R1282_([0-9]+)_", "", sampleID)
sampleID = sub("[\\_]+[\\(]+Clariom_S_Human+[\\)]+[\\.]+CEL", "", sampleID)

Patient.code <- sub(".*?(\\d+).*", "\\1", sampleID)
Visit <- sub(".*?(V.*)","\\1", sampleID)

fAnn <- fData(Remicade_eset)
colnames(fAnn) <- c("featureID", "ENTREZID", "SYMBOL", "Description", "dataType")

pAnn <- data.frame(sampleID=sampleID, Patient.code=Patient.code, Visit=Visit, row.names=colnames(Remicade_eset))


Remicade_eset <- ExpressionSet(exprs(Remicade_eset), phenoData = AnnotatedDataFrame(pAnn)
                               , featureData = AnnotatedDataFrame(fAnn)
                               , annotation = 'ClariomS.db')

#Add batch classification to phenotypic data

batch.classification <- read.csv(paste(dataDir, "/GX/CEL.files.list.batch.csv", sep=""),stringsAsFactors = F)

pData(Remicade_eset)$batch <- batch.classification[match(row.names(pData(Remicade_eset)), batch.classification$sample),'batch']

Remicade_eset.GX <- Remicade_eset
rm(Remicade_eset)
gc()


#filter control probes

Remicade_eset.GX <- Remicade_eset.GX[!is.na(fData(Remicade_eset.GX)$ENTREZID),]

#Add patient classification according to response status
R <- c("20", "22", "24", "28", "29", "31", "33", "35", "39", "40", "42", "44", "46", "47", "48")
NR <-c("23", "26", "30", "32", "36", "37")
IM <- c("21", "27", "38")

pData(Remicade_eset.GX)$Patient.code <- as.character(pData(Remicade_eset.GX)$Patient.code)
pData(Remicade_eset.GX)$group <- "NA"
pData(Remicade_eset.GX)$group <- ifelse(pData(Remicade_eset.GX)$Patient.code %in% R, "R", ifelse(pData(Remicade_eset.GX)$Patient.code %in% NR, "NR", ifelse(pData(Remicade_eset.GX)$Patient.code %in% IM, "IM", "")))

#Add clincal data
sample_annot <- read.csv(paste(dataDir, "/GX/sample_annotations.csv", sep=""),stringsAsFactors = F)
pData(Remicade_eset.GX) <- merge(pData(Remicade_eset.GX), sample_annot, by="sampleID", all.x=T)

exprs <- exprs(Remicade_eset.GX)
row.names(exprs) <- paste("GX.", row.names(exprs), sep="")

#creat an expression-set
Remicade_eset.GX <- ExpressionSet(exprs, phenoData = AnnotatedDataFrame(pData(Remicade_eset.GX))
                                  , featureData = AnnotatedDataFrame(fData(Remicade_eset.GX))
                                  , annotation = 'ClariomS.db')



# correct for batch effect using combat

modcombat = model.matrix(~1, data=pData(Remicade_eset.GX))
combat_edata = ComBat(dat=exprs(Remicade_eset.GX), batch=pData(Remicade_eset.GX)$batch, mod=modcombat, par.prior=TRUE)

#test heatmap post batch removal
# filtered according to highly variable genes 
var.gene <- apply(combat_edata, 1, var)
quant.var <- quantile(var.gene, probs = c(0.7,0.8, 0.9))
filtered.combat_edata <- combat_edata[var.gene>quant.var[2],]
test <- exprs(Remicade_eset.GX)
filtered.eset.data <- test[var.gene>quant.var[2],]

filtered.combat_edata_scaled = apply(filtered.combat_edata, 1, scale)
row.names(filtered.combat_edata_scaled) <- pData(Remicade_eset.GX)$sampleID

Heatmap(filtered.combat_edata_scaled, name="log2", row_names_gp = gpar(fontsize = 5), show_row_names = TRUE, show_column_names = FALSE)+
  Heatmap(pData(Remicade_eset.GX)$group, name = "Response",col=structure(c("green", "red", "grey"), names=c("R","NR", "IM")),width = unit(5, "mm"))+
  Heatmap(pData(Remicade_eset.GX)$batch, name = "batch",col=structure(c("red", "orange", "yellow", "darkred"), names=c("1","2", "3", "4")),width = unit(5, "mm"))+
  Heatmap(pData(Remicade_eset.GX)$Visit, name = "Visit",col=structure(c("cyan", "blue", "purple", "pink"), names=c("V1","V2", "V3", "V4")),width = unit(5, "mm"))

#prepare updated expression set with corrected values, and without outliers

Remicade_eset.GX.post.batch.correction <- Remicade_eset.GX
exprs(Remicade_eset.GX.post.batch.correction) <- combat_edata

row.names(pData(Remicade_eset.GX.post.batch.correction)) <- make.names(as.character(pData(Remicade_eset.GX.post.batch.correction)$sampleID), unique=TRUE)
row.names(pData(Remicade_eset.GX.post.batch.correction)) <- sub("X", "", row.names(pData(Remicade_eset.GX.post.batch.correction)))

colnames(exprs(Remicade_eset.GX.post.batch.correction)) <- row.names(pData(Remicade_eset.GX.post.batch.correction))


exprs <- exprs(Remicade_eset.GX.post.batch.correction)
row.names(exprs) <- paste("GX.", row.names(exprs), sep="")
fData(Remicade_eset.GX.post.batch.correction)$featureID <- paste("GX.", fData(Remicade_eset.GX.post.batch.correction)$featureID, sep="")
row.names(fData(Remicade_eset.GX.post.batch.correction)) <- fData(Remicade_eset.GX.post.batch.correction)$featureID
Remicade_eset.GX.post.batch.correction <- ExpressionSet(exprs, phenoData = AnnotatedDataFrame(pData(Remicade_eset.GX.post.batch.correction))
                                                        , featureData = AnnotatedDataFrame(fData(Remicade_eset.GX.post.batch.correction))
                                                        , annotation = 'ClariomS.db')



######CyTOF processing #####

# load total percentage of granulocytes and PBMCs
tot_pct <- read.csv(paste(dataDir, "/CyTOF/IFX time-course 2016 normed all patients _CR_Exported_Stats_Jul-09-2018_01-59-AM.csv", sep=""),stringsAsFactors = F)

#Run CitrusCodeLymph.R in the CyTOF/CyTOF/FCS CyTOF files/Lymphocytes/0.01 directory
#load data of citrus output 
#load lymphocytes
load(paste(dataDir, "/CyTOF/FCS CyTOF files/Lymphocytes/0.01/citrusOutFinal.RData", sep=""))
lymph.abund.citrus <- citrusRes$General$abundances$foldFeatureSet$allFeatures
lymph.marker.citrus <- citrusRes$General$medians$foldFeatureSet$allFeatures
# cluster annotation
cell.annot_pbmc <- read.csv(paste(dataDir, "/CyTOF/Citrus annotations.csv", sep=""),stringsAsFactors = F)
marker.annot.pbmc <- read.csv(paste(dataDir, "/CyTOF/marker_annotation_lymph.2018.csv", sep=""),stringsAsFactors = F)

feature.data.lymph.marker <- data.frame(full.name=colnames(lymph.marker.citrus), stringsAsFactors = FALSE)
feature.data.lymph.marker$dataSource <- rep("PBMC", nrow(feature.data.lymph.marker))
feature.data.lymph.marker$dataType <- rep("CyTOF-PBMC", nrow(feature.data.lymph.marker))
feature.data.lymph.marker$clusterID <- sub(".*?(\\d+).*", "\\1", feature.data.lymph.marker$full.name)

#marker annotation

feature.data.lymph.marker$channelID0 <- sub(".*?(\\(.*\\)).*","\\1", feature.data.lymph.marker$full.name)
feature.data.lymph.marker$channelID0 <- paste(feature.data.lymph.marker$channelID0, "Di", sep="")
feature.data.lymph.marker$marker <-  marker.annot.pbmc[match(feature.data.lymph.marker$channelID0, marker.annot.pbmc$Marker),'marker']
feature.data.lymph.marker$featureID <- paste("PB", feature.data.lymph.marker$clusterID, feature.data.lymph.marker$marker, sep=".")

citrusRes.clustering.merged <- citrusRes$foldClustering$allClustering$clustering$merge

for( i in seq(nrow(cell.annot_pbmc))) {
  print(i)
  if(is.na(cell.annot_pbmc$Start[i])) {next}
  if(is.na(cell.annot_pbmc$End[i])) {
    res <- citrus::citrus.getClusterDecendants(cell.annot_pbmc$Start[i],citrusRes.clustering.merged)
    res <- as.character(res)
  } 
  else {
    all <- citrus::citrus.getClusterDecendants(cell.annot_pbmc$Start[i],citrusRes.clustering.merged)
    to.reduced <- citrus::citrus.getClusterDecendants(cell.annot_pbmc$End[i],citrusRes.clustering.merged)
    all <- as.character(all)
    to.reduced <- as.character(to.reduced)
    while(is.na(cell.annot_pbmc$Start[i+1])) {
      to.reduced2 <- citrus::citrus.getClusterDecendants(cell.annot_pbmc$End[i+1],citrusRes.clustering.merged)
      to.reduced <- c(to.reduced,to.reduced2)
      i <- i+1
    }
    res <- all[!(all %in% to.reduced)]}
  final.results[[i]] <- res[res %in% citrusRes$General$abundances$foldFeatureSet$allLargeEnoughClusters]
}

clusters.Decendants <- final.results

clusters.Decendants.largeEnough <- lapply(clusters.Decendants, function(l) {data.frame(l, stringsAsFactors = FALSE)})

map.clusters.to.cellType.df <- rbindlist(clusters.Decendants.largeEnough, idcol="index")

cell.annot_pbmc$index <- seq(1:nrow(cell.annot_pbmc))
map.clusters.to.cellType.df$cellType <- cell.annot_pbmc[match(map.clusters.to.cellType.df$index, cell.annot_pbmc$index),'Annotation']
map.clusters.to.cellType.df$Start <- cell.annot_pbmc[match(map.clusters.to.cellType.df$index, cell.annot_pbmc$index),'Start']

map.clusters.to.cellType.df.merged <-aggregate(Start ~ l , data=map.clusters.to.cellType.df, min)
map.clusters.to.cellType.df.merged2 <- map.clusters.to.cellType.df.merged
map.clusters.to.cellType.df.merged2$cellType <- cell.annot_pbmc[match(map.clusters.to.cellType.df.merged$Start, cell.annot_pbmc$Start),'Annotation']
map.clusters.to.cellType.df.merged2 <- map.clusters.to.cellType.df.merged2[!map.clusters.to.cellType.df.merged2$l %in% cell.annot_pbmc$Start,]

map.clusters.to.cellType.df.merged.include.start <- data.frame(l=cell.annot_pbmc$Start, Start=cell.annot_pbmc$Start, cellType=cell.annot_pbmc$Annotation, stringsAsFactors = FALSE)
map.clusters.to.cellType.df.final <- rbind(map.clusters.to.cellType.df.merged2, map.clusters.to.cellType.df.merged.include.start)

feature.data.lymph.marker$cellType <- map.clusters.to.cellType.df.final[match(feature.data.lymph.marker$clusterID,map.clusters.to.cellType.df.final$l),'cellType']
feature.data.lymph.marker$displayName <- paste(feature.data.lymph.marker$cellType, feature.data.lymph.marker$marker, sep="-")
feature.data.lymph.marker$is.median <- NULL
feature.data.lymph.marker$featureID <- paste("PB", feature.data.lymph.marker$clusterID, feature.data.lymph.marker$marker, sep=".")

#For abundance marker

feature.data.lymph.abundance <- data.frame(full.name=colnames(lymph.abund.citrus), stringsAsFactors = FALSE)

feature.data.lymph.abundance$dataSource <- rep("PBMC", nrow(feature.data.lymph.abundance))
feature.data.lymph.abundance$dataType <- rep("CyTOF-PBMC", nrow(feature.data.lymph.abundance))
feature.data.lymph.abundance$clusterID <- sub(".*?(\\d+).*", "\\1", feature.data.lymph.abundance$full.name)
feature.data.lymph.abundance$cellType <- map.clusters.to.cellType.df.final[match(feature.data.lymph.abundance$clusterID,map.clusters.to.cellType.df.final$l),'cellType']
feature.data.lymph.abundance$displayName <- paste(feature.data.lymph.abundance$cellType,"Abundance", sep="-")
feature.data.lymph.abundance$featureID <- paste("PB", feature.data.lymph.abundance$clusterID, sep=".")

#prepare pData for the pb.cytof data
phenotypic.data.lymph <- data.frame(sample.full.name=row.names(lymph.abund.citrus), stringsAsFactors = FALSE)
phenotypic.data.lymph$sampleID <- sub(".*?(HR-\\d+-V\\d+).*", "\\1" , phenotypic.data.lymph$sample.full.name)  
phenotypic.data.lymph$Patient.code <-   sub(".*?(\\d+).*", "\\1", phenotypic.data.lymph$sampleID)
phenotypic.data.lymph$Visit <-   sub(".*?(V\\d+).*", "\\1", phenotypic.data.lymph$sampleID)
row.names(phenotypic.data.lymph) <- phenotypic.data.lymph$sampleID

batch.lymph.data <- read.csv(paste(dataDir, "/CyTOF/phenotypic.data.lymph.batch.csv", sep=""), stringsAsFactors = FALSE)
phenotypic.data.lymph$batch <- batch.lymph.data[match(phenotypic.data.lymph$sampleID, batch.lymph.data$sampleID),'batch']
phenotypic.data.lymph$group <- pData(Remicade_eset.GX.post.batch.correction)[match(phenotypic.data.lymph$sampleID, pData(Remicade_eset.GX.post.batch.correction)$sampleID),'group']
phenotypic.data.lymph$group <- "NA"
phenotypic.data.lymph$group <- ifelse(phenotypic.data.lymph$Patient.code %in% R, "R", ifelse(phenotypic.data.lymph$Patient.code %in% NR, "NR", ifelse(phenotypic.data.lymph$Patient.code %in% IM, "IM", "")))

#combine markers and abundance data
row.names(lymph.abund.citrus) <-phenotypic.data.lymph$sampleID 
colnames(lymph.abund.citrus) <- feature.data.lymph.abundance$featureID
row.names(lymph.marker.citrus) <-phenotypic.data.lymph$sampleID 
colnames(lymph.marker.citrus) <- feature.data.lymph.marker$featureID

test <- cbind(lymph.abund.citrus, lymph.marker.citrus)
test <- data.frame(t(test), stringsAsFactors = FALSE)

feature.pb <- rbind.fill(feature.data.lymph.abundance, feature.data.lymph.marker)
row.names(feature.pb) <- feature.pb$featureID
colnames(test) <- gsub("\\.", "-", colnames(test))

Remicade_eset.PB <- ExpressionSet(as.matrix(test), phenoData = AnnotatedDataFrame(phenotypic.data.lymph)
                                  , featureData = AnnotatedDataFrame(feature.pb)
                                  , annotation = 'CyTOF')


#Granulocytes

#Run CitrusCodeGran2.R in the /CyTOF/FCS CyTOF files/Grans/0.02 directory
load(paste(dataDir, "/CyTOF/FCS CyTOF files/Grans/0.02/citrusOutFinal.RData", sep=""))

Gran.abund.citrus <- citrusRes$General$abundances$foldFeatureSet$allFeatures
Gran.marker.citrus <- citrusRes$General$medians$foldFeatureSet$allFeatures
# cluster annotation

marker.annot.gran <- read.csv(paste(dataDir, "/CyTOF/marker_annotation_grans.2018.csv", sep=""),stringsAsFactors = F)


# cluster annotation
feature.data.grans.marker <- data.frame(full.name=colnames(Gran.marker.citrus), stringsAsFactors = FALSE)
feature.data.grans.marker$dataSource <- rep("GRAN", nrow(feature.data.grans.marker))
feature.data.grans.marker$dataType <- rep("CyTOF-GRAN", nrow(feature.data.grans.marker))
feature.data.grans.marker$clusterID <- sub(".*?(\\d+).*", "\\1", feature.data.grans.marker$full.name)

#marker annotation

feature.data.grans.marker$channelID0 <- sub(".*?(\\(.*\\)).*","\\1", feature.data.grans.marker$full.name)
feature.data.grans.marker$channelID0 <- paste(feature.data.grans.marker$channelID0, "Di", sep="")
feature.data.grans.marker$marker <-  marker.annot.gran[match(feature.data.grans.marker$channelID0, marker.annot.gran$Marker),'Label']
feature.data.grans.marker$featureID <- paste("GR", feature.data.grans.marker$clusterID, feature.data.grans.marker$marker, sep=".")

feature.data.grans.marker$cellType <- rep("Granulocytes", nrow(feature.data.grans.marker))
feature.data.grans.marker$displayName <- paste(feature.data.grans.marker$cellType, feature.data.grans.marker$marker, sep="-")

#For abundance marker

feature.data.grans.abundance <- data.frame(full.name=colnames(Gran.abund.citrus), stringsAsFactors = FALSE)

feature.data.grans.abundance$dataSource <- rep("GRAN", nrow(feature.data.grans.abundance))
feature.data.grans.abundance$dataType <- rep("CyTOF-GRAN", nrow(feature.data.grans.abundance))
feature.data.grans.abundance$clusterID <- sub(".*?(\\d+).*", "\\1", feature.data.grans.abundance$full.name)
feature.data.grans.abundance$cellType <- rep("Granulocytes", nrow(feature.data.grans.abundance))
feature.data.grans.abundance$displayName <- paste(feature.data.grans.abundance$cellType,"Abundance", sep="-")
feature.data.grans.abundance$featureID <- paste("GR", feature.data.grans.abundance$clusterID, sep=".")

#preapare pData for the pb.cytof data
phenotypic.data.grans <- data.frame(sample.full.name=row.names(Gran.abund.citrus), stringsAsFactors = FALSE)
phenotypic.data.grans$sampleID <- sub(".*?(HR-\\d+-V\\d+).*", "\\1" , phenotypic.data.grans$sample.full.name)  
phenotypic.data.grans$Patient.code <-   sub(".*?(\\d+).*", "\\1", phenotypic.data.grans$sampleID)
phenotypic.data.grans$Visit <-   sub(".*?(V\\d+).*", "\\1", phenotypic.data.grans$sampleID)
row.names(phenotypic.data.grans) <- phenotypic.data.grans$sampleID

batch.grans.data <- read.csv(paste(dataDir, "/CyTOF/phenotypic.data.lymph.batch.csv", sep=""),stringsAsFactors = F)

batch.grans.data$batch <- batch.grans.data[match(phenotypic.data.grans$sampleID, batch.grans.data$sampleID),'batch']
batch.grans.data$group <- pData(Remicade_eset.GX.post.batch.correction)[match(phenotypic.data.grans$sampleID, pData(Remicade_eset.GX.post.batch.correction)$sampleID),'group']
batch.grans.data$group <- "NA"
batch.grans.data$group <- ifelse(phenotypic.data.grans$Patient.code %in% R, "R", ifelse(phenotypic.data.grans$Patient.code %in% NR, "NR", ifelse(phenotypic.data.grans$Patient.code %in% IM, "IM", "")))


#combine markers and abundance data
row.names(Gran.abund.citrus) <-phenotypic.data.grans$sampleID 
colnames(Gran.abund.citrus) <- feature.data.grans.abundance$featureID
row.names(Gran.marker.citrus) <-phenotypic.data.grans$sampleID 
colnames(Gran.marker.citrus) <- feature.data.grans.marker$featureID

test <- cbind(Gran.abund.citrus, Gran.marker.citrus)
test <- data.frame(t(test), stringsAsFactors = FALSE)

feature.gr <- rbind.fill(feature.data.grans.abundance, feature.data.grans.marker)
row.names(feature.gr) <- feature.gr$featureID
colnames(test) <- gsub("\\.", "-", colnames(test))

Remicade_eset.GR <- ExpressionSet(as.matrix(test), phenoData = AnnotatedDataFrame(phenotypic.data.grans)
                                  , featureData = AnnotatedDataFrame(feature.gr)
                                  , annotation = 'CyTOF')

#filter clusters
#filter clusters and markers from Granulocyte eset

#exclude descendant clusters of 1079994 and clusters 1079950, 1079981 
citrusRes.clustering.merged <- citrusRes$foldClustering$allClustering$clustering$merge
Decendants.1079994<- citrus::citrus.getClusterDecendants(1079994,citrusRes.clustering.merged)

clusters.Decendants.largeEnough <- Decendants.1079994[Decendants.1079994 %in% citrusRes$General$abundances$foldFeatureSet$allLargeEnoughClusters]
clusters.Decendants.largeEnough <- as.character(clusters.Decendants.largeEnough)
to.exclude.clusters <- c(clusters.Decendants.largeEnough, "1079950", "1079981")
index.to.exclude <- grepl(paste(as.character(to.exclude.clusters), collapse="|"),fData(Remicade_eset.GR)$featureID)
Remicade_eset.GR.final <- Remicade_eset.GR[!index.to.exclude,]

#filter markers

Remicade_eset.GR.final <- Remicade_eset.GR.final[fData(Remicade_eset.GR.final)$marker %in% marker.annot.gran[marker.annot.gran$To_keep %in% "Y",'Label']| is.na(fData(Remicade_eset.GR.final)$marker),]

pData(Remicade_eset.GR.final)$batch <- batch.lymph.data[match(pData(Remicade_eset.GR.final)$sampleID, batch.lymph.data$sampleID),'batch']
pData(Remicade_eset.GR.final)$group <- pData(Remicade_eset.GX.post.batch.correction)[match(pData(Remicade_eset.GR.final)$sampleID, pData(Remicade_eset.GX.post.batch.correction)$sampleID),'group']
saveRDS(Remicade_eset.GR.final, "Remicade_eset.GR.final.rds")


#combine PB and GR expression sets

Remicade_eset.PB.for.merge <- Remicade_eset.PB
pData(Remicade_eset.PB.for.merge)$sample.full.name <- NULL

Remicade_eset.GR.for.merge <- Remicade_eset.GR.final
pData(Remicade_eset.GR.for.merge)$sample.full.name <- NULL


Remicade_eset.CyTOF <- combine(Remicade_eset.PB.for.merge, Remicade_eset.GR.for.merge)

#batch correction for CyTOF data
#before batch correction

exprs <- t(exprs(Remicade_eset.CyTOF))

var.gene <- apply(exprs, 1, var)
quant.var <- quantile(var.gene, probs = c(0.7,0.8, 0.9))
filtered.exprs <- exprs[,var.gene>quant.var[2]]

filtered.exprs.scaled <- apply(filtered.exprs, 2, scale)
row.names(filtered.exprs.scaled) <- row.names(filtered.exprs)

Heatmap(filtered.exprs.scaled, name="abundance/arcSin(marker.int)", row_names_gp = gpar(fontsize = 5), show_row_names = TRUE, show_column_names = FALSE)+
  Heatmap(pData(Remicade_eset.CyTOF)$batch, name = "batch",col=structure(c("red", "pink", "orange", "yellow", "darkred"), names=c("1","2", "3", "4", "5")),width = unit(5, "mm"))+
  Heatmap(pData(Remicade_eset.CyTOF)$group, name = "Response",col=structure(c("green", "red", "grey"), names=c("R","NR", "IM")),width = unit(5, "mm"))+
  Heatmap(pData(Remicade_eset.CyTOF)$Visit, name = "Visit",col=structure(c("cyan", "blue", "purple", "pink"), names=c("V1","V2", "V3", "V4")),width = unit(5, "mm"))

modcombat = model.matrix(~1, data=pData(Remicade_eset.CyTOF))
combat_edata = ComBat(dat=exprs(Remicade_eset.CyTOF), batch=pData(Remicade_eset.CyTOF)$batch, mod=modcombat, par.prior=TRUE)
combat_edata <- t(combat_edata)
filtered.combat_edata <- combat_edata[,var.gene>quant.var[2]]

combat_edata_scaled = apply(filtered.combat_edata, 2, scale)
row.names(combat_edata_scaled) <- pData(Remicade_eset.CyTOF)$sampleID


#after batch correction
Heatmap(combat_edata_scaled, name="log2", row_names_gp = gpar(fontsize = 5), show_row_names = TRUE, show_column_names = FALSE)+
  Heatmap(pData(Remicade_eset.CyTOF)$batch, name = "batch",col=structure(c("red", "pink", "orange", "yellow", "darkred"), names=c("1","2", "3", "4", "5")),width = unit(5, "mm"))+
  Heatmap(pData(Remicade_eset.CyTOF)$group, name = "Response",col=structure(c("green", "red", "grey"), names=c("R","NR", "IM")),width = unit(5, "mm"))+
  Heatmap(pData(Remicade_eset.CyTOF)$Visit, name = "Visit",col=structure(c("cyan", "blue", "purple", "pink"), names=c("V1","V2", "V3", "V4")),width = unit(5, "mm"))

#prepare updated expression set with corrected values, and without outliers

Remicade_eset.CyTOF.post.batch.correction <- Remicade_eset.CyTOF
combat_edata_t <- t(combat_edata)
exprs(Remicade_eset.CyTOF.post.batch.correction) <- combat_edata_t
row.names(pData(Remicade_eset.CyTOF.post.batch.correction)) <- make.names(as.character(pData(Remicade_eset.CyTOF.post.batch.correction)$sampleID), unique=TRUE)
row.names(pData(Remicade_eset.CyTOF.post.batch.correction)) <- sub("X", "", row.names(pData(Remicade_eset.CyTOF.post.batch.correction)))
row.names(pData(Remicade_eset.CyTOF.post.batch.correction)) <- gsub("/.", "-", row.names(pData(Remicade_eset.CyTOF.post.batch.correction)))

#Add total granulocytes and lymphocytes to expressionset
exprs <- exprs(Remicade_eset.CyTOF.post.batch.correction)
exprs_t <- t(exprs)
PB.TOTAL <- tot_pct[match(row.names(exprs_t), tot_pct$sample.name),'Lymphocytes']
GR.TOTAL <- tot_pct[match(row.names(exprs_t), tot_pct$sample.name),'Granulocytes']
exprs.final <- rbind(exprs, PB.TOTAL, GR.TOTAL)

features <- fData(Remicade_eset.CyTOF.post.batch.correction)
feature.PB.GR <- data.frame(dataSource=c("PBMC", "GRAN"), featureID=c("PB.TOTAL", "GR.TOTAL"), clusterID=c("TOTAL","TOTAL"), cellType=c("CC-CyTOF-PBMC", "CC-CyTOF-GRAN"), displayName=c("CC-CyTOF-PBMC - Abundance", "CC-CyTOF-GRAN - Abundance"), stringsAsFactors = FALSE)
features.final <- rbind.fill(features, feature.PB.GR)
row.names(features.final) <- c(row.names(features), row.names(feature.PB.GR))
row.names(pData(Remicade_eset.CyTOF.post.batch.correction)) <- pData(Remicade_eset.CyTOF.post.batch.correction)$sampleID

Remicade_eset.CyTOF.post.batch.correction.final <- ExpressionSet(as.matrix(exprs.final), phenoData = AnnotatedDataFrame(pData(Remicade_eset.CyTOF.post.batch.correction))
                                                                 , featureData = AnnotatedDataFrame(features.final)
                                                                 , annotation = 'CyTOF')




#####process Luminex pre-processing######

dataDir_LU <- paste(dataDir, "/Luminex files", sep="")
setwd(dataDir_LU)

#Data processing

MIN_BEAD_NUM = 50

#function for Luminex data processing
Lu.preprocess.function <- function(file.name.data, file.name.out, plate.design) {
  
  #Read in MFI luminex data
  wb.data <- loadWorkbook(file.name.data)
  wb.out <- loadWorkbook(file.name.out)
  
  lx=list()
  
  lx$medians = readWorksheet(wb.data,'Median',rownames=1,check.names=TRUE,header=TRUE)
  lx$counts = readWorksheet(wb.data,'Count',rownames=1,check.names=TRUE,header=TRUE)
  lx$CV = readWorksheet(wb.data,'CV',rownames=1,check.names=TRUE,header=TRUE)
  relDataCols=c(2:46)
  #Read in plate-patient mapping
  wb.design <- loadWorkbook(plate.design)
  pdesign1 = readWorksheet(wb.design,'Sheet1',header=TRUE,endRow = 9,endCol = 13,rownames=1)
  rows = rownames(pdesign1)
  cols = colnames(pdesign1)
  sampleMap = data.frame(platePos = rownames(lx$medians),sampleInWell = NA,TotalEvents = 0,highCVflag=0,lowCountFlag=0,row.names=1)
  for (rw in 1:length(rows)) {
    for (cl in 1:length(cols)) {
      sampleMap[paste0(rows[rw],sub('X','',cols[cl])),'sampleInWell'] = pdesign1[rw,cl]
    }
  }
  sampleMap[rownames(lx$medians),'TotalEvents'] = lx$medians[,'Total.Events']
  
  
  #Initate dump matrices
  cytData.sample = list();
  mes = c('medians','counts','CV')
  for (m in mes) {
    cytData.sample[[m]] = matrix(NA,length(sampleMap[,'sampleInWell']),length(relDataCols),dimnames=list(sort(sampleMap[,'sampleInWell']),colnames(lx[[m]])[relDataCols]))
  }
  
  #Reogrnaize sample format for easy manual inspection
  sampleTypes = sort(unique(sampleMap[,'sampleInWell']))
  for (s in sampleTypes) {
    relWells = rownames(sampleMap)[sampleMap[,'sampleInWell'] %in% s]
    for(m in mes) {
      cytData.sample[[m]][which(rownames(cytData.sample[[m]]) == s),] = data.matrix(lx[[m]][relWells,relDataCols])
      rownames(cytData.sample[[m]])[which(rownames(cytData.sample[[m]]) == s)] = paste(rownames(cytData.sample[[m]])[which(rownames(cytData.sample[[m]]) == s)],relWells,sep='-')
    }
  }
  for(m in mes) {
    createSheet(wb.out,name=m)
    writeWorksheet(wb.out,cytData.sample[[m]][sort(rownames(cytData.sample[[m]])),],sheet=m,rownames='SampleInfo')
  }              
  
  #Cerate Avg Medians (pre-filtering)
  createSheet(wb.out,name='AvgMedians')
  createSheet(wb.out,name='AvgMedianCV')
  cytData.mean = list()
  cytData.mean$medians = matrix(NA,length(unique(sampleMap[,'sampleInWell'])),length(relDataCols),dimnames=list(unique(sampleMap[,'sampleInWell']),colnames(lx$medians)[relDataCols]))
  cytData.mean$CV = matrix(NA,length(unique(sampleMap[,'sampleInWell'])),length(relDataCols),dimnames=list(unique(sampleMap[,'sampleInWell']),colnames(lx$medians)[relDataCols]))
  
  for (s in sampleTypes) {
    relWells = rownames(sampleMap)[sampleMap[,'sampleInWell'] %in% s]
    ids = paste(s,relWells,sep='-')
    if(length(relWells) > 1) {
      cytData.mean$medians[s,] = apply(cytData.sample$medians[ids,],2,mean,na.rm=TRUE)
      cytData.mean$CV[s,] = 100*(apply(cytData.sample$medians[ids,],2,sd,na.rm=TRUE)/cytData.mean$medians[s,])
    } else { #There is only one replicate of the sample
      cytData.mean$medians[s,] = cytData.sample$medians[ids,]
      cytData.mean$CV[s,] = NA
    }
  }
  writeWorksheet(wb.out,cytData.mean$medians[sort(rownames(cytData.mean$medians)),],sheet='AvgMedians',rownames='SampleInfo')
  writeWorksheet(wb.out,cytData.mean$CV[sort(rownames(cytData.mean$medians)),],sheet='AvgMedianCV',rownames='SampleInfo')
  
  #Pre-filtering of samples based on CV or counts and make matrix of means
  createSheet(wb.out,name='filterSampleMed')
  cytData.sample.filter = cytData.sample
  for (s in sampleTypes) {
    relWells = rownames(sampleMap)[sampleMap[,'sampleInWell'] %in% s]
    ids = paste(s,relWells,sep='-')
    #Global number of beads measured in the well must be high enough (suggests clog otherwise)
    totalEventFlag  = sampleMap[relWells,'TotalEvents'] < MIN_BEAD_NUM*dim(cytData.sample$medians)[2]
    sampleMap[relWells,'lowCountFlag'] = totalEventFlag
    if(sum(totalEventFlag) > 0) {
      cytData.sample.filter$medians[ids[which(totalEventFlag == TRUE)],] = NA
    }
    #Cytokine specific filter
    for (i in ids) {
      for(cyt in colnames(cytData.sample$medians)) {
        print(i)
        print(cyt)
        #DO NOT omit if one replicate has acceptable bead count and CV's for both replicates are less than 25%.
        if(!is.na(cytData.mean$CV[s,cyt])) {
          if(cytData.mean$CV[s,cyt] <= 25 & sum(cytData.sample$count[ids,cyt] >= 20) >= 1) {
            next;
          }} else {
            if(cytData.sample$count[i,cyt] < 20) {
              # Omit low bead count in wells and/or targets-counts< 30 per target either entire well or only a target cytokine omit or flag.
              cytData.sample.filter$medians[i,cyt] = NA
            } else { #You have a high bead count but CV is out the roof
              if(cytData.mean$CV[s,cyt] > 40 | is.na(cytData.mean$CV[s,cyt])) {
                cytData.sample.filter$medians[i,cyt] = NA
              }
            }
          }
      }
    }
  }
  writeWorksheet(wb.out,cytData.sample.filter$medians[sort(rownames(cytData.sample.filter$medians)),],sheet='filterSampleMed',rownames='SampleInfo')
  
  
  createSheet(wb.out,name='filterAvgMed')
  cytData.mean$filteredMedians = matrix(NA,length(unique(sampleMap[,'sampleInWell'])),length(relDataCols),dimnames=list(unique(sampleMap[,'sampleInWell']),colnames(lx$medians)[relDataCols]))
  
  for (s in sampleTypes) {
    relWells = rownames(sampleMap)[sampleMap[,'sampleInWell'] %in% s]
    ids = paste(s,relWells,sep='-')
    cytData.mean$filteredMedians[s,] = apply(cytData.sample.filter$medians[ids,],2,mean,na.rm=TRUE)
  }
  writeWorksheet(wb.out,cytData.mean$filteredMedians[sort(rownames(cytData.mean$filteredMedians)),],sheet='filterAvgMed',rownames='SampleInfo')
  
  saveWorkbook(wb.out)
  
}

#for plate 1  
file.name.data <- "NOV9TH2016PLATE1ELINA_20161109_061154.xlsx"
file.name.out <- "NOV9TH2016PLATE1ELINA_20161109_061154_preprocessed.xlsx"
plate.design <- "Luminex plate1.xlsx"

Lu.preprocess.function(file.name.data, file.name.out, plate.design)

#for plate 2

file.name.data <- "NOV9TH2016PLATE2ELINA_20161109_075729.xlsx"
file.name.out <- "NOV9TH2016PLATE2ELINA_20161109_075729_preprocessed.xlsx"
plate.design <- "Luminex plate2.xlsx"

Lu.preprocess.function(file.name.data, file.name.out, plate.design)


#calculate net MFI by reducing blank values

LU.plate1 <- loadWorkbook("NOV9TH2016PLATE1ELINA_20161109_061154_preprocessed.xlsx")
LU.plate1.filterAVGMed <- readWorksheet(LU.plate1,'filterAvgMed',rownames=1,check.names=TRUE,header=TRUE)

#reduce first row values (blank) from each matrix row
LU.plate1.filterAVGMed.Net <- do.call('rbind', lapply(seq(nrow(LU.plate1.filterAVGMed)), function(i) {
  data.frame(LU.plate1.filterAVGMed[i,], stringsAsFactors = FALSE)-data.frame(LU.plate1.filterAVGMed[1,],stringsAsFactors = FALSE)
}))

LU.plate2 <- loadWorkbook("NOV9TH2016PLATE2ELINA_20161109_075729_preprocessed.xlsx")
LU.plate2.filterAVGMed <- readWorksheet(LU.plate2,'filterAvgMed',rownames=1,check.names=TRUE,header=TRUE)

#reduce first row values (blank) from each matrix row
LU.plate2.filterAVGMed.Net <- do.call('rbind', lapply(seq(nrow(LU.plate2.filterAVGMed)), function(i) {
  data.frame(LU.plate2.filterAVGMed[i,], stringsAsFactors = FALSE)-data.frame(LU.plate2.filterAVGMed[1,],stringsAsFactors = FALSE)
}))

#before combining Lu plates data, prepare feature and phenotipic data for expression set and batch correction

#phenotypic data
wb.design1 <- loadWorkbook("Luminex plate1.xlsx")
pdesign1 = readWorksheet(wb.design1,'Sheet1',header=TRUE,endRow = 9,endCol = 13,rownames=1)
plate1.samples <- as.vector(as.matrix(pdesign1))

wb.design2 <- loadWorkbook("Luminex plate2.xlsx")
pdesign2 = readWorksheet(wb.design2,'Sheet1',header=TRUE,endRow = 9,endCol = 13,rownames=1)
plate2.samples <- as.vector(as.matrix(pdesign2))

plate1.samples.full <- sub(".*?(HR)","HR-", plate1.samples)
plate1.samples.full <- sub("(HR-\\d+)(V\\d+.*)","\\1-\\2", plate1.samples.full)
plate1.samples.full <- gsub(" ", "-", plate1.samples.full)


plate2.samples.full <- sub(".*?(HR)","HR-", plate2.samples)
plate2.samples.full <- sub("(HR-\\d+)(V\\d+.*)","\\1-\\2", plate2.samples.full)

plate2.samples.full[plate2.samples.full %in% c("Teva 53", "Teva 55", "Blank")] <- c("Teva-53-1", "Teva-55-1", "Blank-1")
plate2.samples.full <- gsub(" ", "-", plate2.samples.full)

phenotypic.data.LU <- data.frame(sample.full.name=c(row.names(LU.plate1.filterAVGMed.Net), row.names(LU.plate2.filterAVGMed.Net)), stringsAsFactors = FALSE)
phenotypic.data.LU$sampleID <-sub(".*?(HR)","HR-", phenotypic.data.LU$sample.full.name)
phenotypic.data.LU$sampleID <- sub("(HR-\\d+)(V\\d+.*)","\\1-\\2", phenotypic.data.LU$sampleID)

phenotypic.data.LU$Patient.code <-   sub(".*?(\\d+).*", "\\1", phenotypic.data.LU$sampleID)
phenotypic.data.LU$Visit <-   sub(".*?(V\\d+).*", "\\1", phenotypic.data.LU$sampleID)
row.names(phenotypic.data.LU) <- make.names(phenotypic.data.LU$sampleID, unique = TRUE, allow_ = TRUE)
row.names(phenotypic.data.LU) <- gsub("\\.", "-", row.names(phenotypic.data.LU))
phenotypic.data.LU$batch <- ifelse(row.names(phenotypic.data.LU) %in% plate1.samples.full, "1", ifelse(row.names(phenotypic.data.LU) %in% plate2.samples.full, "2", "NA"))

#feature data data
feature.data.LU <- read.csv("Luminex analytes.csv", stringsAsFactors = FALSE)
feature.data.LU$featureID <- gsub(" ", "\\.", feature.data.LU$featureID)
row.names(feature.data.LU) <- feature.data.LU$featureID

#exprs data
LU.filterAVGMed.Net.combined <- rbind(LU.plate1.filterAVGMed.Net, LU.plate2.filterAVGMed.Net)
row.names(LU.filterAVGMed.Net.combined) <- gsub(" ", "\\.", row.names(LU.filterAVGMed.Net.combined))
row.names(LU.filterAVGMed.Net.combined)[row.names(LU.filterAVGMed.Net.combined) %in% c("Blank1", "Teva.531", "Teva.551")] <- c("Blank.1", "Teva.53.1", "Teva.55.1")
row.names(LU.filterAVGMed.Net.combined) <- sub(".*?(HR)","HR-", row.names(LU.filterAVGMed.Net.combined))
row.names(LU.filterAVGMed.Net.combined) <- sub("(HR-\\d+)(V\\d+.*)","\\1-\\2", row.names(LU.filterAVGMed.Net.combined))
row.names(LU.filterAVGMed.Net.combined) <- gsub("\\.", "-", row.names(LU.filterAVGMed.Net.combined))
colnames(LU.filterAVGMed.Net.combined) <- paste("LU.", colnames(LU.filterAVGMed.Net.combined), sep="")

LU.filterAVGMed.Net.combined <- t(LU.filterAVGMed.Net.combined)

#generate Luminex expression set
LU.filterAVGMed.Net.combined <- LU.filterAVGMed.Net.combined[,!grepl("Standard|Teva|Blank", colnames(LU.filterAVGMed.Net.combined))]
phenotypic.data.LU <- phenotypic.data.LU[!grepl("Standard|Teva|Blank", row.names(phenotypic.data.LU)),]
Remicade_eset.LU <- ExpressionSet(LU.filterAVGMed.Net.combined, phenoData = AnnotatedDataFrame(phenotypic.data.LU)
                                  , featureData = AnnotatedDataFrame(feature.data.LU)
                                  , annotation = 'Luminex')


#batch correction for Lu data
#before batch correction

Remicade_eset.LU.final <- Remicade_eset.LU[,!pData(Remicade_eset.LU)$sampleID %in% c("Standard 1", "Standard 2", "Standard 3", "Standard 4", "Standard 5", "Standard 6", "Blank")]
exprs <- t(exprs(Remicade_eset.LU.final))
exprs.scaled <- apply(exprs, 2, scale)
row.names(exprs.scaled) <- row.names(exprs)

Heatmap(exprs.scaled, name="log2", row_names_gp = gpar(fontsize = 5), show_row_names = TRUE, show_column_names = FALSE)+
  Heatmap(pData(Remicade_eset.LU.final)$batch, name = "batch",col=structure(c("red", "orange", "yellow", "darkred"), names=c("1","2", "3", "4")),width = unit(5, "mm"))

modcombat = model.matrix(~1, data=pData(Remicade_eset.LU.final))
combat_edata = ComBat(dat=exprs(Remicade_eset.LU.final), batch=pData(Remicade_eset.LU.final)$batch, mod=modcombat, par.prior=TRUE)

combat_edata_scaled = apply(combat_edata, 1, scale)
row.names(combat_edata_scaled) <- pData(Remicade_eset.LU.final)$sampleID


#after batch correction
Heatmap(combat_edata_scaled, name="log2", row_names_gp = gpar(fontsize = 5), show_row_names = TRUE, show_column_names = FALSE)+
  Heatmap(pData(Remicade_eset.LU.final)$batch, name = "batch",col=structure(c("red", "orange", "yellow", "darkred"), names=c("1","2", "3", "4")),width = unit(5, "mm"))
#Heatmap(pData(Remicade_eset.LU)$Visit, name = "Visit",col=structure(c("cyan", "blue", "purple", "pink"), names=c("V1","V2", "V3", "V4")),width = unit(5, "mm"))

#prepare updated expression set with corrected values, and without outliers

Remicade_eset.LU.post.batch.correction <- Remicade_eset.LU.final
exprs(Remicade_eset.LU.post.batch.correction) <- combat_edata
row.names(pData(Remicade_eset.GX.post.batch.correction)) <- make.names(as.character(pData(Remicade_eset.GX.post.batch.correction)$sampleID), unique=TRUE)
row.names(pData(Remicade_eset.GX.post.batch.correction)) <- sub("X", "", row.names(pData(Remicade_eset.GX.post.batch.correction)))

#remove controls
Remicade_eset.LU.post.batch.correction <- Remicade_eset.LU.post.batch.correction[,!pData(Remicade_eset.LU.post.batch.correction)$sampleID %in%  c("Teva 53", "Teva 55")]

#remove samples that do not have gx data
Remicade_eset.LU.post.batch.correction <- Remicade_eset.LU.post.batch.correction[,pData(Remicade_eset.LU.post.batch.correction)$Patient.code %in% unique(pData(Remicade_eset.GX.post.batch.correction)$Patient.code)]

pData(Remicade_eset.LU.post.batch.correction)$sample.full.name <- NULL



#####Generate adjusted GX (cell-centered)####

setwd(dataDir)

#Define major blood celltypes for the adjustment based on CyTOF
major <- c(GRAN = 'GR.TOTAL'
           , CD8 = 'PB.1079991', CD4 = 'PB.1079995', CD19 = 'PB.1079966', CD14 = 'PB.1079983', NK = c('PB.1079974', 'PB.1079942'))

P <- exprs(Remicade_eset.CyTOF.post.batch.correction.final)[major, ]
rownames(P) <- names(major)

P[-1L, ] <- sweep(P[-1L, ], 2L, exprs(Remicade_eset.CyTOF.post.batch.correction.final)['PB.TOTAL', ], '*')

P <- rbind(P, NK=P['NK1', ]+P['NK2', ])
P <- P[!row.names(P) %in% c('NK1','NK2'),]

#adjust data by major celltypes
#test for colinearity

P.cor <-cor(t(P),t(P), method="spearman")
Heatmap(P.cor, name="spearman.cor", row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12), show_row_names = TRUE, show_column_names = TRUE)

# correct gene expression within each visit separately

# for V1
P.sub.V1 <- P[,grep("V1", colnames(P))]
Remicade_eset.GX.post.batch.correction.V1 <- exprs(Remicade_eset.GX.post.batch.correction)[,grep("V1", colnames(exprs(Remicade_eset.GX.post.batch.correction)))]

#using cellmix#
AG.V1 <- ged(Remicade_eset.GX.post.batch.correction.V1, P.sub.V1 ,method = 'correction')
AG.V1 <- log2(exprs(AG.V1))

#for V2
P.sub.V2 <- P[,grep("V2", colnames(P))]
Remicade_eset.GX.post.batch.correction.V2 <- exprs(Remicade_eset.GX.post.batch.correction)[,grep("V2", colnames(exprs(Remicade_eset.GX.post.batch.correction)))]

AG.V2 <- ged(Remicade_eset.GX.post.batch.correction.V2, P.sub.V2 ,method = 'correction')
AG.V2 <- log2(exprs(AG.V2))

#for V3
P.sub.V3 <- P[,grep("V3", colnames(P))]
Remicade_eset.GX.post.batch.correction.V3 <- exprs(Remicade_eset.GX.post.batch.correction)[,grep("V3", colnames(exprs(Remicade_eset.GX.post.batch.correction)))]

AG.V3 <- ged(Remicade_eset.GX.post.batch.correction.V3, P.sub.V3 ,method = 'correction')
AG.V3 <- log2(exprs(AG.V3))

#combine all visits
AG <- cbind(AG.V1, AG.V2, AG.V3)

#incorporate to expressionset
gxeset_adjust <- Remicade_eset.GX.post.batch.correction
exprs(gxeset_adjust) <- AG[featureNames(gxeset_adjust), sampleNames(gxeset_adjust)]
fData(gxeset_adjust)$dataType <- 'AGX'
# hack to enable rbinding gx and agx features
featureNames(gxeset_adjust) <- paste0('_', featureNames(gxeset_adjust))
validObject(gxeset_adjust)
i <- which(apply(exprs(gxeset_adjust)<=0, 1L, any))
table(fData(gxeset_adjust)[i, ]$isImmune)

# Hack to remove extra underscore
featureNames(gxeset_adjust) <- gsub('_', '', featureNames(gxeset_adjust))
fData(gxeset_adjust)$featureID <- gsub('_', '', fData(gxeset_adjust)$featureID)
fData(gxeset_adjust)$featureID <- gsub('GX', 'AG', fData(gxeset_adjust)$featureID)
row.names(gxeset_adjust) <- gsub('GX', 'AG', row.names(fData(gxeset_adjust)))

Remicade_eset.AGX.post.batch.correction<- gxeset_adjust


#####combine all expression sets for all dataTypes#####

# combine all data types
#combine.expression.data
exprs.data <- rbind.fill(data.frame(exprs(Remicade_eset.GX.post.batch.correction), stringsAsFactors = FALSE), 
                         data.frame(exprs(Remicade_eset.AGX.post.batch.correction), stringsAsFactors = FALSE),
                         data.frame(exprs(Remicade_eset.CyTOF.post.batch.correction.final), stringsAsFactors = FALSE),
                         data.frame(exprs(Remicade_eset.LU.post.batch.correction), stringsAsFactors = FALSE))


row.names(exprs.data) <- c(row.names(exprs(Remicade_eset.GX.post.batch.correction)),
                           row.names(exprs(Remicade_eset.AGX.post.batch.correction)), 
                           row.names(exprs(Remicade_eset.CyTOF.post.batch.correction.final)), 
                           row.names(exprs(Remicade_eset.LU.post.batch.correction)))

exprs.data <- as.matrix(exprs.data) 
colnames(exprs.data) <- gsub("\\.", "-", colnames(exprs.data)) 

#combine fData
feature.data <- rbind.fill(data.frame(fData(Remicade_eset.GX.post.batch.correction), stringsAsFactors = FALSE), 
                           data.frame(fData(Remicade_eset.AGX.post.batch.correction), stringsAsFactors = FALSE),
                           data.frame(fData(Remicade_eset.CyTOF.post.batch.correction.final), stringsAsFactors = FALSE),
                           data.frame(fData(Remicade_eset.LU.post.batch.correction), stringsAsFactors = FALSE))

row.names(feature.data) <- feature.data$featureID

esetALL <- ExpressionSet(exprs.data, phenoData = AnnotatedDataFrame(pData(Remicade_eset.GX.post.batch.correction))
                         , featureData = AnnotatedDataFrame(feature.data)
                         , annotation = 'Analysis.july.2018')

#####Generate Fold-change expression set of V2 and V3 relative to V1#####

esetALL.FC <- substractReference(esetALL, 'Patient.code', 'Visit', ref = 'V1', op = '-')


#####Differential expression analysis#####

#using mixed effects linear models

#Filter eset
#filter out features for which have NA in more than 20% of the samples
NA.perc = apply(exprs(esetALL), 1, function(x) sum(is.na(x))/length(x)*100)
filtered.esetALL<-esetALL[!NA.perc>20,]

#For the remaining NAs put average value of the feature
exprs(filtered.esetALL)[Features.with.NA,] <- which(apply(exprs(filtered.esetALL), 1, function(x) any(is.na(x))))
test <- do.call('rbind',lapply(Features.with.NA, function(row.index) {
  exprs(filtered.esetALL)[row.index,is.na(exprs(filtered.esetALL)[row.index,])] <- mean(exprs(filtered.esetALL)[row.index,], na.rm=T)
  return(exprs(filtered.esetALL)[row.index,])
}))

featureNames(filtered.esetALL) <- gsub("-", "", featureNames(filtered.esetALL))
featureNames(filtered.esetALL)[grepl(".CD25", featureNames(filtered.esetALL))] <- as.character(do.call('rbind',lapply(featureNames(filtered.esetALL)[grepl(".CD25", featureNames(filtered.esetALL))], function(x) gsub("^\\s+|\\s+$", "", x))))
fData(filtered.esetALL)$featureID <- featureNames(filtered.esetALL)

filtered.esetALL.R <- filtered.esetALL[,pData(filtered.esetALL)$group %in% "R"]
filtered.esetALL.NR <- filtered.esetALL[,pData(filtered.esetALL)$group %in% c("IM", "NR")]

saveRDS(filtered.esetALL, "filtered.esetALL.rds")
saveRDS(filtered.esetALL.R, "filtered.esetALL.R.rds")
saveRDS(filtered.esetALL.NR, "filtered.esetALL.NR.rds")

#Identification of features that changed over time using mixed effect linear models

lmer <- differential.over.time.func("filtered.esetALL.R")
lmer.NR <- differential.over.time.func("filtered.esetALL.NR")

saveRDS(lmer, "lmer.rds")
saveRDS(lmer.NR, "lmer.NR.rds")

#Calculate permutation based p-value
#permute time point within a patient

#For responders
#Use permutation identifier for shuffeling between V1 and V2

#For V1V2
filtered.esetALL.R.V1V2<- filtered.esetALL.R[,pData(filtered.esetALL.R)$Visit %in% c("V1", "V2")]
CTRL <- how(within=Within("free"),
            blocks=factor(pData(filtered.esetALL.R.V1V2)$Patient.code),
            complete=F,maxperm=1000,
            observed = TRUE)

perm.V1V2 <- shuffleSet(nrow(t(exprs(filtered.esetALL.R.V1V2))), nset = 1000, control = CTRL)

#For V1V3
filtered.esetALL.R.V1V3 <- filtered.esetALL.R[,pData(filtered.esetALL.R)$Visit %in% c("V1", "V3")]
Var.vec <- apply(exprs(filtered.esetALL.R.V1V3), 1, var)
filtered.esetALL.R.V1V3.sub <- filtered.esetALL.R.V1V3[!Var.vec==0,]

saveRDS(filtered.esetALL.R.V1V3, version=2, "filtered.esetALL.R.V1V3.rds")
CTRL <- how(within=Within("free"),
            blocks=factor(pData(filtered.esetALL.R.V1V3)$Patient.code),
            complete=F,maxperm=1000,
            observed = TRUE)

perm.V1V3 <- shuffleSet(nrow(t(exprs(filtered.esetALL.R.V1V3))), nset = 1000, control = CTRL)


#For non-responders
#For non-responders V1V2
filtered.esetALL.NR.V1V2 <- filtered.esetALL.NR[,pData(filtered.esetALL.NR)$Visit %in% c("V1", "V2")]
CTRL <- how(within=Within("free"),
            blocks=factor(pData(filtered.esetALL.NR.V1V2)$Patient.code),
            complete=F,maxperm=1000,
            observed = TRUE)

perm.NR.V1V2 <- shuffleSet(nrow(t(exprs(filtered.esetALL.NR.V1V2))), nset = 1000, control = CTRL)


#For non-responders V1V3
filtered.esetALL.NR.V1V3 <- filtered.esetALL.NR[,pData(filtered.esetALL.NR)$Visit %in% c("V1", "V3")]
CTRL <- how(within=Within("free"),
            blocks=factor(pData(filtered.esetALL.NR.V1V3)$Patient.code),
            complete=F,maxperm=1000,
            observed = TRUE)

perm.NR.V1V3 <- shuffleSet(nrow(t(exprs(filtered.esetALL.NR.V1V3))), nset = 1000, control = CTRL)


#generate shuffled expressen sets by visits and calculate permuted coefficients

#run the perm.func to generate permuted coefficients
combined.perm.res.R.V1V2 <- perm.func("perm.V1V2","filtered.esetALL.R.V1V2")
combined.perm.res.R.V1V3 <- perm.func("perm.V1V3","filtered.esetALL.R.V1V3")
combined.perm.res.NR.V1V2 <- perm.func("perm.NR.V1V2", "filtered.esetALL.NR.V1V2")
combined.perm.res.NR.V1V3 <- perm.func("perm.NR.V1V3", "filtered.esetALL.NR.V1V23")

#Calculate p-value and FDR based on the permuted coefficient values

lmer<-readRDS("lmer.rds")
combined.perm.res<-readRDS("combined.perm.res.V1V2.rds")

P.FDR.perm.R.V2 <- permutation.derived.FDR.func(lmer, "combined.perm.res.filtered.esetALL.R.V1V2", "V2")
P.FDR.perm.R.V3 <- permutation.derived.FDR.func(lmer, "combined.perm.res.filtered.esetALL.R.V1V3", "V3")
P.FDR.perm.NR.V2 <- permutation.derived.FDR.func(lmer.NR, "combined.perm.res.filtered.esetALL.NR.V1V2", "V2")
P.FDR.perm.NR.V3 <- permutation.derived.FDR.func(lmer.NR, "combined.perm.res.filtered.esetALL.NR.V1V3", "V3")

Perm.Pvalue.V2 <- P.FDR.perm.R.V2[[1]]
Perm.Pvalue.V3 <- P.FDR.perm.R.V3[[1]]
Perm.Pvalue.V2.NR <- P.FDR.perm.NR.V2[[1]]
Perm.Pvalue.V3.NR <- P.FDR.perm.NR.V3[[1]]

saveRDS("Perm.Pvalue.V2.rds")
saveRDS("Perm.Pvalue.V3.rds")
saveRDS("Perm.Pvalue.V2.NR.rds")
saveRDS("Perm.Pvalue.V3.NR.rds")

#Test changes over time between responders and non-responders
#using an equal sample-size comparison, by subsampling 
#Supp. Fig2

#Generate 200 subsamples of responders
R.patients <- as.character(unique(pData(filtered.esetALL)[pData(filtered.esetALL)$group %in% "R", "Patient.code"]))
downsampling.group.size <- length(unique(pData(filtered.esetALL)[!pData(filtered.esetALL)$group %in% "R", "Patient.code"]))
nrSamples <- 200

res.samp.R.200 <- sapply(1:nrSamples, simplify=T, function(k) {
  samp <- R.patients[sample(1:length(downsampling.group.size), length(downsampling.group.size), replace=F)]
})

saveRDS(res.samp.R.200, paste(dataDir, "/res.samp.R.200.rds", sep=""))
samp.file.name <- "res.samp.R.200"
eset.file.name <- "filtered.esetALL"

combined.subsample.res <- subsampling.func(samp.file.name, eset.file.name)

#read subsampling results from the cluster
combined.subsample.res<- readRDS(paste(dataDir, "/combined.subsample.res.rds", sep=""))
#set BH adjusted p-value for each dataType separately

lmer.BH.R.subsampling <- lapply(combined.subsample.res, function(m) {
  m$featureID <- row.names(m)
  m$dataTypes <- sub("\\..+", "", m$featureID)
  m$V2.BH.FDR<-c()
  m$V3.BH.FDR<-c()
  m$V2V3.BH.FDR <- c()
  dataTypes <- unique(sub("\\..+", "", m$featureID))
  
  lmer.BH<-do.call('rbind', lapply(dataTypes, function(cur.dt) {
    cur.lmer <- m[grep(paste0(cur.dt), m$featureID),]
    cur.lmer$V2.BH.FDR[!is.na(cur.lmer$V2.Pvalue)]<-p.adjust(cur.lmer$V2.Pvalue[!is.na(cur.lmer$V2.Pvalue)], method = "BH", n = length(cur.lmer$V2.Pvalue))
    cur.lmer$V3.BH.FDR[!is.na(cur.lmer$V3.Pvalue)]<-p.adjust(cur.lmer$V3.Pvalue[!is.na(cur.lmer$V3.Pvalue)], method = "BH", n = length(cur.lmer$V3.Pvalue))
    cur.lmer$V2V3.BH.FDR[!is.na(cur.lmer$V3.2.Pvalue)]<-p.adjust(cur.lmer$V3.2.Pvalue[!is.na(cur.lmer$V3.2.Pvalue)], method = "BH", n = length(cur.lmer$V3.2.Pvalue))
    
    return(cur.lmer)
  }))
  
  return(lmer.BH)
})

#plot per data type separately

dataTypes <-  unique(sub("\\..+", "", lmer.BH.R.subsampling[[1]]$featureID))

count.for.thresholds.res <- lapply(dataTypes, function(cur.dt) {
  count.for.thresholds.list.all.subsampling.cur.data <- lapply(lmer.BH.R.subsampling, function (m) {
    cur.lmer.BH <- m[grep(paste0(cur.dt),m$featureID),]
    cur.lmer.BH.ordered<-cur.lmer.BH[order(cur.lmer.BH$V2.BH.FDR, decreasing=FALSE),]
    count.for.thresholds <- ldply(seq(from = 0.01, to = 1, by = 0.01), function(cur.P) c(cur.P, sum(cur.lmer.BH.ordered$V2.BH.FDR[!is.na(cur.lmer.BH.ordered$V2.BH.FDR)] <= cur.P)))
    colnames(count.for.thresholds) <- c("BH.FDR", "Count")
    return(count.for.thresholds)})
  count.for.thresholds.long <- rbindlist(count.for.thresholds.list.all.subsampling.cur.data, idcol=TRUE)
  return(count.for.thresholds.long)
})


count.for.thresholds.res.V3 <- lapply(dataTypes, function(cur.dt) {
  count.for.thresholds.list.all.subsampling.cur.data <- lapply(lmer.BH.R.subsampling, function (m) {
    cur.lmer.BH <- m[grep(paste0(cur.dt),m$featureID),]
    cur.lmer.BH.ordered<-cur.lmer.BH[order(cur.lmer.BH$V3.BH.FDR, decreasing=FALSE),]
    count.for.thresholds <- ldply(seq(from = 0.01, to = 1, by = 0.01), function(cur.P) c(cur.P, sum(cur.lmer.BH.ordered$V3.BH.FDR[!is.na(cur.lmer.BH.ordered$V3.BH.FDR)] <= cur.P)))
    colnames(count.for.thresholds) <- c("BH.FDR", "Count")
    return(count.for.thresholds)})
  count.for.thresholds.long <- rbindlist(count.for.thresholds.list.all.subsampling.cur.data, idcol=TRUE)
  return(count.for.thresholds.long)
})

count.for.thresholds.res.V3.2 <- lapply(dataTypes, function(cur.dt) {
  count.for.thresholds.list.all.subsampling.cur.data <- lapply(lmer.BH.R.subsampling, function (m) {
    cur.lmer.BH <- m[grep(paste0(cur.dt),m$featureID),]
    cur.lmer.BH.ordered<-cur.lmer.BH[order(cur.lmer.BH$V2V3.BH.FDR, decreasing=FALSE),]
    count.for.thresholds <- ldply(seq(from = 0.01, to = 1, by = 0.01), function(cur.P) c(cur.P, sum(cur.lmer.BH.ordered$V2V3.BH.FDR[!is.na(cur.lmer.BH.ordered$V2V3.BH.FDR)] <= cur.P)))
    colnames(count.for.thresholds) <- c("BH.FDR", "Count")
    return(count.for.thresholds)})
  count.for.thresholds.long <- rbindlist(count.for.thresholds.list.all.subsampling.cur.data, idcol=TRUE)
  return(count.for.thresholds.long)
})


count.for.thresholds.res <- lapply(count.for.thresholds.res, function(m) {
  m$".id" <- factor(paste("subsample",m$".id",sep="."))
  colnames(m) <- c("subsample.num", "BH.FDR", "Count")
  return(m)
})

names(count.for.thresholds.res) <-dataTypes 

FDR_plots.V2 <- lapply(seq(length(count.for.thresholds.res)), function(i) {
  ggplot(count.for.thresholds.res[[i]], aes(x = as.numeric(Count), y = as.numeric(BH.FDR), colour=subsample.num))+ geom_line() + scale_x_log10() + ylab("V2.BH.FDR.threshold")+ xlab("Number of features")+
    theme_bw()+
    ggtitle(paste("data type", names(count.for.thresholds.res)[i]))+
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          strip.text = element_text(size=14), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),                                                                  
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"), 
          axis.ticks.length=unit(.15, "cm"), 
          axis.text.x = element_text(size=14), 
          axis.text.y = element_text(size=14),
          legend.position="none") })

grid.arrange(FDR_plots.V2[[1]], FDR_plots.V2[[2]],  FDR_plots.V2[[3]],FDR_plots.V2[[4]],FDR_plots.V2[[5]], ncol=5) # Write the grid.arrange in the file


count.for.thresholds.res.V3 <- lapply(count.for.thresholds.res.V3, function(m) {
  m$".id" <- factor(paste("subsample",m$".id",sep="."))
  colnames(m) <- c("subsample.num", "BH.FDR", "Count")
  return(m)
})

names(count.for.thresholds.res.V3) <-dataTypes 

FDR_plots.V3 <- lapply(seq(length(count.for.thresholds.res.V3)), function(i) {
  ggplot(count.for.thresholds.res.V3[[i]], aes(x = as.numeric(Count), y = as.numeric(BH.FDR), colour=subsample.num))+ geom_line() + scale_x_log10() + ylab("V3.BH.FDR.threshold")+ xlab("Number of features")+
    theme_bw()+
    ggtitle(paste("data type", names(count.for.thresholds.res.V3)[i]))+
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          strip.text = element_text(size=14), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),                                                                  
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"), 
          axis.ticks.length=unit(.15, "cm"), 
          axis.text.x = element_text(size=14), 
          axis.text.y = element_text(size=14),
          legend.position="none") })

grid.arrange(FDR_plots.V3[[1]], FDR_plots.V3[[2]],  FDR_plots.V3[[3]],FDR_plots.V3[[4]],FDR_plots.V3[[5]], ncol=5) # Write the grid.arrange in the file



count.for.thresholds.res.V3.2 <- lapply(count.for.thresholds.res.V3.2, function(m) {
  m$".id" <- factor(paste("subsample",m$".id",sep="."))
  colnames(m) <- c("subsample.num", "BH.FDR", "Count")
  return(m)
})

names(count.for.thresholds.res.V3.2) <-dataTypes 

FDR_plots.V3.2 <- lapply(seq(length(count.for.thresholds.res.V3.2)), function(i) {
  ggplot(count.for.thresholds.res.V3.2[[i]], aes(x = as.numeric(Count), y = as.numeric(BH.FDR), colour=subsample.num))+ geom_line() + scale_x_log10() + ylab("V3.2.BH.FDR.threshold")+ xlab("Number of features")+
    theme_bw()+
    ggtitle(paste("data type", names(count.for.thresholds.res.V3.2)[i]))+
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          strip.text = element_text(size=14), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),                                                                  
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"), 
          axis.ticks.length=unit(.15, "cm"), 
          axis.text.x = element_text(size=14), 
          axis.text.y = element_text(size=14),
          legend.position="none") })

grid.arrange(FDR_plots.V3.2[[1]], FDR_plots.V3.2[[2]],  FDR_plots.V3.2[[3]],FDR_plots.V3.2[[4]],FDR_plots.V3.2[[5]], ncol=5) # Write the grid.arrange in the file

#For NR
lmer.NR.BH <- readRDS("lmer.NR.BH.rds")
count.for.thresholds.NR.V2<-lapply(dataTypes, function(cur.dt) {
  cur.lmer.BH <- lmer.NR.BH[grep(paste0(cur.dt),lmer.NR.BH$featureID),]
  cur.lmer.BH.ordered<-cur.lmer.BH[order(cur.lmer.BH$V2.BH.FDR, decreasing=FALSE),]
  count.for.thresholds <- ldply(seq(from = 0.01, to = 1, by = 0.01), function(cur.P) c(cur.P, sum(cur.lmer.BH.ordered$V2.BH.FDR[!is.na(cur.lmer.BH.ordered$V2.BH.FDR)] <= cur.P)))
  colnames(count.for.thresholds) <- c("BH.FDR", "Count")
  return(count.for.thresholds)
})


count.for.thresholds.NR.V3<-lapply(dataTypes, function(cur.dt) {
  cur.lmer.BH <- lmer.NR.BH[grep(paste0(cur.dt),lmer.NR.BH$featureID),]
  cur.lmer.BH.ordered<-cur.lmer.BH[order(cur.lmer.BH$V3.BH.FDR, decreasing=FALSE),]
  count.for.thresholds <- ldply(seq(from = 0.01, to = 1, by = 0.01), function(cur.P) c(cur.P, sum(cur.lmer.BH.ordered$V3.BH.FDR[!is.na(cur.lmer.BH.ordered$V3.BH.FDR)] <= cur.P)))
  colnames(count.for.thresholds) <- c("BH.FDR", "Count")
  return(count.for.thresholds)
})


count.for.thresholds.NR.V2V3<-lapply(dataTypes, function(cur.dt) {
  cur.lmer.BH <- lmer.NR.BH[grep(paste0(cur.dt),lmer.NR.BH$featureID),]
  cur.lmer.BH.ordered<-cur.lmer.BH[order(cur.lmer.BH$V2V3.BH.FDR, decreasing=FALSE),]
  count.for.thresholds <- ldply(seq(from = 0.01, to = 1, by = 0.01), function(cur.P) c(cur.P, sum(cur.lmer.BH.ordered$V2V3.BH.FDR[!is.na(cur.lmer.BH.ordered$V2V3.BH.FDR)] <= cur.P)))
  colnames(count.for.thresholds) <- c("BH.FDR", "Count")
  return(count.for.thresholds)
})




count.for.threshold.R.NR.function <- function (count.for.thresholds.res, count.for.thresholds.NR) {
  plot.count.for.threshold.R.NR <- lapply(1:5, function(i) {
    m <- count.for.thresholds.res[[i]]
    data_wide.count.for.thresholds.res <- dcast(m, BH.FDR ~ subsample.num, value.var="Count")
    data_wide.count.for.thresholds.res$mean <- apply(data_wide.count.for.thresholds.res, 1, function(x) {mean(x[!names(x) %in% "BH.FDR"])})
    SEM <- apply(data_wide.count.for.thresholds.res, 1, function(x) {sd(x[!names(x) %in% "BH.FDR"])/sqrt(length(x[!names(x) %in% "BH.FDR"]))})
    data_wide.count.for.thresholds.res$count.min <- data_wide.count.for.thresholds.res$mean-SEM
    data_wide.count.for.thresholds.res$count.max <-data_wide.count.for.thresholds.res$mean + SEM
    data_wide.count.for.thresholds.res$mean <- data_wide.count.for.thresholds.res$mean+1
    data_wide.count.for.thresholds.res$count.min <- data_wide.count.for.thresholds.res$count.min+1
    data_wide.count.for.thresholds.res$count.max <- data_wide.count.for.thresholds.res$count.max+1
    data_wide.count.for.thresholds.res$count.min <- ifelse(data_wide.count.for.thresholds.res$count.min<1, 1, data_wide.count.for.thresholds.res$count.min)
    data_wide.count.for.thresholds.res$count.max <- ifelse(data_wide.count.for.thresholds.res$count.max<1, 1, data_wide.count.for.thresholds.res$count.max)
    m.NR <- count.for.thresholds.NR[[i]]
    m.NR$Count <- m.NR$Count+1
    combined.count <- merge(data_wide.count.for.thresholds.res, m.NR, by="BH.FDR")
    colnames(combined.count)[ncol(combined.count)] <- "NR.count"
    combined.count <- data.frame(combined.count, stringsAsFactors = F)
    combined.count.long <- reshape2::melt(combined.count[,colnames(combined.count) %in% c("BH.FDR", "mean", "NR.count")], id.vars="BH.FDR")
    combined.count.long$count.min <- c(data_wide.count.for.thresholds.res$count.min, rep(0, nrow(combined.count)))
    combined.count.long$count.max <- c(data_wide.count.for.thresholds.res$count.max, rep(0, nrow(combined.count)))
    ggplot(combined.count.long,aes(BH.FDR, value, color=variable))+
      coord_flip()+
      geom_line(size = 0.7, alpha = 0.7)+
      scale_y_continuous(trans="log10")+
      scale_colour_manual(values=c("cyan3", "red"))+
      ggtitle(paste("data type", dataTypes[i]))+
      geom_ribbon(aes(ymin=count.min, ymax=count.max, x=BH.FDR), alpha = 0.2, fill = "cyan3", linetype=0)+
      theme_bw()+
      theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
            legend.position="none",
            strip.text = element_text(size=10), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black"),
            axis.ticks.length=unit(.15, "cm"),
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10))
  })
  grid.arrange(plot.count.for.threshold.R.NR[[1]], plot.count.for.threshold.R.NR[[2]],  plot.count.for.threshold.R.NR[[3]],plot.count.for.threshold.R.NR[[4]],plot.count.for.threshold.R.NR[[5]], ncol=5) # Write the grid.arrange in the file
}

#plot count results
count.for.threshold.R.NR.function(count.for.thresholds.res, count.for.thresholds.NR.V2)
count.for.threshold.R.NR.function(count.for.thresholds.res.V3, count.for.thresholds.NR.V3)
count.for.threshold.R.NR.function(count.for.thresholds.res.V3.2, count.for.thresholds.NR.V2V3)

