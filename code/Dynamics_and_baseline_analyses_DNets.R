#---
#title: "Dynamics and baseline analysis DNets"
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
library(ggbeeswarm))
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
library(pROC)
library(ROCR)

options(stringsAsFactors = F)

code_dir = "./"
source(paste0(code_dir,'Disruption_Networks_functions.R'))

#define the data directory path:
dataDir = '...'
setwd(dataDir)

# Read data

filtered.esetALL <- readRDS(paste(dataDir, "/filtered.esetALL.rds", sep=""))
Perm.Pvalue.V2 <- readRDS(paste(dataDir, "/Perm.Pvalue.V2.rds", sep=""))
Perm.Pvalue.V3 <- readRDS(paste(dataDir, "/Perm.Pvalue.V3.rds", sep=""))
Perm.Pvalue.V2.NR <- readRDS(paste(dataDir, "/Perm.Pvalue.V2.NR.rds", sep=""))
Perm.Pvalue.V3.NR <- readRDS(paste(dataDir, "/Perm.Pvalue.V3.NR.rds", sep=""))
lmer <- readRDS(paste(dataDir, "/lmer.rds", sep=""))
lmer.NR <- readRDS(paste(dataDir, "/lmer.NR.rds", sep=""))
CRP.FCV.df <- read.csv(paste(dataDir, "/GX/CRP.FCV.df.csv", sep=""),stringsAsFactors = F)


#####Inflammatory axis#####
#Figure 1

#use the external axis for R-NR dynamics description
#GEO94648

#Process external data

eset.geo94648 <- getGEO("GSE94648", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]

#test if the data is in log scale

is_logscale <- function (x) 
{
  ex <-  x
  qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 
                                  1), na.rm = T))
  LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 
                              0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 
                                       2)
  !LogC
}

ifelse(is_logscale(exprs(eset.geo94648)), print("eset is in log scale"), print("eset is not log scaled"))

#test normalization (equal medians and quantiles)

test <- data.frame(exprs(eset.geo94648), stringsAsFactors=F)
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

#Add annotatio to probe IDs

eset.geo94648 <- annotateEset(eset.geo94648, hgu133plus2.db)
fData(eset.geo94648)$PROBEID <- as.character(fData(eset.geo94648)$PROBEID)
fData(eset.geo94648)$ENTREZID <- U133plus2.anno.table[match(featureNames(eset.geo94648),U133plus2.anno.table$V1),'V2']
fData(eset.geo94648)$ENTREZID <- as.character(fData(eset.geo94648)$ENTREZID)
fData(eset.geo94648)$SYMBOL <- NULL
fData(eset.geo94648)$GENENAME <- NULL

Description <- select(org.Hs.eg.db, as.character(fData(eset.geo94648)$ENTREZID),keytype = "ENTREZID",  'GENENAME')
fData(eset.geo94648)$Description <- Description[match(fData(eset.geo94648)$ENTREZID, Description$ENTREZID),'GENENAME']
SYMBOL <- select(org.Hs.eg.db, as.character(fData(eset.geo94648)$ENTREZID),keytype = "ENTREZID",  'SYMBOL')
fData(eset.geo94648)$SYMBOL <- SYMBOL[match(fData(eset.geo94648)$ENTREZID, SYMBOL$ENTREZID),'SYMBOL']

eset.geo94648 <- eset.geo94648[!is.na(fData(eset.geo94648)$ENTREZID),]


#Filter and impute data

#Filter features with NA values in more than 20% of the samples 
NA.perc = apply(exprs(eset.geo94648), 1, function(x) sum(is.na(x))/length(x)*100)
filtered.eset.geo94648<-eset.geo94648[!NA.perc>20,]

#For the remaining NAs put average value of the feature

Features.with.NA = which(apply(exprs(filtered.eset.geo94648), 1, function(x) any(is.na(x))))
exprs(filtered.eset.geo94648)[Features.with.NA,] <- do.call('rbind',lapply(Features.with.NA, function(row.index) {
  exprs(filtered.eset.geo94648)[row.index,is.na(exprs(filtered.eset.geo94648)[row.index,])] <- mean(exprs(filtered.eset.geo94648)[row.index,], na.rm=T)
  return(exprs(filtered.eset.geo94648)[row.index,])
}))


#Test DEGs with limma
#Compare active-control separately in CD and UC

#in UC
  filtered.eset.geo94648.Active.control.UC <- filtered.eset.geo94648[,as.character(pData(filtered.eset.geo94648)$characteristics_ch1.2) %in% c("endoscopic_activity: Active", "endoscopic_activity: n.a.") & as.character(pData(filtered.eset.geo94648)$characteristics_ch1.1) %in% c("case_phenotype: Colitis", "case_phenotype: Control")]
  design.Active.control.UC <- model.matrix(~0+ as.character(pData(filtered.eset.geo94648.Active.control.UC)$characteristics_ch1.2))
  colnames(design.Active.control.UC) = c("Active","Control")

  #fit linear model

  data.fit<-lmFit(exprs(filtered.eset.geo94648.Active.control.UC), design.Active.control.UC)
  contrast.matrix = makeContrasts(Active-Control,levels=design.Active.control.UC)
  data.fit.con = contrasts.fit(data.fit,contrast.matrix)
  data.fit.eb = eBayes(data.fit.con)
  limma_results.df.Active.control.UC = topTable(data.fit.eb, coef=1, n=length(data.fit.eb$coefficients), adjust="none")
  limma_results.df.Active.control.UC$featureID <- row.names(limma_results.df.Active.control.UC)
  limma_results.df.Active.control.UC$BH.FDR <- p.adjust(limma_results.df.Active.control.UC$P.Value, method = "BH", n = length(limma_results.df.Active.control.UC$P.Value))
  limma_results.df.sig.Active.control.UC<-limma_results.df.Active.control.UC[limma_results.df.Active.control.UC$BH.FDR<0.1,]

#in CD

  filtered.eset.geo94648.Active.control.CD <- filtered.eset.geo94648[,as.character(pData(filtered.eset.geo94648)$characteristics_ch1.2) %in% c("endoscopic_activity: Active", "endoscopic_activity: n.a.") & as.character(pData(filtered.eset.geo94648)$characteristics_ch1.1) %in% c("case_phenotype: Crohn", "case_phenotype: Control")]
  design.Active.control.CD <- model.matrix(~0+ as.character(pData(filtered.eset.geo94648.Active.control.CD)$characteristics_ch1.2))
  colnames(design.Active.control.CD) = c("Active","Control")

  #fit linear model

  data.fit<-lmFit(exprs(filtered.eset.geo94648.Active.control.CD), design.Active.control.CD)
  contrast.matrix = makeContrasts(Active-Control,levels=design.Active.control.CD)
  data.fit.con = contrasts.fit(data.fit,contrast.matrix)
  data.fit.eb = eBayes(data.fit.con)
  limma_results.df.Active.control.CD = topTable(data.fit.eb, coef=1, n=length(data.fit.eb$coefficients), adjust="none")
  limma_results.df.Active.control.CD$featureID <- row.names(limma_results.df.Active.control.CD)
  limma_results.df.Active.control.CD$BH.FDR <- p.adjust(limma_results.df.Active.control.CD$P.Value, method = "BH", n = length(limma_results.df.Active.control.CD$P.Value))
  limma_results.df.sig.Active.control.CD<-limma_results.df.Active.control.CD[limma_results.df.Active.control.CD$BH.FDR<0.1,]
  

DE.tot.final <- rbind.fill(limma_results.df.sig.Active.control.UC, limma_results.df.sig.Active.control.CD)
DE.tot.final <- DE.tot.final[!duplicated(DE.tot.final$featureID), ]
DE.tot.final$ENTREZID <- fData(filtered.eset.geo94648)[match(DE.tot.final$featureID, row.names(fData(filtered.eset.geo94648))),'ENTREZID']
DE.tot.final$ENTREZID <- as.character(DE.tot.final$ENTREZID)  
DE.tot.final$SYMBOL <- fData(filtered.eset.geo94648)[match(DE.tot.final$featureID, row.names(fData(filtered.eset.geo94648))),'SYMBOL']


#plot PCA

filtered.eset.geo94648.DE <-filtered.eset.geo94648[row.names(fData(filtered.eset.geo94648)) %in% unique(DE.tot.final$featureID),]

all.dataType.pca.geo94648.DE<- prcomp(t(exprs(filtered.eset.geo94648.DE)),
                                         center = TRUE,
                                         scale = TRUE) 


pc_dat.geo94648.DE<- data.frame(all.dataType.pca.geo94648.DE$x[,1:20],
                                   Disease=as.factor(as.character(pData(filtered.eset.geo94648.DE)$characteristics_ch1.1)),
                                   Disease.activity=as.factor(pData(filtered.eset.geo94648.DE)$characteristics_ch1.2))


percentage <- round(all.dataType.pca.geo94648.DE$sdev^2 / sum(all.dataType.pca.geo94648.DE$sdev^2) * 100, 2)
pc_dat.geo94648.DE$label <- row.names(pc_dat.geo94648.DE)

#test which delta PC predict gradient of Active,Inactive and
#control group (ordinal response variable)

dat <- pc_dat.geo94648.DE[,!colnames(pc_dat.geo94648.DE) %in% c("Disease", "Disease.activity", "label")]

x <- as.matrix(dat)
y <- pc_dat.geo94648.DE$Disease.activity
names(y) <- row.names(x)
y <- factor(y, levels=c("endoscopic_activity: Active", "endoscopic_activity: Inactive", "endoscopic_activity: n.a."))            

#glmnet is by

fit <- glmnetcr(x,y, alpha=1, nlambda = 100, pmax=1000, maxit=1000)
BIC.step <- select.glmnetcr(fit, which = "BIC")
coefficients<-as.matrix(coef(fit, s = BIC.step))
coefficients[[2]]

#According to ordinal lasso PC1 and PC4 were selected

ggplot(pc_dat.geo94648.DE,aes(x=PC1, y=PC4, col=Disease.activity, shape=Disease)) + 
  geom_point(size=3) + 
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC1 (", percentage[1], "%) ", sep="")) +ylab(paste("PC4 (", percentage[4], "%) ", sep=""))

#Project our in-house cohort data on external PCA based the inflammatory axis 

filtered.eset.DE <- filtered.esetALL[fData(filtered.esetALL)$ENTREZID %in% DE.tot.final$ENTREZID & fData(filtered.esetALL)$dataType %in% "GX",]
exprs.filtered.eset.DE <- t(exprs(filtered.eset.DE))
fData(filtered.eset.DE)$PROBID <- fData(filtered.eset.geo94648.DE)[match(fData(filtered.eset.DE)$ENTREZID, fData(filtered.eset.geo94648.DE)$ENTREZID),'PROBEID']
exprs.filtered.eset.DE.t <- data.frame(t(exprs.filtered.eset.DE), stringsAsFactors = F)
exprs.filtered.eset.DE.t$PROBEID <- fData(filtered.eset.DE)[match(row.names(exprs.filtered.eset.DE.t), fData(filtered.eset.DE)$featureID),'PROBID']
exprs.filtered.eset.DE.t$meanEXP <- apply(exprs.filtered.eset.DE.t[,!colnames(exprs.filtered.eset.DE.t) %in% "PROBEID"], 1, mean)
exprs.filtered.eset.DE.t.sub <- aggregate(exprs.filtered.eset.DE.t, by=list(exprs.filtered.eset.DE.t$PROBEID) , FUN=max)
exprs.filtered.eset.DE.t.sub$meanEXP <- NULL
row.names(exprs.filtered.eset.DE.t.sub) <- exprs.filtered.eset.DE.t.sub$PROBEID
exprs.filtered.eset.DE.t.sub <- exprs.filtered.eset.DE.t.sub[,!names(exprs.filtered.eset.DE.t.sub) %in% c("PROBEID", "Group.1", "meanEXP")]
exprs.filtered.eset.DE.sub <- data.frame(t(exprs.filtered.eset.DE.t.sub), stringsAsFactors = F)
colnames(exprs.filtered.eset.DE.sub) <- sub("^X", "", colnames(exprs.filtered.eset.DE.sub))

#test only genes that exists in both data sets
filtered.eset.geo94648.DE.final <- filtered.eset.geo94648.DE[row.names(fData(filtered.eset.geo94648.DE)) %in% colnames(exprs.filtered.eset.DE.sub),]
exprs.filtered.eset.geo94648.DE.final <- data.frame(t(exprs(filtered.eset.geo94648.DE.final)), stringsAsFactors = F)

all.dataType.pca.geo94648.DE<- prcomp(exprs.filtered.eset.geo94648.DE.final,
                                      scale = TRUE,
                                      retx=TRUE) 

#plot the PCA of external data

pc_dat.geo94648.DE<- data.frame(all.dataType.pca.geo94648.DE$x[,1:20],
                                Disease=as.factor(as.character(pData(filtered.eset.geo94648.DE.final)$characteristics_ch1.1)),
                                Disease.activity=as.factor(pData(filtered.eset.geo94648.DE.final)$characteristics_ch1.2))

percentage <- round(all.dataType.pca.geo94648.DE$sdev^2 / sum(all.dataType.pca.geo94648.DE$sdev^2) * 100, 2)

ggplot(pc_dat.geo94648.DE,aes(x=PC1, y=PC4, col=Disease.activity, shape=Disease)) + 
  geom_point(size=3) + 
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC1 (", percentage[1], "%) ", sep="")) +ylab(paste("PC4 (", percentage[4], "%) ", sep=""))


#Plot selected PCs by ordinal lasso and variance exlained

PC.selected.df <- data.frame(PC=names(coefficients[[2]]), coef=abs(coefficients[[2]]), stringsAsFactors = F)
PC.selected.df <- PC.selected.df[1:20,]
PC.selected.df$Variance.explained <-percentage[1:20]
PC.selected.df$combined <- (scale(PC.selected.df$coef)+scale(PC.selected.df$Variance.explained))/2

#Calculate coordinates of the external inflammatory axis

x1 <- mean(pc_dat.geo94648.DE[pc_dat.geo94648.DE$Disease.activity %in% "endoscopic_activity: Active",'PC1'])
x2 <-  mean(pc_dat.geo94648.DE[pc_dat.geo94648.DE$Disease.activity %in% "endoscopic_activity: n.a.",'PC1'])
y1 <- mean(pc_dat.geo94648.DE[pc_dat.geo94648.DE$Disease.activity %in% "endoscopic_activity: Active",'PC4'])
y2 <- mean(pc_dat.geo94648.DE[pc_dat.geo94648.DE$Disease.activity %in% "endoscopic_activity: n.a.",'PC4'])

external.health.axis <- c(x2-x1, y2-y1)
m <- (y2-y1)/(x2-x1)

intercept <- y1-(m*x1)


#Plot PCA based on the informative PCs
#Fig 1b

ggplot(pc_dat.geo94648.DE,aes(x=PC1, y=PC4, col=Disease.activity, shape=Disease)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2),
               size=0.8,
               arrow = arrow(length = unit(0.4, "cm"), type="closed"),col='darkgrey', linetype=2)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC1 (", percentage[1], "%) ", sep="")) +ylab(paste("PC4 (", percentage[4], "%) ", sep=""))


#now project the new data on the loadings of the external data
c.fun<-function(df, center, scale) {
  return((df-center)/scale )
}

centeredData<-apply(exprs.filtered.eset.DE.sub, 1, FUN=c.fun, all.dataType.pca.geo94648.DE$center, all.dataType.pca.geo94648.DE$scale  )

pcs<-t(all.dataType.pca.geo94648.DE$rotation) %*% centeredData
pcs.t <- data.frame(t(pcs), stringsAsFactors = F)

pc_dat <- pcs.t[,1:20]

#project the new data PCA
row.names(pc_dat) <- gsub("\\.", "-", row.names(pc_dat))

pc_dat$group <- pData(filtered.eset.DE)[match(row.names(pc_dat), pData(filtered.eset.DE)$sampleID),'group']
pc_dat$group <- ifelse(pc_dat$group %in% "IM", "NR", pc_dat$group)
pc_dat$Patient.code <- pData(filtered.eset.DE)[match(row.names(pc_dat), pData(filtered.eset.DE)$sampleID),'Patient.code']
pc_dat$Visit <- pData(filtered.eset.DE)[match(row.names(pc_dat), pData(filtered.eset.DE)$sampleID),'Visit']


ggplot(pc_dat,aes(x=PC1, y=PC4, col=group, group=Patient.code, abel = Patient.code)) + 
  geom_point() + 
  geom_path(arrow=arrow(length = unit(0.3, "cm"), type="closed"))  +
  geom_text(aes(label = Patient.code) ,hjust=0, vjust=0)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC1 (", percentage[1], "%) ", sep="")) +ylab(paste("PC4 (", percentage[4], "%) ", sep=""))

pc_dat.V1V3 <- pc_dat[pc_dat$Visit %in% c("V1", "V3"),]
pc_dat.V1V2 <- pc_dat[pc_dat$Visit %in% c("V1", "V2"),]

ggplot(pc_dat.V1V3,aes(x=PC1, y=PC4, col=group, label=Patient.code, group=Patient.code)) + 
  geom_point() + 
  geom_path(arrow=arrow(length = unit(0.3, "cm"), type="closed"))  +
  geom_segment(aes(x=x3, xend=x4, y=y3, yend=y4),
               size=0.8,
               arrow = arrow(length = unit(0.4, "cm"), type="closed"),col='darkgrey', linetype=2)+

  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))


ggplot(pc_dat.V1V2,aes(x=PC1, y=PC4, col=group, label=Patient.code, group=Patient.code)) + 
  geom_point() + 
  geom_path(arrow=arrow(length = unit(0.3, "cm"), type="closed"))  +
  geom_segment(aes(x=x3, xend=x4, y=y3, yend=y4),
               size=0.8,
               arrow = arrow(length = unit(0.4, "cm"), type="closed"),col='darkgrey', linetype=2)+
  geom_text(aes(label = Patient.code) ,hjust=0, vjust=0)+
  
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  
  xlab(paste("PC1 (", percentage[1], "%) ", sep="")) +ylab(paste("PC4 (", percentage[4], "%) ", sep=""))


#Calculate projection of V1V2 and V1V3 
#for each patient on the external axis

pc_coordinates.on.external.health.axis <- pc_dat
pc_coordinates.on.external.health.axis$projection <- do.call('rbind', lapply(seq(nrow(pc_dat)), function(i) {
  dot.vector <- c(pc_dat[i,'PC1'],pc_dat[i,'PC4'])
  projection <- external.health.axis%*%dot.vector/(sqrt(external.health.axis[1]^2+external.health.axis[2]^2))
  return(projection)
}))

pc_coordinates.on.external.health.axis.rearranged <- dcast(pc_coordinates.on.external.health.axis, Patient.code +group ~ Visit, value.var="projection")
pc_coordinates.on.external.health.axis.rearranged$V1V3.projection <- pc_coordinates.on.external.health.axis.rearranged$V3-pc_coordinates.on.external.health.axis.rearranged$V1
pc_coordinates.on.external.health.axis.rearranged$V1V2.projection <- pc_coordinates.on.external.health.axis.rearranged$V2-pc_coordinates.on.external.health.axis.rearranged$V1
pc_coordinates.on.external.health.axis.rearranged$V2V3.projection <- pc_coordinates.on.external.health.axis.rearranged$V3-pc_coordinates.on.external.health.axis.rearranged$V2

pc_coordinates.on.external.health.axis.rearranged.long <- melt(pc_coordinates.on.external.health.axis.rearranged[,!colnames(pc_coordinates.on.external.health.axis.rearranged) %in% c("V1", "V2", "V3")])
pc_coordinates.on.external.health.axis.rearranged.long$class <- paste(pc_coordinates.on.external.health.axis.rearranged.long$group, pc_coordinates.on.external.health.axis.rearranged.long$variable, sep=".")
pc_coordinates.on.external.health.axis.rearranged.long$class <- sub(".projection", "", pc_coordinates.on.external.health.axis.rearranged.long$class)
pc_coordinates.on.external.health.axis.rearranged.long$class <- factor(pc_coordinates.on.external.health.axis.rearranged.long$class, levels=c("R.V1V2", "R.V2V3", "R.V1V3", "NR.V1V2", "NR.V2V3", "NR.V1V3"))
pc_coordinates.on.external.health.axis.rearranged.long$group <- factor(pc_coordinates.on.external.health.axis.rearranged.long$group, levels=c("R", "NR"))

pc_coordinates.on.external.health.axis.rearranged.long$variable <- factor(pc_coordinates.on.external.health.axis.rearranged.long$variable, levels=c("V1V2.projection", "V2V3.projection", "V1V3.projection"))

#plot combined projection of v1v2, v2v3 and v1v3 in a single plot
ggplot(pc_coordinates.on.external.health.axis.rearranged.long[!pc_coordinates.on.external.health.axis.rearranged.long$Patient.code %in% c("37"),], aes(x = group, y = as.numeric(value), fill=group, colour=factor(group))) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  ylab("PCA based projection on health axis")+
  scale_fill_manual(values=c("NR"="red", "R"="cyan3"))+
  scale_color_manual(values=c("NR"="black", "R"="black"))+
  facet_wrap(~variable, scales = "fix")+
  geom_hline(yintercept=0, linetype="dashed", color="darkgrey")+
  theme_bw()+theme(strip.text = element_text(size=14), panel.grid.major = element_blank(),
                   axis.text.x = element_text(size=10),
                   axis.text.y = element_text(size=10),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "greater"), size=4, comparisons=list(c(1,2)), label.y = 30)


#Permutation test for distance differences

pc_coordinates.on.external.health.axis.rearranged.long <- pc_coordinates.on.external.health.axis.rearranged.long[!pc_coordinates.on.external.health.axis.rearranged.long$Patient.code %in% "37",]
perm.test.projected.distance <- do.call('rbind', lapply(unique(pc_coordinates.on.external.health.axis.rearranged.long$variable), function(cur.facet) {
nsim <- 100000
res <- numeric(nsim) 
dat <- pc_coordinates.on.external.health.axis.rearranged.long[pc_coordinates.on.external.health.axis.rearranged.long$variable %in% cur.facet,]
for (i in 1:nsim) {
  perm <- sample(nrow(dat))
  bdat <- transform(dat,perm=dat[perm,])
  res[i] <- mean(bdat[bdat$group=="R","perm.value"])-
    mean(bdat[bdat$group=="NR","perm.value"])
}
obs <- mean(dat[dat$group=="R","value"])-
  mean(dat[dat$group=="NR","value"])
res <- c(res,obs)
hist(res,col=muted("red"),las=1,main="")
abline(v=obs,col="red")

p <- mean(res>obs)

return(data.frame(variable=cur.facet, group1="R", group2="NR", p=p, method="permutation test", stringsAsFactors = F))
}))

perm.test.projected.distance$group1 <- c("R.V1V3", "R.V1V2", "R.V2V3")
perm.test.projected.distance$group2 <- c("NR.V1V3", "NR.V1V2", "NR.V2V3")

#For extraction of xy position 
projection.df_stat <- pc_coordinates.on.external.health.axis.rearranged.long  %>% 
  pairwise_wilcox_test(value ~ class, comparisons=list(c("R.V1V2","NR.V1V2"), c("R.V2V3","NR.V2V3"), c("R.V1V3","NR.V1V3"))) %>% 
  add_xy_position() 
projection.df_stat$xmin <- 1
projection.df_stat$xmax <- 2

projection.df_stat.sub.final <- merge(projection.df_stat, perm.test.projected.distance, by=c("group1", "group2"))
projection.df_stat.sub.final$p.y <- round(projection.df_stat.sub.final$p.y, 4)

#Plot projected distance  - Fig 1c
#for each patient on the external axis

ggplot(pc_coordinates.on.external.health.axis.rearranged.long , aes(x = class, y = as.numeric(value), colour="black")) +
  geom_boxplot(alpha=1, outlier.shape = NA, aes(fill=group)) +
  scale_color_manual(values=c("black"="black"))+
  scale_fill_manual(values=c("R"="cyan3", "NR"="Indianred"))+
  facet_wrap(~variable, scale="free", ncol=5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  ylab("PCA based projection on health axis")+
  theme_bw()+
  geom_hline(yintercept=0, linetype="dashed", color="darkgrey", size=0.1)+
  theme(strip.text = element_text(size=12), panel.grid.major = element_blank(),
        axis.text.y=element_text( size=10),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border =element_blank(), 
        axis.ticks.length=unit(.15, "cm"))+
  ggpubr::stat_pvalue_manual(data = projection.df_stat.sub.final, label = "p.y", size=4)


#Plot correlation between V1V2 and V2V3
#for responders and non-responders separately
#Fig 1d


pc_coordinates.on.external.health.axis.rearranged$group <- factor(pc_coordinates.on.external.health.axis.rearranged$group, levels=c("R", "NR"))

ggplot(pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",],aes(x=V1V2.projection, y=V2V3.projection, fill=group)) + 
  geom_point(size=3, shape=21)+theme_bw()+
  scale_fill_manual(values=c("R"="cyan3", "NR"="indianred1"))+ 
  scale_color_manual(values=c("black"))+
  geom_smooth(data=pc_coordinates.on.external.health.axis.rearranged[pc_coordinates.on.external.health.axis.rearranged$group %in% "R",], method='lm', se=T, aes(color="black"))+
  theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))

corr.test(pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",]$V1V2.projection[pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",]$group %in% "R"], pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",]$V2V3.projection[pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",]$group %in% "R"], method="spearman")$r
corr.test(pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",]$V1V2.projection[pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",]$group %in% "R"], pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",]$V2V3.projection[pc_coordinates.on.external.health.axis.rearranged[!pc_coordinates.on.external.health.axis.rearranged$Patient.code %in% "37",]$group %in% "R"], method="spearman")$p

corr.test(pc_coordinates.on.external.health.axis.rearranged$V1V2.projection[pc_coordinates.on.external.health.axis.rearranged$group %in% "NR"], pc_coordinates.on.external.health.axis.rearranged$V2V3.projection[pc_coordinates.on.external.health.axis.rearranged$group %in% "NR"], method="spearman")$r
corr.test(pc_coordinates.on.external.health.axis.rearranged$V1V2.projection[pc_coordinates.on.external.health.axis.rearranged$group %in% "NR"], pc_coordinates.on.external.health.axis.rearranged$V2V3.projection[pc_coordinates.on.external.health.axis.rearranged$group %in% "NR"], method="spearman")$p


#draw patient individual progress
#on the inflammatory axis
#Fig 1b , right

pc_coordinates.on.external.health.axis.rearranged.V1V2 <- pc_coordinates.on.external.health.axis.rearranged[,colnames(pc_coordinates.on.external.health.axis.rearranged) %in% c('Patient.code', "group", "V1", "V2")]
pc_coordinates.on.external.health.axis.rearranged.V2V3 <- pc_coordinates.on.external.health.axis.rearranged[,colnames(pc_coordinates.on.external.health.axis.rearranged) %in% c('Patient.code', "group", "V2", "V3")]
pc_coordinates.on.external.health.axis.rearranged.V2V3$class <- paste(pc_coordinates.on.external.health.axis.rearranged.V2V3$Patient.code, "V2V3", sep=".")
pc_coordinates.on.external.health.axis.rearranged.V1V2$class <- paste(pc_coordinates.on.external.health.axis.rearranged.V1V2$Patient.code, "V1V2", sep=".")


colnames(pc_coordinates.on.external.health.axis.rearranged.V2V3) <- c("Patient.code", "group", "start", "end", "class")
colnames(pc_coordinates.on.external.health.axis.rearranged.V1V2) <- c("Patient.code", "group", "start", "end", "class")

pc_coordinates.on.external.health.axis.rearranged.combined <- rbind(pc_coordinates.on.external.health.axis.rearranged.V1V2, pc_coordinates.on.external.health.axis.rearranged.V2V3)
pc_coordinates.on.external.health.axis.rearranged.combined$Absolute.baseline <- pc_coordinates.on.external.health.axis.rearranged[match(pc_coordinates.on.external.health.axis.rearranged.combined$Patient.code, pc_coordinates.on.external.health.axis.rearranged$Patient.code),'V1']
pc_coordinates.on.external.health.axis.rearranged.combined <- pc_coordinates.on.external.health.axis.rearranged.combined[order(pc_coordinates.on.external.health.axis.rearranged.combined$Absolute.baseline, pc_coordinates.on.external.health.axis.rearranged.combined$Patient.code, pc_coordinates.on.external.health.axis.rearranged.combined$class),]
pc_coordinates.on.external.health.axis.rearranged.combined$class <- factor(pc_coordinates.on.external.health.axis.rearranged.combined$class, levels=pc_coordinates.on.external.health.axis.rearranged.combined$class)

ggplot(pc_coordinates.on.external.health.axis.rearranged.combined,
       aes(y = class, col=group)) +
  labs(x = "Location on health axis", y = "Patient") +
  scale_color_manual(values=c("R"="cyan3", "NR"="indianred1"))+
  geom_segment(aes(x = pc_coordinates.on.external.health.axis.rearranged.combined$start,
                   y = class,
                   xend = pc_coordinates.on.external.health.axis.rearranged.combined$end,
                   yend = class),
               size = 1, arrow=arrow(length = unit(0.3, "cm"), type="closed"))+
  geom_hline(data=pc_coordinates.on.external.health.axis.rearranged.combined, yintercept=seq(0.5,nrow(pc_coordinates.on.external.health.axis.rearranged.combined), by=2))+
  scale_y_discrete(limits = rev(levels(pc_coordinates.on.external.health.axis.rearranged.combined$class)))+
  theme_bw()+
  theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))


pc_coordinates.on.external.health.axis.rearranged.V1V3 <- pc_coordinates.on.external.health.axis.rearranged[,colnames(pc_coordinates.on.external.health.axis.rearranged) %in% c('Patient.code', "group", "V1", "V2")]
colnames(pc_coordinates.on.external.health.axis.rearranged.V1V3) <- c("Patient.code", "group", "start", "end")
pc_coordinates.on.external.health.axis.rearranged.V1V3$delta <- pc_coordinates.on.external.health.axis.rearranged.V1V3$end-pc_coordinates.on.external.health.axis.rearranged.V1V3$start
pc_coordinates.on.external.health.axis.rearranged.V1V3$Absolute.baseline <- pc_coordinates.on.external.health.axis.rearranged[match(pc_coordinates.on.external.health.axis.rearranged.V1V3$Patient.code, pc_coordinates.on.external.health.axis.rearranged$Patient.code),'V1']
pc_coordinates.on.external.health.axis.rearranged.V1V3 <- pc_coordinates.on.external.health.axis.rearranged.V1V3[!pc_coordinates.on.external.health.axis.rearranged.V1V3$Patient.code %in% "37",]
pc_coordinates.on.external.health.axis.rearranged.V1V3 <- pc_coordinates.on.external.health.axis.rearranged.V1V3[order(pc_coordinates.on.external.health.axis.rearranged.V1V3$group, pc_coordinates.on.external.health.axis.rearranged.V1V3$delta),]
pc_coordinates.on.external.health.axis.rearranged.V1V3$Patient.code <- factor(pc_coordinates.on.external.health.axis.rearranged.V1V3$Patient.code, levels=pc_coordinates.on.external.health.axis.rearranged.V1V3$Patient.code)


ggplot(pc_coordinates.on.external.health.axis.rearranged.V1V3,
       aes(y = Patient.code, col=group)) +
  labs(x = "Location on health axis W2-Baseline", y = "Patient") +
  scale_color_manual(values=c("R"="cyan3", "NR"="indianred1"))+
  geom_segment(aes(x = pc_coordinates.on.external.health.axis.rearranged.V1V3$start,
                   y = Patient.code,
                   xend = pc_coordinates.on.external.health.axis.rearranged.V1V3$end,
                   yend = Patient.code),
               size = 0.5, arrow=arrow(length = unit(0.3, "cm"), type="closed"))+
  scale_y_discrete(limits = levels(pc_coordinates.on.external.health.axis.rearranged.V1V3$Patient.code))+
  theme_bw()+
  theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))



######Chracterization of normal infliximab dynamics (in responders)####
#Figure 2


#Summary of features changing over time
#Figure 2b

#extract significant features from responders and non-responders using permutation with FDR threshold of 0.15
#thresholds for FDR <0.15

Perm.Pvalue.R.V2.thrFDR0.15<-subset(Perm.Pvalue.V2, Perm.Pvalue.V2$perm.pvalue<0.010 & grepl("AG", Perm.Pvalue.V2$featureID)|
                                      Perm.Pvalue.V2$perm.pvalue<0.015 & grepl("PB", Perm.Pvalue.V2$featureID)|
                                      Perm.Pvalue.V2$perm.pvalue<0 & grepl("GR", Perm.Pvalue.V2$featureID)|
                                      Perm.Pvalue.V2$perm.pvalue<0 & grepl("LU", Perm.Pvalue.V2$featureID))



Perm.Pvalue.R.V3.thrFDR0.15<-subset(Perm.Pvalue.V3, Perm.Pvalue.V3$perm.pvalue<=0.001 & grepl("AG", Perm.Pvalue.V3$featureID)|
                                      Perm.Pvalue.V3$perm.pvalue<=0.01 & grepl("PB", Perm.Pvalue.V3$featureID)|
                                      Perm.Pvalue.V3$perm.pvalue<=0.01 & grepl("GR", Perm.Pvalue.V3$featureID)|
                                      Perm.Pvalue.V3$perm.pvalue<0 & grepl("LU", Perm.Pvalue.V3$featureID))


Perm.Pvalue.NR.V2.thrFDR0.15 <- subset(Perm.Pvalue.V2.NR, Perm.Pvalue.V2.NR$perm.pvalue<=0.001 & grepl("AG", Perm.Pvalue.V2.NR$featureID)|
                                         Perm.Pvalue.V2.NR$perm.pvalue<0 & grepl("PB", Perm.Pvalue.V2.NR$featureID)|
                                         Perm.Pvalue.V2.NR$perm.pvalue<0 & grepl("GR", Perm.Pvalue.V2.NR$featureID)|
                                         Perm.Pvalue.V2.NR$perm.pvalue<0 & grepl("LU", Perm.Pvalue.V2.NR$featureID))


Perm.Pvalue.NR.V3.thrFDR0.15 <- subset(Perm.Pvalue.V3.NR, Perm.Pvalue.V3.NR$perm.pvalue<0 & grepl("AG", Perm.Pvalue.V3.NR$featureID)|
                                         Perm.Pvalue.V3.NR$perm.pvalue<0 & grepl("PB", Perm.Pvalue.V3.NR$featureID)|
                                         Perm.Pvalue.V3.NR$perm.pvalue<0 & grepl("GR", Perm.Pvalue.V3.NR$featureID)|
                                         Perm.Pvalue.V3.NR$perm.pvalue<0 & grepl("LU", Perm.Pvalue.V3.NR$featureID))

#Add Venn diagram comparing overlapping features

VennDiagram.res <- VennDiagram::get.venn.partitions(list(Perm.Pvalue.R.V2.thrFDR0.15=Perm.Pvalue.R.V2.thrFDR0.15$featureID, Perm.Pvalue.R.V3.thrFDR0.15=Perm.Pvalue.R.V3.thrFDR0.15$featureID, Perm.Pvalue.NR.V2.thrFDR0.15=Perm.Pvalue.NR.V2.thrFDR0.15$featureID))
grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(list(Perm.Pvalue.R.V2.thrFDR0.15=Perm.Pvalue.R.V2.thrFDR0.15$featureID, Perm.Pvalue.R.V3.thrFDR0.15=Perm.Pvalue.R.V3.thrFDR0.15$featureID, Perm.Pvalue.NR.V2.thrFDR0.15=Perm.Pvalue.NR.V2.thrFDR0.15$featureID), NULL))

#for response network select only responders features

changesV2 <- data.frame(featureID=unique(c(Perm.Pvalue.R.V2.thrFDR0.15$featureID)), stringsAsFactors = F)
changesV3 <- data.frame(featureID=unique(c(Perm.Pvalue.R.V3.thrFDR0.15$featureID)), stringsAsFactors = F)

changesV2$SYMBOL <- fData(filtered.esetALL)[match(changesV2$featureID, fData(filtered.esetALL)$featureID),'SYMBOL']
changesV3$SYMBOL <- fData(filtered.esetALL)[match(changesV3$featureID, fData(filtered.esetALL)$featureID),'SYMBOL']

#Prepare supplementary table of permutation results of selected dynamic features
#For responders
lmer.sub <- lmer[lmer$featureID %in% unique(c(changesV2$featureID,changesV3$featureID)),]
Perm.selected.features <- merge(lmer.sub, Perm.Pvalue.V2, by="featureID", all.x=T)
Perm.selected.features <- merge(Perm.selected.features, Perm.Pvalue.V3, by="featureID", all.x=T)
Perm.selected.features$SYMBOL <- fData(filtered.esetALL)[match(Perm.selected.features$featureID, fData(filtered.esetALL)$featureID),'displayName.final']

#For non-responders

lmer.NR.sub <- lmer.NR[lmer.NR$featureID %in% unique(c(Perm.Pvalue.NR.V2.thrFDR0.15$featureID,Perm.Pvalue.NR.V3.thrFDR0.15$featureID)),]
Perm.selected.features.NR <- merge(lmer.NR.sub, Perm.Pvalue.V2.NR, by="featureID", all.x=T)
Perm.selected.features.NR <- merge(Perm.selected.features.NR, Perm.Pvalue.V3.NR, by="featureID", all.x=T)
Perm.selected.features.NR$SYMBOL <- fData(filtered.esetALL)[match(Perm.selected.features.NR$featureID, fData(filtered.esetALL)$featureID),'displayName.final']



######network propagation######

#add nodes for network propagation using stringdb

linker.genes.final.featureIDs.V2 <- network.propagation.function(changesV2, "V2", filtered.esetALL)
linker.genes.final.featureIDs.V3 <- network.propagation.function(changesV3, "V3", filtered.esetALL)

#select linkers if they are significantly enriched in the core list by a threshold of FDR<0.01 
linker.genes.final.featureIDs.V2 <- linker.genes.final.featureIDs.V2[linker.genes.final.featureIDs.V2$BH.FDR<0.01,]
linker.genes.final.featureIDs.V3 <- linker.genes.final.featureIDs.V3[linker.genes.final.featureIDs.V3$BH.FDR<0.01,]

#Add featureID to linkers

filtered.esetALL.AG <- filtered.esetALL[fData(filtered.esetALL)$dataType %in% "AGX",]
linker.genes.final.featureIDs.V2$featureID <- fData(filtered.esetALL.AG)[match(linker.genes.final.featureIDs.V2$ENTRZID, fData(filtered.esetALL.AG)$ENTREZID),'featureID']
linker.genes.final.featureIDs.V3$featureID <- fData(filtered.esetALL.AG)[match(linker.genes.final.featureIDs.V3$ENTRZID, fData(filtered.esetALL.AG)$ENTREZID),'featureID']

#combine the linker genes to changes final

changesV2.final <- unique(c(changesV2$featureID, linker.genes.final.featureIDs.V2$featureID))
changesV3.final <- unique(c(changesV3$featureID, linker.genes.final.featureIDs.V3$featureID))

changes.final.all <- unique(c(changesV2.final, changesV3.final))

#Supplementary table of linker genes

linker.genes.Supp.table <- merge(linker.genes.final.featureIDs.V2, linker.genes.final.featureIDs.V3, by="linker.name", all=T)
linker.genes.Supp.table <- linker.genes.Supp.table[!(linker.genes.Supp.table$featureID.x %in% unique(c(changesV2$featureID, changesV3$featureID))|linker.genes.Supp.table$featureID.y %in% unique(c(changesV2$featureID, changesV3$featureID))),]

linker.genes.Supp.table$combined.enrichment.p <- ifelse(!is.na(linker.genes.Supp.table$enrichment.p.x) & !is.na(linker.genes.Supp.table$enrichment.p.y), 
(linker.genes.Supp.table$enrichment.p.x+linker.genes.Supp.table$enrichment.p.y)/2, 
ifelse(!is.na(linker.genes.Supp.table$enrichment.p.x) & is.na(linker.genes.Supp.table$enrichment.p.y), linker.genes.Supp.table$enrichment.p.x, linker.genes.Supp.table$enrichment.p.y))

linker.genes.Supp.table$combined.BH.FDR <- ifelse(!is.na(linker.genes.Supp.table$BH.FDR.x) & !is.na(linker.genes.Supp.table$BH.FDR.y), 
                                                        (linker.genes.Supp.table$BH.FDR.x+linker.genes.Supp.table$BH.FDR.y)/2, 
                                                        ifelse(!is.na(linker.genes.Supp.table$BH.FDR.x) & is.na(linker.genes.Supp.table$BH.FDR.y), linker.genes.Supp.table$BH.FDR.x, linker.genes.Supp.table$BH.FDR.y))



#######Test changes in cell abundances#######
#PCA for abundances only
#Plot Fig 2a


# compute log FC
filtered.esetALL.abundance<-filtered.esetALL[grep("Abundance|CyTOF-PBMC|CyTOF-GRAN|GR.TOTAL|PB.TOTAL", fData(filtered.esetALL)$displayName.final),]
filtered.esetALL.abundance.FC <- substractReference(filtered.esetALL.abundance, 'Patient.code', 'Visit', ref = 'V1', op = '/')
exprs(filtered.esetALL.abundance.FC)[is.infinite(exprs(filtered.esetALL.abundance.FC)) | is.nan(exprs(filtered.esetALL.abundance.FC)) ] <- NA

#craete additional column symbabun with abundance or gene symbol 
#for use later to aggregate rows for network visualization
fData(filtered.esetALL.abundance.FC)$symbabun<-fData(filtered.esetALL.abundance.FC)$SYMBOL
fData(filtered.esetALL.abundance.FC)$symbabun[grep('CC',fData(filtered.esetALL.abundance.FC)$dataType)]<-fData(filtered.esetALL.abundance.FC)$displayName[grep('CC',fData(filtered.esetALL.abundance.FC)$dataType)]
fData(filtered.esetALL.abundance.FC)$displayName.final<-as.character(fData(filtered.esetALL.abundance.FC)$displayName)
fData(filtered.esetALL.abundance.FC)$displayName.final <- ifelse(is.na(fData(filtered.esetALL.abundance.FC)$displayName.final), fData(filtered.esetALL.abundance.FC)$symbabun, fData(filtered.esetALL.abundance.FC)$displayName.final)

filtered.esetALL.abundance.FC <- filtered.esetALL.abundance.FC[,!pData(filtered.esetALL.abundance.FC)$Patient.code %in% "37"]

dat_abundance <- data.frame(exprs(filtered.esetALL.abundance.FC), stringsAsFactors = FALSE)
dat_abundance$displayName.final <- fData(filtered.esetALL.abundance.FC)[match(row.names(dat_abundance),fData(filtered.esetALL.abundance.FC)$featureID),'displayName.final']
dat_abundance$clust <- fData(filtered.esetALL.abundance.FC)[match(row.names(dat_abundance),fData(filtered.esetALL.abundance.FC)$featureID),'clusterID']
dat_abundance[dat_abundance$displayName.final %in% c("CyTOF-PBMC","CyTOF-GRAN"),'clust'] <- 1080000

#aggregate cells by the major canonical cellular frequencies  
dat_abundance.merged <-  dat_abundance %>% 
  group_by(displayName.final) %>% 
  dplyr::slice(which.max(clust))

dat_abundance.merged <- data.frame(dat_abundance.merged, stringsAsFactors = F)
row.names(dat_abundance.merged) <- dat_abundance.merged$displayName.final

dat_abundance.merged.for.pca <- dat_abundance.merged[,!colnames(dat_abundance.merged) %in% c("displayName.final", "clust")]
dat_abundance.merged.for.pca <- dat_abundance.merged.for.pca[rowSums(dat_abundance.merged.for.pca)!=0,]

dat_abundance.merged.for.pca <- dat_abundance.merged.for.pca[!row.names(dat_abundance.merged.for.pca) %in% c("Granulocytes-Abundance", "T cells-Abundance"),]
dat_abundance.merged.for.pca.sub <- dat_abundance.merged.for.pca[!row.names(dat_abundance.merged.for.pca) %in% c("CyTOF-PBMC"),]
row.names(dat_abundance.merged.for.pca.sub)[row.names(dat_abundance.merged.for.pca.sub) %in% 'CyTOF-GRAN'] <- "Granulocyte-Abundance"


all.dataType.pca.Abundance <- prcomp(t(dat_abundance.merged.for.pca.sub),
                                     center = TRUE,
                                     scale. = TRUE) 

groups <- as.factor(dat_abundance$group)

pc_dat_abundance<- data.frame(PC1 = all.dataType.pca.Abundance$x[,1], PC2= all.dataType.pca.Abundance$x[,2], group=as.factor(pData(filtered.esetALL.abundance.FC)$group), patient=as.factor(pData(filtered.esetALL.abundance.FC)$Patient.code), visit=as.factor(pData(filtered.esetALL.abundance.FC)$Visit))

pc_dat_abundance$group <- ifelse(as.character(pc_dat_abundance$group) %in% "IM", "NR", as.character(pc_dat_abundance$group))

percentage <- round(all.dataType.pca.Abundance$sdev^2 / sum(all.dataType.pca.Abundance$sdev^2) * 100, 2)

pc_dat_abundance$class <- paste(pc_dat_abundance$group, pc_dat_abundance$visit, sep=".")

centroids <- aggregate(cbind(PC1,PC2)~group,pc_dat_abundance,mean)
centroids$label <- 1

ggplot(pc_dat_abundance,aes(x=PC1, y=PC2, col=group)) + 
  geom_point() + 
  geom_path(aes(group=patient), arrow=arrow(length = unit(0.3, "cm"), type="closed"))  +
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border =element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC1 (", percentage[1], "%)", sep="")) +ylab(paste("PC2 (", percentage[2], "%)", sep=""))+
  stat_ellipse(geom = "polygon", type="euclid", aes(group=group, fill = group), alpha = 0.2)+
  geom_point(data=centroids, size=5, shape=17)


#Compute PERMANOVA for differemce in composition
dat_abundance.merged.for.pca.sub_t <- data.frame(t(dat_abundance.merged.for.pca.sub), stringsAsFactors = F)

dat_abundance.merged.for.pca.sub_t$group <- pData(filtered.esetALL)[match(gsub("\\.", "-", row.names(dat_abundance.merged.for.pca.sub_t)), pData(filtered.esetALL)$sampleID),'group']
dat_abundance.merged.for.pca.sub_t$group <- ifelse(dat_abundance.merged.for.pca.sub_t$group %in% "IM", "NR", dat_abundance.merged.for.pca.sub_t$group)
dat_abundance.merged.for.pca.sub_t$Visit <- pData(filtered.esetALL)[match(gsub("\\.", "-", row.names(dat_abundance.merged.for.pca.sub_t)), pData(filtered.esetALL)$sampleID),'Visit']

res.permanova <- adonis(dat_abundance.merged.for.pca.sub_t[dat_abundance.merged.for.pca.sub_t$Visit %in% "V2",!colnames(dat_abundance.merged.for.pca.sub_t) %in% c("Visit", "group") ]~dat_abundance.merged.for.pca.sub_t[dat_abundance.merged.for.pca.sub_t$Visit %in% "V2",'group'],  method="euclidean", perm=10000) 
p.value <- res.permanova$aov.tab$`Pr(>F)`[1]

res.permanova <- adonis(dat_abundance.merged.for.pca.sub_t[,!colnames(dat_abundance.merged.for.pca.sub_t) %in% c("Visit", "group") ]~dat_abundance.merged.for.pca.sub_t[,'group'],  method="euclidean", perm=10000) 
p.value <- res.permanova$aov.tab$`Pr(>F)`[1]


#PCA loadings (PC2)
#Plot Supp. Fig 1


loadings.PC2 <- data.frame(all.dataType.pca.Abundance$rotation[, colnames(all.dataType.pca.Abundance$rotation) %in% "PC2"], stringsAsFactors = F)
loadings.PC2$cell <- row.names(loadings.PC2)
names(loadings.PC2) <- c("PC2", "Cell")
loadings.PC2 <- loadings.PC2[order(loadings.PC2$PC2, decreasing=T),]
loadings.PC2$Cell <- factor(loadings.PC2$Cell, levels=c(loadings.PC2$Cell))

ggplot(data=loadings.PC2, aes(x=Cell, y=PC2, fill=muted("blue"))) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c(muted("blue")))+
  theme_bw()+
  theme(strip.text = element_text(size=14), panel.grid.major = element_blank(),
        axis.text.y=element_text( size=12),
        axis.text.x=element_text( size=12, angle =90),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border =element_blank(), 
        axis.ticks.length=unit(.15, "cm"))
  
#Test correlation of cell frequency fold-changes with CRP-fold changes
#Plot Supp. Fig 1

dat.merged.abundance.final.t <- data.frame(t(dat_abundance.merged.for.pca.sub), stringsAsFactors = F)
dat.merged.abundance.final.t.V2 <- dat.merged.abundance.final.t[grepl("V2", row.names(dat.merged.abundance.final.t)),]
dat.merged.abundance.final.t.V2$Patient.code <- pData(filtered.esetALL)[match(row.names(dat.merged.abundance.final.t.V2), gsub("-", "\\.", pData(filtered.esetALL)$sampleID)),'Patient.code'] 

dat.merged.abundance.final.t.V2.res <- merge(dat.merged.abundance.final.t.V2, CRP.FC.df, by="Patient.code")
M2 <- dat.merged.abundance.final.t.V2.res[,!colnames(dat.merged.abundance.final.t.V2.res) %in% c("CRP.FCV1V2", "CRP.FCV1V3", "CRP.FCV2V3", "Patient.code")]
M2 <- apply(M2,2,as.numeric)
M<-corr.test(dat.merged.abundance.final.t.V2.res$CRP.FCV1V2, M2, method="spearman")

dat.merged.abundance.final.t.V3 <- dat.merged.abundance.final.t[grepl("V3", row.names(dat.merged.abundance.final.t)),]
dat.merged.abundance.final.t.V3$Patient.code <- pData(filtered.esetALL)[match(row.names(dat.merged.abundance.final.t.V3), gsub("-", "\\.", pData(filtered.esetALL)$sampleID)),'Patient.code'] 
dat.merged.abundance.final.t.V3.res <- merge(dat.merged.abundance.final.t.V3, CRP.FC.df, by="Patient.code")
M3 <- dat.merged.abundance.final.t.V3.res[,!colnames(dat.merged.abundance.final.t.V3.res) %in% c("CRP.FCV1V2", "CRP.FCV1V3", "CRP.FCV2V3", "Patient.code")]
M3 <- apply(M3,2,as.numeric)

M<-corr.test(c(dat.merged.abundance.final.t.V2.res$CRP.FCV1V2, dat.merged.abundance.final.t.V3.res$CRP.FCV1V3), rbind(M2,M3), method="spearman", adjust="none")


M.long <- data.frame(r=t(M$r), p=t(M$p), stringsAsFactors = F)
M.long$cell <- row.names(M.long)
M.long$label <- ifelse(M.long$p<0.05, "*", "") 

M.long <- M.long[order(M.long$r),]
M.long$cell <- factor(M.long$cell, levels=M.long$cell)

ggplot(M.long,aes(x=cell, y=r, fill=r)) + 
  coord_flip()+
  ylab("correlation of cell proportion and CRP fold changes")+
  geom_bar(stat = 'identity') + 
  scale_colour_gradient2(low = muted("red"), mid = "grey",high = muted("blue"))+
  scale_fill_gradient2(low = muted("red"), mid = "grey",high = muted("blue"))+
  geom_text(data = M.long, label = M.long$label, y= ifelse(M.long$r>0, M.long$r+0.03, M.long$r-0.03), size=6)+
  theme_bw()+
  theme(strip.text = element_text(size=14), panel.grid.major = element_blank(),
        axis.text.y=element_text( size=12),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border =element_blank(), 
        axis.ticks.length=unit(.15, "cm"))


#Show specific scatterplot for 
#monocytes correlation with CRP
#Plot Fig 2a right

dat.merged.abundance.final.t.V2.res.sub <- dat.merged.abundance.final.t.V2.res[,colnames(dat.merged.abundance.final.t.V2.res) %in% c("Monocytes.Abundance", "CRP.FCV1V2", "Patient.code")]
colnames(dat.merged.abundance.final.t.V2.res.sub) <- c("Patient.code", "Monocytes.Abundance", "CRP.FC")
dat.merged.abundance.final.t.V2.res.sub$Visit <- "V1V2"

dat.merged.abundance.final.t.V3.res.sub <- dat.merged.abundance.final.t.V3.res[,colnames(dat.merged.abundance.final.t.V3.res) %in% c("Monocytes.Abundance", "CRP.FCV1V3", "Patient.code")]
colnames(dat.merged.abundance.final.t.V3.res.sub) <- c("Patient.code", "Monocytes.Abundance", "CRP.FC")
dat.merged.abundance.final.t.V3.res.sub$Visit <- "V1V3"


dat.merged.abundance.final.t.V2.V3.res.sub <- rbind(dat.merged.abundance.final.t.V2.res.sub, dat.merged.abundance.final.t.V3.res.sub)
dat.merged.abundance.final.t.V2.V3.res.sub$group <- pData(filtered.esetALL)[match(dat.merged.abundance.final.t.V2.V3.res.sub$Patient.code, pData(filtered.esetALL)$Patient.code),'group']
dat.merged.abundance.final.t.V2.V3.res.sub$group <- ifelse(dat.merged.abundance.final.t.V2.V3.res.sub$group %in% "IM", "NR", dat.merged.abundance.final.t.V2.V3.res.sub$group)
dat.merged.abundance.final.t.V2.V3.res.sub$class <- paste(dat.merged.abundance.final.t.V2.V3.res.sub$group, dat.merged.abundance.final.t.V2.V3.res.sub$Visit, sep=".")
dat.merged.abundance.final.t.V2.V3.res.sub$Monocytes.Abundance.log <- log10(dat.merged.abundance.final.t.V2.V3.res.sub$Monocytes.Abundance)


ggplot(dat.merged.abundance.final.t.V2.V3.res.sub, aes(x=Monocytes.Abundance.log, y=CRP.FC)) + 
  ylab("CRP FC")+
  xlab("Monocyte abundance FC (log10)")+
  geom_point(aes(color = class), size = 4) + 
  geom_point(shape = 1, color = "black", size = 4) +
  scale_colour_manual(values=c("R.V1V2"="cyan4","R.V1V3"="cyan2", "NR.V1V2"="#F8766D","NR.V1V3"="rosybrown1"))+
  theme_bw()+
  theme(strip.text = element_text(size=14), panel.grid.major = element_blank(),
        axis.text.y=element_text( size=12),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border =element_blank(), 
        axis.ticks.length=unit(.15, "cm"))+
  stat_smooth(method = "lm", fullrange = TRUE, color="black", se=F)


#Compute paired changes of absolute frequency within each response group
#Plot Supp. Fig 1a

fData(filtered.esetALL.abundance)$symbabun<-fData(filtered.esetALL.abundance)$SYMBOL
fData(filtered.esetALL.abundance)$symbabun[grep('CC',fData(filtered.esetALL.abundance)$dataType)]<-fData(filtered.esetALL.abundance)$displayName[grep('CC',fData(filtered.esetALL.abundance)$dataType)]
fData(filtered.esetALL.abundance)$displayName.final <- ifelse(is.na(fData(filtered.esetALL.abundance)$displayName.final), fData(filtered.esetALL.abundance)$symbabun, fData(filtered.esetALL.abundance)$displayName.final)
filtered.esetALL.abundance <- filtered.esetALL.abundance[,!pData(filtered.esetALL.abundance)$Patient.code %in% "37"]
dat_abundance.abs <- data.frame(exprs(filtered.esetALL.abundance), stringsAsFactors = FALSE)
dat_abundance.abs$displayName.final <- fData(filtered.esetALL.abundance)[match(row.names(dat_abundance.abs),fData(filtered.esetALL.abundance)$featureID),'displayName.final']
dat_abundance.abs$clust <- fData(filtered.esetALL.abundance)[match(row.names(dat_abundance.abs),fData(filtered.esetALL.abundance)$featureID),'clusterID']
dat_abundance.abs[dat_abundance.abs$displayName.final %in% c("CyTOF-PBMC","CyTOF-GRAN"),'clust'] <- 1080000

dat_abundance.abs.merged <-  dat_abundance.abs %>%
  group_by(displayName.final) %>%
  dplyr::slice(which.max(clust))
  
dat_abundance.abs.merged <- data.frame(dat_abundance.abs.merged, stringsAsFactors = F)
row.names(dat_abundance.abs.merged) <- dat_abundance.abs.merged$displayName.final
dat_abundance.merged.for.pca <- dat_abundance.merged[,!colnames(dat_abundance.merged) %in% c("displayName.final", "clust")]
dat_abundance.abs.merged.for.pca <- dat_abundance.abs.merged[,!colnames(dat_abundance.abs.merged) %in% c("displayName.final", "clust")]
dat_abundance.abs.merged.for.pca <- dat_abundance.abs.merged.for.pca[rowSums(dat_abundance.abs.merged.for.pca)!=0,]  

dat_abundance.abs.merged.for.pca <- dat_abundance.abs.merged.for.pca[!row.names(dat_abundance.abs.merged.for.pca) %in% c("Granulocytes-Abundance", "T cells-Abundance"),]
dat_abundance.abs.merged.for.pca.sub <- dat_abundance.abs.merged.for.pca[!row.names(dat_abundance.abs.merged.for.pca) %in% c("CyTOF-PBMC"),]

dat_abundance.abs.merged.for.pca[row.names(dat_abundance.abs.merged.for.pca) %in% "CyTOF-GRAN",] <- dat_abundance.abs.merged.for.pca[row.names(dat_abundance.abs.merged.for.pca) %in% "CyTOF-GRAN",]/100
dat_abundance.abs.merged.for.pca[row.names(dat_abundance.abs.merged.for.pca) %in% "CyTOF-PBMC",] <- dat_abundance.abs.merged.for.pca[row.names(dat_abundance.abs.merged.for.pca) %in% "CyTOF-PBMC",]/100

dat_abundance.abs.merged.for.pca$cell <- row.names(dat_abundance.abs.merged.for.pca)
dat_abundance.merged.abs.long <- reshape2::melt(dat_abundance.abs.merged.for.pca[,!colnames(dat_abundance.abs.merged.for.pca) %in% "clust"], id.vars = c("cell"))

dat_abundance.merged.abs.long$variable <- gsub("\\.", "-",dat_abundance.merged.abs.long$variable )
dat_abundance.merged.abs.long$group <- pData(filtered.esetALL)[match(dat_abundance.merged.abs.long$variable, pData(filtered.esetALL)$sampleID),'group']
dat_abundance.merged.abs.long$group <- ifelse(dat_abundance.merged.abs.long$group %in% "IM", "NR", dat_abundance.merged.abs.long$group)
dat_abundance.merged.abs.long$Visit <- pData(filtered.esetALL)[match(dat_abundance.merged.abs.long$variable, pData(filtered.esetALL)$sampleID),'Visit']
dat_abundance.merged.abs.long$class <- paste(dat_abundance.merged.abs.long$group, dat_abundance.merged.abs.long$Visit, sep=".")

dat_abundance.merged.abs.long$class <- factor(dat_abundance.merged.abs.long$class, levels=c("R.V1", "R.V2", "R.V3", "NR.V1", "NR.V2", "NR.V3"))

dat_abundance.merged.abs.long$cell <- factor(dat_abundance.merged.abs.long$cell, levels=c("CyTOF-GRAN", "CyTOF-PBMC", "Monocytes-Abundance", "Tregs-Abundance", "Na?ve CD4 Tcell-Abundance", "CD4+central memory Tcells-Abundance", "gdT cells (CD8)-Abundance", "Na?ve CD8 Tcell-Abundance", "CD8+ effector memory T cells-Abundance", "gdT cells (CD4)-Abundance", "B cells -Abundance", "Anti-inflammatory monocytes-Abundance", "Plasmacytoid dendritic cells-Abundance", "NK cells-Abundance", "CD4+ effector memory T cells-Abundance", "CD4+ T cells-Abundance", "CD8+ T cells-Abundance", "CD8+ effector T cells-Abundance"))

dat_abundance.merged.abs.long$class2 <- ifelse(dat_abundance.merged.abs.long$Visit %in% "V1", paste(dat_abundance.merged.abs.long$group, "D0", sep="."), ifelse(dat_abundance.merged.abs.long$Visit %in% "V2", paste(dat_abundance.merged.abs.long$group, "W2", sep="."), paste(dat_abundance.merged.abs.long$group, "W14", sep=".")))

dat_abundance.merged.abs.long$class2 <- factor(dat_abundance.merged.abs.long$class2, levels=c("R.D0", "R.W2", "R.W14", "NR.D0", "NR.W2", "NR.W14"))
ggplot(dat_abundance.merged.abs.long, aes(x = class2, y = as.numeric(value), fill=group, colour="black")) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  scale_color_manual(values="black")+
  facet_wrap(~cell, scale="free", ncol=3)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  ylab("Cell frequency")+
  theme_bw()+
  theme(strip.text = element_text(size=12), panel.grid.major = element_blank(),
        axis.text.y=element_text( size=10),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border =element_blank(), 
        axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox",size=3.5, paired=T, comparisons=list(c(1,2),c(1,3), c(4,5), c(4,6)))

#Plot boxplot of monocytes abundance changes
#Fig 2a, center

ggplot(dat_abundance.merged.abs.long[dat_abundance.merged.abs.long$cell %in% "Monocytes-Abundance",], aes(x = class, y = as.numeric(value), fill=group, colour="black")) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  scale_color_manual(values="black")+
  facet_wrap(~cell, scale="free", ncol=5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size=7))+
  ylab("Cell frequency")+
  theme_bw()+
  theme(strip.text = element_text(size=12), panel.grid.major = element_blank(),
        axis.text.y=element_text( size=10),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border =element_blank(), 
        axis.ticks.length=unit(.15, "cm"))+
  stat_compare_means(method = "wilcox",size=3.5, paired=T, comparisons=list(c(1,2),c(1,3), c(4,5), c(4,6)))


#######Identify significant edges based on two reference networks of Fold change in responders######

filtered.esetALL.FC.no.abundance  <- substractReference(filtered.esetALL[!grepl("Abundance", fData(filtered.esetALL)$displayName),], 'Patient.code', 'Visit', ref = 'V1', op = '-')

#combine fold change values of abundance with log2.FC of expression 

combined.expression.set <- rbind.fill(data.frame(exprs(filtered.esetALL.FC.no.abundance), stringsAsFactors = F), data.frame(exprs(filtered.esetALL.abundance.FC), stringsAsFactors = F))
row.names(combined.expression.set) <- c(row.names(exprs(filtered.esetALL.FC.no.abundance)), row.names(exprs(filtered.esetALL.abundance.FC)))
colnames(combined.expression.set) <- gsub("\\.", "-", colnames(combined.expression.set))
combined.fData <- rbind.fill(fData(filtered.esetALL.FC.no.abundance), fData(filtered.esetALL.abundance.FC))
row.names(combined.fData) <- combined.fData$featureID

combined.expression.set.mat <- as.matrix(combined.expression.set)
row.names(combined.expression.set.mat) <- row.names(combined.expression.set)

filtered.esetALL.FC.final <- ExpressionSet(combined.expression.set.mat, phenoData = AnnotatedDataFrame(pData(filtered.esetALL.FC.no.abundance))
                                           , featureData = AnnotatedDataFrame(combined.fData)
                                           , annotation = 'ClariomS.db')

filtered.esetALL.FC.final.R <- filtered.esetALL.FC.final[,pData(filtered.esetALL.FC.final)$group %in% "R"]


#Construct co-expression matrix based on Spearman's r

#Define significant edges

cor.edges.df.FC.V1V2 <- sig.edge.func(filtered.esetALL.FC.final.R, "V2", changesV2.final)
cor.edges.df.FC.V1V3 <- sig.edge.func(filtered.esetALL.FC.final.R,"V3", changesV3.final) 

cor.edges.df.FC.V1V2$Visit <- "FCV1V2"
cor.edges.df.FC.V1V3$Visit <- "FCV1V3"

cor.edges.df.FC.combined <- rbind(cor.edges.df.FC.V1V2, cor.edges.df.FC.V1V3)

cor.edges.df.FC.V1V2$Var1displayName <- fData(filtered.esetALL.FC.final)[match(as.character(cor.edges.df.FC.V1V2$Var1), as.character(fData(filtered.esetALL.FC.final)$featureID)),'displayName.final'] 
cor.edges.df.FC.V1V2$Var2displayName <- fData(filtered.esetALL.FC.final)[match(as.character(cor.edges.df.FC.V1V2$Var2), as.character(fData(filtered.esetALL.FC.final)$featureID)),'displayName.final'] 

cor.edges.df.FC.V1V3$Var1displayName <- fData(filtered.esetALL.FC.final)[match(as.character(cor.edges.df.FC.V1V3$Var1), as.character(fData(filtered.esetALL.FC.final)$featureID)),'displayName.final'] 
cor.edges.df.FC.V1V3$Var2displayName <- fData(filtered.esetALL.FC.final)[match(as.character(cor.edges.df.FC.V1V3$Var2), as.character(fData(filtered.esetALL.FC.final)$featureID)),'displayName.final'] 

#group feature edges by name
cor.edges.df.FC.V1V2.aggregated <- cor.edges.df.FC.V1V2 %>%
  group_by(Var1displayName, Var2displayName) %>%
  dplyr::slice(which.max(cor))


cor.edges.df.FC.V1V3.aggregated <- cor.edges.df.FC.V1V3 %>%
  group_by(Var1displayName, Var2displayName) %>%
  dplyr::slice(which.max(cor))

cor.edges.df.FC.V1V2.aggregated <- as.data.frame(cor.edges.df.FC.V1V2.aggregated, stringsAsFactors = F)
cor.edges.df.FC.V1V3.aggregated <- as.data.frame(cor.edges.df.FC.V1V3.aggregated, stringsAsFactors = F)


########Calculate 'Disruption Network'##############
#use responders' V1V2 FC network as a backbone
#calculate NR-based drops for each edge in each sample
#use responders' drops as NULL distribution values to determine significance
#use cell centered expression for genes

filtered.esetALL.FC.final.R.changes.V1V2 <- filtered.esetALL.FC.final.R[fData(filtered.esetALL.FC.final.R)$featureID %in% changesV2.final & !fData(filtered.esetALL.FC.final.R)$dataType %in% "GX" ,pData(filtered.esetALL.FC.final.R)$Visit %in% "V2"]
filtered.esetALL.FC.final.NR.changes.V1V2 <- filtered.esetALL.FC.final.NR[fData(filtered.esetALL.FC.final.NR)$featureID %in% changesV2.final & !fData(filtered.esetALL.FC.final.NR)$dataType %in% "GX",pData(filtered.esetALL.FC.final.NR)$Visit %in% "V2"]

sig.features.dataType <- filtered.esetALL.FC.final.R.changes.V1V2[fData(filtered.esetALL.FC.final.R.changes.V1V2)$featureID %in% changesV2.final & !fData(filtered.esetALL.FC.final.R.changes.V1V2)$dataType %in% "GX",]
dataType<-unique(fData(sig.features.dataType)$dataType)

cellTypes.V2 =unique(subset(fData(sig.features.dataType), grepl('Abundance|CyTOF-PBMC|CyTOF-GRAN', fData(sig.features.dataType)$displayName))$featureID)
cellTypes.V2<-cellTypes.V2[!is.na(cellTypes.V2)]

filtered.esetALL.FC.final.R.changes.V1V3 <- filtered.esetALL.FC.final.R[fData(filtered.esetALL.FC.final.R)$featureID %in% changesV3.final & !fData(filtered.esetALL.FC.final.R)$dataType %in% "GX" ,pData(filtered.esetALL.FC.final.R)$Visit %in% "V3"]
filtered.esetALL.FC.final.NR.changes.V1V3 <- filtered.esetALL.FC.final.NR[fData(filtered.esetALL.FC.final.NR)$featureID %in% changesV3.final & !fData(filtered.esetALL.FC.final.NR)$dataType %in% "GX",pData(filtered.esetALL.FC.final.NR)$Visit %in% "V3"]

sig.features.dataType <- filtered.esetALL.FC.final.R.changes.V1V3[fData(filtered.esetALL.FC.final.R.changes.V1V3)$featureID %in% changesV3.final & !fData(filtered.esetALL.FC.final.R.changes.V1V3)$dataType %in% "GX",]
dataType<-unique(fData(sig.features.dataType)$dataType)

cellTypes.V3 =unique(subset(fData(sig.features.dataType), grepl('Abundance|CyTOF-PBMC|CyTOF-GRAN', fData(sig.features.dataType)$displayName))$featureID)
cellTypes.V3<-cellTypes.V3[!is.na(cellTypes.V3)]

#Generated empirical null distribution of dropouts 
#('normal response' dropouts)

#calculate normal drops of V2FC

cor.responders.drop.fun(changesV2.final, "V2", filtered.esetALL.FC.final.R.changes.V1V2, cellTypes.V2)

list.rds.files <- paste("cor.responders.drop", "V2", seq(length(unique(pData(filtered.esetALL.FC.final.R.changes.V1V2)$Patient.code))), "rds", sep=".")
df.l.responders.drop.CI.FC.V2<-do.call(cbind, lapply(list.rds.files, function(cur.file) {
  x <- readRDS(cur.file)
  drop.dir <- x[, colnames(x) %in% "drop"]*x[,colnames(x) %in% "dir"]
  rm(x)
  gc()
  return(drop.dir)
}))

#add the type column

df.l.responders.drop.CI.FC.V2<-data.frame(df.l.responders.drop.CI.FC.V2, stringsAsFactors=FALSE)
df.l.responders.drop.CI.FC.V2.1 <- readRDS(paste("cor.responders.drop", "V2", "1", "rds", sep="."))
df.l.responders.drop.CI.FC.V2$Var1<-df.l.responders.drop.CI.FC.V2.1[,'Var1']
df.l.responders.drop.CI.FC.V2$Var2<-df.l.responders.drop.CI.FC.V2.1[,'Var2']
df.l.responders.drop.CI.FC.V2$type<-df.l.responders.drop.CI.FC.V2.1[,'type']
colnames(df.l.responders.drop.CI.FC.V2)[!colnames(df.l.responders.drop.CI.FC.V2) %in% c("Var1", "Var2", "type")] <- row.names(t(exprs(filtered.esetALL.FC.final.R[,pData(filtered.esetALL.FC.final.R)$Visit %in% "V2"])))
df.l.responders.drop.CI.FC.V2[is.na(df.l.responders.drop.CI.FC.V2)] <- 0

#calculate normal drops of V3FC

cor.responders.drop.fun(changesV3.final, "V3", filtered.esetALL.FC.final.R.changes.V1V3, cellTypes.V3)
list.rds.files <- paste("cor.responders.drop", "V3", seq(length(unique(pData(filtered.esetALL.FC.final.R.changes.V1V2)$Patient.code))), "rds", sep=".")
df.l.responders.drop.CI.FC.V3<-do.call(cbind, lapply(list.rds.files, function(cur.file) {
  x <- readRDS(cur.file)
  drop.dir <- x[, colnames(x) %in% "drop"]*x[,colnames(x) %in% "dir"]
  rm(x)
  gc()
  return(drop.dir)
}))

df.l.responders.drop.CI.FC.V3<-data.frame(df.l.responders.drop.CI.FC.V3, stringsAsFactors=FALSE)
df.l.responders.drop.CI.FC.V3.1 <- readRDS(paste("cor.responders.drop", "V3", "1", "rds", sep="."))
df.l.responders.drop.CI.FC.V3$Var1<-df.l.responders.drop.CI.FC.V3.1[,'Var1']
df.l.responders.drop.CI.FC.V3$Var2<-df.l.responders.drop.CI.FC.V3.1[,'Var2']
df.l.responders.drop.CI.FC.V3$type<-df.l.responders.drop.CI.FC.V3.1[,'type']

colnames(df.l.responders.drop.CI.FC.V3)[!colnames(df.l.responders.drop.CI.FC.V3) %in% c("Var1", "Var2", "type")] <- row.names(t(exprs(filtered.esetALL.FC.final.R[,pData(filtered.esetALL.FC.final.R)$Visit %in% "V3"])))
df.l.responders.drop.CI.FC.V3[is.na(df.l.responders.drop.CI.FC.V3)] <- 0

#Test outliers in responders drop matrix

df.l.responders.drop.CI.FC.V2$Var1 <- as.character(df.l.responders.drop.CI.FC.V2$Var1)
df.l.responders.drop.CI.FC.V2$Var2 <- as.character(df.l.responders.drop.CI.FC.V2$Var2)
cor.edges.df.FC.V1V2.aggregated$Var1 <- as.character(cor.edges.df.FC.V1V2.aggregated$Var1)
cor.edges.df.FC.V1V2.aggregated$Var2 <- as.character(cor.edges.df.FC.V1V2.aggregated$Var2)

df.l.responders.drop.CI.FC.sig.cor.V2 <- merge(cor.edges.df.FC.V1V2.aggregated[,colnames(cor.edges.df.FC.V1V2.aggregated) %in% c("Var1", "Var2")], df.l.responders.drop.CI.FC.V2, by=c("Var1", "Var2"), all.x=T)

df.l.responders.drop.CI.FC.V3$Var1 <- as.character(df.l.responders.drop.CI.FC.V3$Var1)
df.l.responders.drop.CI.FC.V3$Var2 <- as.character(df.l.responders.drop.CI.FC.V3$Var2)
cor.edges.df.FC.V1V3.aggregated$Var1 <- as.character(cor.edges.df.FC.V1V3.aggregated$Var1)
cor.edges.df.FC.V1V3.aggregated$Var2 <- as.character(cor.edges.df.FC.V1V3.aggregated$Var2)

df.l.responders.drop.CI.FC.sig.cor.V3 <- merge(cor.edges.df.FC.V1V3.aggregated[,colnames(cor.edges.df.FC.V1V3.aggregated) %in% c("Var1", "Var2")], df.l.responders.drop.CI.FC.V3, by=c("Var1", "Var2"), all.x=T)


#For V2
inf.vec <- apply(df.l.responders.drop.CI.FC.sig.cor.V2[,!colnames(df.l.responders.drop.CI.FC.sig.cor.V2) %in% c("Var1", "Var2", "type")], 1, function(x) any(is.infinite(x)))
df.l.responders.drop.CI.sig.cor.for.pca <- df.l.responders.drop.CI.FC.sig.cor.V2[!inf.vec,]

var.vec <- apply(df.l.responders.drop.CI.sig.cor.for.pca[,!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")], 1, function(x) var(x))
df.l.responders.drop.CI.sig.cor.for.pca <- df.l.responders.drop.CI.sig.cor.for.pca[!var.vec==0,]

normal.drop.pca<- prcomp(t(df.l.responders.drop.CI.sig.cor.for.pca[,!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")]),
                         center = T,
                         scale. = T) 


pc_dat.normal.drop<- data.frame(PC1 = normal.drop.pca$x[,1], PC2= normal.drop.pca$x[,2], PC3= normal.drop.pca$x[,3],
                                group=pData(filtered.esetALL.FC.final)[match(colnames(df.l.responders.drop.CI.sig.cor.for.pca)[!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")], pData(filtered.esetALL.FC.final)$sampleID),'group'],
                                stringsAsFactors = F)
pc_dat.normal.drop$label <- row.names(pc_dat.normal.drop)

percentage <- round(normal.drop.pca$sdev^2 / sum(normal.drop.pca$sdev^2) * 100, 2)

ggplot(pc_dat.normal.drop,aes(x=PC1, y=PC2, col=group, lable=label)) + 
  geom_point(size=3) + 
  scale_color_manual(values=c("cyan3"))+
  geom_text(aes(label=label),hjust=0.5, vjust=-1, size=5)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC1 (", percentage[1], "%) ", sep="")) +ylab(paste("PC2 (", percentage[2], "%) ", sep=""))

ggplot(pc_dat.normal.drop,aes(x=PC2, y=PC3, col=group, lable=label)) + 
  geom_point(size=3) + 
  scale_color_manual(values=c("cyan3"))+
  geom_text(aes(label=label),hjust=0.5, vjust=-1, size=5)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC2 (", percentage[2], "%) ", sep="")) +ylab(paste("PC3 (", percentage[3], "%) ", sep=""))


anno.df = data.frame(group = pData(filtered.esetALL.FC.final)[match(colnames(df.l.responders.drop.CI.sig.cor.for.pca[,!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")]), pData(filtered.esetALL.FC.final)$sampleID),'group'], stringsAsFactors = F)
anno.df$group <- as.character(anno.df$group)

ha <-   HeatmapAnnotation(df = anno.df,
                          col = list(group =c("R" = "cyan3")))


Ht1 <- Heatmap(df.l.responders.drop.CI.sig.cor.for.pca[,!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")], name="Normal drop FC.V2", 
        row_names_gp = gpar(fontsize = 5), show_row_names = FALSE, show_column_names = TRUE, top_annotation = ha)


#For V3
inf.vec <- apply(df.l.responders.drop.CI.FC.sig.cor.V3[,!colnames(df.l.responders.drop.CI.FC.sig.cor.V3) %in% c("Var1", "Var2", "type")], 1, function(x) any(is.infinite(x)))
df.l.responders.drop.CI.sig.cor.for.pca <- df.l.responders.drop.CI.FC.sig.cor.V3[!inf.vec,]

var.vec <- apply(df.l.responders.drop.CI.sig.cor.for.pca[,!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")], 1, function(x) var(x))
df.l.responders.drop.CI.sig.cor.for.pca <- df.l.responders.drop.CI.sig.cor.for.pca[!var.vec==0,]

normal.drop.pca<- prcomp(t(df.l.responders.drop.CI.sig.cor.for.pca[,!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")]),
                         center = T,
                         scale. = T) 


pc_dat.normal.drop<- data.frame(PC1 = normal.drop.pca$x[,1], PC2= normal.drop.pca$x[,2], PC3= normal.drop.pca$x[,3],
                                #group=as.factor(ifelse(as.character(pData(filtered.combined.eset.updated.adj.FC)$final.group) %in% "IM", "Non-responder", pData(filtered.combined.eset.updated.adj.FC)$final.group)), 
                                group=pData(filtered.esetALL.FC.final)[match(colnames(df.l.responders.drop.CI.sig.cor.for.pca)[!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")], pData(filtered.esetALL.FC.final)$sampleID),'group'],
                                stringsAsFactors = F)
pc_dat.normal.drop$label <- row.names(pc_dat.normal.drop)

percentage <- round(normal.drop.pca$sdev^2 / sum(normal.drop.pca$sdev^2) * 100, 2)

ggplot(pc_dat.normal.drop,aes(x=PC1, y=PC2, col=group, lable=label)) + 
  geom_point(size=3) + 
  scale_color_manual(values=c("cyan3"))+
  geom_text(aes(label=label),hjust=0.5, vjust=1, size=5)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC1 (", percentage[1], "%) ", sep="")) +ylab(paste("PC2 (", percentage[2], "%) ", sep=""))

ggplot(pc_dat.normal.drop,aes(x=PC2, y=PC3, col=group, lable=label)) + 
  geom_point(size=3) + 
  scale_color_manual(values=c("cyan3"))+
  geom_text(aes(label=label),hjust=0.5, vjust=1, size=5)+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  xlab(paste("PC2 (", percentage[2], "%) ", sep="")) +ylab(paste("PC3 (", percentage[3], "%) ", sep=""))


anno.df = data.frame(group = pData(filtered.esetALL.FC.final)[match(colnames(df.l.responders.drop.CI.sig.cor.for.pca[,!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")]), pData(filtered.esetALL.FC.final)$sampleID),'group'], stringsAsFactors = F)
anno.df$group <- as.character(anno.df$group)

ha <-   HeatmapAnnotation(df = anno.df,
                          col = list(group =c("R" = "cyan3")))


Ht2 <- Heatmap(df.l.responders.drop.CI.sig.cor.for.pca[,!colnames(df.l.responders.drop.CI.sig.cor.for.pca) %in% c("Var1", "Var2", "type")], name="Normal drop FC.V3", 
        row_names_gp = gpar(fontsize = 5), show_row_names = FALSE, show_column_names = TRUE, top_annotation = ha)


#Calculate non-responders drops

#add each iteration one NUR and recalculate correlation matrix. 
#calculate the r drop for each pair of features as a result of the NUR patient addiation. 

#combine responders and non responders for calculating drop

#now the whole responders correlation network is used as cc_n and cc_n1 is the correlation 
#after addition of one non-responder

cor.disruption.fun(changesV2.final, "V2", filtered.esetALL.FC.final.R.changes.V1V2, filtered.esetALL.FC.final.NR.changes.V1V2, "V2", "NR", cellTypes.V2)

list.rds.files <- paste("cor.NR.drop", "V2", seq(length(unique(pData(filtered.esetALL.FC.final.NR.changes.V1V2)$Patient.code))), "rds", sep=".")
df.l.FC.V2.NUR.final.to.sig.edges<-lapply(list.rds.files, function(cur.file) {
  x <- readRDS(cur.file)
  x$Var1 <- as.character(x$Var1)
  x$Var2 <- as.character(x$Var2)
  cor.edges.df.FC.V1V2.aggregated$Var1 <- as.character(cor.edges.df.FC.V1V2.aggregated$Var1)
  cor.edges.df.FC.V1V2.aggregated$Var2 <- as.character(cor.edges.df.FC.V1V2.aggregated$Var2)
  tmp <- merge(x, cor.edges.df.FC.V1V2.aggregated[,colnames(cor.edges.df.FC.V1V2.aggregated) %in% c("Var1", "Var2")], by.x=c("Var1", "Var2"), by.y=c("Var1", "Var2"), all.y=T)
  return(tmp)
})


cor.disruption.fun(changesV3.final, "V3", filtered.esetALL.FC.final.R.changes.V1V3, filtered.esetALL.FC.final.NR.changes.V1V3, "V3", "NR", cellTypes.V3)


list.rds.files <- paste("cor.NR.drop","V3", seq(length(unique(pData(filtered.esetALL.FC.final.NR.changes.V1V3)$Patient.code))), "rds", sep=".")
df.l.FC.V3.NUR.final.to.sig.edges<-lapply(list.rds.files, function(cur.file) {
  x <- readRDS(cur.file)
  x$Var1 <- as.character(x$Var1)
  x$Var2 <- as.character(x$Var2)
  cor.edges.df.FC.V1V3.aggregated$Var1 <- as.character(cor.edges.df.FC.V1V3.aggregated$Var1)
  cor.edges.df.FC.V1V3.aggregated$Var2 <- as.character(cor.edges.df.FC.V1V3.aggregated$Var2)
  tmp <- merge(x, cor.edges.df.FC.V1V3.aggregated[,colnames(cor.edges.df.FC.V1V3.aggregated) %in% c("Var1", "Var2")], by.x=c("Var1", "Var2"), by.y=c("Var1", "Var2"), all.y=T)
  return(tmp)
})

#Determine NR-drop significance for each edge in each sample by calculation of left-tail percentile, 
#within the null distribution of the normal dropouts for each edge

df.l.FC.V2.NR.percentile.final <- percentile.and.FDR.calc.function(df.l.FC.V2.NUR.final.to.sig.edges, df.l.responders.drop.CI.FC.sig.cor.V2)
df.l.FC.V3.NR.percentile.final <- percentile.and.FDR.calc.function(df.l.FC.V3.NUR.final.to.sig.edges, df.l.responders.drop.CI.FC.sig.cor.V3)

names(df.l.FC.V2.NR.percentile.final) <- row.names(t(exprs(filtered.esetALL.FC.final.NR[,pData(filtered.esetALL.FC.final.NR)$Visit %in% "V2"])))
names(df.l.FC.V3.NR.percentile.final) <- row.names(t(exprs(filtered.esetALL.FC.final.NR[,pData(filtered.esetALL.FC.final.NR)$Visit %in% "V3"])))

#Generate final drop matrix 

drop.matrix.aggregation.function <- function(df.l.NR.percentile.final, sig.threshold, abs.drop.intensity.threshold) {
  
  #calculate disruption of features by drop if they have disruption percentile<sig.threshold, direction of disruption<0, and correlation FDR<abs.drop.intensity.threshold; otherwize disruption=0
  df.l.NUR.sub.percentile.disruption<-sapply(seq(length(df.l.NR.percentile.final)), simplify=FALSE, function(m) {
    disruption<- ((df.l.NR.percentile.final[[m]]$disrup.FDR < sig.threshold) & (abs(df.l.NR.percentile.final[[m]]$drop.dir)>abs.drop.intensity.threshold) & (df.l.NR.percentile.final[[m]]$drop.dir<0))*df.l.NR.percentile.final[[m]]$drop.dir
    df.l.NR.percentile.final[[m]]$disruption <- disruption
    return(df.l.NR.percentile.final[[m]])
  })
  
  df.l.NUR.sub.percentile.disruption <- sapply(df.l.NUR.sub.percentile.disruption, simplify=FALSE, function(m) {
    sub <- m[,colnames(m) %in% c("Var1", "Var2", "disruption")]
    return(sub)
  })
  
  #generate disruption matrix
  
  drop_matrix.pre <- rbindlist(df.l.NUR.sub.percentile.disruption, idcol=TRUE)
  colnames(drop_matrix.pre)[1] <- "id"
  drop_matrix <- dcast(drop_matrix.pre, Var1 + Var2 ~ id, value.var="disruption")
  
  colnames(drop_matrix)[!colnames(drop_matrix) %in% c("Var1", "Var2")] <- names(df.l.NR.percentile.final)
  
  names <- sapply(nrow(drop_matrix), function (m) {paste(drop_matrix$Var1,"-", drop_matrix$Var2)})
  row.names(drop_matrix)<-make.names(names, unique = TRUE)
  
  drop_matrix[drop_matrix=="NA"] <- "0"
  drop_matrix <- data.frame(drop_matrix, stringsAsFactors=FALSE)
  
  drop_matrix$Var1 <- as.character(drop_matrix$Var1)
  drop_matrix$Var2 <- as.character(drop_matrix$Var2)
  
  return(drop_matrix)
}


drop_matrix.V2.NUR <- drop.matrix.aggregation.function(df.l.FC.V2.NR.percentile.final, 0.1, 0.4)
drop_matrix.V3.NUR <- drop.matrix.aggregation.function(df.l.FC.V3.NR.percentile.final, 0.1, 0.4)

#combine drop matrix

abs.drop.intensity.threshold <- 0.4
colnames(df.l.responders.drop.CI.FC.V2) <- gsub("-", "\\.", colnames(df.l.responders.drop.CI.FC.V2))
drop_matrix.R.NR.V2<- merge(drop_matrix.V2.NUR,df.l.responders.drop.CI.FC.V2[,!colnames(df.l.responders.drop.CI.FC.V2) %in% "type"], by=c("Var1", "Var2"))
drop_matrix.R.NR.V2[,!colnames(drop_matrix.R.NR.V2) %in% c("Var1", "Var2")][drop_matrix.R.NR.V2[,!colnames(drop_matrix.R.NR.V2) %in% c("Var1", "Var2")]>(-abs.drop.intensity.threshold)] <- 0


colnames(df.l.responders.drop.CI.FC.V3) <- gsub("-", "\\.", colnames(df.l.responders.drop.CI.FC.V3))
drop_matrix.R.NR.V3<- merge(drop_matrix.V3.NUR,df.l.responders.drop.CI.FC.V3[,!colnames(df.l.responders.drop.CI.FC.V3) %in% "type"], by=c("Var1", "Var2"))
drop_matrix.R.NR.V3[,!colnames(drop_matrix.R.NR.V3) %in% c("Var1", "Var2")][drop_matrix.R.NR.V3[,!colnames(drop_matrix.R.NR.V3) %in% c("Var1", "Var2")]>(-abs.drop.intensity.threshold)] <- 0


#Identify pathways that fits the network structure

#download gene-sets and combine them to one list
# 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/'

geneSets.with.Go <- do.call('append', list(readList("h.all.v7.0.symbols.gmt"), readList("c2.cp.v7.0.symbols.gmt")))
geneSets.with.Go <- do.call('append', list(geneSets.with.Go, readList("c5.bp.v7.0.symbols.gmt")))

nodes.in.all.pathways <- do.call('rbind', lapply(geneSets.with.Go, function(x) data.frame(x, stringsAsFactors = F)))
nodes.in.all.pathways$pathway <- as.character(row.names(nodes.in.all.pathways))
nodes.in.all.pathways$pathway <- gsub("[.][^.]+$", "", nodes.in.all.pathways$pathway)

#Determine pathway sizes
pathway.size <- data.frame(table(nodes.in.all.pathways$pathway), stringsAsFactors = F)

#filter out too large pathways (above 400 genes)

geneSets.with.Go.sub <- geneSets.with.Go[names(geneSets.with.Go) %in% as.character(pathway.size[pathway.size$Freq<400,'Var1'])]
nodes.in.all.pathways.sub <- nodes.in.all.pathways[nodes.in.all.pathways$pathway %in% as.character(pathway.size[pathway.size$Freq<400,'Var1']),]

#map each edge node to a pathway
drop_matrix.R.NR.V2.for.classification <- drop_matrix.R.NR.V2
drop_matrix.R.NR.V2.for.classification$Var1SYMBOL <- fData(filtered.esetALL)[match(drop_matrix.R.NR.V2.for.classification$Var1, fData(filtered.esetALL)$featureID),'SYMBOL']
drop_matrix.R.NR.V2.for.classification$Var2SYMBOL <- fData(filtered.esetALL)[match(drop_matrix.R.NR.V2.for.classification$Var2, fData(filtered.esetALL)$featureID),'SYMBOL']

drop_matrix.R.NR.V2.for.classification$Var1.class <- lapply(seq(length(drop_matrix.R.NR.V2.for.classification$Var1SYMBOL)), function(cur.symb) {
  nodes.in.all.pathways.sub[which(nodes.in.all.pathways.sub$x %in% drop_matrix.R.NR.V2.for.classification$Var1SYMBOL[cur.symb]),'pathway']
})

drop_matrix.R.NR.V2.for.classification$Var2.class <- lapply(seq(length(drop_matrix.R.NR.V2.for.classification$Var2SYMBOL)), function(cur.symb) {
  nodes.in.all.pathways.sub[which(nodes.in.all.pathways.sub$x %in% drop_matrix.R.NR.V2.for.classification$Var2SYMBOL[cur.symb]),'pathway']
})

#identify edges which both nodes are annotated, only one node annotated and both nodes not-annotated 
#In case of full annotation, return pathway intersection, i.e. the pathways with agreement of both nodes of  the edge
drop_matrix.R.NR.V2.for.classification$edge.annotation <- do.call('rbind', lapply(seq(nrow(drop_matrix.R.NR.V2.for.classification)), function(i) {
  ifelse(length(drop_matrix.R.NR.V2.for.classification$Var1.class[i][[1]])==0 &  length(drop_matrix.R.NR.V2.for.classification$Var2.class[i][[1]])==0, "both.not.annotated",
         ifelse(length(drop_matrix.R.NR.V2.for.classification$Var1.class[i][[1]])==0 & length(drop_matrix.R.NR.V2.for.classification$Var2.class[i][[1]])>0|length(drop_matrix.R.NR.V2.for.classification$Var1.class[i][[1]])>0 & length(drop_matrix.R.NR.V2.for.classification$Var2.class[i][[1]])==0, "only.one.annotated",
                list(intersect(drop_matrix.R.NR.V2.for.classification$Var1.class[i][[1]], drop_matrix.R.NR.V2.for.classification$Var2.class[i][[1]]))))
}))

#count edge by type (pathway when both nodes are annotated, edges of only one node annotated, and both.not.annotated)

res.summary.edge.classification <- data.frame(table(unlist(drop_matrix.R.NR.V2.for.classification$edge.annotation)), stringsAsFactors = F)
colnames(res.summary.edge.classification) <- c("pathway", "count")
geneSets.Go.sub.filtered.names <- res.summary.edge.classification[res.summary.edge.classification$count>5, 'pathway']
geneSets.Go.sub.filtered.names <- geneSets.Go.sub.filtered.names[!geneSets.Go.sub.filtered.names %in% c("only.one.annotated", "both.not.annotated")]

net.nodes.by.pathways <- do.call('rbind', lapply(geneSets.Go.sub.filtered.names, function(cur.path) {
  submat <- drop_matrix.R.NR.V2.for.classification[grepl(cur.path, drop_matrix.R.NR.V2.for.classification$edge.annotation),]
  node.list.in.module.and.net.edge <- unique(c(submat$Var1, submat$Var2))
  return(data.frame(path=cur.path, featureID=node.list.in.module.and.net.edge, stringsAsFactors = F))
})) 

net.nodes.by.pathways$SYMBOL <- fData(filtered.esetALL)[match(net.nodes.by.pathways$featureID, fData(filtered.esetALL)$featureID),'SYMBOL']

net.nodes.by.pathways.df <- data.frame(table(net.nodes.by.pathways$path), stringsAsFactors = F)
colnames(net.nodes.by.pathways.df) <- c("path", "size")

net.nodes.by.pathways <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% as.character(net.nodes.by.pathways.df[net.nodes.by.pathways.df$size>=5,'path']),]

#Supplementary table of network dropouts and pathway annotation 
drop_matrix.R.NR.V2.SuppTable <- drop_matrix.R.NR.V2.for.classification
drop_matrix.R.NR.V2.SuppTable$edgeAnnot <- do.call('rbind', lapply(drop_matrix.R.NR.V2.for.classification$edge.annotation, function(i) str_c(i,collapse='|') ))
drop_matrix.R.NR.V2.SuppTable$edge.annotation <- NULL
drop_matrix.R.NR.V2.SuppTable$Var1.class <- NULL
drop_matrix.R.NR.V2.SuppTable$Var2.class <- NULL

drop.matrix.with.pathways.V2 <- merge(cor.edges.df.FC.V1V2.aggregated, drop_matrix.R.NR.V2.SuppTable, by=c("Var1", "Var2"))

write.csv(drop.matrix.with.pathways.V2, "corr.net.responders.V2.SuppTable.csv")


#Normal anti-TNF response in responders
#Global enrichment using GSEA

R.V2.for.GSEA <- lmer[!lmer$dataTypes %in% "GX",]
R.V2.for.GSEA$SYMBOL <- fData(esetALL)[match(R.V2.for.GSEA$featureID, fData(esetALL)$featureID),'SYMBOL']
R.V2.for.GSEA <- R.V2.for.GSEA[,colnames(R.V2.for.GSEA) %in% c("SYMBOL", "V2.Estimate", "V3.Estimate", "V2.Pvalue", "V3.Pvalue")]

R.V2.for.GSEA.aggregated <- aggregate(R.V2.for.GSEA[,!colnames(R.V2.for.GSEA) %in% "SYMBOL"], list(SYMBOL= R.V2.for.GSEA$SYMBOL), mean)

R.V2.for.GSEA.aggregated$logp <- -log10(R.V2.for.GSEA.aggregated$V2.Pvalue)*sign(R.V2.for.GSEA.aggregated$V2.Estimate)
R.V2.for.GSEA.aggregated$logp.V3 <- -log10(R.V2.for.GSEA.aggregated$V3.Pvalue)*sign(R.V2.for.GSEA.aggregated$V3.Estimate)

pval <- 0.3
GO_file<- "h.all.v7.0.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)

gene_list.V2.R.net <- R.V2.for.GSEA.aggregated$logp
names(gene_list.V2.R.net) <- R.V2.for.GSEA.aggregated$SYMBOL

gene_list.V2.R.net <- gene_list.V2.R.net[order(gene_list.V2.R.net, decreasing=T)]

fgRes.hallmark.R.net.V2 <- fgsea::fgsea(pathways = myGO, 
                                        stats = gene_list.V2.R.net,
                                        minSize=10,
                                        maxSize=400,
                                        nperm=1000) %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < !!pval)

#for canonical pathways

GO_file<- "c2.cp.v7.0.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)

fgRes.canonical.R.net.V2 <- fgsea::fgsea(pathways = myGO, 
                                         stats = gene_list.V2.R.net,
                                         minSize=10,
                                         maxSize=400,
                                         nperm=1000) %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < !!pval)

#for GO terms
GO_file<- "c5.bp.v7.0.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)

fgRes.GO.R.net.V2 <- fgsea::fgsea(pathways = myGO, 
                                  stats = gene_list.V2.R.net,
                                  minSize=10,
                                  maxSize=400,
                                  nperm=1000) %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < !!pval)


fgRes.combined.R.V2 <- rbind(fgRes.hallmark.R.net.V2, fgRes.canonical.R.net.V2, fgRes.GO.R.net.V2)

fgRes.combined.R.V2.sig.and.in.network <- fgRes.combined.R.V2[fgRes.combined.R.V2$padj<0.15,]
fgRes.combined.R.V2.sig.and.in.network <- fgRes.combined.R.V2.sig.and.in.network[fgRes.combined.R.V2.sig.and.in.network$pathway %in% geneSets.Go.sub.filtered.names,]


#Test for responders V3

pval <- 0.3
GO_file<- "h.all.v7.0.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)

gene_list.V3.R.net <- R.V2.for.GSEA.aggregated$logp.V3
names(gene_list.V3.R.net) <- R.V2.for.GSEA.aggregated$SYMBOL

gene_list.V3.R.net <- gene_list.V3.R.net[order(gene_list.V3.R.net, decreasing=T)]

fgRes.hallmark.R.net.V3 <- fgsea::fgsea(pathways = myGO, 
                                        stats = gene_list.V3.R.net,
                                        minSize=10,
                                        maxSize=600,
                                        nperm=1000) %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < !!pval)

#for canonical pathways

GO_file<- "c2.cp.v7.0.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)

fgRes.canonical.R.net.V3 <- fgsea::fgsea(pathways = myGO, 
                                         stats = gene_list.V3.R.net,
                                         minSize=10,
                                         maxSize=600,
                                         nperm=1000) %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < !!pval)

#for GO terms

GO_file<- "c5.bp.v7.0.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)

fgRes.GO.R.net.V3 <- fgsea::fgsea(pathways = myGO, 
                                  stats = gene_list.V3.R.net,
                                  minSize=10,
                                  maxSize=600,
                                  nperm=1000) %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < !!pval)


fgRes.combined.R.V3 <- rbind(fgRes.hallmark.R.net.V3, fgRes.canonical.R.net.V3, fgRes.GO.R.net.V3)

fgRes.combined.R.V3.sig.and.in.network <- fgRes.combined.R.V3[fgRes.combined.R.V3$padj<0.15,]
fgRes.combined.R.V3.sig.and.in.network <- fgRes.combined.R.V3.sig.and.in.network[fgRes.combined.R.V3.sig.and.in.network$pathway %in% geneSets.Go.sub.filtered.names,]

Global.enriched.dynamic.path.list <- unique(rbind(fgRes.combined.R.V2.sig.and.in.network, fgRes.combined.R.V3.sig.and.in.network)$pathway)


#Test that pathways are connected in the network topology

edge <- drop_matrix.R.NR.V2.for.classification
edge <- edge[do.call('rbind', lapply(edge$edge.annotation, function(i) length(i)[[1]]))>0,]
edge <- edge[!do.call('rbind',lapply(edge$edge.annotation, function(i) grepl("only.one.annotated|both.not.annotated", i[[1]]))),]
node<-unique(c(edge$Var1SYMBOL, edge$Var2SYMBOL))
d <- data.frame(Var1SYMBOL=edge$Var1SYMBOL,Var2SYMBOL=edge$Var2SYMBOL, stringsAsFactors = F)
net <- graph.data.frame(d, vertices=node, directed=FALSE)


Net.enrichment.pathways.density <- do.call('rbind', lapply(unique(net.nodes.by.pathways$path), function(cur.path) {
  print(cur.path)
  node.list.in.module.and.net.edge <- unique(net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID'])
  node.list.in.module.and.net.edge.SYMBOL <- unique(net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'SYMBOL'])
  net.sub <- induced_subgraph(net, node.list.in.module.and.net.edge.SYMBOL)
  subnet.diameter <- diameter(net.sub, unconnected=T) #sub-network diameter (D), which uses the maximum length of all shortest paths between any two connected nodes
  DG <- edge_density(net.sub)  
  CE <- centr_degree(net.sub)$centralization/centr_degree_tmax(net.sub)
  clu <- components(net.sub)
  num.of.conected.components <- clu$no
  max.clust <- max(clu$csize)
  max.clust.perc <- max(clu$csize)/length(clu$membership)*100
  return(data.frame(path=cur.path,  net.diameter=subnet.diameter, DG=DG, CE=CE, num.of.conected.components=num.of.conected.components,max.clust=max.clust, max.clust.perc=max.clust.perc, stringsAsFactors = F))
}))

#Remove unconnected pathways in the net and pathways in which the largest cluster is less than 50%
connected.path.list <- unique(net.nodes.by.pathways$path)[unique(net.nodes.by.pathways$path) %in% Net.enrichment.pathways.density[Net.enrichment.pathways.density$max.clust>=10 & Net.enrichment.pathways.density$max.clust.perc>50,'path']]


#Compare density significance by comparison to random networks of the same size
random <- lapply(seq(200), function(m) {
  sample_gnm(length(node), nrow(edge), directed=FALSE, loops=FALSE)
})

random <- lapply(random, function(cur.net) {
  V(cur.net)$name <- names(V(net))
  return(cur.net)
})

Net.enrichment.pathways.random.res <- lapply(1:100, function(i) {
  print(paste("cur.random", i, sep=" "))
  Net.enrichment.pathways.random <- do.call('rbind', lapply(connected.path.list, function(cur.path) {
    cur.random.net <- random[[i]]
    submat <- drop_matrix.R.NR.V2.for.classification.all.test.GO[grepl(cur.path, drop_matrix.R.NR.V2.for.classification.all.test.GO$edge.annotation),]
    node.list.in.module.and.net.edge.SYMBOL <- unique(c(submat$Var1SYMBOL, submat$Var2SYMBOL))
    net.sub <- induced_subgraph(cur.random.net, node.list.in.module.and.net.edge.SYMBOL)
    subnet.diameter <- diameter(net.sub) #sub-network diameter (D), which uses the maximum length of all shortest paths between any two connected nodes
    # graph density (GD), a score that can also be defined as the local clustering coefficient. formula: 2E/n(n-1)
    DG <- edge_density(net.sub)  
    CE <- centr_degree(net.sub)$centralization/centr_degree_tmax(net.sub)
    return(data.frame(path=cur.path, net.diameter=subnet.diameter, DG=DG, CE=CE, stringsAsFactors = F))
  }))
  return(Net.enrichment.pathways.random)
})


#For DG:
DG.perm <- do.call('cbind', lapply(Net.enrichment.pathways.random.res, function(cur.m) {
  cur.m[,colnames(cur.m) %in% "DG"]
}))
colnames(DG.perm) <- paste("random.net", seq(1:100), sep=".")

Net.enrichment.pathways.sub.sig <- Net.enrichment.pathways.pathway.density[Net.enrichment.pathways.pathway.density$path %in% connected.path.list,]
Net.enrichment.pathways.sub.sig.random.DG.combined <- cbind(Net.enrichment.pathways.sub.sig[,colnames(Net.enrichment.pathways.sub.sig) %in% c("path", "DG")], DG.perm)
Net.enrichment.pathways.sub.sig.random.DG.combined$p <- apply(Net.enrichment.pathways.sub.sig.random.DG.combined, 1, function(x) {
  sum(as.numeric(as.character(x[!names(x) %in% c("path", "DG", "p", "FDR")]))>as.numeric(as.character(x[names(x) %in% c("DG")])))/200
})


#For net diameter:

diameter.perm <- do.call('cbind', lapply(Net.enrichment.pathways.random.res, function(cur.m) {
  cur.m[,colnames(cur.m) %in% "net.diameter"]
}))

colnames(diameter.perm) <- paste("random.net", seq(1:100), sep=".")

Net.enrichment.pathways.sub.sig.random.diameter.combined <- cbind(Net.enrichment.pathways.sub.sig[,colnames(Net.enrichment.pathways.sub.sig) %in% c("path", "net.diameter")], diameter.perm)
Net.enrichment.pathways.sub.sig.random.diameter.combined$p <- apply(Net.enrichment.pathways.sub.sig.random.diameter.combined, 1, function(x) {
  sum(as.numeric(as.character(x[!names(x) %in% c("path", "net.diameter", "p", "FDR")]))<as.numeric(as.character(x[names(x) %in% c("net.diameter")])))/200
})


#Test overlapping by jaccard index and exclude redundant pathways
net.nodes.by.pathways$path <- as.character(net.nodes.by.pathways$path)
net.nodes.by.pathways.sub <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% connected.path.list, ]

path.comb <- combinations(length(unique(net.nodes.by.pathways.sub$path)),2,unique(net.nodes.by.pathways.sub$path))
jaccrad.normal.similarity.for.merging.only.net <- do.call('rbind', lapply(1:nrow(path.comb), function(i) {
  nodes.normal.list1 <-  net.nodes.by.pathways.sub[net.nodes.by.pathways.sub$path %in% path.comb[i,1], 'SYMBOL']
  nodes.normal.list2 <-  net.nodes.by.pathways.sub[net.nodes.by.pathways.sub$path %in% path.comb[i,2], 'SYMBOL']
  intersect.var <- intersect(nodes.normal.list1, nodes.normal.list2)
  union.var <- union(nodes.normal.list1, nodes.normal.list2)
  Jaccard.similarity <- length(intersect.var)/length(union.var)
  return(data.frame(cur.path.normal1=path.comb[i,1], cur.path.normal2=path.comb[i,2], Jaccard.similarity=Jaccard.similarity, length.path1=length(nodes.normal.list1), length.path2=length(nodes.normal.list2), intersect=length(intersect.var), union=length(union.var), stringsAsFactors = F))
}))

jaccrad.normal.similarity.for.merging.only.net <- jaccrad.normal.similarity.for.merging.only.net[!jaccrad.normal.similarity.for.merging.only.net$cur.path.normal1==jaccrad.normal.similarity.for.merging.only.net$cur.path.normal2,]
list.of.redundant.path <- jaccrad.normal.similarity.for.merging.only.net[jaccrad.normal.similarity.for.merging.only.net$Jaccard.similarity>0.75,]

list.of.redundant.path.res <- apply(list.of.redundant.path, 1, function(x) {
  ifelse(x['length.path1']<x['length.path2'], x['cur.path.normal1'],
         ifelse(x['length.path1']==x['length.path2'],"", 
                ifelse(x['length.path1']>x['length.path2'], x['cur.path.normal2'], NA)))
})

list.of.redundant.path.res <- unique(list.of.redundant.path.res)
list.of.redundant.path.res <- list.of.redundant.path.res[!list.of.redundant.path.res %in% ""]

connected.path.list.final <- connected.path.list[!connected.path.list %in% list.of.redundant.path.res]
connected.path.list.final <- as.character(connected.path.list.final)

Global.enriched.dynamic.path.list.final <- Global.enriched.dynamic.path.list[Global.enriched.dynamic.path.list %in% connected.path.list.final]

#subset the gene set

geneSets.with.Go.sub <- geneSets.with.Go.sub[names(geneSets.with.Go.sub) %in% as.character(connected.path.list.final)]


#Test significant chage over time for pathway score
GSEA.clust.FC.summary.V1V2.wilcoxon.nodes <-  do.call('rbind', lapply(Global.enriched.dynamic.path.list.final, function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  esetMat.V1 <- filtered.esetALL.R[fData(filtered.esetALL.R)$featureID %in% node.list.in.module.and.net.edge, as.character(pData(filtered.esetALL.R)$Visit) %in% "V1"]
  expMat.V1 <- data.frame(exprs(esetMat.V1), stringsAsFactors = F)
  expMat.V1 <- data.frame(t(expMat.V1), stringsAsFactors = F)
  expMat.V1$time <- "V1"
  esetMat.V2 <- filtered.esetALL.R[fData(filtered.esetALL.R)$featureID %in% node.list.in.module.and.net.edge , as.character(pData(filtered.esetALL.R)$Visit) %in% "V2"]
  expMat.V2 <- data.frame(exprs(esetMat.V2), stringsAsFactors = F)
  expMat.V2 <- data.frame(t(expMat.V2), stringsAsFactors = F)
  expMat.V2$time <- "V2"
  combined.expMat <-rbind(expMat.V1, expMat.V2) 
  combined.expMat.scaled <- combined.expMat
  combined.expMat.scaled[,!colnames(combined.expMat.scaled) %in% "time"] <- apply(combined.expMat[,!colnames(combined.expMat) %in% "time"], 2, scale)
  row.names(combined.expMat.scaled) <- row.names(combined.expMat)
  combined.expMat.scaled$time <- factor(combined.expMat.scaled$time)
  print(cur.path)
  combined.expMat.scaled$module.score <- apply(combined.expMat.scaled[,!colnames(combined.expMat.scaled) %in% "time"], 1,  function(x) mean(x, na.rm=T))
  p.value=wilcox.test(module.score ~ time, paired=T, data=combined.expMat.scaled)$p.value 
  #res.permanova <- adonis(as.matrix(combined.expMat[,!colnames(combined.expMat) %in% "time"]) ~combined.expMat$time, data = combined.expMat, method="euclidean") 
  #p.value <- res.permanova$aov.tab$`Pr(>F)`[1]
  df <- data.frame(clust=cur.path, p.valueV2=p.value,stringsAsFactors = F)
  return(df)
}))



GSEA.clust.FC.summary.V1V3.wilcoxon.node <-  do.call('rbind', lapply(Global.enriched.dynamic.path.list.final, function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  esetMat.V1 <- filtered.esetALL.R[fData(filtered.esetALL.R)$featureID %in% node.list.in.module.and.net.edge, as.character(pData(filtered.esetALL.R)$Visit) %in% "V1"]
  expMat.V1 <- data.frame(exprs(esetMat.V1), stringsAsFactors = F)
  expMat.V1 <- data.frame(t(expMat.V1), stringsAsFactors = F)
  expMat.V1$time <- "V1"
  esetMat.V3 <- filtered.esetALL.R[fData(filtered.esetALL.R)$featureID %in% node.list.in.module.and.net.edge , as.character(pData(filtered.esetALL.R)$Visit) %in% "V3"]
  expMat.V3 <- data.frame(exprs(esetMat.V3), stringsAsFactors = F)
  expMat.V3 <- data.frame(t(expMat.V3), stringsAsFactors = F)
  expMat.V3$time <- "V3"
  combined.expMat <-rbind(expMat.V1, expMat.V3) 
  combined.expMat.scaled <- combined.expMat
  combined.expMat.scaled[,!colnames(combined.expMat.scaled) %in% "time"] <- apply(combined.expMat[,!colnames(combined.expMat) %in% "time"], 2, scale)
  row.names(combined.expMat.scaled) <- row.names(combined.expMat)
  combined.expMat.scaled$time <- factor(combined.expMat.scaled$time)
  print(cur.path)
  combined.expMat.scaled$module.score <- apply(combined.expMat.scaled[,!colnames(combined.expMat.scaled) %in% "time"], 1,  function(x) mean(x, na.rm=T))
  p.value=wilcox.test(module.score ~ time, paired=T, data=combined.expMat.scaled)$p.value 
  #res.permanova <- adonis(as.matrix(combined.expMat[,!colnames(combined.expMat) %in% "time"]) ~combined.expMat$time, data = combined.expMat, method="euclidean") 
  #p.value <- res.permanova$aov.tab$`Pr(>F)`[1]
  df <- data.frame(clust=cur.path, p.valueV3=p.value,stringsAsFactors = F)
  return(df)
}))



GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$BHV2 <- p.adjust(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$p.valueV2, method="BH")
GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$labelV2 <- ifelse(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$BHV2<0.005, "***", ifelse(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$BHV2<0.01, "**", ifelse(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$BHV2<0.05, "*", "")))
GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$p.valueV3 <- GSEA.clust.FC.summary.V1V3.wilcoxon.node[match(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$clust, GSEA.clust.FC.summary.V1V3.wilcoxon.node$clust),'p.valueV3']
GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$BHV3 <- p.adjust(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$p.valueV3, method="BH")
GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$labelV3 <- ifelse(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$BHV3<0.005, "***", ifelse(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$BHV3<0.01, "**", ifelse(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$BHV3<0.05, "*", "")))
GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$module.size <- net.nodes.by.pathways.df[match(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$clust,net.nodes.by.pathways.df$path),'size']


filtered.esetALL.FC.final.R <- filtered.esetALL.FC.final[,pData(filtered.esetALL.FC.final)$group %in% "R"]

GSEA.clust.FC.summary.V1V2 <- do.call('rbind', lapply(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$clust, function(cur.path) {
  node.list <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  esetMat.V2 <- filtered.esetALL.FC.final.R[fData(filtered.esetALL.FC.final.R)$featureID %in% node.list , as.character(pData(filtered.esetALL.FC.final.R)$Visit) %in% "V2"]
  esetMat.V2 <- data.frame(exprs(esetMat.V2), stringsAsFactors = F)
  expMat.V2.median <- apply(esetMat.V2, 2, function(x) median(x, na.rm=T))
  esetMat.V3 <- filtered.esetALL.FC.final.R[fData(filtered.esetALL.FC.final.R)$featureID %in% node.list , as.character(pData(filtered.esetALL.FC.final.R)$Visit) %in% "V3"]
  expMat.V3 <- data.frame(exprs(esetMat.V3), stringsAsFactors = F)
  expMat.V3.median <- apply(expMat.V3, 2, function(x) median(x, na.rm=T))
  df <- data.frame(clust=cur.path, Avg.medianV2=mean(expMat.V2.median), Avg.medianV3=mean(expMat.V3.median),SEM.V2=sd(expMat.V2.median, na.rm=T)/(sqrt(sum(!is.na(expMat.V2.median)))), SEM.V3=sd(expMat.V3.median, na.rm=T)/(sqrt(sum(!is.na(expMat.V3.median)))) ,stringsAsFactors = F)
  return(df)
}))

GSEA.clust.FC.summary.V1V2.final <- merge(GSEA.clust.FC.summary.V1V2, GSEA.clust.FC.summary.V1V2.wilcoxon.nodes, by="clust")

plot(GSEA.clust.FC.summary.V1V2.final$module.size, abs(GSEA.clust.FC.summary.V1V2.final$Avg.medianV2), method="spearman")  
corr.test(GSEA.clust.FC.summary.V1V2.final$module.size, abs(GSEA.clust.FC.summary.V1V2.final$Avg.medianV2), method="spearman")$p

#adjust median FC for module size

model1 <- lm(GSEA.clust.FC.summary.V1V2.final$Avg.medianV2 ~ GSEA.clust.FC.summary.V1V2.final$module.size ,data=GSEA.clust.FC.summary.V1V2.final, na.action="na.exclude")
Avg.medianV2_adj <- resid(model1) + coef(model1)[1]
GSEA.clust.FC.summary.V1V2.final$Avg.medianV2_adj <- Avg.medianV2_adj


model1 <- lm(GSEA.clust.FC.summary.V1V2.final$Avg.medianV3 ~ GSEA.clust.FC.summary.V1V2.final$module.size ,data=GSEA.clust.FC.summary.V1V2.final, na.action="na.exclude")
Avg.medianV3_adj <- resid(model1) + coef(model1)[1]
GSEA.clust.FC.summary.V1V2.final$Avg.medianV3_adj <- Avg.medianV3_adj


#Comparing the number of enriched dynamic pathways between 
#early and late response periodes in responders
#Supp. Fig 3

#plot per data type separately

GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$log.p.valueV2 <- -log10(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$p.valueV2)
GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$log.p.valueV3 <- -log10(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes$p.valueV3)


ggplot(GSEA.clust.FC.summary.V1V2.wilcoxon.nodes, aes(x=log.p.valueV2, y=log.p.valueV3)) + 
  ylab("Differential pathway score V1-V3 in responders \n (-log10(p))")+
  xlab("Differential pathway score V1-V2 in responders \n (-log10(p))")+
  geom_point(position = position_quasirandom()) + 
  theme_bw()+
 # geom_abline(slope=1, intercept=0)+
  geom_hline(yintercept=1.6, linetype="dashed")+
  geom_vline(xintercept=1.6, linetype="dashed")+
  #ylim(0,4)+xlim(0,4)+
  theme(strip.text = element_text(size=14), panel.grid.major = element_blank(),
        axis.text.y=element_text( size=12),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border =element_blank(), 
        axis.ticks.length=unit(.15, "cm"))


#Add scatterplot for pathway visualization: FC vs betweeness
#Fig 2c

edge<-drop_matrix.R.NR.V2.for.classification
edge <- edge[do.call('rbind', lapply(edge$edge.annotation, function(i) length(i)[[1]]))>0,]
edge <- edge[!do.call('rbind',lapply(edge$edge.annotation, function(i) grepl("only.one.annotated|both.not.annotated", i[[1]]))),]
node<-unique(c(edge$Var1, edge$Var2))

d <- data.frame(Var1=edge$Var1,Var2=edge$Var2, stringsAsFactors = F)
net <- graph.data.frame(d, vertices=node, directed=FALSE)

betweeness<- igraph::betweenness(net, normalized = T)
betweeness.df <- data.frame(featureID=names(betweeness), betweeness=betweeness, stringsAsFactors = F )
betweeness.df$SYMBOL <- fData(filtered.esetALL)[match(betweeness.df$featureID, fData(filtered.esetALL)$featureID),'SYMBOL']

#calculate for each module the average betweeness

betweeness.path <-  do.call('rbind', lapply(geneSets.Go.sub.filtered.names[!geneSets.Go.sub.filtered.names %in% c("only.one.annotated", "both.not.annotated")], function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  betweeness.df.sub <- betweeness.df[betweeness.df$featureID %in% node.list.in.module.and.net.edge,]
  df <- data.frame(clust=cur.path, betweeness=quantile(betweeness.df.sub$betweeness, 0.8), stringsAsFactors = F)
  return(df)
}))


GSEA.clust.FC.summary.V1V2.final$betweeness <- betweeness.path[match(GSEA.clust.FC.summary.V1V2.final$clust, betweeness.path$clust),'betweeness']

#calculate for each module the average degree

degree<- igraph::degree(net, normalized = T)
degree.df <- data.frame(featureID=names(degree), degree=degree, stringsAsFactors = F )
degree.df$SYMBOL <- fData(filtered.esetALL)[match(degree.df$featureID, fData(filtered.esetALL)$featureID),'SYMBOL']

degree.path <-  do.call('rbind', lapply(GSEA.clust.FC.summary.V1V2.final$clust, function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  degree.df.sub <- degree.df[degree.df$featureID %in% node.list.in.module.and.net.edge,]
  degree.df.sub <- degree.df.sub[order(degree.df.sub$degree, decreasing=T),]
  df <- data.frame(clust=cur.path, degree=quantile(degree.df.sub$degree, 0.8), stringsAsFactors = F)
  return(df)
}))

GSEA.clust.FC.summary.V1V2.final$degree <- degree.path[match(GSEA.clust.FC.summary.V1V2.final$clust, degree.path$clust),'degree']


#adjust betweenness and degree for module size

model1 <- lm(GSEA.clust.FC.summary.V1V2.final$betweeness ~ GSEA.clust.FC.summary.V1V2.final$module.size ,data=GSEA.clust.FC.summary.V1V2.final, na.action="na.exclude")
betweeness_adj <- resid(model1) + coef(model1)[1]
GSEA.clust.FC.summary.V1V2.final$betweeness_adj <- betweeness_adj


model1 <- lm(GSEA.clust.FC.summary.V1V2.final$degree ~ GSEA.clust.FC.summary.V1V2.final$module.size ,data=GSEA.clust.FC.summary.V1V2.final, na.action="na.exclude")
degree_adj <- resid(model1) + coef(model1)[1]
GSEA.clust.FC.summary.V1V2.final$degree_adj <- degree_adj

GSEA.clust.FC.summary.V1V2.final$order <- scale(GSEA.clust.FC.summary.V1V2.final$degree)+scale(GSEA.clust.FC.summary.V1V2.final$betweeness)
GSEA.clust.FC.summary.V1V2.final <- GSEA.clust.FC.summary.V1V2.final[order(GSEA.clust.FC.summary.V1V2.final$order, decreasing=T),]
labeled.path <- GSEA.clust.FC.summary.V1V2.final[1:35,'clust']
GSEA.clust.FC.summary.V1V2.final$label <- ifelse(GSEA.clust.FC.summary.V1V2.final$clust %in% labeled.path, GSEA.clust.FC.summary.V1V2.final$clust, "")

#GSEA.clust.FC.summary.V1V2.final <- read.csv("GSEA.clust.FC.summary.V1V2.final.csv", stringsAsFactors = F)
GSEA.clust.FC.summary.V1V2.final$logpV2 <- (-log10(GSEA.clust.FC.summary.V1V2.final$p.valueV2))
GSEA.clust.FC.summary.V1V2.final$colorV2 <- ifelse(GSEA.clust.FC.summary.V1V2.final$BHV2>0.05 , 0, GSEA.clust.FC.summary.V1V2.final$Avg.medianV2_adj)                
GSEA.clust.FC.summary.V1V2.final$abs_Avg.medianV2_adj <- abs(GSEA.clust.FC.summary.V1V2.final$Avg.medianV2_adj)

GSEA.clust.FC.summary.V1V2.final <- read.csv("GSEA.clust.FC.summary.V1V2.final.csv", stringsAsFactors = F)


ggplot(GSEA.clust.FC.summary.V1V2.final, aes(x=betweeness, y=degree, fill=colorV2, label=label))+
  geom_point(aes(size=5, shape="1", color="black"))+
  scale_color_manual(values="black")+
  scale_fill_gradient2(low =muted("red"), mid = "white", high=muted("blue"), breaks=c(-0.2,0,0.2),  labels=c("-0.2","0", "0.2"), midpoint = 0, limits = c(-0.4,0.4))+
  scale_shape_manual(values=c(21))+
  xlim(0,0.023)+
  # scale_size(range = c(0.5,10))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   #panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  geom_text_repel(data=GSEA.clust.FC.summary.V1V2.final, aes(label=label), size=4, 
                  direction = "y", 
                  #hjust = -0.5, 
                  vjust=0,
                  max.iter=10000,
                  #force = 40,
                  nudge_x = 0.05-GSEA.clust.FC.summary.V1V2.final$betweeness,
                  #nudge_y      = 40,
                  #direction    = "both",n
                  segment.color = "grey50")


pdf(file="plot2D.pdf", width=10, height=8)
g1
dev.off()

GSEA.clust.FC.summary.V1V2.final.for.Supp.table <- GSEA.clust.FC.summary.V1V2.final
GSEA.clust.FC.summary.V1V2.final.for.Supp.table <- merge(GSEA.clust.FC.summary.V1V2.final.for.Supp.table, Net.enrichment.pathways.density, all.x=T, by.x="clust", by.y="path")
GSEA.clust.FC.summary.V1V2.final.for.Supp.table <- merge(GSEA.clust.FC.summary.V1V2.final.for.Supp.table, fgRes.combined.R.V2.sig.and.in.network, all.x=T, by.x="clust", by.y="pathway")
GSEA.clust.FC.summary.V1V2.final.for.Supp.table <- merge(GSEA.clust.FC.summary.V1V2.final.for.Supp.table, fgRes.combined.R.V3.sig.and.in.network, all.x=T, by.x="clust", by.y="pathway")
GSEA.clust.FC.summary.V1V2.final.for.Supp.table$leadingEdge.x <- do.call('rbind', lapply(GSEA.clust.FC.summary.V1V2.final.for.Supp.table$leadingEdge.x, function(x) paste(x, collapse="|")))
GSEA.clust.FC.summary.V1V2.final.for.Supp.table$leadingEdge.y <- do.call('rbind', lapply(GSEA.clust.FC.summary.V1V2.final.for.Supp.table$leadingEdge.y, function(x) paste(x, collapse="|")))



#Plot Supp. Fig 3b
#Sort pathways by fold change

GSEA.clust.FC.summary.V1V2.final.sub <- GSEA.clust.FC.summary.V1V2.final[order(GSEA.clust.FC.summary.V1V2.final$Avg.medianV2_adj, GSEA.clust.FC.summary.V1V2.final$degree, GSEA.clust.FC.summary.V1V2.final$betweeness),]

GSEA.clust.FC.summary.V1V2.final.sub$clust <- factor(GSEA.clust.FC.summary.V1V2.final.sub$clust, levels=rev(GSEA.clust.FC.summary.V1V2.final.sub$clust))

GSEA.clust.FC.summary.V1V2.final.sub$colorV2.binary <- ifelse(GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV2_adj<0, "downregulated", "upregulated")
GSEA.clust.FC.summary.V1V2.final.sub$colorV3.binary <- ifelse(GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV3_adj<0, "downregulated", "upregulated")

GSEA.clust.FC.summary.V1V2.final.sub <- GSEA.clust.FC.summary.V1V2.final.sub[!GSEA.clust.FC.summary.V1V2.final.sub$label %in% "",]

GSEA.clust.FC.summary.V1V2.final.sub$logpV3 <- (-log10(GSEA.clust.FC.summary.V1V2.final.sub$p.valueV3))
GSEA.clust.FC.summary.V1V2.final.sub$logpV2 <- (-log10(GSEA.clust.FC.summary.V1V2.final.sub$p.valueV2))

GSEA.clust.FC.summary.V1V2.final.sub <- GSEA.clust.FC.summary.V1V2.final.sub[order(GSEA.clust.FC.summary.V1V2.final.sub$logpV2*sign(GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV2_adj), decreasing=T),]
GSEA.clust.FC.summary.V1V2.final.sub$clust <- factor(GSEA.clust.FC.summary.V1V2.final.sub$clust, levels=GSEA.clust.FC.summary.V1V2.final.sub$clust)
GSEA.clust.FC.summary.V1V2.final.sub <- rbind(head(GSEA.clust.FC.summary.V1V2.final.sub, n=60), tail(GSEA.clust.FC.summary.V1V2.final.sub, n=10))

GSEA.clust.FC.summary.V1V2.final.sub <- GSEA.clust.FC.summary.V1V2.final.sub[order(GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV2_adj),]
GSEA.clust.FC.summary.V1V2.final.sub$clust <- factor(GSEA.clust.FC.summary.V1V2.final.sub$clust, levels=GSEA.clust.FC.summary.V1V2.final.sub$clust)

ggplot(GSEA.clust.FC.summary.V1V2.final.sub, aes(x=clust, y=Avg.medianV2_adj, fill=colorV2.binary)) + 
  geom_bar(position=position_dodge(), stat="identity") + 
  xlab("Pathway")+ylab("log2(FC.V2)")+
  coord_flip()+
  ylim(-0.6, 0.6)+
  scale_fill_manual(values = c("downregulated"=muted("red"), "upregulated"= muted("blue")))+
  geom_text(data = GSEA.clust.FC.summary.V1V2.final.sub, label = GSEA.clust.FC.summary.V1V2.final.sub$labelV2, y= ifelse(GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV2_adj>0, GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV2_adj+0.1, GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV2_adj-0.1))+
  geom_errorbar(aes(ymin=GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV2_adj-GSEA.clust.FC.summary.V1V2.final.sub$SEM.V2, ymax=GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV2_adj+GSEA.clust.FC.summary.V1V2.final.sub$SEM.V2),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")+
  
  ggplot(GSEA.clust.FC.summary.V1V2.final.sub, aes(x=clust, y=Avg.medianV3_adj, fill=colorV3.binary)) + 
  geom_bar(position=position_dodge(), stat="identity") + xlab("Pathway")+ylab("log2(FC.V3)")+
  coord_flip()+
   ylim(-0.6, 0.6)+
    scale_fill_manual(values = c("downregulated"=muted("red"), "upregulated"= muted("blue")))+
  geom_text(data = GSEA.clust.FC.summary.V1V2.final.sub, label = GSEA.clust.FC.summary.V1V2.final.sub$labelV3, y= ifelse(GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV3_adj>0, GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV3_adj+0.07, GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV3_adj-0.07))+
  geom_errorbar(aes(ymin=GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV3_adj-GSEA.clust.FC.summary.V1V2.final.sub$SEM.V3, ymax=GSEA.clust.FC.summary.V1V2.final.sub$Avg.medianV3_adj+GSEA.clust.FC.summary.V1V2.final.sub$SEM.V3),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                  # axis.text.y = element_text(size=8),
                   axis.text.y = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))


#Test addition of CRP correlation parameter
# Plot Supp. Fig 3d

FC.summary.correlation.with.CRP <- do.call('rbind', lapply(GSEA.clust.FC.summary.V1V2.final$clust, function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  esetMat.V1 <- filtered.esetALL[fData(filtered.esetALL)$featureID %in% node.list.in.module.and.net.edge, as.character(pData(filtered.esetALL)$Visit) %in% "V1" & pData(filtered.esetALL)$group %in% c("R", "IM", "NR")]
  expMat.V1 <- data.frame(exprs(esetMat.V1), stringsAsFactors = F)
  expMat.V1 <- data.frame(t(expMat.V1), stringsAsFactors = F)
  expMat.V1$time <- "V1"
  esetMat.V2 <- filtered.esetALL[fData(filtered.esetALL)$featureID %in% node.list.in.module.and.net.edge , as.character(pData(filtered.esetALL)$Visit) %in% "V2" & pData(filtered.esetALL)$group %in% c("R", "IM", "NR")]
  expMat.V2 <- data.frame(exprs(esetMat.V2), stringsAsFactors = F)
  expMat.V2 <- data.frame(t(expMat.V2), stringsAsFactors = F)
  expMat.V2$time <- "V2"
  esetMat.V3 <- filtered.esetALL[fData(filtered.esetALL)$featureID %in% node.list.in.module.and.net.edge , as.character(pData(filtered.esetALL)$Visit) %in% "V3" & pData(filtered.esetALL)$group %in% c("R", "IM", "NR")]
  expMat.V3 <- data.frame(exprs(esetMat.V3), stringsAsFactors = F)
  expMat.V3 <- data.frame(t(expMat.V3), stringsAsFactors = F)
  expMat.V3$time <- "V3"
  # combined.expMat <-rbind(expMat.V1, expMat.V2) 
  
  combined.expMat <-rbind(expMat.V1, expMat.V2, expMat.V3) 
  combined.expMat.scaled <- combined.expMat
  combined.expMat.scaled[,!colnames(combined.expMat.scaled) %in% "time"] <- apply(combined.expMat[,!colnames(combined.expMat) %in% "time"], 2, scale)
  row.names(combined.expMat.scaled) <- row.names(combined.expMat)
  combined.expMat.scaled$time <- factor(combined.expMat.scaled$time)
  combined.expMat.scaled$module.score <- apply(combined.expMat.scaled[,!colnames(combined.expMat.scaled) %in% "time"], 1,  function(x) mean(x, na.rm=T))
  
  row.names(combined.expMat.scaled) <- gsub("\\.", "-", row.names(combined.expMat.scaled))
  combined.expMat.scaled$CRP <- pData(filtered.esetALL)[match(row.names(combined.expMat.scaled), pData(filtered.esetALL)$sampleID),'CRP']
  
  cor.res.V2=corr.test(combined.expMat.scaled$module.score[combined.expMat.scaled$time %in% c("V1", "V2")], combined.expMat.scaled$CRP[combined.expMat.scaled$time %in% c("V1", "V2")], method="spearman")
  cor.res.V3=corr.test(combined.expMat.scaled$module.score[combined.expMat.scaled$time %in% c("V1", "V3")], combined.expMat.scaled$CRP[combined.expMat.scaled$time %in% c("V1", "V3")], method="spearman")
  cor.res.all <- corr.test(combined.expMat.scaled$module.score, combined.expMat.scaled$CRP, method="spearman")
  
  return(data.frame(path=cur.path, r=cor.res.all$r, p=cor.res.all$p, r.V2=cor.res.V2$r, p.V2=cor.res.V2$p,r.V3=cor.res.V3$r, p.V3=cor.res.V3$p, stringsAsFactors=F))
}))


GSEA.clust.FC.summary.V1V2.final.with.CRP <- merge(GSEA.clust.FC.summary.V1V2.final, FC.summary.correlation.with.CRP, by.x="clust", by.y="path")

GSEA.clust.FC.summary.V1V2.final.with.CRP <- GSEA.clust.FC.summary.V1V2.final.with.CRP[GSEA.clust.FC.summary.V1V2.final.with.CRP$clust %in% GSEA.clust.FC.summary.V1V2.final.sub$clust,]
GSEA.clust.FC.summary.V1V2.final.with.CRP$r.FDR <- p.adjust(GSEA.clust.FC.summary.V1V2.final.with.CRP$p, method="BH")
GSEA.clust.FC.summary.V1V2.final.with.CRP$color.r <- ifelse(GSEA.clust.FC.summary.V1V2.final.with.CRP$r<0 & GSEA.clust.FC.summary.V1V2.final.with.CRP$r.FDR<0.05, "negative", 
                                                            ifelse(GSEA.clust.FC.summary.V1V2.final.with.CRP$r>0 & GSEA.clust.FC.summary.V1V2.final.with.CRP$r.FDR<0.05, "positive", "NS")) 

GSEA.clust.FC.summary.V1V2.final.with.CRP <- GSEA.clust.FC.summary.V1V2.final.with.CRP[order(GSEA.clust.FC.summary.V1V2.final.with.CRP$r),]
GSEA.clust.FC.summary.V1V2.final.with.CRP$clust <- factor(GSEA.clust.FC.summary.V1V2.final.with.CRP$clust, levels=GSEA.clust.FC.summary.V1V2.final.with.CRP$clust)

ggplot(GSEA.clust.FC.summary.V1V2.final.with.CRP, aes(x=clust, y=r, fill=color.r)) + 
  geom_bar(position=position_dodge(), stat="identity") + 
  xlab("Pathway")+ylab("Spearman's r")+
  coord_flip()+
  scale_fill_manual(values = c("negative"=muted("red"), "positive"= muted("blue"), "NS"="lightgrey"))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   axis.text.y = element_text(size=8),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"), 
                   legend.position = "none")



#Cell contribution for pathway fold change
#Test contribution by linear model and beta coefficient measurement
#Plot Supp. Fig 3c 

major <- c(GRAN = 'GR.TOTAL', CD8 = 'PB.1079991', CD4 = 'PB.1079995', CD19 = 'PB.1079966', CD14 = 'PB.1079983', NK = c('PB.1079974', 'PB.1079942'))
P <- exprs(filtered.esetALL.FC.final)[major, ]
rownames(P) <- names(major)
P <- rbind(P, NK=P['NK1', ]+P['NK2', ])
P <- P[!row.names(P) %in% c('NK1','NK2'),]
P_t <- data.frame(t(P), stringsAsFactors = F)
row.names(P_t) <- gsub("-", "\\.", row.names(P_t))
P_t$sampleID <- row.names(P_t)


P <- t(P_t)
P <- P[!row.names(P) %in% "sampleID",]
P.beta <- apply(P[!row.names(P) %in% "sampleID",], 1, function(x) scale(as.numeric(x), center = TRUE, scale = TRUE))
P.beta <- t(P.beta)
colnames(P.beta) <- colnames(P)
P <- P.beta

if(identical(gsub("-", "\\.", colnames(exprs(filtered.esetALL.FC.final[fData(filtered.esetALL.FC.final)$dataType %in% "GX",]))), colnames(P))) {

fit.res.for.beta <- lsfit(t(P), 2^t(exprs(filtered.esetALL.FC.final[fData(filtered.esetALL.FC.final)$dataType %in% "GX",])), wt = NULL, intercept = TRUE,
                          yname = row.names(exprs(filtered.esetALL.FC.final[fData(filtered.esetALL.FC.final)$dataType %in% "GX",])))
} else {print("not identical")}

fit.res.bet <-  data.frame(t(fit.res.for.beta$coefficients), stringsAsFactors = F)
fit.res.bet$Intercept <- NULL 
fit.res.bet$SYMBOL <- fData(filtered.esetALL)[match(row.names(fit.res.bet), fData(filtered.esetALL)$featureID),'SYMBOL']
fit.res.bet$featureID <- row.names(fit.res.bet)

FC.V1V2.beta.of.cell.contributions <- do.call('rbind', lapply(as.character(GSEA.clust.FC.summary.V1V2.final.sub$clust), function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  fit.res.bet.sub <- fit.res.bet[fit.res.bet$featureID %in% gsub("^AG", "GX", node.list.in.module.and.net.edge),]
    beta.module <- apply(fit.res.bet.sub[,!colnames(fit.res.bet.sub) %in% c("SYMBOL","featureID")], 2, mean)
  beta.module <- data.frame(beta.module, stringsAsFactors = F)
  return(data.frame(path=cur.path, beta.module=beta.module, cells=row.names(beta.module), stringsAsFactors = F))
}))

FC.V1V2.beta.of.cell.contributions.sub <- FC.V1V2.beta.of.cell.contributions[FC.V1V2.beta.of.cell.contributions$path %in% GSEA.clust.FC.summary.V1V2.final.sub$clust,]
FC.V1V2.beta.of.cell.contributions.sub$path <-factor(FC.V1V2.beta.of.cell.contributions.sub$path, levels=levels(GSEA.clust.FC.summary.V1V2.final.sub$clust)) 


FC.V1V2.beta.of.cell.contributions.sub.wide <- reshape2::dcast(FC.V1V2.beta.of.cell.contributions.sub, path ~ cells, value.var="beta.module")
row.names(FC.V1V2.beta.of.cell.contributions.sub.wide) <- FC.V1V2.beta.of.cell.contributions.sub.wide$path
FC.V1V2.beta.of.cell.contributions.sub.wide$path <- NULL 


ht <- Heatmap(FC.V1V2.beta.of.cell.contributions.sub.wide,
              rect_gp = gpar(col = "white", lwd = 2), row_title = NULL,
              col=c( "#084594", "#4292C6", "#9ECAE1","#DEEBF7" ,"#f0e3ab", "#f49d77", "#ce5f5f", "#b23030"),
              name="cell contribution (beta)", row_names_gp = gpar(fontsize = 8), 
              column_names_gp = gpar(fontsize = 7), show_row_names = TRUE, show_column_names = T,
              column_names_rot = 90, column_names_max_height = unit(10,"cm"))

draw(ht, padding = unit(c(2, 2, 2, 90), "mm"), heatmap_legend_side = "left") #bottom, left, top, right paddings



#Calculate disruption in the geneset level

patients <- colnames(drop_matrix.R.NR.V2)[!colnames(drop_matrix.R.NR.V2) %in% c("Var1", "Var2")]
patients <- patients[!patients %in% "HR.37.V2"]
GSEA.Module.disruption.measurment.results <- GSEA.Module.disruption.measurment.function(connected.path.list.final, drop_matrix.R.NR.V2, "V2")

GSEA.Module.disruption.measurment.results <- GSEA.Module.disruption.measurment.results[complete.cases(GSEA.Module.disruption.measurment.results),]

GSEA.Module.disruption.measurment.results$patient <- gsub("\\.", "-", GSEA.Module.disruption.measurment.results$patient)
GSEA.Module.disruption.measurment.results$group <- pData(filtered.esetALL)[match(GSEA.Module.disruption.measurment.results$patient, pData(filtered.esetALL)$sampleID),'group']
GSEA.Module.disruption.measurment.results$group <- ifelse(GSEA.Module.disruption.measurment.results$group %in% "IM", "NR", GSEA.Module.disruption.measurment.results$group)


#For disrupted edges

GSEA.Module.disruption.measurment.results.wide <- reshape2::dcast(GSEA.Module.disruption.measurment.results, patient+group ~ cluster, value.var="percentage.disrupt.in.module")

GSEA.Module.disruption.measurment.results.wide.final <- data.frame(t(GSEA.Module.disruption.measurment.results.wide), stringsAsFactors = FALSE)
colnames(GSEA.Module.disruption.measurment.results.wide.final) <- GSEA.Module.disruption.measurment.results.wide.final[row.names(GSEA.Module.disruption.measurment.results.wide.final) %in% "group",]
GSEA.Module.disruption.measurment.results.wide.tmp <- GSEA.Module.disruption.measurment.results.wide.final[-(1:2),]

GSEA.Module.disruption.measurment.results.wide.final <- apply(GSEA.Module.disruption.measurment.results.wide.tmp,2, as.numeric)
row.names(GSEA.Module.disruption.measurment.results.wide.final) <- row.names(GSEA.Module.disruption.measurment.results.wide.tmp)

edge.disruption.range <- lapply(seq(0, 1, by=0.1), function(cur.percentile) {
  cur.threshold <- quantile(GSEA.Module.disruption.measurment.results.wide.final, cur.percentile)
  for.num.disrupted.NR.edge <- apply(GSEA.Module.disruption.measurment.results.wide.final[,grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/8*100)
  for.num.disrupted.R.edge <- apply(GSEA.Module.disruption.measurment.results.wide.final[,!grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/15*100)
  return(data.frame(percentile=cur.percentile, num.disrupted.NR.by.edge=for.num.disrupted.NR.edge, num.disrupted.R.by.edge=for.num.disrupted.R.edge, stringsAsFactors = F))
})


#percentage of disruption by intensity

GSEA.Module.disruption.measurment.results.wide <- reshape2::dcast(GSEA.Module.disruption.measurment.results, patient+group ~ cluster, value.var="mean.drop.intensity")

GSEA.Module.disruption.measurment.results.wide.final <- data.frame(t(GSEA.Module.disruption.measurment.results.wide), stringsAsFactors = FALSE)
colnames(GSEA.Module.disruption.measurment.results.wide.final) <- GSEA.Module.disruption.measurment.results.wide.final[row.names(GSEA.Module.disruption.measurment.results.wide.final) %in% "group",]
GSEA.Module.disruption.measurment.results.wide.tmp <- GSEA.Module.disruption.measurment.results.wide.final[-(1:2),]

GSEA.Module.disruption.measurment.results.wide.final <- apply(GSEA.Module.disruption.measurment.results.wide.tmp,2, as.numeric)
row.names(GSEA.Module.disruption.measurment.results.wide.final) <- row.names(GSEA.Module.disruption.measurment.results.wide.tmp)

intensity.disruption.range <- lapply(seq(0, 1, by=0.1), function(cur.percentile) {
  cur.threshold <- quantile(GSEA.Module.disruption.measurment.results.wide.final, cur.percentile)
  for.num.disrupted.NR.edge <- apply(GSEA.Module.disruption.measurment.results.wide.final[,grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/8*100)
  for.num.disrupted.R.edge <- apply(GSEA.Module.disruption.measurment.results.wide.final[,!grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/15*100)
  return(data.frame(percentile=cur.percentile, num.disrupted.NR.by.edge=for.num.disrupted.NR.edge, num.disrupted.R.by.edge=for.num.disrupted.R.edge, stringsAsFactors = F))
})

#percentage of disrupted nodes

GSEA.Module.disruption.measurment.results.wide <- reshape2::dcast(GSEA.Module.disruption.measurment.results, patient+group ~ cluster, value.var="percentage.of.disrupted.nodes")
GSEA.Module.disruption.measurment.results.wide.final <- data.frame(t(GSEA.Module.disruption.measurment.results.wide), stringsAsFactors = FALSE)
colnames(GSEA.Module.disruption.measurment.results.wide.final) <- GSEA.Module.disruption.measurment.results.wide.final[row.names(GSEA.Module.disruption.measurment.results.wide.final) %in% "group",]
GSEA.Module.disruption.measurment.results.wide.tmp <- GSEA.Module.disruption.measurment.results.wide.final[-(1:2),]

GSEA.Module.disruption.measurment.results.wide.final <- apply(GSEA.Module.disruption.measurment.results.wide.tmp,2, as.numeric)
row.names(GSEA.Module.disruption.measurment.results.wide.final) <- row.names(GSEA.Module.disruption.measurment.results.wide.tmp)

node.disruption.range <- lapply(seq(0, 1, by=0.1), function(cur.percentile) {
  cur.threshold <- quantile(GSEA.Module.disruption.measurment.results.wide.final, cur.percentile)
  for.num.disrupted.NR.edge <- apply(GSEA.Module.disruption.measurment.results.wide.final[,grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/8*100)
  for.num.disrupted.R.edge <- apply(GSEA.Module.disruption.measurment.results.wide.final[,!grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/15*100)
  return(data.frame(percentile=cur.percentile, num.disrupted.NR.by.edge=for.num.disrupted.NR.edge, num.disrupted.R.by.edge=for.num.disrupted.R.edge, stringsAsFactors = F))
})


#Test agreement between parameters across a dynamic range of percentiles 
dynamic.range.venn.diagram <- do.call('rbind', lapply(seq(length(edge.disruption.range)), function(cur.percentile) {
  cur.edge.disruption.range <- edge.disruption.range[[cur.percentile]]
  cur.edge.disruption.range.sub <- cur.edge.disruption.range[cur.edge.disruption.range$num.disrupted.NR.by.edge>60  ,]
  
  cur.intensity.disruption.range <- intensity.disruption.range[[cur.percentile]]
  cur.intensity.disruption.range.sub <- cur.intensity.disruption.range[cur.intensity.disruption.range$num.disrupted.NR.by.edge>60,]
  
  cur.node.disruption.range <- node.disruption.range[[cur.percentile]]
  cur.node.disruption.range.sub <- cur.node.disruption.range[cur.node.disruption.range$num.disrupted.NR.by.edge>60,]
  
  first.intersect <- intersect(row.names(cur.edge.disruption.range.sub), row.names(cur.intensity.disruption.range.sub))
  second.intersect <- intersect(first.intersect, row.names(cur.node.disruption.range.sub))
  count.df <- data.frame(percentile=unique(cur.edge.disruption.range$percentile)[1],  edge=nrow(cur.edge.disruption.range.sub), intensity=nrow(cur.intensity.disruption.range.sub), node=nrow(cur.node.disruption.range.sub), shared= length(second.intersect), stringsAsFactors = F)
  return(count.df)
}))


dynamic.range.venn.diagram.long <- reshape2::melt(dynamic.range.venn.diagram, id.vars=c("percentile"))

ggplot(dynamic.range.venn.diagram.long, aes(x=percentile, y=value, color=variable))+
  geom_line()+
  ylab("# disrupted pathways")+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                               axis.text.y = element_text(size=8),
                               panel.grid.minor = element_blank(),                                                                  
                               strip.background = element_blank(),
                               panel.border = element_blank(), 
                               axis.ticks.length=unit(.15, "cm"), 
                               legend.position = "left")



dynamic.range.venn.diagram.shared.pathways <- do.call('rbind', lapply(seq(length(intensity.disruption.range)), function(cur.percentile) {
  cur.edge.disruption.range <- edge.disruption.range[[cur.percentile]]
  cur.edge.disruption.range.sub <- cur.edge.disruption.range[cur.edge.disruption.range$num.disrupted.NR.by.edge>60 ,]
  cur.intensity.disruption.range <- intensity.disruption.range[[cur.percentile]]
  cur.intensity.disruption.range.sub <- cur.intensity.disruption.range[cur.intensity.disruption.range$num.disrupted.NR.by.edge>60 ,]
  cur.node.disruption.range <- node.disruption.range[[cur.percentile]]
 cur.node.disruption.range.sub <- cur.node.disruption.range[cur.node.disruption.range$num.disrupted.NR.by.edge>60 ,]
  first.intersect <- intersect(row.names(cur.edge.disruption.range.sub), row.names(cur.intensity.disruption.range.sub))
  second.intersect <- intersect(first.intersect, row.names(cur.node.disruption.range.sub))
  percentile <- seq(0,1, by=0.1)[cur.percentile]
  count.df <- data.frame(percentile=rep(percentile, length(second.intersect)), shared.pathways=second.intersect, stringsAsFactors = F)
  return(count.df)
}))

#select the pathways with agreement of all 3 parameters in the highest percentile

selected.percentile <- 0.8
dynamic.range.venn.diagram.shared.pathways$percentile <- as.numeric(as.character(dynamic.range.venn.diagram.shared.pathways$percentile))
dynamic.range.venn.diagram.shared.pathways$shared.pathways <- as.character(dynamic.range.venn.diagram.shared.pathways$shared.pathways)
shared.disrupted.pathways <- dynamic.range.venn.diagram.shared.pathways[dynamic.range.venn.diagram.shared.pathways$percentile==selected.percentile,'shared.pathways']


#Plot a scatterplot of representative highly disrupted edge
#Plot Supp. Fig 4a

#test manually correlation between paired genes
#IL5 and KAZN as examples

filtered.esetALL.FC.final.R.V2 <-  filtered.esetALL.FC.final.R[,pData(filtered.esetALL.FC.final.R)$Visit %in% "V2"]
filtered.esetALL.FC.final.NR.V2 <- filtered.esetALL.FC.final.NR[,pData(filtered.esetALL.FC.final.NR)$Visit %in% "V2"] 

genex <- "AG.TC0100006954.hg.1"
geney <- "AG.TC0500012019.hg.1"
genexName <- fData(filtered.esetALL)[match(genex, fData(filtered.esetALL)$featureID),'displayName']
geneyName <- fData(filtered.esetALL)[match(geney, fData(filtered.esetALL)$featureID),'displayName']


df.R <- data.frame(x=t(exprs(filtered.esetALL.FC.final.R.V2[fData(filtered.esetALL.FC.final.R.V2)$featureID %in% genex,])),
                    y=t(exprs(filtered.esetALL.FC.final.R.V2[fData(filtered.esetALL.FC.final.R.V2)$featureID %in% geney,])),
                    stringsAsFactors = F)

ggplot(df.R, aes(x=df.R[,colnames(df.R) %in% genex], y=df.R[,colnames(df.R) %in% geney]))+
  geom_point()+ geom_smooth(method = lm, se = FALSE)+
  theme_bw()

corr.test(df.R[,paste(genex)], df.R[,paste(geney)], method='spearman')$r

#put all regression lines in one plot
NR.samples.regression <- do.call('rbind', lapply(colnames(filtered.esetALL.FC.final.NR.V2), function(cur.patient) {
  additional.sample <- t((exprs(filtered.esetALL.FC.final.NR.V2[fData(filtered.esetALL.FC.final.NR.V2)$featureID %in% c(genex, geney), pData(filtered.esetALL.FC.final.NR.V2)$sampleID %in% cur.patient])))
  colnames(additional.sample) <- c(paste(genex), paste(geney))
  df.disrupted <- rbind(df.R, additional.sample)
  df.disrupted$patient <- rep(cur.patient, nrow( df.disrupted))
  df.disrupted$colour <- c(rep("grey75", dim(filtered.esetALL.FC.final.R.V2)[2]), "red")
  return(df.disrupted)
}))


Normal.samples.regression <- do.call('rbind', lapply(colnames(filtered.esetALL.FC.final.R.V2), function(cur.patient) {
  df.disrupted <- df.R[!row.names(df.R) %in% cur.patient,]
  df.disrupted$patient <- rep(cur.patient, nrow( df.disrupted))
  df.disrupted$colour <- "grey75"
  return(df.disrupted)
}))


samples.regression.combined <- rbind(Normal.samples.regression, NR.samples.regression)
samples.regression.combined$group <- pData(filtered.esetALL.FC.final)[match(samples.regression.combined$patient, pData(filtered.esetALL.FC.final)$sampleID),'group']
samples.regression.combined$group <- ifelse(samples.regression.combined$group %in% "IM", "NR", samples.regression.combined$group)
line.color.for.each.patient <- unique(samples.regression.combined$patient)
samples.regression.combined$colour <- factor(samples.regression.combined$colour)
line.color.for.each.patient.final <- ifelse(line.color.for.each.patient %in% pData(filtered.esetALL.FC.final)[grepl("^R", pData(filtered.esetALL.FC.final)$group), 'sampleID'], "grey30", 
                                            ifelse(line.color.for.each.patient %in% pData(filtered.esetALL.FC.final)[grepl("IM|NR", pData(filtered.esetALL.FC.final)$group), 'sampleID'], "red", "cyan3"))


#plot only normal.samples (responders)

Normal.samples.regression$group <- pData(filtered.esetALL.FC.final)[match(Normal.samples.regression$patient, pData(filtered.esetALL.FC.final)$sampleID),'group']
line.color.for.each.patient <- unique(Normal.samples.regression$patient)
Normal.samples.regression$colour <- factor(Normal.samples.regression$colour)
line.color.for.each.patient.final <- ifelse(line.color.for.each.patient %in% pData(filtered.esetALL.FC.final)[grepl("^R", pData(filtered.esetALL.FC.final)$group), 'sampleID'], "grey30", 
                                            ifelse(line.color.for.each.patient %in% pData(filtered.esetALL.FC.final)[grepl("NR|IM", pData(filtered.esetALL.FC.final)$group), 'sampleID'], "red", "cyan3"))


colnames(Normal.samples.regression)[1:2] <- fData(filtered.esetALL)[match(colnames(Normal.samples.regression)[1:2], fData(filtered.esetALL)$featureID),'SYMBOL']

ggplot(Normal.samples.regression, aes(x=Normal.samples.regression[,colnames(Normal.samples.regression) %in% colnames(Normal.samples.regression)[1]], y=Normal.samples.regression[,colnames(Normal.samples.regression) %in% colnames(Normal.samples.regression)[2]]))+
  geom_point(aes(fill=colour), pch=21, size=3)+ 
  geom_smooth(method = lm, se = F, aes(group=patient, colour="grey75", fill="grey75"), size=0.2, alpha=0.1)+
  scale_fill_manual(values=c("grey75"="cyan3"))+
  scale_color_manual(values=c("grey75"="cyan3"))+
  theme_bw()+
  ylim(-1.2,0.6)+
  xlab(colnames(Normal.samples.regression)[1])+ylab(colnames(Normal.samples.regression)[2])+
  theme(axis.text.x = element_text( hjust = 1),strip.text = element_text(size=14), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))


#calculate combined regression (responders and non-responders)

samples.regression.combined <- rbind(Normal.samples.regression, NR.samples.regression)
samples.regression.combined$group <- pData(filtered.esetALL.FC.final)[match(samples.regression.combined$patient, pData(filtered.esetALL.FC.final)$sampleID),'group']
samples.regression.combined$group <- ifelse(samples.regression.combined$group %in% "IM", "NR", samples.regression.combined$group)
line.color.for.each.patient <- unique(samples.regression.combined$patient)
samples.regression.combined$colour <- factor(samples.regression.combined$colour)
line.color.for.each.patient.final <- ifelse(line.color.for.each.patient %in% pData(filtered.esetALL.FC.final)[grepl("^R", pData(filtered.esetALL.FC.final)$group), 'sampleID'], "grey30", 
                                            ifelse(line.color.for.each.patient %in% pData(filtered.esetALL.FC.final)[grepl("IM|NR", pData(filtered.esetALL.FC.final)$group), 'sampleID'], "red", "cyan3"))

colnames(samples.regression.combined)[1:2] <- fData(filtered.esetALL)[match(colnames(samples.regression.combined)[1:2], fData(filtered.esetALL)$featureID),'SYMBOL']
ggplot(samples.regression.combined, aes(x=samples.regression.combined[,colnames(samples.regression.combined) %in% colnames(samples.regression.combined)[1]], y=samples.regression.combined[,colnames(samples.regression.combined) %in% colnames(samples.regression.combined)[2]]))+
  geom_point(aes(fill=colour), pch=21, size=3)+ 
  geom_smooth(data=samples.regression.combined[samples.regression.combined$group %in% "NR",],
              aes(x=samples.regression.combined[samples.regression.combined$group %in% "NR",colnames(samples.regression.combined) %in% colnames(samples.regression.combined)[1]], 
                  y=samples.regression.combined[samples.regression.combined$group %in% "NR",colnames(samples.regression.combined) %in% colnames(samples.regression.combined)[2]], 
                  group=patient, color="red", fill="red"),
              method = lm, se = T, size=0.2, alpha=0.1)+
  geom_smooth(data=samples.regression.combined[samples.regression.combined$group %in% "R",],
              aes(x=samples.regression.combined[samples.regression.combined$group %in% "R",colnames(samples.regression.combined) %in% colnames(samples.regression.combined)[1]], 
                  y=samples.regression.combined[samples.regression.combined$group %in% "R",colnames(samples.regression.combined) %in% colnames(samples.regression.combined)[2]], 
                  group=group, fill="grey75", color="grey75"),
              method = lm, se = T, size=1.5, alpha=0.7)+
  scale_fill_manual(values=c("R"="black", "grey75"="cyan3","red"="red"))+
  scale_colour_manual(values=c("R"="black", "grey75"="cyan3","red"="red"))+
  theme_bw()+
  #xlim(-0.5,0.4)+ylim(-1.2,0.6)+
  xlab(colnames(samples.regression.combined)[1])+ylab(colnames(samples.regression.combined)[2])+
  theme(axis.text.x = element_text( hjust = 1),strip.text = element_text(size=14), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))



#Test relationship between wilcoxon FC p value between R and NR and disruption measurments
#Plot Supp. Fig 4b

disruption.count.V2 <- disruption.count.fun(cor.edges.df.FC.V1V2.aggregated, drop_matrix.V2.NUR, "V2")
disruption.count.V3 <- disruption.count.fun(cor.edges.df.FC.V1V3.aggregated, drop_matrix.V3.NUR, "V3")


#Plot of V2 for percentage of disrupted individuals

filtered.esetALL.FC.final.R.V2.noGX <- filtered.esetALL.FC.final.R[!fData(filtered.esetALL.FC.final.R)$dataType %in% "GX",pData(filtered.esetALL.FC.final.R)$Visit %in% "V2"]
filtered.esetALL.FC.final.NR.V2.noGX <- filtered.esetALL.FC.final.NR[!fData(filtered.esetALL.FC.final.NR)$dataType %in% "GX",pData(filtered.esetALL.FC.final.NR)$Visit %in% "V2"]

wilcox.p.summary.V2 <- do.call('rbind', lapply(seq(nrow(exprs(filtered.esetALL.FC.final.R.V2.noGX))), function(i)  {
  res.p <- wilcox.test(exprs(filtered.esetALL.FC.final.R.V2.noGX)[i,], exprs(filtered.esetALL.FC.final.NR.V2.noGX)[i,])$p.value
  df <- data.frame(featureID=row.names(exprs(filtered.esetALL.FC.final.R.V2.noGX))[i], wilcox.p=res.p, stringsAsFactors = F)
  return(df)
}))


disruption.count.V2$wilcox.p <- wilcox.p.summary.V2[match(disruption.count.V2$node, wilcox.p.summary.V2$featureID),'wilcox.p']
disruption.count.V2$dataType <- fData(filtered.esetALL)[match(disruption.count.V2$node, fData(filtered.esetALL)$featureID),'dataType']

disruption.count.V2.final <- do.call('rbind', lapply(unique(disruption.count.V2$dataType), function(cur.dt) {
  m <- disruption.count.V2[disruption.count.V2$dataType %in% cur.dt,]
  m$wilcox.FDR <- p.adjust(m$wilcox.p, method = "BH", n = length(m$wilcox.p))
  return(m)
}))


disruption.count.V2.final$percentage.max.patient.count.for.disruption.V2 <- disruption.count.V2$max.patient.count.for.disruption.V2/9*100
disruption.count.V2.final$sig.color <- ifelse(disruption.count.V2.final$percentage.max.patient.count.for.disruption.V2>=(6/9*100) & disruption.count.V2.final$wilcox.p<=0.1, "darkseagreen", 
                                              ifelse(disruption.count.V2.final$percentage.max.patient.count.for.disruption.V2>=(6/9*100) & disruption.count.V2.final$wilcox.p>0.1, "indianred3", 
                                                     ifelse(disruption.count.V2.final$percentage.max.patient.count.for.disruption.V2<(6/9*100) & disruption.count.V2.final$wilcox.p<=0.1, "darkorange3", "white")))

disruption.count.V2.final$sig.color <- ifelse(disruption.count.V2.final$sig.color %in% "darkorange3" & disruption.count.V2.final$wilcox.FDR<0.1, "orange", ifelse(disruption.count.V2.final$sig.color %in% "darkorange3" & disruption.count.V2.final$wilcox.FDR>=0.1, "bisque",  disruption.count.V2.final$sig.color))

disruption.count.V2.final$log.wilcoxon.p <- -log10(disruption.count.V2.final$wilcox.p)

#Plot of V2 for mean drop intensity
disruption.count.V2.final$sig.color <- ifelse(disruption.count.V2.final$mean.drop.intensity.V2<=quantile(disruption.count.V2$mean.drop.intensity.V2, 0.1) & disruption.count.V2.final$wilcox.FDR<=0.15, "darkseagreen", 
                                              ifelse(disruption.count.V2.final$mean.drop.intensity.V2<=quantile(disruption.count.V2$mean.drop.intensity.V2, 0.1) & disruption.count.V2.final$wilcox.FDR>0.15, "indianred3", 
                                                     ifelse(disruption.count.V2.final$mean.drop.intensity.V2>=quantile(disruption.count.V2$mean.drop.intensity.V2, 0.1) & disruption.count.V2.final$wilcox.FDR<=0.15, "darkorange3", "grey")))


ggplot(disruption.count.V2.final, aes(x=mean.drop.intensity.V2, y=log.wilcoxon.p, fill=sig.color, colour=sig.color))+
  geom_point(size=2)+
  scale_color_manual(values=c("darkorange3"="orange","bisque"="bisque", "grey"="grey", "indianred3"="indianred3","darkseagreen"="darkseagreen"))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_blank(), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("-log10(p.value)")+xlab("mean drop intensity")+
  geom_hline(yintercept=-log10(0.0001), color = "grey")+
  geom_vline(xintercept=quantile(disruption.count.V2$mean.drop.intensity.V2, 0.1), color = "grey")


disruption.count.V2.final$disrupted.edge.sig.color <- ifelse(disruption.count.V2.final$disrupted.edge.ratio.V2>=quantile(disruption.count.V2$disrupted.edge.ratio.V2, 0.9) & disruption.count.V2.final$wilcox.FDR<=0.1, "darkseagreen", 
                                              ifelse(disruption.count.V2.final$disrupted.edge.ratio.V2>=quantile(disruption.count.V2$disrupted.edge.ratio.V2, 0.9) & disruption.count.V2.final$wilcox.FDR>0.1, "indianred3", 
                                                     ifelse(disruption.count.V2.final$disrupted.edge.ratio.V2>=quantile(disruption.count.V2$disrupted.edge.ratio.V2, 0.9) & disruption.count.V2.final$wilcox.FDR<=0.1, "darkorange3", "grey")))


ggplot(disruption.count.V2.final, aes(x=disrupted.edge.ratio.V2, y=log.wilcoxon.p, colour=disrupted.edge.sig.color))+
  geom_point(size=2)+
  scale_color_manual(values=c("orange"="orange","bisque"="bisque", "grey"="grey", "indianred3"="indianred3","darkseagreen"="darkseagreen"))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("-log10(p.value)")+xlab("Disrupted edge ratio")+
  geom_hline(yintercept=-log10(0.0001), color = "grey")+
  geom_vline(xintercept=quantile(disruption.count.V2$disrupted.edge.ratio.V2, 0.9), color = "grey")


disruption.count.V2.final$sig.color.both <- ifelse(disruption.count.V2.final$mean.drop.intensity.V2<=quantile(disruption.count.V2$mean.drop.intensity.V2, 0.1) & disruption.count.V2.final$disrupted.edge.ratio.V2>=quantile(disruption.count.V2$disrupted.edge.ratio.V2, 0.9), "darkseagreen", 
                                              ifelse(disruption.count.V2.final$mean.drop.intensity.V2<=quantile(disruption.count.V2$mean.drop.intensity.V2, 0.1) & disruption.count.V2.final$disrupted.edge.ratio.V2<quantile(disruption.count.V2$disrupted.edge.ratio.V2, 0.9), "indianred3", 
                                                     ifelse(disruption.count.V2.final$mean.drop.intensity.V2>=quantile(disruption.count.V2$mean.drop.intensity.V2, 0.1) & disruption.count.V2.final$disrupted.edge.ratio.V2>=quantile(disruption.count.V2$disrupted.edge.ratio.V2, 0.9), "darkorange3", "grey")))


disruption.count.V2.final$sig.wilcox <- ifelse(disruption.count.V2.final$wilcox.FDR<0.1, "black", disruption.count.V2.final$sig.color.both)

ggplot(disruption.count.V2.final, aes(x=disrupted.edge.ratio.V2, y=mean.drop.intensity.V2, fill=sig.color.both, color=sig.wilcox,  stroke = 2))+
  geom_point(shape=21, alpha=0.8, aes(size=sig.wilcox))+
  scale_color_manual(values=c("darkorange3"="orange","bisque"="bisque", "grey"="grey", "indianred3"="indianred3","darkseagreen"="darkseagreen", "black"="black"))+
  scale_fill_manual(values=c("darkorange3"="orange","bisque"="bisque", "grey"="grey", "indianred3"="indianred3","darkseagreen"="darkseagreen"))+
 # scale_size_manual(values=c("black"=8, "grey"=2))+
  scale_size_manual(values=c("darkorange3"=2,"bisque"=2, "grey"=2, "indianred3"=2,"darkseagreen"=2, "black"=5))+
  theme_bw()+theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),                                                                  
                   strip.background = element_blank(),
                   panel.border = element_rect(colour = "black"), 
                   axis.ticks.length=unit(.15, "cm"))+
  ylab("Mean drop intensity)")+xlab("Disrupted edge ratio")+
  geom_vline(xintercept=quantile(disruption.count.V2$disrupted.edge.ratio.V2, 0.9), color = "grey")+
  geom_hline(yintercept=quantile(disruption.count.V2$mean.drop.intensity.V2, 0.1), color = "grey")


disruption.count.V2.final$displayName <- fData(filtered.esetALL)[match(disruption.count.V2.final$node, fData(filtered.esetALL)$featureID),'displayName.final']


#plot disruption by shared disrupted modules according to max percentile (0.8)
#Figure 3c & Supp. Fig 4c

#by edge
GSEA.Module.disruption.measurment.results.wide <- reshape2::dcast(GSEA.Module.disruption.measurment.results, patient+group ~ cluster, value.var="percentage.disrupt.in.module")
GSEA.Module.disruption.measurment.results.wide <- GSEA.Module.disruption.measurment.results.wide[,colnames(GSEA.Module.disruption.measurment.results.wide) %in% c("patient", "group", connected.path.list.final)]

anno.df = data.frame(group = pData(filtered.esetALL)[match(GSEA.Module.disruption.measurment.results.wide$patient, pData(filtered.esetALL)$sampleID),'group'],
                     stringsAsFactors = F)

anno.df$group <- ifelse(anno.df$group %in% "IM", "NR", anno.df$group)
anno.df$group <- as.character(anno.df$group)

row.names(GSEA.Module.disruption.measurment.results.wide) <- paste(GSEA.Module.disruption.measurment.results.wide$group,GSEA.Module.disruption.measurment.results.wide$patient, sep=".")
GSEA.Module.disruption.measurment.results.wide$patient <- NULL
GSEA.Module.disruption.measurment.results.wide$group <- NULL


GSEA.Module.disruption.measurment.results.wide <- data.frame(t(GSEA.Module.disruption.measurment.results.wide), stringsAsFactors = FALSE)
GSEA.Module.disruption.measurment.results.wide.final <- apply(GSEA.Module.disruption.measurment.results.wide,2, as.numeric)
row.names(GSEA.Module.disruption.measurment.results.wide.final) <- row.names(GSEA.Module.disruption.measurment.results.wide)
cur.threshold <- quantile(GSEA.Module.disruption.measurment.results.wide.final, selected.percentile)
GSEA.Module.disruption.measurment.results.wide.final <- GSEA.Module.disruption.measurment.results.wide.final[row.names(GSEA.Module.disruption.measurment.results.wide.final) %in% dynamic.range.venn.diagram.shared.pathways[dynamic.range.venn.diagram.shared.pathways$percentile==selected.percentile,'shared.pathways'],]

ha <-   HeatmapAnnotation(df =anno.df ,
                          col = list(group =c("R" = "cyan3", "NR" = "red")))


Heatmap(GSEA.Module.disruption.measurment.results.wide.final, name="percentage of disrupted edges", row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, show_column_names = TRUE, top_annotation = ha) 

for.num.disrupted.NR.edge <- apply(GSEA.Module.disruption.measurment.results.wide.final[,grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/8*100)
for.num.disrupted.R.edge <- apply(GSEA.Module.disruption.measurment.results.wide.final[,!grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/15*100)

ha1 = rowAnnotation(patients.count = anno_lines(cbind(for.num.disrupted.NR.edge, for.num.disrupted.R.edge),  
                                                gp = gpar(col = c("red", "cyan3", fontsize = 10)), add_points = TRUE, pt_gp = gpar(col = c("red", "cyan3")), pch = c(16)), show_annotation_name = FALSE)

col_fun =c("grey", "#f0e3ab", "#f49d77", "#ce5f5f", "#b23030")

ht.edge <- Heatmap(GSEA.Module.disruption.measurment.results.wide.final,
                   rect_gp = gpar(col = "black", lwd = 2), row_title = NULL,
                   name="percentage of disrupted edges", row_names_gp = gpar(fontsize = 8), 
                   top_annotation = ha,
                   column_names_gp = gpar(fontsize = 10), show_row_names = TRUE, show_column_names = T,
                   col =col_fun,
                   column_names_max_height = unit(10,"cm"), right_annotation = ha1,
                   clustering_method_row="average", clustering_distance_row="euclidean",
                   clustering_method_columns="average", clustering_distance_columns = "euclidean")

#by mean drop intensity

GSEA.Module.disruption.measurment.results.wide <- reshape2::dcast(GSEA.Module.disruption.measurment.results, patient+group ~ cluster, value.var="mean.drop.intensity")
GSEA.Module.disruption.measurment.results.wide <- GSEA.Module.disruption.measurment.results.wide[,colnames(GSEA.Module.disruption.measurment.results.wide) %in% c("patient", "group", connected.path.list.final)]

anno.df = data.frame(group = pData(filtered.esetALL)[match(GSEA.Module.disruption.measurment.results.wide$patient, pData(filtered.esetALL)$sampleID),'group'],
                     stringsAsFactors = F)

anno.df$group <- ifelse(anno.df$group %in% "IM", "NR", anno.df$group)
anno.df$group <- as.character(anno.df$group)

row.names(GSEA.Module.disruption.measurment.results.wide) <- paste(GSEA.Module.disruption.measurment.results.wide$group,GSEA.Module.disruption.measurment.results.wide$patient, sep=".")
GSEA.Module.disruption.measurment.results.wide$patient <- NULL
GSEA.Module.disruption.measurment.results.wide$group <- NULL


GSEA.Module.disruption.measurment.results.wide <- data.frame(t(GSEA.Module.disruption.measurment.results.wide), stringsAsFactors = FALSE)
GSEA.Module.disruption.measurment.results.wide.final <- apply(GSEA.Module.disruption.measurment.results.wide,2, as.numeric)
row.names(GSEA.Module.disruption.measurment.results.wide.final) <- row.names(GSEA.Module.disruption.measurment.results.wide)
cur.threshold <- quantile(GSEA.Module.disruption.measurment.results.wide.final, selected.percentile)
GSEA.Module.disruption.measurment.results.wide.final <- GSEA.Module.disruption.measurment.results.wide.final[row.names(GSEA.Module.disruption.measurment.results.wide.final) %in% dynamic.range.venn.diagram.shared.pathways[dynamic.range.venn.diagram.shared.pathways$percentile %in% selected.percentile,'shared.pathways' ],]

ha <-   HeatmapAnnotation(df =anno.df ,col = list(group =c("R" = "cyan3", "NR" = "red")))


Heatmap(GSEA.Module.disruption.measurment.results.wide.final, name="percentage of disrupted edges", row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, show_column_names = TRUE, top_annotation = ha) 

for.num.disrupted.NR <- apply(GSEA.Module.disruption.measurment.results.wide.final[,grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/8*100)
for.num.disrupted.R <- apply(GSEA.Module.disruption.measurment.results.wide.final[,!grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/15*100)

ha1 = rowAnnotation(patients.count = anno_lines(cbind(for.num.disrupted.NR, for.num.disrupted.R),  
                                                gp = gpar(col = c("red", "cyan3", fontsize = 10)), add_points = TRUE, pt_gp = gpar(col = c("red", "cyan3")), pch = c(16)), show_annotation_name = FALSE)


col_fun =c("grey", "#f0e3ab", "#f49d77", "#ce5f5f", "#b23030")

ht.mean.drop <- Heatmap(GSEA.Module.disruption.measurment.results.wide.final,
                        rect_gp = gpar(col = "black", lwd = 2), row_title = NULL,
                        name="mean drop intensity", row_names_gp = gpar(fontsize = 10), top_annotation = ha,
                        column_names_gp = gpar(fontsize = 10), show_row_names = TRUE, show_column_names = T,
                        #column_names_rot = 30, 
                        col =col_fun,
                        column_names_max_height = unit(8,"cm"), right_annotation = ha1,
                        clustering_method_row="average", clustering_distance_row="euclidean",
                        clustering_method_columns="average", clustering_distance_columns = "euclidean")


##By node

GSEA.Module.disruption.measurment.results.wide <- reshape2::dcast(GSEA.Module.disruption.measurment.results, patient+group ~ cluster, value.var="percentage.of.disrupted.nodes")
GSEA.Module.disruption.measurment.results.wide <- GSEA.Module.disruption.measurment.results.wide[,colnames(GSEA.Module.disruption.measurment.results.wide) %in% c("patient", "group", connected.path.list.final)]

anno.df = data.frame(group = pData(filtered.esetALL)[match(GSEA.Module.disruption.measurment.results.wide$patient, pData(filtered.esetALL)$sampleID),'group'],
                     stringsAsFactors = F)

anno.df$group <- ifelse(anno.df$group %in% "IM", "NR", anno.df$group)
anno.df$group <- as.character(anno.df$group)

row.names(GSEA.Module.disruption.measurment.results.wide) <- paste(GSEA.Module.disruption.measurment.results.wide$group,GSEA.Module.disruption.measurment.results.wide$patient, sep=".")
GSEA.Module.disruption.measurment.results.wide$patient <- NULL
GSEA.Module.disruption.measurment.results.wide$group <- NULL


GSEA.Module.disruption.measurment.results.wide <- data.frame(t(GSEA.Module.disruption.measurment.results.wide), stringsAsFactors = FALSE)
GSEA.Module.disruption.measurment.results.wide.final <- apply(GSEA.Module.disruption.measurment.results.wide,2, as.numeric)
row.names(GSEA.Module.disruption.measurment.results.wide.final) <- row.names(GSEA.Module.disruption.measurment.results.wide)
cur.threshold <- quantile(GSEA.Module.disruption.measurment.results.wide.final, selected.percentile)

GSEA.Module.disruption.measurment.results.wide.final <- GSEA.Module.disruption.measurment.results.wide.final[row.names(GSEA.Module.disruption.measurment.results.wide.final) %in% dynamic.range.venn.diagram.shared.pathways[dynamic.range.venn.diagram.shared.pathways$percentile %in% selected.percentile,'shared.pathways' ],]

ha <-   HeatmapAnnotation(df =anno.df ,col = list(group =c("R" = "cyan3", "NR" = "red")))


Heatmap(GSEA.Module.disruption.measurment.results.wide.final, name="percentage of disrupted edges", row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, show_column_names = TRUE, top_annotation = ha) 

for.num.disrupted.NR <- apply(GSEA.Module.disruption.measurment.results.wide.final[,grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/8*100)
for.num.disrupted.R <- apply(GSEA.Module.disruption.measurment.results.wide.final[,!grepl("NR", colnames(GSEA.Module.disruption.measurment.results.wide.final))], 1, function(x) sum(x>cur.threshold)/15*100)

ha1 = rowAnnotation(patients.count = anno_lines(cbind(for.num.disrupted.NR, for.num.disrupted.R),  
                                                gp = gpar(col = c("red", "cyan3", fontsize = 10)), add_points = TRUE, pt_gp = gpar(col = c("red", "cyan3")), pch = c(16)), show_annotation_name = FALSE)


ht.mean.node <- Heatmap(GSEA.Module.disruption.measurment.results.wide.final,
                        rect_gp = gpar(col = "black", lwd = 2), row_title = NULL,
                        name="Percentage of disrupted nodes", row_names_gp = gpar(fontsize = 10), top_annotation = ha,
                        column_names_gp = gpar(fontsize = 10), show_row_names = TRUE, show_column_names = T,
                        col =col_fun,
                        column_names_max_height = unit(8,"cm"), right_annotation = ha1,
                        clustering_method_row="average", clustering_distance_row="euclidean",
                        clustering_method_columns="average", clustering_distance_columns = "euclidean")

draw(ht.edge+ht.mean.drop+ht.mean.node, padding = unit(c(2, 2, 2, 100), "mm"), heatmap_legend_side = "bottom", annotation_legend_side="bottom") #bottom, left, top, right paddings


#Extract disrupted meta-module sub-network based on the disrupted pathways
#Generate cytoscape object for network visualization
#Plot Fig 3e


FCV1V2.disrupted.pathways <- shared.disrupted.pathways
drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted <- drop_matrix.R.NR.V2.for.classification[grepl(paste(FCV1V2.disrupted.pathways, collapse="|"), drop_matrix.R.NR.V2.for.classification$edge.annotation),]
drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted$Var1displayName <- fData(filtered.esetALL)[match(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted$Var1, fData(filtered.esetALL)$featureID),'displayName.final']
drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted$Var2displayName <- fData(filtered.esetALL)[match(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted$Var2, fData(filtered.esetALL)$featureID),'displayName.final']

FCV1V2.disrupted.nodes <- unique(c(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted$Var1, drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted$Var2))
SYMBOL <- fData(filtered.esetALL)[match(FCV1V2.disrupted.nodes, fData(filtered.esetALL)$featureID),'SYMBOL']

filtered.esetALL.FC.final.gx <- filtered.esetALL.FC.final[fData(filtered.esetALL.FC.final)$dataType %in% "GX" & fData(filtered.esetALL.FC.final)$featureID %in% sub("AG", "GX", unique(c(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted$Var1, drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted$Var2))),]

filtered.esetALL.abundance<-filtered.esetALL[grep("Abundance|CyTOF-PBMC|CyTOF-GRAN|GR.TOTAL|PB.TOTAL", fData(filtered.esetALL)$displayName.final),]
filtered.esetALL.abundance.FC <- substractReference(filtered.esetALL.abundance, 'Patient.code', 'Visit', ref = 'V1', op = '/')
exprs(filtered.esetALL.abundance.FC)[is.infinite(exprs(filtered.esetALL.abundance.FC)) | is.nan(exprs(filtered.esetALL.abundance.FC)) ] <- NA

#craete additional column symbabun with abundance or gene symbol 
#for use later to aggregate rows for network visualization
fData(filtered.esetALL.abundance.FC)$symbabun<-fData(filtered.esetALL.abundance.FC)$SYMBOL
fData(filtered.esetALL.abundance.FC)$symbabun[grep('CC',fData(filtered.esetALL.abundance.FC)$dataType)]<-fData(filtered.esetALL.abundance.FC)$displayName[grep('CC',fData(filtered.esetALL.abundance.FC)$dataType)]
fData(filtered.esetALL.abundance.FC)$displayName.final<-as.character(fData(filtered.esetALL.abundance.FC)$displayName)
fData(filtered.esetALL.abundance.FC)$displayName.final <- ifelse(is.na(fData(filtered.esetALL.abundance.FC)$displayName.final), fData(filtered.esetALL.abundance.FC)$symbabun, fData(filtered.esetALL.abundance.FC)$displayName.final)

dat_abundance <- data.frame(exprs(filtered.esetALL.abundance.FC), stringsAsFactors = FALSE)
dat_abundance$displayName.final <- fData(filtered.esetALL.abundance.FC)[match(row.names(dat_abundance),fData(filtered.esetALL.abundance.FC)$featureID),'displayName.final']
dat_abundance$clust <- fData(filtered.esetALL.abundance.FC)[match(row.names(dat_abundance),fData(filtered.esetALL.abundance.FC)$featureID),'clusterID']
dat_abundance[dat_abundance$displayName.final %in% c("CyTOF-PBMC","CyTOF-GRAN"),'clust'] <- 1080000

dat_abundance.merged <-  dat_abundance %>% 
  group_by(displayName.final) %>% 
  dplyr::slice(which.max(clust))

dat_abundance.merged <- data.frame(dat_abundance.merged, stringsAsFactors = F)
row.names(dat_abundance.merged) <- dat_abundance.merged$displayName.final

dat_abundance.merged <- dat_abundance.merged[,!colnames(dat_abundance.merged) %in% c("displayName.final", "clust")]
dat_abundance.merged <- dat_abundance.merged[rowSums(dat_abundance.merged)!=0,]

dat_abundance.merged <- dat_abundance.merged[!row.names(dat_abundance.merged) %in% c("Granulocytes-Abundance", "T cells-Abundance"),]
dat_abundance.merged <- dat_abundance.merged[!row.names(dat_abundance.merged) %in% c("CyTOF-PBMC"),]
row.names(dat_abundance.merged)[row.names(dat_abundance.merged) %in% 'CyTOF-GRAN'] <- "Granulocyte-Abundance"

filtered.esetALL.FC.final.gx.exprs <- data.frame(exprs(filtered.esetALL.FC.final.gx), stringsAsFactors = F)
if(identical(colnames(filtered.esetALL.FC.final.gx.exprs), colnames(dat_abundance.merged))) {
  cor.res <- corr.test(t(filtered.esetALL.FC.final.gx.exprs), t(dat_abundance.merged), method="spearman")
  cor.res.r.long <- reshape2::melt(cor.res$r)
  cor.res.p.long <- reshape2::melt(cor.res$p)
  cor.res.df <- merge(cor.res.r.long, cor.res.p.long, by=c("Var1", "Var2"))
  colnames(cor.res.df) <- c("Var1", "Var2", "cor", "cor.p")
  cor.res.df$Var1displayName <- fData(filtered.esetALL)[match(cor.res.df$Var1, fData(filtered.esetALL)$featureID),'displayName.final']
  cor.res.df$Var2displayName <- cor.res.df$Var2
  cor.res.df.sub <- cor.res.df[abs(cor.res.df$cor)>0.5,]
  cor.res.df.sub$FDR <- p.adjust(cor.res.df.sub$cor.p, method="BH")
  cor.res.df.sub <- cor.res.df.sub[cor.res.df.sub$FDR<0.05,]
} else { print("error: not identical")}

cor.res.df.sub$Var1displayName <- as.character(cor.res.df.sub$Var1displayName)
cor.res.df.sub$Var2displayName <- as.character(cor.res.df.sub$Var2displayName)

drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final <- rbind.fill(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted, cor.res.df.sub)

patient.NR.V2 <- pData(filtered.esetALL)[pData(filtered.esetALL)$Visit %in% "V2" & !pData(filtered.esetALL)$group %in% "R", 'sampleID']
patient.NR.V2 <- gsub("-", "\\.", patient.NR.V2)

for.disruption.count.by.edge <- drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final[,colnames(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final) %in% c("Var1", "Var2", patient.NR.V2)]
for.disruption.count.by.edge$mean.drop <- apply(for.disruption.count.by.edge[,!colnames(for.disruption.count.by.edge) %in% c("Var1", "Var2")], 1, function(x) {mean(x)})
for.disruption.count.by.edge$num.patients <- apply(for.disruption.count.by.edge[,!colnames(for.disruption.count.by.edge) %in% c("Var1", "Var2")], 1, function(x) { sum(x<0)/length(x)*100})

drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final <- merge(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final, for.disruption.count.by.edge[, colnames(for.disruption.count.by.edge) %in% c("Var1", "Var2", "mean.drop", "num.patients")], by=c("Var1", "Var2"), all.x=T)

drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final[drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final$edge.annotation %in% "NULL",'edge.annotation'] <- "Abundance"

drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape <- do.call('rbind', lapply(c(as.character(shared.disrupted.pathways), "Abundance"), function(cur.path) {
  print(cur.path)
  sub <- drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final[grepl(cur.path, drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final$edge.annotation),]
  sub$path.map <- cur.path
  return(sub)
}))

colnames(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape)

edge<-drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape
node<-unique(c(edge$Var1displayName, edge$Var2displayName))

d <- data.frame(Var1=edge$Var1displayName,Var2=edge$Var2displayName, stringsAsFactors = F)
net <- graph.data.frame(d, vertices=node, directed=FALSE)

betweeness<- betweenness(net, normalize=T)
betweeness <- data.frame(feature=names(betweeness), betweeness=betweeness, stringsAsFactors = F)

degree<- igraph::degree(net, normalize=T)
degree <- data.frame(feature=names(degree), degree=degree, stringsAsFactors = F)

feature.betweenness.degree <- merge(degree, betweeness, by='feature')
feature.betweenness.degree$combined <- feature.betweenness.degree$degree+feature.betweenness.degree$betweeness

drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$degree.Var1 <- feature.betweenness.degree[match(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var1displayName, feature.betweenness.degree$feature),'degree']
drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$degree.Var2 <- feature.betweenness.degree[match(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var2displayName, feature.betweenness.degree$feature),'degree']
drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$betweenness.Var1 <- feature.betweenness.degree[match(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var1displayName, feature.betweenness.degree$feature),'betweeness']
drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$betweenness.Var2 <- feature.betweenness.degree[match(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var2displayName, feature.betweenness.degree$feature),'betweeness']

#Test enrichment

res.enrich.kegg.FCV1V2.disrupted <- res.enrich.res.function(unique(c(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var1, drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var2)), "KEGG", 1)
res.enrich.kegg.FCV1V2.disrupted.top <- head(res.enrich.kegg.FCV1V2.disrupted, n=5)
res.enrich.kegg.FCV1V2.disrupted.top$set <- "KEGG"

res.enrich.go.FCV1V2.disrupted <- res.enrich.res.function(unique(c(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var1, drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var2)), "GO", 0.05)
res.enrich.go.FCV1V2.disrupted.top <- head(res.enrich.go.FCV1V2.disrupted, n=5)
res.enrich.go.FCV1V2.disrupted.top$set <- "GO"

res.enrich.Reactome.FCV1V2.disrupted <- res.enrich.res.function(unique(c(drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var1, drop_matrix.R.NR.V2.for.classification.FCV1V2.disrupted.final.for.cytoscape$Var2)), "ReactomePA", 1)
res.enrich.Reactome.FCV1V2.disrupted.top <- head(res.enrich.Reactome.FCV1V2.disrupted, n=5)
res.enrich.Reactome.FCV1V2.disrupted.top$set <- "Reactome"

res.enrich.FCV1V2.disrupted.combined <- rbind(res.enrich.kegg.FCV1V2.disrupted.top, res.enrich.go.FCV1V2.disrupted.top, res.enrich.Reactome.FCV1V2.disrupted.top)

res.enrich.FCV1V2.disrupted.combined[res.enrich.FCV1V2.disrupted.combined$ID %in% "GO:0110053",'Description'] <- "regulation of actin filament organization"


res.enrich.FCV1V2.disrupted.combined$log.qValue <- -log10(res.enrich.FCV1V2.disrupted.combined$qvalue)
res.enrich.FCV1V2.disrupted.combined$set <- factor(res.enrich.FCV1V2.disrupted.combined$set, levels=c("Reactome","KEGG", "GO"))
res.enrich.FCV1V2.disrupted.combined <- res.enrich.FCV1V2.disrupted.combined[order(res.enrich.FCV1V2.disrupted.combined$set,res.enrich.FCV1V2.disrupted.combined$log.qValue, decreasing=T),]
res.enrich.FCV1V2.disrupted.combined$Description <- factor(res.enrich.FCV1V2.disrupted.combined$Description, levels=rev(res.enrich.FCV1V2.disrupted.combined$Description))

ggplot(res.enrich.FCV1V2.disrupted.combined,aes(x=Description, y=log.qValue, fill=set)) + 
  ylab("-log10 (q-value)")+
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values=c("GO" = "#b23030", "KEGG" = "#f49d77","Reactome" = "#f0e3ab"))+
  theme_bw()+
  theme(strip.text = element_text(size=14), panel.grid.major = element_blank(),
        axis.text.y=element_text( size=12),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border =element_blank(), 
        axis.ticks.length=unit(.15, "cm"))+
  coord_flip()



#Compare core gene centrality with equally-sized node-set of random combinations

feature.betweenness.degree.no.abundance <- feature.betweenness.degree[!grepl("Abundance|CyTOF-GRAN", feature.betweenness.degree$feature),]

real.mean.betweenness <- mean(feature.betweenness.degree.no.abundance[feature.betweenness.degree.no.abundance$feature %in% c("RAC1", "PAK1", "HCK"), 'combined'])

perm.centrality <- do.call('rbind', lapply(1:10000, function(i) {
  samp <- feature.betweenness.degree.no.abundance[sample(1:length(feature.betweenness.degree.no.abundance$combined), 4, replace=T),'combined']
  return(mean(samp))
}))

mean(perm.centrality>real.mean.betweenness)  

#percentile of monocytes in terms of degree and betweenness

sum(feature.betweenness.degree$degree <= feature.betweenness.degree$degree[feature.betweenness.degree$feature %in% "Monocytes-Abundance"])/length(feature.betweenness.degree$degree)
sum(feature.betweenness.degree$betweeness <= feature.betweenness.degree$betweeness[feature.betweenness.degree$feature %in% "Monocytes-Abundance"])/length(feature.betweenness.degree$betweeness)



#Test degree and betweeness of the different nodes and their 
#relationship to disruption
#Plot Fig. 3d

edge<-cor.edges.df.FC.combined.aggregated
node<-unique(c(edge$Var1, edge$Var2))

d <- data.frame(Var1=edge$Var1,Var2=edge$Var2, stringsAsFactors = F)
net <- graph.data.frame(d, vertices=node, directed=FALSE)

betweeness<- betweenness(net, normalize=T)

#node list of the shared disrupted modules

node.list.df <-  do.call('rbind', lapply(seq(length(shared.disrupted.pathways)), function(cur.path) {
  node.list <-net.nodes.by.pathways[net.nodes.by.pathways$path %in% shared.disrupted.pathways[cur.path], 'featureID']
  node.list <- unique(node.list)
  return(data.frame(path=shared.disrupted.pathways[cur.path], node=node.list, stringsAsFactors = F))
}))

betweeness <- data.frame(featureID=names(betweeness), betweeness=betweeness, stringsAsFactors = F)
betweeness$betweeness.mod <- 15-(-log10(ifelse(betweeness$betweeness==0, 0.00000000000001, betweeness$betweeness)))
betweeness$disruption.status <- ifelse(betweeness$featureID %in% node.list.df$node, "Nodes in disrupted modules", "Other")
betweeness$disruption.status <- factor(betweeness$disruption.status)

ggplot(betweeness, aes(x = betweeness))+
  stat_density(aes(fill = disruption.status, color="black"), alpha = 0.4, position="identity", trim=T) +
  facet_zoom(xlim = c(0, 0.01))+
  scale_color_manual(values=c("black"="black", "Other"="cyan3", "Nodes in disrupted modules"="indianred"))+
  geom_vline(data=betweeness[betweeness$disruption.status %in% "Nodes in disrupted modules",], aes(xintercept = median(betweeness), color=disruption.status), linetype = "dashed", size=1)+
  geom_vline(data=betweeness[betweeness$disruption.status %in% "Other",], aes(xintercept = median(betweeness), color = disruption.status), linetype = "dashed", size=1)+
  theme_bw()+
  theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))


#Test significance using permutation test:

set.seed(102) 
nsim <- 10000
res <- numeric(nsim) 
for (i in 1:nsim) {
  perm <- sample(nrow(betweeness))
  bdat <- transform(betweeness,betweeness=betweeness[perm])
  res[i] <- mean(bdat[bdat$disruption.status=="Nodes in disrupted modules","betweeness"])-
    mean(bdat[bdat$disruption.status=="Other","betweeness"])
}
obs <- mean(betweeness[betweeness$disruption.status=="Nodes in disrupted modules","betweeness"])-
  mean(betweeness[betweeness$disruption.status=="Other","betweeness"])
res <- c(res,obs)
hist(res,col=muted("red"),las=1,main="")
abline(v=obs,col="red")

mean(res>=obs)         

#For degree:

degree<- igraph::degree(net)

degree <- data.frame(featureID=names(degree), degree=degree, stringsAsFactors = F)
degree$disruption.status <- ifelse(degree$featureID %in% node.list.df$node, "Nodes in disrupted modules", "Other")

ggplot(degree, aes(x = degree))+
  geom_density(aes(fill = disruption.status), alpha = 0.4) +
  #xlim(0,0.01)+
  geom_vline(data=degree[degree$disruption.status %in% "Nodes in disrupted modules",], aes(xintercept = median(degree), color=disruption.status), linetype = "dashed", size=1)+
  geom_vline(data=degree[degree$disruption.status %in% "Other",], aes(xintercept = median(degree), color = disruption.status), linetype = "dashed", size=1)+
  theme_bw()+
  theme(strip.text = element_text(size=9), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        panel.grid.minor = element_blank(),                                                                  
        strip.background = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.length=unit(.15, "cm"))

wilcox.test(degree[degree$disruption.status %in% "Other",'degree'], degree[degree$disruption.status %in% "Nodes in disrupted modules",'degree'])

##Test using permutation test:

set.seed(102) 
nsim <- 10000
res <- numeric(nsim) 
for (i in 1:nsim) {
  perm <- sample(nrow(degree))
  bdat <- transform(degree,degree=degree[perm])
  res[i] <- mean(bdat[bdat$disruption.status=="Nodes in disrupted modules","degree"])-
    mean(bdat[bdat$disruption.status=="Other","degree"])
}
obs <- mean(degree[degree$disruption.status=="Nodes in disrupted modules","degree"])-
  mean(bdat[bdat$disruption.status=="Other","degree"])
res <- c(res,obs)
hist(res,col=muted("red"),breaks=100, xlim=c(-11,100))
abline(v=obs,col="red")

mean(res>=obs)    



#Compare differences of disrupted modules at baseline
#Plot Fig. 4

Disrupted.clust.FC.summary.V1.wilcoxon <-  do.call('rbind', lapply(seq(length(shared.disrupted.pathways)), function(cur.path) {
  
  print(shared.disrupted.pathways[cur.path])
  node.list <-net.nodes.by.pathways[net.nodes.by.pathways$path %in% shared.disrupted.pathways[cur.path], 'featureID']
  node.list <- unique(node.list)
  esetMat.eset <- filtered.esetALL[fData(filtered.esetALL)$featureID %in% node.list, !pData(filtered.esetALL)$sampleID %in% "HR-37-V1"]
  expMat <- data.frame(exprs(esetMat.eset), stringsAsFactors = F)
  expMat <- data.frame(t(expMat), stringsAsFactors = F)
  expMat[,!colnames(expMat) %in% "Visit"] <- apply(expMat[,!colnames(expMat) %in% "Visit"], 2, scale)
  row.names(expMat) <- gsub("\\.", "-", row.names(expMat))
  expMat$Visit <- pData(filtered.esetALL)[match(row.names(expMat), pData(filtered.esetALL)$sampleID),'Visit']
  expMat$group <- pData(filtered.esetALL)[match(row.names(expMat), pData(filtered.esetALL)$sampleID),'group']
  expMat$group <- ifelse(expMat$group %in% "IM", "NR", expMat$group)
  print(shared.disrupted.pathways[cur.path])
  expMat$module.score <- apply(expMat[,!colnames(expMat) %in% c("group", "Visit")], 1,  function(x) mean(x, na.rm=T))
  expMat.V1 <- expMat[expMat$Visit %in% "V1",]
  obs.diff <- mean(expMat.V1[expMat.V1$group=="R","module.score"])-
    mean(expMat.V1[expMat.V1$group=="NR","module.score"])
  
  mean.module.score.R <- mean(expMat.V1[expMat.V1$group %in% "R", 'module.score'], na.rm=T)
  mean.module.score.NR <- mean(expMat.V1[expMat.V1$group %in% "NR", 'module.score'], na.rm=T)
  SEM.module.score.R <- sd(expMat.V1[expMat.V1$group %in% "R", 'module.score'], na.rm=T)/sqrt(length(expMat.V1[expMat.V1$group %in% "R", 'module.score']))
  SEM.module.score.NR <- sd(expMat.V1[expMat.V1$group %in% "NR", 'module.score'], na.rm=T)/sqrt(length(expMat.V1[expMat.V1$group %in% "NR", 'module.score']))
  p.value <- wilcox.test(module.score ~ group, paired=F, data=expMat.V1)$p.value 
  p.permanova <- adonis(expMat.V1[,!colnames(expMat.V1) %in% c("group", "Visit", "module.score")] ~ expMat.V1$group, method = "euclidean",
                        permutations = 100000)
  p.permanova.final <- p.permanova$aov.tab$`Pr(>F)`[1]
  df <- data.frame(clust=shared.disrupted.pathways[cur.path], p.valueV1=p.value, obs.diff=obs.diff, 
                   mean.module.score.R=mean.module.score.R, mean.module.score.NR=mean.module.score.NR, SEM.module.score.R=SEM.module.score.R, SEM.module.score.NR=SEM.module.score.NR, 
                   p.permanova=p.permanova.final, stringsAsFactors = F)
  return(df)
}))

Disrupted.clust.FC.summary.V1.wilcoxon$FDR <- p.adjust(Disrupted.clust.FC.summary.V1.wilcoxon$p.permanova, method="BH")



#To prioritize pathways at baseline, test correlation
#to CRP and network centrality

disrupted.pathways.correlation.with.CRP.FC <- do.call('rbind', lapply(shared.disrupted.pathways, function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  esetMat.V2 <- filtered.esetALL.FC.final[fData(filtered.esetALL.FC.final)$featureID %in% node.list.in.module.and.net.edge, as.character(pData(filtered.esetALL.FC.final)$Visit) %in% "V2" & pData(filtered.esetALL.FC.final)$group %in% c("R", "IM", "NR") & !pData(filtered.esetALL.FC.final)$Patient.code %in% "37"]
  esetMat.V2 <- data.frame(exprs(esetMat.V2), stringsAsFactors = F)
  esetMat.V2 <- data.frame(t(esetMat.V2), stringsAsFactors = F)
  esetMat.V2$time <- "V2"
  esetMat.V3 <- filtered.esetALL.FC.final[fData(filtered.esetALL.FC.final)$featureID %in% node.list.in.module.and.net.edge, as.character(pData(filtered.esetALL.FC.final)$Visit) %in% "V3" & pData(filtered.esetALL.FC.final)$group %in% c("R", "IM", "NR") &  !pData(filtered.esetALL.FC.final)$Patient.code %in% "37"]
  esetMat.V3 <- data.frame(exprs(esetMat.V3), stringsAsFactors = F)
  esetMat.V3 <- data.frame(t(esetMat.V3), stringsAsFactors = F)
  esetMat.V3$time <- "V3"
  
  combined.expMat <-esetMat.V3 
  combined.expMat$time <- factor(combined.expMat$time)
  combined.expMat$module.score <- apply(combined.expMat[,!colnames(combined.expMat) %in% "time"], 1,  function(x) mean(x, na.rm=T))
  
  row.names(combined.expMat) <- gsub("\\.", "-", row.names(combined.expMat))
  combined.expMat$CRP.FC <- ifelse(combined.expMat$time %in% "V2",  
                                pData(filtered.esetALL)[match(row.names(combined.expMat), pData(filtered.esetALL)$sampleID),'CRP.FCV1V2'], 
                                pData(filtered.esetALL)[match(row.names(combined.expMat), pData(filtered.esetALL)$sampleID),'CRP.FCV1V3'])
  
  cor.res.all <- corr.test(combined.expMat$module.score, combined.expMat$CRP.FC, method="pearson")
  
  return(data.frame(path=cur.path, r=cor.res.all$r, p=cor.res.all$p,stringsAsFactors=F))
}))



disrupted.pathways.correlation.with.CRP.FC$fdr.r <- p.adjust(disrupted.pathways.correlation.with.CRP.FC$p, method="BH")
Disrupted.clust.FC.summary.V1.wilcoxon.final <- merge(Disrupted.clust.FC.summary.V1.wilcoxon.eigen.genes.full.PC1, disrupted.pathways.correlation.with.CRP.FC, by.x="clust", by.y="path")

Disrupted.clust.FC.summary.V1.wilcoxon.final$logPv1 <- -log10(Disrupted.clust.FC.summary.V1.wilcoxon.final$p.value)

Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final <- merge(Disrupted.clust.FC.summary.V1.wilcoxon, disrupted.pathways.correlation.with.CRP.FC, by.x="clust", by.y="path")


#Test centrality

degree<- igraph::degree(net, normalized = T)
degree.df <- data.frame(featureID=names(degree), degree=degree, stringsAsFactors = F )
degree.df$SYMBOL <- fData(filtered.esetALL)[match(degree.df$featureID, fData(filtered.esetALL)$featureID),'SYMBOL']
degree.path <-  do.call('rbind', lapply(shared.disrupted.pathways, function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  degree.df.sub <- degree.df[degree.df$featureID %in% node.list.in.module.and.net.edge,]
  degree.df.sub <- degree.df.sub[order(degree.df.sub$degree, decreasing=T),]
  df <- data.frame(clust=cur.path, degree=quantile(degree.df.sub$degree, 0.8), stringsAsFactors = F)
  return(df)
}))


betweenness<- igraph::betweenness(net, normalized = T)
betweenness.df <- data.frame(featureID=names(betweenness), betweenness=betweenness, stringsAsFactors = F )
betweenness.df$SYMBOL <- fData(filtered.esetALL)[match(betweenness.df$featureID, fData(filtered.esetALL)$featureID),'SYMBOL']

betweenness.path <-  do.call('rbind', lapply(shared.disrupted.pathways, function(cur.path) {
  node.list.in.module.and.net.edge <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.path, 'featureID']
  betweenness.df.sub <- betweenness.df[betweenness.df$featureID %in% node.list.in.module.and.net.edge,]
  betweenness.df.sub <- betweenness.df.sub[order(betweenness.df.sub$betweenness, decreasing=T),]
  df <- data.frame(clust=cur.path, betweenness=quantile(betweenness.df.sub$betweenness, 0.8), stringsAsFactors = F)
  return(df)
}))

Disrupted.clust.FC.summary.V1.wilcoxon.final <- merge(Disrupted.clust.FC.summary.V1.wilcoxon.final, degree.path, by.x="clust", by.y="clust")
Disrupted.clust.FC.summary.V1.wilcoxon.final$label <- ifelse(Disrupted.clust.FC.summary.V1.wilcoxon.final$clust %in% c("GO_POSITIVE_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION", "GO_POSITIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION", "GO_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION"), Disrupted.clust.FC.summary.V1.wilcoxon.final$clust, "")
Disrupted.clust.FC.summary.V1.wilcoxon.final$label <- ifelse(Disrupted.clust.FC.summary.V1.wilcoxon.final$label %in% "GO_POSITIVE_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION", "Positive regulation of fiber organization", 
                                                             ifelse(Disrupted.clust.FC.summary.V1.wilcoxon.final$label %in% "GO_POSITIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION", "Positive regulation of cytoskeleton organization",
                                                                    ifelse(Disrupted.clust.FC.summary.V1.wilcoxon.final$label %in% "GO_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION", "Regulation of fiber organization","" )))

Disrupted.clust.FC.summary.V1.wilcoxon.final$FDR <- p.adjust(Disrupted.clust.FC.summary.V1.wilcoxon.final$p.value, method="BH")

Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final <- merge(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final, degree.path, by.x="clust", by.y="clust")
Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final <- merge(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final, betweenness.path, by.x="clust", by.y="clust")

Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$log10FDR.permanova <- -log10(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$FDR)
Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final <- Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final[order(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$log10FDR.permanova ,Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$betweenness, decreasing=F),]
Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$clust <- factor(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$clust, levels=Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$clust)

Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$label <- ifelse(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$clust %in% "GO_POSITIVE_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION", "Positive regulation of fiber organization", 
                                                             ifelse(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$clust %in% "GO_POSITIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION", "Positive regulation of cytoskeleton organization",
                                                                    ifelse(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$clust %in% "GO_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION", "Regulation of fiber organization","" )))


bold.labels <- ifelse(levels(Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$clust) %in% c("GO_POSITIVE_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION", "GO_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION", "GO_REGULATION_OF_PROTEIN_COMPLEX_ASSEMBLY", "GO_POSITIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION"), yes = "bold", no = "plain")

ggplot(data=Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final) + 
geom_bar(mapping = aes(x = clust, y = log10FDR.permanova ,  fill = betweenness), stat = "identity") + 
  scale_fill_gradientn(limits = c(0,0.003), colours=c("grey", "#f0e3ab", "#b23030"), breaks=c(0,0.001, 0.002, 0.003), labels=c("0","0.001", "0.002", "0.003"))+
  geom_line(mapping = aes(x = clust, y = r*2.5/0.5, group=1)) + 
  geom_point(mapping = aes(x = clust, y = r*2.5/0.5), size = 3, shape = 21, fill = "white") + 
  geom_hline(yintercept = 1, linetype="dashed", color="black")+
  scale_y_continuous(
    name = expression("NPMANOVA based FDR (-log10) "), 
    sec.axis = sec_axis(~ . * 0.5 / 2.5 , name = "CRP cor (FC)"), 
    limits = c(0, 2.5)) + 
  coord_flip()+
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(face = bold.labels)
  )


#Test baseline differences of the disrupted pathways in the node level
# Plot Supp. Fig 5a

V1.disrupted.pathways <- Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final[Disrupted.clust.FC.summary.V1.wilcoxon.permanova.final$FDR<0.1,'clust']
V1.module.disrupted.nodes <- unique(net.nodes.by.pathways[net.nodes.by.pathways$path %in% as.character(V1.disrupted.pathways), 'featureID'])
V1.module.disrupted.nodes <- as.character(V1.module.disrupted.nodes)
V1.module.disrupted.nodes.exp <- filtered.esetALL[fData(filtered.esetALL)$featureID %in% V1.module.disrupted.nodes, pData(filtered.esetALL)$Visit %in% "V1"]
V1.module.disrupted.nodes.exp.mat <- data.frame(t(exprs(V1.module.disrupted.nodes.exp)), stringsAsFactors = F)
V1.module.disrupted.nodes.exp.mat$group <- pData(filtered.esetALL)[match(row.names(V1.module.disrupted.nodes.exp.mat), pData(filtered.esetALL)$sampleID),'group']
V1.module.disrupted.nodes.exp.mat$group <- ifelse(V1.module.disrupted.nodes.exp.mat$group %in% "IM", "NR", V1.module.disrupted.nodes.exp.mat$group)
V1.module.disrupted.nodes.exp.mat$group <- factor(V1.module.disrupted.nodes.exp.mat$group, levels=c("R", "NR"))

V1.module.disrupted.nodes.exp.mat.long <- reshape2::melt(V1.module.disrupted.nodes.exp.mat, id.vars="group")

V1.module.disrupted.nodes.wilcoxon <- do.call('rbind', lapply(unique(V1.module.disrupted.nodes), function(cur.gene) {
  submat <- V1.module.disrupted.nodes.exp.mat.long[V1.module.disrupted.nodes.exp.mat.long$variable %in% cur.gene,]
  p <- wilcox.test(value~ group, data=submat)$p.value
  return(data.frame(gene=cur.gene, p=p, stringsAsFactors=F))
}))

V1.module.disrupted.nodes.wilcoxon$SYMBOL <- fData(filtered.esetALL)[match(V1.module.disrupted.nodes.wilcoxon$gene, fData(filtered.esetALL)$featureID),'SYMBOL']
V1.module.disrupted.nodes.wilcoxon$SYMBOL <- fData(filtered.esetALL)[match(V1.module.disrupted.nodes.wilcoxon$gene, fData(filtered.esetALL)$featureID),'displayName.final']

V1.module.disrupted.nodes.wilcoxon$FDR <- p.adjust(V1.module.disrupted.nodes.wilcoxon$p, method="BH")

saveRDS(V1.module.disrupted.nodes.wilcoxon, "V1.module.disrupted.nodes.wilcoxon.rds")

esetMat <- filtered.esetALL[fData(filtered.esetALL)$featureID %in% V1.module.disrupted.nodes.wilcoxon[V1.module.disrupted.nodes.wilcoxon$FDR<0.2,'gene'], pData(filtered.esetALL)$Visit %in% "V1"]

esetMat.exprs <- data.frame(t(exprs(esetMat)), stringsAsFactors = F)
esetMat.exprs <- esetMat.exprs[!row.names(esetMat.exprs) %in% 'featureID',] 
esetMat.exprs$group <- pData(filtered.esetALL)[match(row.names(esetMat.exprs), pData(filtered.esetALL)$sampleID),'group']
esetMat.exprs$group <- ifelse( esetMat.exprs$group %in% "IM", "NR",  esetMat.exprs$group)
esetMat.exprs$sampleID <- row.names(esetMat.exprs)


anno.df = data.frame(group = esetMat.exprs$group,
                     stringsAsFactors = F)


ha <-   rowAnnotation(df =anno.df ,
                      col = list(group =c("R" = "cyan3", "NR" = "red")))


colnames(esetMat.exprs)[!colnames(esetMat.exprs) %in% c("sampleID", "group")] <- fData(filtered.esetALL)[match(colnames(esetMat.exprs)[!colnames(esetMat.exprs) %in% c("sampleID", "group")], fData(filtered.esetALL)$featureID),'SYMBOL'] 
esetMat.exprs.scaled <- esetMat.exprs
esetMat.exprs.scaled[,!colnames(esetMat.exprs.scaled) %in% c("sampleID", "group")] <- apply(esetMat.exprs.scaled[,!colnames(esetMat.exprs.scaled) %in% c("sampleID", "group")] , 2, scale)


Heatmap(as.matrix(esetMat.exprs.scaled[,!colnames(esetMat.exprs.scaled) %in% c("sampleID", "group")]), 
        name="Baseline log2 expression", 
        rect_gp = gpar(col = "black", lwd = 2), row_title = NULL,
        col=c( "#084594", "#4292C6", "#9ECAE1","#DEEBF7" ,col_fun),
        row_names_gp = gpar(fontsize = 8), 
        cluster_rows = T,
        right_annotation = ha,
        column_names_gp = gpar(fontsize = 8), show_row_names = TRUE, show_column_names = T,
        #column_names_rot = 30, 
        column_names_max_height = unit(10,"cm"), 
        clustering_method_row="ward.D2", clustering_distance_row="euclidean")



#Test correlation of the fiber organization module unadjusted expression with cell abundances
#Supp. fig 5b

nodes.in.relevant.pathways.sub <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% "GO_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION",]
esetMat <- filtered.esetALL[fData(filtered.esetALL)$featureID %in% gsub("^AG", "GX", nodes.in.relevant.pathways.sub$featureID) & fData(filtered.esetALL)$dataType %in% "GX", ]

esetMat.exprs <- apply(exprs(esetMat),2, as.numeric)
row.names(esetMat.exprs) <- row.names(exprs(esetMat))
esetMat.exprs <- data.frame(esetMat.exprs, stringsAsFactors = F)
esetMat.exprs$featureID <- row.names(esetMat.exprs)

esetMat.disrupted.genes <- esetMat.exprs
esetMat.abundance <- dat.merged.for.pca.not.FC.sub

#do not include redundant cell type
esetMat.abundance <- esetMat.abundance[!row.names(esetMat.abundance) %in% c("T cells-Abundance", "Anti-inflammatory monocytes-Abundance"),]
esetMat.abundance <- esetMat.abundance[,!colnames(esetMat.abundance) %in% c("displayName.final", "clust")]
esetMat.disrupted.genes <- esetMat.disrupted.genes[,!colnames(esetMat.disrupted.genes) %in% "featureID"]

esetMat.disrupted.genes <- data.frame(t(esetMat.disrupted.genes), stringsAsFactors = F)
esetMat.abundance <- data.frame(t(esetMat.abundance), stringsAsFactors = F)

esetMat.disrupted.genes <- esetMat.disrupted.genes[order(row.names(esetMat.disrupted.genes)),]
esetMat.abundance <- esetMat.abundance[order(row.names(esetMat.abundance)),]

esetMat.disrupted.genes$group <- pData(filtered.esetALL)[match(gsub("\\.", "-", row.names(esetMat.disrupted.genes)), pData(filtered.esetALL)$sampleID),'group']
esetMat.abundance$group <- pData(filtered.esetALL)[match(gsub("\\.", "-", row.names(esetMat.abundance)), pData(filtered.esetALL)$sampleID),'group']

esetMat.disrupted.genes.sub <- esetMat.disrupted.genes[esetMat.disrupted.genes$group %in% "R",]
esetMat.abundance.sub <- esetMat.abundance[esetMat.abundance$group %in% "R",]

if(identical(row.names(esetMat.disrupted.genes.sub), row.names(esetMat.abundance.sub))) {
  
  cor.cells.disrupted.genes <- corr.test(esetMat.abundance.sub[,!colnames(esetMat.abundance.sub) %in% "group"], esetMat.disrupted.genes.sub[,!colnames(esetMat.disrupted.genes.sub) %in% "group"], method="spearman", adjust="none")
} else {
  "not identical row.names"
}

cor.cells.disrupted.genes.r.long <- reshape2::melt(cor.cells.disrupted.genes$r)
cor.cells.disrupted.genes.p.long <- reshape2::melt(cor.cells.disrupted.genes$p)
cor.cells.disrupted.genes.df <- cbind(cor.cells.disrupted.genes.r.long, cor.cells.disrupted.genes.p.long$value)
colnames(cor.cells.disrupted.genes.df) <- c("Var1", "Var2", "cor", "adjusted.p")
cor.cells.disrupted.genes.df <- cor.cells.disrupted.genes.df[!is.na(cor.cells.disrupted.genes.df$cor),]
cor.cells.disrupted.genes.df <- cor.cells.disrupted.genes.df[!cor.cells.disrupted.genes.df$cor %in% "NA",]
cor.cells.disrupted.genes.df$displayName.final <- fData(filtered.esetALL)[match(cor.cells.disrupted.genes.df$Var2, fData(filtered.esetALL)$featureID), 'displayName.final']
cor.cells.disrupted.genes.df$cor <- as.numeric(as.character(cor.cells.disrupted.genes.df$cor))
cor.cells.disrupted.genes.df$adjusted.p <- as.numeric(as.character(cor.cells.disrupted.genes.df$adjusted.p))
cor.cells.disrupted.genes.df.sub <- cor.cells.disrupted.genes.df[abs(cor.cells.disrupted.genes.df$cor)>=0.5 & cor.cells.disrupted.genes.df$adjusted.p<0.05,]

#plot correlation matrix

colnames(cor.cells.disrupted.genes$r) <- fData(filtered.esetALL)[match(colnames(cor.cells.disrupted.genes$r), fData(filtered.esetALL)$featureID),'displayName.final']
colnames(cor.cells.disrupted.genes$p) <- fData(filtered.esetALL)[match(colnames(cor.cells.disrupted.genes$p), fData(filtered.esetALL)$featureID),'displayName.final']

corrplot(cor.cells.disrupted.genes$r, method="number")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor.cells.disrupted.genes$r, method="color", col=col(200),  
         order="original", 
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = cor.cells.disrupted.genes$p, sig.level = 0.01, insig = "blank" )


