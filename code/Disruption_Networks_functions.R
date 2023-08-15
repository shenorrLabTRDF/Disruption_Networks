
#Disruption Networks functions

#Generate expression set of Fold Change values using substractReference function

substractReference <- function(object, id, y, ref = NULL, pheno = NULL, op = "-") 
{
  if (dim(object)[2] > 0) {
    if (isString(id)) 
      id <- object[[id]]
    if (isString(y)) 
      y <- object[[y]]
  }
  stopifnot(length(y) == ncol(object))
  stopifnot(length(id) == ncol(object))
  stopifnot(!any(is.na(id)))
  stopifnot(is.null(ref) || isString(ref))
  if (!is.factor(y)) 
    y <- as.factor(y)
  if (!is.null(ref)) 
    y <- relevel(y, ref = ref)
  ref <- levels(y)[1L]
  iref <- which(y %in% ref)
  iref <- setNames(iref, id[iref])
  stopifnot(!anyDuplicated(names(iref)))
  with_ref <- setdiff(which(id %in% names(iref)), iref)
  op <- match.fun(op)
  x <- exprs(object)
  x0 <- matrix(NA, nrow(x), ncol(x), dimnames = dimnames(x))
  lapply(with_ref, function(j) {
    j0 <- iref[id[j]]
    x0[, j] <<- op(x[, j], x[, j0])
    if (!is.null(pheno)) {
      lapply(pheno, function(v) {
        p <- pData(object)[[v]]
        p[j] <- op(p[j], p[j0])
        pData(object)[[v]] <<- p
      })
    }
  })
  object <- object[, with_ref, drop = FALSE]
  exprs(object) <- x0[, with_ref, drop = FALSE]
  object
}


#Compute differential analysis using mixed effect linear models

lsmeansStats <- function(fit, model_stats = NULL){
  
  if( is.null(model_stats) ){
    message("* Computing summary statistics")
    model_stats <- sapply0(fit, lsmeans::lsmeans, pairwise ~ Visit, adjust = 'tukey')
  }
  
  message("* Computing contrasts statistics")
  fit_stats <- lapply(model_stats, function(x){
    stat <- summary(x)
    
    # means
    st <- stat$lsmeans
    rn <- as.character(st[[1L]])
    st <- as.matrix(st[-1L])
    rownames(st) <- paste0(rn, '_0')
    colnames(st)[1L] <- 'estimate'
    m <- st
    
    # contrasts
    st <- stat$contrasts
    rn <- as.character(st[[1L]])
    st <- as.matrix(st[-1L])
    rownames(st) <- rn
    contr <- st
    
    res <- rbind.fill.matrix(m, st)
    rownames(res) <- c(rownames(m), rownames(st))
    res
  })
  
  class(fit_stats) <- 'cytnf_lsmeans'
  fit_stats
}

sapply0 <- function(...) sapply(..., simplify = FALSE)

fitLMER <- function (formula, data, variables = NULL, pattern = NULL, scale = FALSE, 
          fit.only = FALSE, envir = parent.frame()) 
{
  pe <- envir
  if (is(data, "ExpressionSet")) {
    variables <- t(exprs(data))
    data <- pData(data)
  }
  if (!is.null(variables)) 
    vars <- variables
  else if (!is.null(pattern)) 
    vars <- grep(pattern, colnames(data))
  else stop("Must provide either argument `variables` or `pattern`")
  if (!is.null(dim(vars))) {
    data <- cbind(data, vars)
    vars <- colnames(vars)
  }
  else if (is.integer(vars)) 
    vars <- colnames(data)[vars]
  model <- paste0(as.character(formula), collapse = " ")
  message("* Fitting model: ", model)
  message("* Variables: ", str_out(vars, total = TRUE))
  message("* Samples: ", str_out(rownames(data), total = TRUE))
  GetLinearR = function(m_data) {
    SigCls <- vars
    if (scale) {
      m_data[, vars] <- scale(m_data[, vars, drop = FALSE])
    }
    fit <- sapply0(SigCls, function(v) {
      try(do.call(lme4::lmer, list(sprintf("%s %s", v, 
                                           model), data = m_data), envir = pe), silent = TRUE)
    })
    i_errs <- which(sapply(fit, is.atomic))
    if (length(i_errs)) {
      print(head(fit[i_errs]))
      fit <- fit[-i_errs]
    }
    if (fit.only) 
      return(fit)
    lsmeansStats(fit)
  }
  GetLinearR(data)
}


#Identification of features that changed over time using mixed effect linear models

differential.over.time.func <- function(eset.file.name) {
  
  eset.file <-readRDS(paste("../data/",  eset.file.name, ".rds", sep=""))
  no_cores<-detectCores()-1
  cl<-makeCluster(no_cores)
  
  clusterEvalQ(cl, {library(reshape);library(dplyr); library(CyTNF);library(plyr);library(Rcpp);library(shenR)})
  
  clusterExport(cl, varlist=c(eset.file.name))
  
  #parallel calculation for each 1000 rows
  res<-parSapply(cl,1:(ceiling(nrow(exprs(eset.file))/1000)-1), simplify=FALSE, function(count) {
    sequence <- seq(from = 0, to = nrow(exprs(eset.file)), by = 1000)+1
    first.seq <- round_any(sequence[count]-1, 1000, f = floor)+1
    cur.seq<-ifelse(count==(ceiling(nrow(exprs(eset.file))/1000)-1),nrow(exprs(eset.file)),sequence[count+1])-1 
    eset.file.sub <- eset.file[first.seq:cur.seq,]
    print(paste("start fitLMER current count:", count, sep=" "))
    cur.fit <- fitLMER(~ Visit + (1|Patient.code), eset.file.sub)
    fit_array<- reshapeFit(cur.fit)
    all <- data.frame(fit_array)
    saveRDS(all, paste("corjobfitLMER.R.",count,".rds"))
  })
  
  stopCluster(cl)
  
  #read and combine results from parallel jobs
  rds <- paste(rep("corjobfitLMER.R. ", ceiling(nrow(exprs(eset.file))/1000)-1), 1:(ceiling(nrow(exprs(eset.file))/1000)-1), rep(" .rds", ceiling(nrow(exprs(eset.file))/1000)-1), sep="")
  lmer <- do.call("rbind", lapply(rds, readRDS))
  saveRDS(paste("lmer", eset.file.name, "rds", sep="."))
  return(lmer)
}

#Calculate permutation based p-value
#generate shuffled expressen sets by visits and calculate permuted coefficients

perm.func <- function(perm.file.name, eset.file.name) {
  perm.file <-readRDS(paste("../data/",  perm.file.name, ".rds", sep=""))
  eset.file <-readRDS(paste("../data/",  eset.file.name, ".rds", sep=""))
  
  lapply(1:(nrow(perm.file)), function(cur.perm) {
    lapply(1:(ceiling(nrow(exprs(eset.file))/1000)-1), function(count) {
      cur.perm.filtered.esetALL <- eset.file
      sampleNames(cur.perm.filtered.esetALL) <- sampleNames(cur.perm.filtered.esetALL)[perm.file[cur.perm,]]
      pData(cur.perm.filtered.esetALL)$Visit <-pData(eset.file)[match(row.names(pData(cur.perm.filtered.esetALL)), pData(eset.file)$sampleID), 'Visit'] 
      pData(cur.perm.filtered.esetALL)$sampleID <-pData(eset.file)[match(row.names(pData(cur.perm.filtered.esetALL)), pData(eset.file)$sampleID), 'sampleID'] 
      pData(cur.perm.filtered.esetALL)$Patient.code <-pData(eset.file)[match(row.names(pData(cur.perm.filtered.esetALL)), pData(eset.file)$sampleID), 'Patient.code'] 
      pData(cur.perm.filtered.esetALL) <- pData(cur.perm.filtered.esetALL)[,colnames(pData(cur.perm.filtered.esetALL)) %in% c("sampleID", "Visit", "Patient.code")]
      
      sequence <- seq(from = 0, to = nrow(exprs(cur.perm.filtered.esetALL)), by = 1000)+1
      first.seq <- round_any(sequence[count]-1, 1000, f = floor)+1
      cur.seq<-ifelse(count==(ceiling(nrow(exprs(eset.file))/1000)-1),nrow(exprs(cur.perm.filtered.esetALL)),sequence[count+1])-1 
      cur.perm.filtered.esetALL.sub <- cur.perm.filtered.esetALL[first.seq:cur.seq,]
      saveRDS(cur.perm.filtered.esetALL.sub, paste("perm",eset.file.name, cur.perm, count,"rds", sep="."))
    })})
  
  
  lapply(1:nrow(perm.file), function(cur.perm) {
    lapply(1:(ceiling(nrow(exprs(eset.file))/1000)-1), function(count) {
      cur.perm.filtered.esetALL <- readRDS(paste("perm",eset.file.name, cur.perm, count,"rds", sep="."))
      cur.fit <- fitLMER(~ Visit + (1|Patient.code), cur.perm.filtered.esetALL)
      fit_array<- reshapeFit(cur.fit)
      all <- data.frame(fit_array)
      saveRDS(all, paste("perm.fitLMER", eset.file.name, cur.perm, count,"rds", sep="."))
    })})
  
  #read and combine results
  perm.res.files <- lapply(1:nrow(perm.file), function(cur.perm) {
    rds.list <- paste("perm.fitLMER", eset.file.name, cur.perm, 1:(ceiling(nrow(exprs(eset.file))/1000)-1),"rds", sep=".")
  })
  
  combined.perm.res <- lapply(perm.res.files, function(cur.rds.perm.list) {
    res.cur.perm <- do.call("rbind", lapply(cur.rds.perm.list, readRDS))
  })
  
  saveRDS(combined.perm.res, paste("combined.perm.res", eset.file.name, "rds", sep="."))
  return(combined.perm.res)
}

#Calculate p-value and FDR based on the permuted coefficient values
permutation.derived.FDR.func <- function(lmer, combined.perm.res.name, Visit) {
  combined.perm.res <-readRDS(paste("../data/",  combined.perm.res.name, ".rds", sep=""))
  permutations.count <- length(combined.perm.res)
  combined.perm.res.beta<-do.call('cbind', lapply(combined.perm.res, function(m) {m[,colnames(m) %in% paste(Visit, "Estimate", sep=".")]}))
  row.names(combined.perm.res.beta)<-row.names(combined.perm.res[[1]])
  colnames(combined.perm.res.beta)<-paste("perm", Visit, "estimate", 1:permutations.count, sep=".")
  
  rm(combined.perm.res)
  gc()
  combined.perm.res.beta<-data.frame(combined.perm.res.beta, stringsAsFactors=F)
  combined.perm.res.beta$featureID<-row.names(combined.perm.res.beta)
  
  merged.estimates <- merge(x =lmer[, c('featureID', paste(Visit,'Estimate', sep="."))], y = combined.perm.res.beta, by = "featureID") 
  
  #for calculation of permutation based p-value 
  # for each row in merged data calculate rate for which fold change of permutations is higher than for real data
  higher.change.rate.real <- apply(merged.estimates, 1, function(x) length(x[3:(permutations.count+2)][abs(as.numeric(x[3:(permutations.count+2)])) > abs(as.numeric(x[2]))])/permutations.count)
  Perm.Pvalue.df <- data.frame(featureID=merged.estimates$featureID, perm.pvalue=higher.change.rate.real, stringsAsFactors = F)
  saveRDS(Perm.Pvalue.df, paste("Perm.Pvalue", combined.perm.res.name, "rds", sep="."))
  
  #for calculation of permutation based FDR
  # for each permutation take all features and compare them to values from other permutations
  # to get a similar rate column as in previous step
  higher.change.rate.permutations <- ldply(3:(permutations.count+2), function(cur.index) {
    apply(merged.estimates, 1, function(cur.feature) length(cur.feature[3:(permutations.count+2)][abs(as.numeric(cur.feature[3:(permutations.count+2)])) > abs(as.numeric(cur.feature[cur.index]))])/(permutations.count-1))
  })
  
  # merge real and permutations change rate comparisons
  higher.change.rate <- rbind(higher.change.rate.real, higher.change.rate.permutations)
  
  rm(higher.change.rate.real, higher.change.rate.permutations)
  gc()
  
  # format data
  higher.change.rate <- data.frame(t(higher.change.rate))
  rownames(higher.change.rate) <- merged.estimates$featureID
  higher.change.rate <- setDT(higher.change.rate, keep.rownames = TRUE)[]
  colnames(higher.change.rate) <- colnames(merged.estimates)
  rm(merged.estimates)
  gc()
  
  # Create a calculated FDR column for each threshold
  dataTypes <- unique(sub("\\..+", "", higher.change.rate$featureID))
  
  FDR.combined <- ldply(dataTypes, function(cur.dt) {
    cur.higher.change.rate <- higher.change.rate[grep(paste0("^",cur.dt), higher.change.rate$featureID),]
    all.permutes.list <- unlist(cur.higher.change.rate[,3:(permutations.count+2)])
    real.list <- unlist(cur.higher.change.rate[,2])
    mean.for.threshold.permutes <- ldply(c(0.000001,0.000005,0.00001, 0.00005,0.0001,0.0005,0.001,0.005,0.01,seq(from = 0.05, to = 1, by = 0.05)), function(cur.P) c(cur.P, length(all.permutes.list[all.permutes.list <= cur.P])/permutations.count))
    count.for.threshold.real <- ldply(c(0.000001,0.000005,0.00001, 0.00005,0.0001,0.0005,0.001,0.005,0.01,seq(from = 0.05, to = 1, by = 0.05)), function(cur.P) c(cur.P, length(real.list[real.list <= cur.P])))
    colnames(mean.for.threshold.permutes) <- c("Prate", "permuteMean")
    colnames(count.for.threshold.real) <- c("Prate", "realCount")
    cur.FDR.combined <- merge(x = count.for.threshold.real, y = mean.for.threshold.permutes)
    cur.FDR.combined$dataType <- cur.dt
    cur.FDR.combined$calcFDR <- cur.FDR.combined$permuteMean/cur.FDR.combined$realCount
    saveRDS(cur.FDR.combined, paste("cur.FDR.combined", cur.dt, "rds", sep="."))
    rm(cur.FDR.combined, cur.higher.change.rate, all.permutes.list, real.list, 
       mean.for.threshold.permutes, count.for.threshold.real)
    gc()
  })
  saveRDS(paste("FDR", combined.perm.res.name, "rds", sep="."))
  return(list(Perm.Pvalue=Perm.Pvalue.df, FDR.combined=FDR.combined))
}

# Subsampling sampling function 

subsampling.func <- function(samp.file.name, eset.file.name) {
  sample.file <-readRDS(paste("../data/",  samp.file.name, ".rds", sep=""))
  eset.file <-readRDS(paste("../data/",  eset.file.name, ".rds", sep=""))
  
  lapply(1:(nrow(sample.file)), function(cur.sample) {
    lapply(1:(ceiling(nrow(exprs(eset.file))/1000)-1), function(count) {
      cur.sample.filtered.esetALL <- eset.file[,sample.file[,cur.sample]]
      sequence <- seq(from = 0, to = nrow(exprs(cur.sample.filtered.esetALL)), by = 1000)+1
      first.seq <- round_any(sequence[count]-1, 1000, f = floor)+1
      cur.seq<-ifelse(count==(ceiling(nrow(exprs(cur.sample.filtered.esetALL))/1000)-1),nrow(exprs(cur.sample.filtered.esetALL)),sequence[count+1])-1 
      cur.sample.filtered.esetALL.sub <- cur.sample.filtered.esetALL[first.seq:cur.seq,]
      saveRDS(cur.sample.filtered.esetALL.sub, paste("sample",eset.file.name, cur.sample, count,"rds", sep="."))
    })})
  
  lapply(1:nrow(sample.file), function(cur.sample) {
    lapply(1:(ceiling(nrow(exprs(eset.file))/1000)-1), function(count) {
      cur.sample.filtered.esetALL <- readRDS(paste("sample",eset.file.name, cur.sample, count,"rds", sep="."))
      cur.fit <- fitLMER(~ Visit + (1|Patient.code), cur.sample.filtered.esetALL)
      fit_array<- reshapeFit(cur.fit)
      all <- data.frame(fit_array)
      saveRDS(all, paste("sample.fitLMER", eset.file.name, cur.sample, count,"rds", sep="."))
    })})
  
  #read and combine results
  sample.res.files <- lapply(1:nrow(sample.file), function(cur.sample) {
    rds.list <- paste("sample.fitLMER", eset.file.name, cur.sample, 1:(ceiling(nrow(exprs(eset.file))/1000)-1),"rds", sep=".")
  })
  
  combined.perm.res <- lapply(sample.res.files, function(cur.rds.sample.list) {
    res.cur.sample <- do.call("rbind", lapply(cur.rds.sample.list, readRDS))
  })
  
  saveRDS(combined.perm.res, paste("combined.subsample.res", "rds", sep="."))
  return(combined.perm.res)
}

#Adjustment based on cell estimates
adj_func <- function(P, eset_exp) {
P_scaled <- apply(P, 1, function(x) scale(x, center = TRUE, scale = TRUE))
fit.res <- lsfit(P_scaled, 2^t(eset_exp), wt = NULL, intercept = TRUE,
                    yname = row.names(eset_exp))
intercept <- fit.res$coefficients[row.names(fit.res$coefficients) %in% "Intercept",]
intercept <- data.frame(intercept)
AG.fitlm <- t(fit.res$residuals) + intercept$intercept 
AG.fitlm <- log2(exprs(AG.fitlm)) #back to log2 scale
}

#Network propagation

network.propagation.function <- function(changes, Visit, esetALL) {
  protein.links <- read.table(gzfile("9606__protein_links.tsv.gz", open = "r"), header=T)
  protein.links.700 <- protein.links[protein.links$combined_score>=700,]
  genes <- unique(na.omit(changes$SYMBOL))
  string_db <- STRINGdb$new(version="10", species=9606, 
                            score_threshold=700, input_directory="~/Shiran/Remicade.analysis.July2018/all.features.together/IFX_paper" )
  sids <- string_db$map(data.frame(genes, stringsAsFactors = FALSE), "genes", removeUnmappedRows = TRUE, takeFirst = TRUE )
  colnames(sids)<-c("ID","STRING_id")
  
  interactors <- protein.links[protein.links$protein1 %in% sids$STRING_id|protein.links$protein2 %in% sids$STRING_id,] 
  interactors <- interactors[interactors$combined_score>=700,]
  
  converter <- read.table("entrez_gene_id.vs.string.v10.28042015.tsv", header=FALSE)
  
  interactors$ENTREZID1 <- converter[match(as.character(interactors$protein1), as.character(converter$V2)),'V1']
  interactors$ENTREZID2 <- converter[match(as.character(interactors$protein2), as.character(converter$V2)),'V1']
  
  library(org.Hs.eg.db)
  library(annotate)
  
  SYMBOL1 <- getSYMBOL(as.character(c(interactors$ENTREZID1[!is.na(interactors$ENTREZID1)], interactors$ENTREZID2[!is.na(interactors$ENTREZID2)])), data='org.Hs.eg')
  symbol.key <- data.frame(ENTREZID=as.character(c(interactors$ENTREZID1[!is.na(interactors$ENTREZID1)], interactors$ENTREZID2[!is.na(interactors$ENTREZID2)])), SYMBOL=SYMBOL1)
  
  interactors$SYMBOL1 <- symbol.key[match(interactors$ENTREZID1, symbol.key$ENTREZID),'SYMBOL']
  interactors$SYMBOL2 <- symbol.key[match(interactors$ENTREZID2, symbol.key$ENTREZID),'SYMBOL']
  
  interactors.final <- interactors[complete.cases(interactors),]
  
  #put the original network genes in SYMBOL1
  i <- as.character(interactors.final$SYMBOL2) %in% genes
  tmp <- interactors.final$SYMBOL2[i]
  interactors.final$SYMBOL2[i] <- interactors.final$SYMBOL1[i]
  interactors.final$SYMBOL1[i] <- tmp
  tmp <- interactors.final$ENTREZID2[i]
  interactors.final$ENTREZID2[i] <- interactors.final$ENTREZID1[i]
  interactors.final$ENTREZID1[i] <- tmp
  
  interactors.final$protein1 <- as.character(interactors.final$protein1)
  interactors.final$protein2 <- as.character(interactors.final$protein2)
  
  tmp <- interactors.final$protein2[i]
  interactors.final$protein2[i] <- interactors.final$protein1[i]
  interactors.final$protein1[i] <- tmp
  
  #test if the linker gene (in symbol2) hubs are enriched in changes
  
  res <- do.call('rbind', lapply(unique(interactors.final$protein2), function (x) {
    protein.links.700$protein1 <- as.character(protein.links.700$protein1)
    protein.links.700$protein2 <- as.character(protein.links.700$protein2)
    linker.gene.hubs <- protein.links.700[protein.links.700$protein1 %in% x|protein.links.700$protein2 %in% x, ]
    if(any((as.character(linker.gene.hubs$protein1) %in% x) ==FALSE)) {
      i <- as.character(linker.gene.hubs$protein2) %in% x
      tmp <- linker.gene.hubs$protein2[i]
      linker.gene.hubs$protein2[i] <- linker.gene.hubs$protein1[i]
      linker.gene.hubs$protein1[i] <- tmp
    }
    
    #perform hypergeometric test
    
    linker.gene.hubs$protein1 <- as.character(linker.gene.hubs$protein1)
    linker.gene.hubs$protein2 <- as.character(linker.gene.hubs$protein2)
    
    hubs.in.changes <-linker.gene.hubs[linker.gene.hubs$protein2 %in% sids$STRING_id,'protein2']
    q <- length(hubs.in.changes)
    changes.genes <- as.character(unique(sids$STRING_id[!is.na(sids$STRING_id)]))
    m <- length(changes.genes)
    gene.not.in.changes <- as.character(unique(fData(esetALL)[!as.character(fData(esetALL)$ENTREZID) %in% changes.genes,'ENTREZID']))
    n <- length(gene.not.in.changes[!is.na(gene.not.in.changes)])
    k <- length(unique(linker.gene.hubs$protein2))
    
    p.val <- ifelse(q==0, 1, 1-phyper(q, m, n, k, lower.tail=TRUE))
    return(data.frame(linker.name=x, enrichment.p=p.val, stringsAsFactors = FALSE))
  }))
  
  
  #res_t <- data.frame(t(res), stringsAsFactors = FALSE)
  res <- res[!duplicated(res$linker.name),]
  res$BH.FDR <- p.adjust(res$enrichment.p, method="BH")
  
  
  linker.genes.final <- res[res$BH.FDR<0.05,]
  
  linker.genes.final$ENTRZID <-  converter[match(as.character(linker.genes.final$linker.name), as.character(converter$V2)),'V1']
  
  linker.genes.final$SYMBOL <- getSYMBOL(as.character(linker.genes.final$ENTRZID), data='org.Hs.eg')
  linker.genes.final.featureIDs <- fData(esetALL)[which(fData(esetALL)$SYMBOL %in% linker.genes.final$SYMBOL),'featureID']
  saveRDS(linker.genes.final, paste(Visit ,"linker.genes.final.rds", sep="."))
  return(linker.genes.final)
}


#Construct co-expression matrix based on Spearman's r
#Define significant edges

sig.edge.func <- function(filtered.esetALL, Visit, changes.final) {
  x <- filtered.esetALL[featureNames(filtered.esetALL) %in% changes.final, pData(filtered.esetALL)$Visit == Visit]
  x <- exprs(x)
  x <- t(x)
  
  cor.Hmisc <- rcorr(x, type="spearman")
  
  rm(x)
  gc()
  
  real_cor <- cor.Hmisc$r
  real_cor[lower.tri(real_cor)]<-"NA"
  
  real_cor.p <- cor.Hmisc$P
  real_cor.p[lower.tri(real_cor.p)]<-NA
  
  cor.edges.list<-list(cor = real_cor, p.value = real_cor.p)
  rm(real_cor, real_cor.p)
  gc()
  
  cor.edges.list <- sapply(names(cor.edges.list), simplify=FALSE, function(n) {
    df <- reshape2::melt(cor.edges.list[[n]], value.name = n)
    df
  })
  
  cor.edges.list<-sapply(cor.edges.list, simplify=FALSE, function(n) {colnames(n)<-c("Var1","Var2","value"); n})
  
  cor.edges.list<-sapply(names(cor.edges.list), simplify=FALSE, function(n) {
    cor.edges.list[[n]][!(is.na(cor.edges.list[[n]][,3])),]})
  
  cor.edges.df<-Reduce(function(...) merge(..., by = c('Var1', 'Var2')), cor.edges.list)
  colnames(cor.edges.df)<-c("Var1","Var2","cor", "cor.p")
  
  rm(cor.edges.list)
  gc()
  
  cor.edges.df$type <- c()
  cor.edges.df[intersect(grep("PB|GR", cor.edges.df$Var1), grep("PB|GR", cor.edges.df$Var2)),'type'] <- "CCCC"
  cor.edges.df[intersect(grep("AG", cor.edges.df$Var1), grep("AG", cor.edges.df$Var2)),'type'] <- "AGAG"
  cor.edges.df[intersect(grep("AG", cor.edges.df$Var1), grep("PB|GR", cor.edges.df$Var2)),'type'] <- "AGCC"
  cor.edges.df[intersect(grep("LU", cor.edges.df$Var1), grep("LU", cor.edges.df$Var2)),'type'] <- "LULU"
  cor.edges.df[intersect(grep("AG", cor.edges.df$Var1), grep("LU", cor.edges.df$Var2)),'type'] <- "AGLU"
  cor.edges.df[intersect(grep("PB|GR", cor.edges.df$Var1), grep("LU", cor.edges.df$Var2)),'type'] <- "CCLU"
  cor.edges.df[intersect(grep("LU", cor.edges.df$Var1), grep("PB|GR", cor.edges.df$Var2)),'type'] <- "CCLU"
  
  
  unique.dt  <- unique(cor.edges.df$type)
  
  
  cor.edges.df$FDR<-c()
  cor.edges.df<-do.call('rbind', lapply(unique.dt, function(cur.dt) {
    cur.cor.edges.df <- cor.edges.df[cor.edges.df$type %in% cur.dt,]
    cur.cor.edges.df$FDR<-p.adjust(cur.cor.edges.df$cor.p, method = "BH", n = length(cur.cor.edges.df$cor.p))
    return(cur.cor.edges.df)
  }))
  
  
  #Filter by FDR and sort
  
  cor.edges.df.sub <-  do.call('rbind', lapply(unique.dt, function(cur.dt) {
    cur.df <- cor.edges.df[cor.edges.df$type %in% cur.dt,]
    ifelse(cur.dt %in% c("AGCC", "AGLU", "CCLU"), return(cur.df[cur.df$FDR<0.1,]), return(cur.df[cur.df$FDR<0.05,]))
  }))
  
  return(cor.edges.df.sub)
}


#Generated empirical null distribution of dropouts 
#('normal response' dropouts)

cor.responders.drop.fun <- function(changes.final, Visit.for.backbone, esetALL.responders, cellTypes) {
  cellTypes <- cellTypes
  assign("cellTypes", cellTypes, envir = .GlobalEnv)
  x <- esetALL.responders[featureNames(esetALL.responders) %in% changes.final, as.character(pData(esetALL.responders)$Visit) %in% Visit.for.backbone]
  x <- exprs(x)
  x <- t(x)
  assign("x", x, envir = .GlobalEnv)
  cc_n.r<-c()
  cc_n1.r <- c()
  cc_n1.r<-cor(x, method="spearman")
  assign("cc_n1.r",cc_n1.r, envir = .GlobalEnv)
  no_cores<-5
  cl<-makeCluster(no_cores, outfile = "debugSubsampling2.txt")
  clusterEvalQ(cl, {library(reshape2);library(plyr);})
  clusterExport(cl, varlist=c("x", "cellTypes", "cc_n1.r", "Visit.for.backbone"), envir = .GlobalEnv)
  res <- parSapply(cl, seq(length(row.names(x))), simplify=F, function(i) { 
    s<-x[!row.names(x) %in% row.names(x)[i],]
    cc_n.r <- cor(s, method="spearman")
    rm(s)
    gc()
    
    delta <- cc_n1.r - cc_n.r
    
    #transformation to Z scores
    
    Z.cc_n.r <- 0.5 * log((1+cc_n.r)/(1-cc_n.r)) 
    Z.cc_n1.r <- 0.5 * log((1+cc_n1.r)/(1-cc_n1.r)) 
    
    diff   <- Z.cc_n1.r - Z.cc_n.r 
    SEdiff <- sqrt( 1/(nrow(x) - 3) + 1/((nrow(x)-1) - 3) ) 
    diff.Z  <- diff/SEdiff 
    
    
    #fisher.z<- function (r1,r2,n1,n2) ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5
    pvalue=(2*(1-pnorm(abs(diff.Z))))/2 #one sided
    pvalue <- matrix(pvalue, nrow = nrow(cc_n.r), ncol = ncol(cc_n.r), 
                     dimnames = dimnames(cc_n.r))
    
    cc_n.r[lower.tri(cc_n.r)]<-NA
    delta[lower.tri(delta)]<-NA
    pvalue[lower.tri(pvalue)]<-NA
    diff.Z[lower.tri(diff.Z)]<-NA
    
    a<-list(cor =cc_n1.r, z.R=Z.cc_n1.r, delta=delta, pvalue=pvalue, diff.Z=diff.Z) 
    
    rm(cc_n.r,delta,pvalue, diff.Z, diff, SEdiff, Z.cc_n.r, Z.cc_n1.r, cc_n1.r, cc_n.r)
    gc()
    
    a <- sapply(names(a), simplify=FALSE, function(n) {
      df <- reshape2::melt(a[[n]], value.name = n)
      return(df)
    })
    
    a<-sapply(a, simplify=FALSE, function(n) {colnames(n)<-c("Var1","Var2","value"); n})
    
    df<-sapply(names(a), simplify=FALSE, function(n) {
      a[[n]][!(is.na(a[[n]][,3])),]})
    rm(a)
    gc()
    
    
    df<-Reduce(function(...) merge(..., by = c('Var1', 'Var2')), df)
    colnames(df)<-c("Var1","Var2","cor", "z.R", "delta", "drop.pvalue", "diff.Z")
    df$dir <- sign(df$cor * df$diff.Z)
    df$drop <- abs(df$diff.Z)
    
    df$type <- c("GG", "GC", "CC")[(df$Var1 %in% cellTypes) + (df$Var2 %in% cellTypes) + 1]
    df <- local({
      df$Var1 <- as.character(df$Var1)
      df$Var2 <- as.character(df$Var2)
      i <- df$type == "GC" & df$Var1 %in% cellTypes
      tmp <- df$Var2[i]
      df$Var2[i] <- df$Var1[i]
      df$Var1[i] <- tmp
      df$Var1 <- as.factor(df$Var1)
      df$Var2 <- as.factor(df$Var2)
      df})
    saveRDS(df, paste("cor.responders.drop",Visit.for.backbone, i, "rds", sep="."))
    rm(df)
    gc()
    #  return(df)
  })
  stopCluster(cl)
}


#Calculate non-responders drops:

#add each iteration one NUR and recalculate correlation matrix. 
#calculate the r drop for each pair of features as a result of the NUR patient addiation. 
#combine responders and non responders for calculating drop
#now the whole responders correlation network is used as cc_n and cc_n1 is the correlation 
#after addition of one non-responder


cor.disruption.fun <- function(changes.final, Visit.for.backbone, esetALL.responders, esetALL.NR, Visit.for.disruption, cellTypes) {
  cellTypes <- cellTypes
  assign("cellTypes", cellTypes, envir = .GlobalEnv)
  x <- esetALL.responders[featureNames(esetALL.responders) %in% changes.final, as.character(pData(esetALL.responders)$Visit) %in% Visit.for.backbone]
  x <- exprs(x)
  x <- t(x)
  assign("x", x, envir = .GlobalEnv)
  
  y<-esetALL.NR[featureNames(esetALL.NR) %in% changes.final, esetALL.NR$Visit == Visit.for.disruption]
  y<-exprs(y)
  assign("y", y, envir = .GlobalEnv)
  
  
  cc_n<-c()
  
  cc_n1.r<-cor(x, method="spearman")
  cc_n.r <- cc_n1.r
  assign("cc_n.r", cc_n.r, envir = .GlobalEnv)
  
  no_cores<-5
  cl<-makeCluster(no_cores, outfile = "debugSubsampling2.txt")
  clusterEvalQ(cl, {library(reshape2);library(plyr);})
  clusterExport(cl, varlist=c("x" ,"y", "cellTypes", "cc_n.r"), envir=.GlobalEnv)
  res <- parSapply(cl, seq(length(row.names(t(y)))), simplify=F, function(i) { 
    s<-rbind(x,t(y)[row.names(t(y)) %in% row.names(t(y))[i],])
    cc_n1.r<-cor(s, method="spearman") #new disrupted correlation matrix with first NUR patient 
    rm(s)
    gc()
    
    delta <- cc_n1.r - cc_n.r
    
    #transformation to Z scores
    
    Z.cc_n.r <- 0.5 * log((1+cc_n.r)/(1-cc_n.r)) 
    Z.cc_n1.r <- 0.5 * log((1+cc_n1.r)/(1-cc_n1.r)) 
    
    diff   <- Z.cc_n1.r - Z.cc_n.r 
    SEdiff <- sqrt( 1/((nrow(x)+1) - 3) + 1/((nrow(x)) - 3) ) 
    diff.Z  <- diff/SEdiff 
    
    #fisher.z<- function (r1,r2,n1,n2) ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5
    pvalue=(2*(1-pnorm(abs(diff.Z))))/2 #one sided
    pvalue <- matrix(pvalue, nrow = nrow(cc_n.r), ncol = ncol(cc_n.r), 
                     dimnames = dimnames(cc_n.r))
    
    cc_n.r[lower.tri(cc_n.r)]<-NA
    delta[lower.tri(delta)]<-NA
    pvalue[lower.tri(pvalue)]<-NA
    diff.Z[lower.tri(diff.Z)]<-NA
    
    a<-list(cor =cc_n1.r, z.R=Z.cc_n1.r, delta=delta, pvalue=pvalue, diff.Z=diff.Z) 
    rm(cc_n.r,delta,pvalue, diff.Z, diff, SEdiff, Z.cc_n.r, Z.cc_n1.r, cc_n1.r, cc_n.r)
    gc()
    
    a <- sapply(names(a), simplify=FALSE, function(n) {
      df <- reshape2::melt(a[[n]], value.name = n)
      df
    })
    
    a<-sapply(a, simplify=FALSE, function(n) {colnames(n)<-c("Var1","Var2","value"); n})
    
    df<-sapply(names(a), simplify=FALSE, function(n) {
      a[[n]][!(is.na(a[[n]][,3])),]})
    rm(a)
    gc()
    
    
    df<-Reduce(function(...) merge(..., by = c('Var1', 'Var2')), df)
    colnames(df)<-c("Var1","Var2","cor", "z.R", "delta", "drop.pvalue", "diff.Z")
    df$dir <- sign(df$cor * df$diff.Z)
    df$drop <- abs(df$diff.Z)
    
    df$type <- c("GG", "GC", "CC")[(df$Var1 %in% cellTypes) + (df$Var2 %in% cellTypes) + 1]
    df <- local({
      df$Var1 <- as.character(df$Var1)
      df$Var2 <- as.character(df$Var2)
      i <- df$type == "GC" & df$Var1 %in% cellTypes
      tmp <- df$Var2[i]
      df$Var2[i] <- df$Var1[i]
      df$Var1[i] <- tmp
      df$Var1 <- as.factor(df$Var1)
      df$Var2 <- as.factor(df$Var2)
      df})
    saveRDS(df, paste("cor.NR.drop", Visit.for.disruption, i, "rds", sep="."))
    rm(df)
    gc()
    
  })
  stopCluster(cl)
}


#Determine NR-drop significance for each edge in each sample by calculation of left-tail percentile, 
#within the null distribution of the normal dropouts for each edge

percentile.and.FDR.calc.function <- function(df.l.disrupted, df.l.refernce.drop.CI) {
  if(identical(df.l.refernce.drop.CI$Var1, df.l.disrupted[[1]]$Var1) &
     identical(df.l.refernce.drop.CI$Var2, df.l.disrupted[[1]]$Var2)) {
    df.l.disrupted.percentile.final <- lapply(df.l.disrupted, function(m) {
      m$drop.dir <- m$drop*m$dir
      for.percentile <- cbind(df.l.refernce.drop.CI[,!colnames(df.l.refernce.drop.CI) %in% c("type", "Var1", "Var2")], m[,colnames(m) %in% c("Var1", "Var2", "drop.dir")])
      for.percentile$percentile <- apply(for.percentile, 1, function(x) {
        sum(as.numeric(x[!names(x) %in% c("Var1", "Var2", "drop.dir", "percentile", "type")]) < as.numeric(x[names(x) %in% "drop.dir"]))/length(x[!names(x) %in% c("Var1", "Var2", "drop.dir", "percentile")]) 
      })
      m$percentile <- for.percentile$percentile
      rm(for.percentile)
      gc()
      return(m)
    })
    names(df.l.disrupted.percentile.final) <- names(df.l.disrupted)
    
    #Add FDR per dataType
    unique.dt <- unique(df.l.disrupted.percentile.final[[1]]$type)
    
    df.l.disrupted.percentile.final <- lapply(df.l.disrupted.percentile.final, function(m) {
      m$FDR<-c()
      tmp<-do.call('rbind', lapply(unique.dt, function(cur.dt) {
        cur.df.l.NUR.percentile <- m[m$type %in% cur.dt,]
        cur.df.l.NUR.percentile$disrup.FDR<-p.adjust(cur.df.l.NUR.percentile$percentile, method = "BH", n = length(cur.df.l.NUR.percentile$percentile))
        return(cur.df.l.NUR.percentile)
      }))
      return(tmp)
    })
  } else {print("not identical")}
}


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



#Calculate disruption in the geneset level

GSEA.Module.disruption.measurment.function <- function(geneSets.sub,  drop_matrix.edited.aggregated, Visit) {
  percentage.of.disruption <- do.call('rbind', lapply(geneSets.sub, function(cur.clust) {
    print(cur.clust)
    nodes.in.clust.featureID <- net.nodes.by.pathways[net.nodes.by.pathways$path %in% cur.clust, "featureID"]
    do.call('rbind', lapply(patients, function (cur.patient) {
      drop.patient <- drop_matrix.edited.aggregated[,colnames(drop_matrix.edited.aggregated) %in% c(paste(cur.patient), "Var1", "Var2")]
      drop.patient.in.module <- drop.patient[drop.patient$Var1 %in% nodes.in.clust.featureID & drop.patient$Var2 %in% nodes.in.clust.featureID,]
      drop.patient.in.module <- drop.patient.in.module[!is.na(drop.patient.in.module[,paste(cur.patient)]),]
      drop.patient.in.module.within.module <- drop.patient.in.module[drop.patient.in.module$Var1 %in% nodes.in.clust.featureID & drop.patient.in.module$Var2 %in% nodes.in.clust.featureID,]
      drop.patient.in.module.count <- drop.patient.in.module[abs(drop.patient.in.module[,cur.patient])>0.4,]
      drop.patient.in.module.count.within.module <- drop.patient.in.module.count[drop.patient.in.module.count$Var1 %in%  nodes.in.clust.featureID & drop.patient.in.module.count$Var2 %in%  nodes.in.clust.featureID,]
      drop.patient.in.module.count.nodes <- unique(c(as.character(drop.patient.in.module.count$Var1), as.character(drop.patient.in.module.count$Var2)))
      drop.patient.in.module.count.nodes <- drop.patient.in.module.count.nodes[drop.patient.in.module.count.nodes %in% nodes.in.clust.featureID]
      percentage.disrupt.in.module <- nrow(drop.patient.in.module.count)/nrow(drop.patient.in.module)*100
      percentage.disrupt.in.module.within.module <- nrow(drop.patient.in.module.count.within.module)/nrow(drop.patient.in.module.within.module)*100
      mean.positive.drop.intensity <- abs(drop.patient.in.module[,colnames(drop.patient.in.module) %in% cur.patient])
      mean.positive.drop.intensity <- mean.positive.drop.intensity[!is.infinite(mean.positive.drop.intensity) & !is.na(mean.positive.drop.intensity)]
      mean.positive.drop.intensity <- mean(mean.positive.drop.intensity)
      mean.drop.intensity.within.module <- abs(drop.patient.in.module.within.module[,cur.patient])
      mean.drop.intensity.within.module <- mean.drop.intensity.within.module[!is.infinite(mean.drop.intensity.within.module) & !is.na(mean.drop.intensity.within.module)]
      mean.drop.intensity.within.module <- mean(mean.drop.intensity.within.module)
      #  #Add constrain of disrupted node- node which is disrupted in at least 50% of its edges in the module
      tmp.node.disrup.perc <- do.call('rbind', lapply(nodes.in.clust.featureID, function(cur.node) {
        percent.disruption.of.node <- nrow(drop.patient.in.module.count[drop.patient.in.module.count$Var1 %in% cur.node|drop.patient.in.module.count$Var2 %in% cur.node,])/nrow(drop.patient.in.module[drop.patient.in.module$Var1 %in% cur.node|drop.patient.in.module$Var2 %in% cur.node,])*100
        return(data.frame(node=cur.node, perc=percent.disruption.of.node, stringsAsFactors = F))
      }))
      percentage.of.disrupted.nodes <- sum(tmp.node.disrup.perc$perc>70, na.rm=T)/length(unique(c(drop.patient.in.module$Var1, drop.patient.in.module$Var2)))*100
      res.enrich.p <- data.frame(cluster=cur.clust, patient=cur.patient, percentage.disrupt.in.module=percentage.disrupt.in.module, percentage.disrupt.in.module.within.module=percentage.disrupt.in.module.within.module, 
                                 mean.drop.intensity=mean.positive.drop.intensity, mean.drop.intensity.within.module=mean.drop.intensity.within.module, percentage.of.disrupted.nodes=percentage.of.disrupted.nodes, 
                                 number.of.edges.in.the.module=nrow(drop.patient.in.module.count.within.module), Total.edge.count=nrow(drop.patient.in.module), stringsAsFactors = F)
      return(res.enrich.p)
    }))
  }))
  return(percentage.of.disruption)
}

#Disruption count per feature 

disruption.count.fun <- function(cor.edges.df.FC.V1V2.aggregated, drop_matrix, Visit) {
  Visit <- Visit
  cur.drop_matrix <- drop_matrix[,colnames(drop_matrix) %in% c("Var1", "Var2")|grepl(Visit, colnames(drop_matrix))]
  disruption.count <- sapply(seq(nrow(cor.edges.df.FC.V1V2.aggregated)), simplify=TRUE, function(m) {
    probMat.index <- intersect(which(as.character(cur.drop_matrix$Var1) %in% as.character(cor.edges.df.FC.V1V2.aggregated$Var1[m])), which(as.character(cur.drop_matrix$Var2) %in% as.character(cor.edges.df.FC.V1V2.aggregated$Var2[m])))
    probMat <- cur.drop_matrix[probMat.index,]
    ProbMat.num <- probMat[,!(colnames(probMat) %in% c("Var1", "Var2"))]
    sum(as.numeric(abs(ProbMat.num))>0.3)
  })
  
  disruption.count <- cbind(as.character(cor.edges.df.FC.V1V2.aggregated$Var1), as.character(cor.edges.df.FC.V1V2.aggregated$Var2), disruption.count)
  disruption.count <- data.frame(disruption.count, stringsAsFactors = FALSE)
  colnames(disruption.count) <- c("Var1", "Var2", "patient.count.for.disruption")
  nodes <- unique(c(disruption.count$Var1, disruption.count$Var2))
  
  disruption.count[is.na(disruption.count)] <- "0"
  
  disruption.count$patient.count.for.disruption <- as.numeric(disruption.count$patient.count.for.disruption)
  
  disruption.count <- do.call('rbind', lapply(nodes, function(m) {
    probmat <- subset(disruption.count, disruption.count$Var1 %in% m | disruption.count$Var2 %in% m)
    backbone.edge.count <- nrow(probmat)
    #disrupted.edge.ratio <- sum(probmat$patient.count.for.disruption>=5)/nrow(probmat)*100 #in at least one patient
    max.ptient.count.for.disruption <- mean(probmat$patient.count.for.disruption)
    probMat.drop.matrix <- cur.drop_matrix[as.character(cur.drop_matrix$Var1) %in% m| as.character(cur.drop_matrix$Var2) %in% m,]
    mean.drop.intensity <- mean(apply(probMat.drop.matrix[,!colnames(probMat.drop.matrix) %in% c("Var1", "Var2")], 2, mean))
    disrupted.edge.ratio <- mean(apply(probMat.drop.matrix[,!colnames(probMat.drop.matrix) %in% c("Var1", "Var2")], 2, function(x) sum(x<=(-0.3))/length(x)))
    df <- data.frame(node=m,  backbone.edge.count= backbone.edge.count, disrupted.edge.ratio=disrupted.edge.ratio, max.patient.count.for.disruption=max.ptient.count.for.disruption, mean.drop.intensity=mean.drop.intensity, stringsAsFactors = FALSE)
    colnames(df) <- c("node", paste("backbone.edge.count", ".", Visit, sep=""), paste("disrupted.edge.ratio", ".", Visit, sep=""), paste("max.patient.count.for.disruption", ".", Visit, sep=""), paste("mean.drop.intensity", ".", Visit, sep=""))
    return(df)
  }))
}

