#!/srv01/technion/shirang/Rscript
message(.libPaths())

rm(list = ls())
library(citrus)


#**************************************************************************************
# IMPORTANT INFORMATION TO LIMIT # OF CORES Running. When this program is run on the cluster we MUST preserve control of the number of threads it uses or said 
# differently - so that it does not steal cores without the job scheduler knowing. There are 3 ways to do so:
# 1) In the call to function citrus.clusterAndMapFolds add the param mc.cores=16
# 2) Set the # of threads in the Rapp call like so: Rclusterpp.setThreads(X);options("mc.cores”=X);. Specifically, I think u should use 1 (ask David to be sure): 
# Rclusterpp.setThreads(1); options("mc.cores"=1);
# 3) Limit through the job scheduler using the P queue like so:
#qsub -q P -l nodes=1:ppn=4 -N VTA_SPBL  runCitrus2class.sh
#**********************************************************************************

# Use this line to limit the number of threads used by clustering
# 
DEFINE_PARAMETERS = TRUE
GENERATE_CLUSTERS = TRUE
COMPUTE_FEATURES_AND_DIFFERENCES = TRUE
PLOT_OUTPUT = TRUE
PLOT_OUTPUT_MAPS=TRUE


#give an example file from all the files, to be used in order to read the column names
#sample<-read.FCS("/storage/md_shenorr/shirang/Remicade.CyTOF.Gran.Fcs.Files/HR-20-V1_cells_found_normalized_Granulocytes.fcs")
#PLOT_NOT_SIG_CLUSTERS = FALSE
# Use this line to limit the number of threads used by clustering

if(DEFINE_PARAMETERS) {
  print('DEFINE_PARAMETERS')
	family = "classification"
	modelTypes = c("sam")
	nFolds = 1
    Rclusterpp.setThreads(8); options("mc.cores"=8);
  #featureTypes = c("medians")
  featureTypes = c("abundances","medians")
	citrusRes = list()
	citrusRes$General = list() #General is the condition name
	cnd = 'General'

	# Change if you want to run from command line
	# dataDirectory = "../"

	#inser comparison name that u choose
	citrusRes$comparison = 'RvsNR'

	#insert the path to the directory on cluster which the FCS filed are saved in
  dataDirectory = "/storage/md_shenorr/shirang/Remicade.CyTOF.Gran.Fcs.Files/"
  #args <- commandArgs(trailingOnly = TRUE)
	citrusRes$minClustSzPrcnt = 0.02
	citrusRes$fileSampleSize = 15000                                            
	print(paste('fileSampleSize = ',citrusRes$fileSampleSize, 'min cluster size',citrusRes$minClustSzPrcnt))
	exp_Type='Grans'
	exp = paste(citrusRes$comparison,'Fsz',(citrusRes$fileSampleSize/1000),'72samples_1206',exp_Type,'Csz',sub('\\.','p',citrusRes$minClustSzPrcnt),sep='')
	citrusRes$folderName = paste('citrusOut',exp,sep='')

  fileList = data.frame(list.files(dataDirectory,pattern='HR-'), stringsAsFactors=FALSE)
  pheno=read.csv("/storage/md_shenorr/shirang/Remicade.CyTOF.Gran.Fcs.Files/Gransfcs.annotation.csv",header=TRUE)
citrusRes$sample.labels = as.factor(pheno[match(fileList[,1],pheno$name),'pheno'])

markersAno=read.csv(file = "/storage/md_shenorr/shirang/Remicade.CyTOF.Gran.Fcs.Files/clustering by phenotyping for grans citrus_ES.csv",header = TRUE, stringsAsFactors=FALSE)
  
clusteringColumns = markersAno[markersAno[[exp_Type]]  %in% c('C','B'),'Marker' ]
  transformColumns = markersAno[markersAno[[exp_Type]]  %in% c('C','B','F'),'Marker' ]
  transformCofactor = 5
  
  medianColumns = markersAno[markersAno[[exp_Type]] %in% c('F','B'),'Marker' ]
	
outputDir=dataDirectory

}


if(GENERATE_CLUSTERS) {
  #**************************************************************************************
  # IMPORTANT - Use the commented out code below when run on cluster (see note on top)
  #citrusRes$foldClustering = citrus.clusterAndMapFolds(citrusRes$combinedFCSSet,clusteringColumns,citrusRes$sample.labels,nFolds,mc.cores=16)
  print('GENERATE_CLUSTERS\n')
  # Read Data
  citrusRes$combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList,citrusRes$fileSampleSize,transformColumns,transformCofactor)
  colnames(citrusRes$combinedFCSSet$fileId) = 'General'
  # Cluster all the data
  citrusRes$foldClustering = citrus.clusterAndMapFolds(citrusRes$combinedFCSSet,clusteringColumns,citrusRes$sample.labels,nFolds)
  print(paste('fileSampleSize = ',citrusRes$fileSampleSize, 'min cluster size',citrusRes$minClustSzPrcnt))
  save(file = paste(outputDir,'/citrusOutClustering.RData',sep=''),citrusRes,family,featureTypes,dataDirectory,outputDir,nFolds)  # save clustering
} else {
  load(file = paste(outputDir,'/citrusOutClustering.RData',sep=''))
}

if(COMPUTE_FEATURES_AND_DIFFERENCES) {
  print('COMPUTE_FEATURES_AND_DIFFERENCES\n')
  citrusRes$General = list() #General is teh condition name
  
  clusterIds = citrus.selectClusters.minimumClusterSize(citrusRes$foldClustering$allClustering,minimumClusterSizePercent = citrusRes$minClustSzPrcnt)
  #for(cnd in citrusRes$General) {
  cnd = 'General'
  
  for(ft in featureTypes) {
    if(ft == 'medians') {
      citrusRes[[cnd]][[ft]]$foldFeatureSet = citrus.calculateFoldFeatureSet(citrusRes$foldClustering,citrusRes$combinedFCSSet,featureType=ft,minimumClusterSizePercent=citrusRes$minClustSzPrcnt,medianColumns=medianColumns)
    } else {
      citrusRes[[cnd]][[ft]]$foldFeatureSet = citrus.calculateFoldFeatureSet(citrusRes$foldClustering,citrusRes$combinedFCSSet,featureType=ft,minimumClusterSizePercent=citrusRes$minClustSzPrcnt) 
    }
    # Build regression models for each model type
    citrusRes[[cnd]][[ft]]$regressionResults = mclapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrusRes[[cnd]][[ft]]$foldFeatureSet,labels=citrusRes$sample.labels,family=family)
    # Plot Results for each model
  }
  #	}  
  save(file = paste(outputDir,'/citrusOutFinal.RData',sep=''),citrusRes,family,featureTypes,dataDirectory,outputDir,nFolds)
} else {
  load(file = paste(outputDir,'/citrusOutFinal.RData',sep=''))
}
# 
if(PLOT_OUTPUT) {
  print('Plot output\n')
  for (conditionName in names(citrusRes$General)) {
    
    cat(paste0("\nPlotting Results for ", conditionName, "\n"))
    conditionOutputDir = paste(outputDir,'/',conditionName,sep='')
    dir.create(conditionOutputDir, showWarnings = F)
    lapply(citrusRes[[cnd]][[conditionName]]$regressionResults,plot,outputDirectory=conditionOutputDir,citrus.foldClustering=citrusRes$foldClustering,citrus.foldFeatureSet=citrusRes[[cnd]][[conditionName]]$foldFeatureSet,citrus.combinedFCSSet=citrusRes$combinedFCSSet)
    
    #     mclapply(citrusRes$General[[conditionName]]$RegressionResults, 
    #              plot, outputDirectory = conditionOutputDir, citrus.foldClustering = citrusRes$foldClustering, 
    #              citrus.foldFeatureSet = citrusRes$General[[conditionName]]$foldFeatureSet, 
    #              citrus.combinedFCSSet = citrusRes$combinedFCSSet, 
    #              family = family, labels = citrusRes$sample.labels)
    #     cat("\n")
  }  
}
citrusRes$foldClustering$allClustering$clusteringColumns<-transformColumns

if(PLOT_OUTPUT_MAPS) {
  hierarchyGraph = citrus.createHierarchyGraph(citrus.clustering = citrusRes$foldClustering$allClustering, 
                                               selectedClusters = citrusRes$General$abundances$foldFeatureSet$allLargeEnoughClusters)
  ######################################################################################################## 
  .getDisplayNames=function(combinedFCSSet,clusteringColumns){
    colLabels = combinedFCSSet$fileChannelNames[[1]][[1]]
    reagentNames = combinedFCSSet$fileReagentNames[[1]][[1]]
    displayNames = colLabels
    displayNames[nchar(reagentNames)>2] = reagentNames[nchar(reagentNames)>2]
    if (all(is.numeric(clusteringColumns))){
      return(displayNames[clusteringColumns])
    } else {
      names(displayNames) = colLabels
      return(as.vector(displayNames[clusteringColumns]))
    }
    
  }
  
  .getClusterMedians = function(clusterId,clusterAssignments,clusterCols,data){
    apply(data[clusterAssignments[[clusterId]],clusterCols],2,median)
  }
  
  ############################################################################################################  
  
  clusterMedians = t(sapply(citrusRes$General$abundances$foldFeatureSet$allLargeEnoughClusters, 
                            .getClusterMedians, clusterAssignments = citrusRes$foldClustering$allClustering$clusterMembership, 
                            data = citrusRes$combinedFCSSet$data,clusterCols = citrusRes$foldClustering$allClustering$clusteringColumns))
  rownames(clusterMedians) = citrusRes$foldFeatureSet$allLargeEnoughClusters
  colnames(clusterMedians) = .getDisplayNames(citrusRes$combinedFCSSet, 
                                              transformColumns)
  
  #clusterMedians=citrusRes$General$medians$foldFeatureSet$allFeatures[1,]
  outputDirectory=  dataDirectory
  citrus.plotClusteringHierarchy(outputFile = file.path(outputDirectory, 
                                                        "markerPlots.pdf"), clusterColors = clusterMedians, 
                                 graph = hierarchyGraph$graph, layout = hierarchyGraph$layout, 
                                 plotSize = hierarchyGraph$plotSize, theme = "white")
  citrus.plotClusteringHierarchy(outputFile = file.path(outputDirectory, 
                                                        "markerPlotsAll.pdf"), clusterColors = clusterMedians, 
                                 graph = hierarchyGraph$graph, layout = hierarchyGraph$layout, 
                                 plotSize = hierarchyGraph$plotSize, theme = "white", 
                                 singlePDF = T, plotClusterIDs = F)
}
}

