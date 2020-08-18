#' Build a deconvolution seed matrix, add the proportional option 
#' @description Use ranger to select features and build a genesInSeed gene matrix
#'
#' @param trainSet  Each row is a gene, and each column is an example of a particular cell type, ie from single cell data 
#' @param genesInSeed  The maximum number of genes in the returned seed matrix (DEFAULT: 200)
#' @param groupSize  The number of groups to break the trainSet into by ADAPTS::scSample (DEFAULT: 30)
#' @param randomize  Set to TRUE randomize the sets selected by ADAPTS::scSample (DEFAULT: TRUE)
#' @param num.trees  The number of trees to be used by ranger (DEFAULT: 1000)
#' @param plotIt  Set to TRUE to plot (DEFAULT: TRUE)
#' @param trainSet.3sam  Optional pre-calculated ADAPTS::scSample(trainSet, groupSize = 3) (DEFAULT: NULL)
#' @param trainSet.30sam  Optional pre-calculated ADAPTS::scSample(trainSet, groupSize=groupSize, randomize=randomize) (DEFAULT: NULL)
#' @param proportional  Set to true to make the training set cell type proportional.  Ignores group size (DEFAULT: FALSE)
#'
#' @export
#' @return A list with condition numbers and gene lists
#' @examples
#' library(ADAPTS)
#' ct1 <- runif(1000, 0, 100)
#' ct2 <- runif(1000, 0, 100)
#' dataMat <- cbind(ct1, ct1, ct1, ct1, ct1, ct1, ct2, ct2, ct2, ct2)
#' rownames(dataMat) <- make.names(rep('gene', nrow(dataMat)), unique=TRUE)
#' noise <- matrix(runif(nrow(dataMat)*ncol(dataMat), -2, 2), nrow = nrow(dataMat), byrow = TRUE)
#' dataMat <- dataMat + noise
#' newSigMat <- buildSeed(trainSet=dataMat)
#' 
buildSeed <- function(trainSet, genesInSeed=200, groupSize=30, randomize=TRUE, num.trees=1000, plotIt=TRUE, trainSet.3sam=NULL, trainSet.30sam=NULL, proportional=FALSE) {
  if(is.null(trainSet.3sam)) {trainSet.3sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = 3, randomize = randomize)}
  if (proportional==TRUE) {
    #colnames(trainSet) <- sub('\\.[0-9]+$', '', colnames(trainSet))
    tsNames <- sub('\\.[0-9]+$', '', colnames(trainSet))
    cellProps <- table(tsNames)
    cellSampleCounts <- 3*ceiling(100*cellProps / sum(cellProps))
    trainList <- list()
    for (curCount in unique(cellSampleCounts)) {
      curClusts <- names(cellSampleCounts)[cellSampleCounts==curCount]
      trainList[[as.character(curCount)]] <- ADAPTS::scSample(RNAcounts = trainSet[, tsNames %in% curClusts], groupSize = curCount, randomize = randomize)
    }
    trainSet.30sam <- do.call(cbind, trainList)
  } else {
    if(is.null(trainSet.30sam)) {trainSet.30sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = groupSize, randomize = randomize)}
  }
  
  clusterIDs <- factor(colnames(trainSet.30sam))
  trainSet.4reg <- t(trainSet.30sam)
  
  rf1 <- ranger::ranger(x=trainSet.4reg, y=clusterIDs, num.trees=num.trees, importance='impurity')
  imp <- ranger::importance(rf1)
  imp <- sort(imp[imp>0], decreasing = TRUE)
  
  
  topGenes <- names(imp)[1:min(genesInSeed, length(imp))]
  #topGenes[!topGenes %in% rownames(trainSet.3sam)]
  
  seedMat <- trainSet.3sam[rownames(trainSet.3sam) %in% topGenes,]
  cellTypes <- sub('\\.[0-9]+$', '', colnames(seedMat))
  seedMat <- t(apply(seedMat, 1, function(x){tapply(x, cellTypes, mean, na.rm=TRUE)}))
  if(plotIt==TRUE) {pheatmap:: pheatmap(seedMat, main='Marker Gene Matrix, train, log2(x+1)') }
  return(seedMat)
}
  
  
#' Load pre-defined clusters
#' @description Load a pre-specified cluster list that can be called by handMetaCluster()
#'
#' @export
#' @return A list of grouped clusters 
#' @examples 
#' handCluster<-loadHandClusters()  
loadHandClusters <- function() {
    metaList <- list()
    metaList[[1]] <- c('Cluster_0')  #DuctalCellType 2 
    metaList[[2]] <- c('Cluster_1', 'Cluster_9', 'Cluster_11')  #DuctalCellType 1 
    metaList[[3]] <- c('Cluster_2', 'Cluster_13', 'Cluster_14')  #Endothelial 
    metaList[[4]] <- c('Cluster_3', 'Cluster_5', 'Cluster_12')  #Stellate_Fibroblast 
    metaList[[5]] <- c('Cluster_4')  #Macrophage 
    metaList[[6]] <- c('Cluster_6')  #TCell
    metaList[[7]] <- c('Cluster_7')  #BCell
    metaList[[8]] <- c('Cluster_8')  #Acinar
    metaList[[9]] <- c('Cluster_10')  #PlasmaCell
    return(metaList)
  }


#' Generate all the signature matrices one time
#' @description  This wrapper is helpful for repetively matrix generation. It generates seed matrix, all-gene matrix, augmneted matrix, shrunk matrix,
#' and all the clustered matrices in one call.
#'
#' @param exprData The gene express data. Each row is a gene, and each column is an example of a particular cell type.
#' @param randomize Set to TRUE randomize the sets selected by ADAPTS::scSample (DEFAULT: TRUE)
#' @param skipShrink Set to TRUE to skip shrinking the signatrure matrix (DEFAULT: TRUE)
#' @param proportional Set to true to make the training set cell type proportional.  Ignores group size (DEFAULT: FALSE)
#' @param handMetaCluster Load in pre-defined meta clusters. Set to NULL to automatically group indistinguishable 
#' cells into same cluster use clustWspillOver(DEFAULT: NULL)
#'
#' @export
#' @return A list of results including prediction accuracy and cell enrichment
#'
#' @examples

scMatrixTest<-function(exprData, randomize = TRUE, skipShrink=FALSE, proportional=FALSE,handMetaCluster=NULL) {
  
  if(randomize==TRUE) {set.seed(Sys.time())}
  resList <- list()
  
  trainTestSet <- ADAPTS::splitSCdata(exprData, numSets=2, randomize = randomize)
  trainSet <- trainTestSet[[1]]
  testSet <-trainTestSet[[2]]
  
  trainSet.30sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = 30, randomize = randomize)
  trainSet.3sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = 3, randomize = randomize)
  
  pseudobulk.test <- data.frame(test=rowSums(testSet))
  pseudobulk.test.counts<-table(sub('\\..*','',colnames(testSet)))
  actFrac.test <- 100 * pseudobulk.test.counts / sum(pseudobulk.test.counts)
  
  
  #seed
  genesInSeed<-100
  seedMat <- buildSeed(trainSet, genesInSeed=genesInSeed, groupSize=30, randomize=TRUE, num.trees=1000, plotIt=TRUE, trainSet.3sam=trainSet.3sam, trainSet.30sam=trainSet.30sam,proportional = proportional)
  resList[['matrix.seed']] <- seedMat
  
  estimates.onTest <- as.data.frame(ADAPTS::estCellPercent.DCQ(seedMat, pseudobulk.test))
  
  colnames(estimates.onTest) <- paste(genesInSeed, 'Marker Genes Seed')
  estimates.onTest$actFrac.test <- round(actFrac.test[rownames(estimates.onTest)],2)
  resList[['estimates.onTest']] <- estimates.onTest
  
  resList[['testAcc.seed']] <- seed2TestAcc <- ADAPTS::calcAcc(estimates=estimates.onTest[,1], reference=estimates.onTest[,2])
  
  
  #All gene
  allGeneSig <- apply(trainSet.3sam, 1, function(x){tapply(x, colnames(trainSet.3sam), mean, na.rm=TRUE)})
  
  estimates.allGene <- as.data.frame(ADAPTS::estCellPercent.DCQ(t(allGeneSig), pseudobulk.test))
  colnames(estimates.allGene)<-'All Gene Matrix'
  
  estimates.onTest<-cbind(estimates.allGene,estimates.onTest)
  resList[['estimates.onTest']] <- estimates.onTest
  
  resList[['testAcc.all']] <- seed2TestAcc <-  ADAPTS::calcAcc(estimates=estimates.onTest[,1],reference=estimates.onTest[,ncol(estimates.onTest)])
  
  
  # Aug
  gList <- ADAPTS::gListFromRF(trainSet=trainSet.30sam)
  
  augTrain <- ADAPTS::AugmentSigMatrix(origMatrix = seedMat, fullData = trainSet.3sam, gList = gList, nGenes = 1:100, newData = trainSet.3sam, plotToPDF = FALSE, pdfDir = '.')
  
  resList[['matrix.aug']] <- augTrain
  
  estimates.augment <- as.data.frame(ADAPTS::estCellPercent.DCQ(augTrain, pseudobulk.test))
  colnames(estimates.augment) <- paste('Augmented Matrix')
  estimates.onTest <- cbind(estimates.augment, estimates.onTest)
  resList[['estimates.onTest']] <- estimates.onTest
  
  resList[['testAcc.aug']] <- seed2TestAcc <- ADAPTS::calcAcc(estimates=estimates.onTest[,1], reference=estimates.onTest[,ncol(estimates.onTest)])
  
  #shrink
  if(skipShrink == FALSE) {
    augTrain.shrink <- ADAPTS::shrinkSigMatrix(sigMatrix=augTrain, numChunks=NULL, verbose=FALSE, plotIt = FALSE, aggressiveMin=TRUE,sigGenesList=NULL, fastStop=TRUE, singleCore=TRUE)
    #pheatmap(augTrain.shrink)
    resList[['matrix.shrink']] <- augTrain.shrink
    
    estimates.shrink <- as.data.frame(ADAPTS::estCellPercent.DCQ(augTrain.shrink, pseudobulk.test))
    colnames(estimates.shrink) <- paste('Shrunk Matrix')
    resList[['estimates.onTest']] <- estimates.onTest <- cbind(estimates.shrink, estimates.onTest)
    
    resList[['testAcc.shrink']] <- seed2TestAcc<- ADAPTS::calcAcc(estimates=estimates.onTest[,1], reference=estimates.onTest[,ncol(estimates.onTest)])}
  
  
  #Clustering
  
  if(!is.null(handMetaCluster)) {
    resList[['allClusters']]  <- loadHandClusters()
  } else {
    varClusts <- ADAPTS::clustWspillOver(sigMatrix = augTrain.shrink, geneExpr = trainSet.3sam)
    resList[['allClusters']]  <- varClusts$allClusters
  }
  
  resList[['allClusters']] 
  
  metaCluster.id <- list()
  for(i in 1:length(resList[['allClusters']])) {
    for (x in resList[['allClusters']][[i]]) {
      metaCluster.id[[x]] <- paste('Meta',i,sep='_')
    }
  }
  metaClust.LUT <- resList[['metaClust.LUT']]  <- unlist(metaCluster.id)
  names(resList[['allClusters']]) <- metaClust.LUT[sapply(resList[['allClusters']], function(x){x[1]})]
  
  #Update 05-19-20: More informative names
  metaNames <- sapply(unique(metaClust.LUT), function(x) { paste(names(metaClust.LUT)[metaClust.LUT == x], collapse='_')})
  metaClust.LUT.MI <- metaClust.LUT
  metaClust.LUT.MI <- metaNames[match(metaClust.LUT.MI, names(metaNames))]
  names(metaClust.LUT.MI) <- names(metaClust.LUT)
  
  metatrainSet<-trainSet
  metatestSet<-testSet
  
  colnames(metatrainSet) <- metaClust.LUT.MI[sub('\\..*','',colnames(metatrainSet))]
  colnames(metatestSet) <- metaClust.LUT.MI[sub('\\..*','',colnames(metatestSet))]
  
  metatrainSet.3sam <- ADAPTS::scSample(RNAcounts = metatrainSet, groupSize = 3, randomize = TRUE)
  metatrainSet.30sam <- ADAPTS::scSample(RNAcounts = metatrainSet, groupSize = 30, randomize = TRUE)
  
  metaclusterIDs <- factor(colnames(metatrainSet.30sam))
  metatrainSet.4reg <- t(metatrainSet.30sam)
  
  metapseudobulk.test <- data.frame(test=rowSums(metatestSet))
  metapseudobulk.test.counts <- table(sub('\\..*','',colnames(metatestSet)))
  meta.actFrac <- 100 * metapseudobulk.test.counts / sum(metapseudobulk.test.counts)
  
  #Metaseed
  genesInSeed<-100
  metaseedMat <-buildSeed(metatrainSet, genesInSeed=genesInSeed, groupSize=30, randomize=TRUE, num.trees=1000, plotIt=FALSE, trainSet.3sam=metatrainSet.3sam, trainSet.30sam=metatrainSet.30sam,proportional = proportional)
  
  resList[['matrix.meta']] <- metaseedMat
  
  estimates.Meta.onTest <- as.data.frame(ADAPTS::estCellPercent.DCQ(metaseedMat, metapseudobulk.test))
  
  colnames(estimates.Meta.onTest) <- paste(genesInSeed, 'Marker Genes Seed')
  estimates.Meta.onTest$actFrac.test <- round(meta.actFrac[rownames(estimates.Meta.onTest)],2)
  
  resList[['estimates.onTest.meta']] <- estimates.Meta.onTest
  
  resList[['testAcc.meta']] <- seed2TestAcc.meta <- ADAPTS::calcAcc(estimates=estimates.Meta.onTest[,1], reference=estimates.Meta.onTest[,2])
  
  
  #meta all gene
  metaallGeneSig <- apply(metatrainSet.3sam, 1, function(x){tapply(x, colnames(metatrainSet.3sam), mean, na.rm=TRUE)})
  metaestimates.allGene <- as.data.frame(ADAPTS::estCellPercent.DCQ(t(metaallGeneSig), metapseudobulk.test))
  colnames(metaestimates.allGene)<-paste('All Gene Meta')
  estimates.Meta.onTest <- cbind(metaestimates.allGene, estimates.Meta.onTest)
  resList[['estimates.onTest.meta']] <- estimates.Meta.onTest
  
  resList[['testAcc.metaAll']] <- seed2TestAcc <- ADAPTS::calcAcc(estimates=estimates.Meta.onTest[,1], reference=estimates.Meta.onTest[,ncol(estimates.Meta.onTest)])
  
  #meta aug
  metagList <- ADAPTS::gListFromRF(trainSet=metatrainSet.30sam)
  sapply(gList,dim)
  
  meta.augTrain <- ADAPTS::AugmentSigMatrix(origMatrix = metaseedMat, fullData = metatrainSet.3sam, gList = metagList, nGenes = 1:100, newData = metatrainSet.3sam, plotToPDF = FALSE, pdfDir = '.')
  
  resList[['matrix.metaAug']] <- meta.augTrain
  
  estimates.Meta.augment <- as.data.frame(ADAPTS::estCellPercent.DCQ(meta.augTrain, metapseudobulk.test))
  colnames(estimates.Meta.augment) <- paste('Augmented Meta')
  resList[['estimates.onTest.meta']] <- estimates.Meta.onTest <- cbind(estimates.Meta.augment, estimates.Meta.onTest)
  
  resList[['testAcc.metaAug']] <- seed2TestAcc.aug.meta <- ADAPTS::calcAcc(estimates=estimates.Meta.onTest[,1], reference=estimates.Meta.onTest[,ncol(estimates.Meta.onTest)])
  
  #meta shrink
  if(skipShrink == FALSE) {
    gc()
    meta.augTrain.shrink <- ADAPTS::shrinkSigMatrix(sigMatrix=meta.augTrain, numChunks=NULL, verbose=FALSE, plotIt = FALSE, aggressiveMin=TRUE, sigGenesList=NULL,fastStop=TRUE, singleCore=TRUE)
    dim(meta.augTrain.shrink)
    #pheatmap(augTrain.shrink)
    resList[['matrix.metaAugShrink']] <- meta.augTrain.shrink
    
    estimates.Meta.shrink <- as.data.frame(ADAPTS::estCellPercent.DCQ(meta.augTrain.shrink, metapseudobulk.test))
    colnames(estimates.Meta.shrink) <- paste('Shrunk Meta')
    resList[['estimates.onTest.meta']] <- estimates.Meta.onTest <- cbind(estimates.Meta.shrink, estimates.Meta.onTest)
    
    resList[['testAcc.metaAugShrink']] <- seed2TestAcc.shrink.meta <- ADAPTS::calcAcc(estimates=estimates.Meta.onTest[,1], reference=estimates.Meta.onTest[,ncol(estimates.Meta.onTest)])
    
  }
  
  return(resList)
}


#' Find out at which iteration the results converge, i.e. the mean results are stable.
#'
#' @param curSeq A sequence of results that generated from each iteration of the loop
#' @param changePer The maximum percentage of change allowed
#' @param winSize  The window size for mean calculation
#'
#' @return The minimum number of iterations needed for the results to converge
#' @export
#'
#' @examples
findConvergenceIter <- function(curSeq, changePer=1, winSize=5) {
  
  #Note, this will remove NAs.  Is that best?  They're caused by bad correlations
  runMean <- sapply(1:length(curSeq), function(x) {sum(curSeq[1:x], na.rm=TRUE)/x})
  #Criteria: running mean has changed less than 5% in the last 5? point
  convIter <- as.numeric(NA)
  if (length(runMean) > winSize) {
    winOffset <- winSize-1
    maxWinChange <- sapply(winSize:length(runMean), function(x) {
      win <- runMean[(x-winOffset):x]
      max(abs((win - mean(win))/mean(win)))
    }) #maxWinChange
    mwcBool <- maxWinChange < changePer/100
    if(any(mwcBool)) {convIter <- winOffset + which(mwcBool)[1]}
  }
  return (convIter)
}



#' A meta analysis for the results from multiple iterations
#' @description Calcualte the mean and the standard deviation of the reaults from all the iterations, and also 
#' test for convergence by % of change with each additional iteration.
#' 
#' @param allResList A list of results generated from all the iterative calls of scMatrixTest
#' @param changePer  The maximum percentage of change allowed for convergence
#'
#' @return The mean and standard deviation of all the results, along with the mumber of iterations needed for the results to converge.
#' @export
#' @examples

meanResults <- function (allResList,changePer=1) {
  testNames <- unique(sub('^.*\\.', '', names(allResList[[1]])))
  testNames <- testNames[!testNames %in% c("onTest", "allClusters", "LUT")]
  compTypes <- names(allResList[[1]][[paste0('testAcc.', testNames[1])]])
  
  allResList <- allResList[!sapply(allResList, function(x){inherits(x,'try-error')})]
  
  compList <- list()
  for (curName in testNames) {
    compList[[curName]] <- list()
    for (curComp in compTypes) {
      compList[[curName]][[curComp]] <- sapply(allResList, function(x){ 
        y <- x[[paste0('testAcc.', curName)]]
        if(!inherits(y, 'try-error')) {return(y[[curComp]])} else {return(as.numeric(NA))}  
        #If only two clusters, set correlation to zero or NA?
      }) #compList[[curName]][[curComp]] <- sapply(allResList, function(x){ 
    } #for (curComp in compTypes) {
  } #for (curName in testNames) {
  
  curColN <- unlist(lapply(compTypes, function(x){c(x, paste0(x,'.sd'))}))
  resMat <- matrix(as.numeric(NA), nrow=length(testNames), ncol=length(curColN),
                   dimnames=list(testNames, curColN))
  convMat <- matrix(as.numeric(NA), nrow=length(testNames), ncol=length(compTypes),
                    dimnames=list(testNames, compTypes))
  for (curName in testNames) {
    for (curComp in compTypes) {
      curMean <- mean(compList[[curName]][[curComp]], na.rm = TRUE)
      curSD <- sd(compList[[curName]][[curComp]], na.rm = TRUE)
      resMat[curName,curComp] <- curMean
      resMat[curName,paste0(curComp,'.sd')] <- curSD
      
      #Also test for convergence by % of change with each additional
      convMat[curName, curComp] <- findConvergenceIter(curSeq=compList[[curName]][[curComp]], changePer=changePer, winSize=5)
      
    } #for (curComp in compTypes) {
  } #for (curName in testNames) {
  colnames(convMat) <- paste0('convIt.',colnames(convMat))
  
  resMat <- as.data.frame(resMat)
  resMat$N <- length(allResList)
  resMat <- cbind(resMat, as.data.frame(convMat))
  
  return(resMat)
}




#' Loop scMatrixtest until convergence
#' @description Iteratively call scMatrixtest numLoops times with the option to fast stop 
#' if correlation, correlation spear, mae and rmse all converge
#' 
#' @param numLoops The number of iterations. Set to null to loop until results converge.
#' @param fastStop Set to TRUE to break the loop when correlation, correlation spear, mae and rmse all converge.
#' @param exprData The single cell matrix
#' @param changePer The maximum percentage of change allowed for convergence#' 
#' @param handMetaCluster Load in pre-defined meta clusters. Set to NULL to automatically group indistinguishable 
#' cells into same cluster use clustWspillOver(DEFAULT: NULL)

#'
#' @return  A list of results generated from all the iterative calls of scMatrixTest
#' @export
#'
#' @examples
loopTillConvergence<-function(numLoops,fastStop,exprData,changePer,handMetaCluster){
  if(is.null(numLoops)){
    fastStop<-TRUE
    numLoops<-1000000
  }
  allResListOut <- list()
  for (i in 1:numLoops) {
    curName <- paste0('res', i)
    allResListOut[[curName]] <- try(scMatrixTest(exprData, randomize = TRUE, proportional=FALSE, handMetaCluster=handMetaCluster))
    if(fastStop==TRUE){
      covtmp<-meanResults(allResList=allResListOut,changePer)[ ,c("convIt.rho.cor", "convIt.spear.rho", "convIt.mae","convIt.rmse")]
      if(all(!is.na(covtmp))) break
    }}
  return(allResListOut)
}
