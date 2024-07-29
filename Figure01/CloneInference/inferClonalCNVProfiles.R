source("~/MB_SCSEQ/Figure01/CloneInference/commonSegmentation.R")

runHMMwithPar <- function(x, e, strength, mu, lambda, nu, kappa, m, eta, gamma, S){
  def.par <- HMMsegment(x, getparam = T)
  if(!missing(e)) def.par$e <- e
  if(!missing(strength)) def.par$strength <- strength
  if(!missing(mu)) def.par$mu <- mu
  if(!missing(lambda)) def.par$lambda <- lambda
  if(!missing(nu)) def.par$nu <- nu
  if(!missing(kappa)) def.par$kappa <- kappa
  if(!missing(m)) def.par$m <- m
  if(!missing(eta)) def.par$eta <- eta
  if(!missing(gamma)) def.par$gamma <- gamma
  if(!missing(S)) def.par$S <- S

  return(HMMsegment(x, param=def.par))
}


runHMMwithParTable <- function(x, params){
  def.par <- HMMsegment(x, getparam = T)
  return(HMMsegment(x, param=params))
}



correctReadcount2 <- function (x, mappability = 0.9, samplesize = 100000, verbose = TRUE) {
  if (length(x$reads) == 0 | length(x$gc) == 0 | length(x$map) == 
      0) {
    stop("Missing one of required columns: reads, gc, map")
  }
  if (verbose) {
    message("Applying filter on data...")
  }
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.05
  range <- quantile(x$reads[x$valid], prob = c(0, 1 - routlier), 
                    na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid], prob = c(doutlier, 1 - 
                                               doutlier), na.rm = TRUE)
  x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] | 
            x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  if (verbose) {
    message("Correcting for GC bias...")
  }
  set <- which(x$ideal)
  select <- sample(set, min(length(set), samplesize))
  rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads/predict(final, x$gc)
  if (verbose) {
    message("Correcting for mappability bias...")
  }
  coutlier <- 0.05
  range <- quantile(x$cor.gc[which(x$valid)], prob = c(0, 
                                                       1 - coutlier), na.rm = TRUE)
  set <- which(x$cor.gc <= range[2] & x$cor.gc > range[1] & x$valid & x$map >= mappability )
  select <- sample(set, min(length(set), samplesize))
  final = approxfun(lowess(x$map[select], x$cor.gc[select]), ties = mean)
  x$cor.map <- x$cor.gc/final(x$map)
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(x)
}


inferClonalCNVs <- function(QDNAseq, HMMCopyCellList, ploidyFilter=0.05, mergeWidth=2, kRange=1:12, 
                            cluster.method=c("ManSegments", "ManBins", "CorBins"), 
                            bootstrap=FALSE, nBoot=49, minCells=5){
  
  cluster.method <- match.arg(cluster.method)
  
  ploidies.to.keep <- as.numeric(names(which((prop.table(table(QDNAseq$ploidy.mod))>ploidyFilter) & table(QDNAseq$ploidy.mod) > 1)))
  
  QDNAseq <- QDNAseq[,which(QDNAseq$ploidy.mod %in% ploidies.to.keep)]

  HMMCopyCellList <- HMMCopyCellList[names(HMMCopyCellList) %in% colnames(QDNAseq)]

  ploidy.res <- list()
  for (cur.ploidy in ploidies.to.keep){

    QDNAseq.cur <- QDNAseq[,QDNAseq$ploidy.mod == cur.ploidy]
    HMMCopyCellList.cur <- HMMCopyCellList[names(HMMCopyCellList) %in% colnames(QDNAseq.cur)]

    switch(cluster.method, ManSegments = {
      cur.man.wardd2 <- clusterMannhattanOnSegments(QDNAseq.cur, mergeWidth)
    }, ManBins = {
      cur.man.wardd2 <- clusterMannhattanOnBins(QDNAseq.cur)
      
    }, CorBins = {
      cur.man.wardd2 <- clusterCorOnBins(QDNAseq.cur)
      
    }
    )
    
    
    ploidy.res[[paste0("Ploidy_", cur.ploidy)]] <- inferCNVforEachCloneOverKRange(QDNAseq.cur, HMMCopyCellList.cur, cur.man.wardd2, kRange = kRange, bootstrap = bootstrap, nBoot, minCells=minCells)
  }
  # tmp.names <- names(ploidy.res[[1]])
  # ploidy.res <- lapply(seq_along(ploidy.res[[1]]), function(x) return(lapply(ploidy.res, `[[`, x)))
  # names(ploidy.res) <-  tmp.names

  return(ploidy.res)
}


clusterMannhattanOnSegments <- function(QDNAseq.cur, mergeWidth){
  cur.segmented <- mapCellsToCommonBreakpoints(QDNAseq.cur, mergeWidth = mergeWidth)
  
  cur.dist <- dist(t(cur.segmented[, c(-1, -2, -3)]), method = "manhattan")
  
  cur.man.wardd2 <- hclust(cur.dist, method = "ward.D2")
  
  return(cur.man.wardd2)
}

clusterMannhattanOnBins <- function(QDNAseq.cur){
  cur.bins <- assayDataElement(QDNAseq.cur, "copynumber")
  
  cur.dist <- dist(t(cur.bins), method = "manhattan")
  
  cur.man.wardd2 <- hclust(cur.dist, method = "ward.D2")
  
  return(cur.man.wardd2)
}

clusterCorOnBins <- function(QDNAseq.cur){
  require(coop)
  cur.bins <- assayDataElement(QDNAseq.cur, "copynumber")
  
  cur.dist <- as.dist(1-coop::pcor(cur.bins, use="complete.obs"))
  
  cur.man.wardd2 <- hclust(cur.dist, method = "ward.D2")
  
  return(cur.man.wardd2)
}

inferCNVforEachCloneOverKRange <- function(QDNAseq.cur, HMMCopyCellList.cur, cellClust, kRange=1:12, bootstrap, nBoot, minCells=5){

    resList <- list()

    for(k in kRange){
      if(k > ncol(QDNAseq.cur)){
        resList[[paste0("k_", k)]] <- list("Cluster_NA" = "Too few cells to cut into this number of clusters")
        next
      }
      cur.clustering <- cutree(cellClust, k=k)

      k.res.list <- list()

      for(clust in unique(cur.clustering)){
        if(sum(cur.clustering == clust)>=minCells){
          print(paste0("Processing Cluster cluster ", clust, " with ", sum(cur.clustering == clust), " cells for ploidy ", unique(QDNAseq.cur$ploidy.mod), " and k = ", k))
          k.res.list[[paste0("Cluster_", clust)]] <- c(inferClonalCNVfromSubclone(QDNAseq.cur[, cur.clustering == clust], HMMCopyCellList.cur[cur.clustering == clust], bootstrap, nBoot), 
                                                       cells=list(colnames(QDNAseq.cur[, cur.clustering == clust])))
        } else {
          k.res.list[[paste0("Cluster_",clust)]] <- c("Too few cells in cluster", cells=list(colnames(QDNAseq.cur[, cur.clustering == clust])))
        }
      }
      resList[[paste0("k_", k)]] <- k.res.list
    }
    return(resList)
}



inferClonalCNVfromSubclone <- function(QDNASeq, HMMCopyCellList, bootstrap, nBoot){#, minCNVStatePrevelance=1e-3){

    stopifnot("All cells don't have the same ploidy in a clone" = 
              length(unique(QDNASeq$ploidy.mod)) == 1)
    
    ## Here we work with ploidy defined as the average copy number of the genome, as opposed to the number of chromosome sets
    ## This is because most of the genome in some of our cells is copy number altered, so the difference between the two is not negligible

    clone.ploidy = mean(QDNASeq$ploidy)
    
    ## Make the pseudobulk object
    pseudobulk <- HMMCopyCellList[[1]][,1:6]
    for(i in seq_along(HMMCopyCellList)){
        if(i == 1) next
        pseudobulk$reads <- pseudobulk$reads + HMMCopyCellList[[i]]$reads
    }

    pseudobulk <- pseudobulk[order(factor(chr, levels = mixedsort(unique(chr))), start)]
    pseudobulk <- correctReadcount2(pseudobulk, verbose = FALSE)

    def.par <- HMMsegment(pseudobulk, getparam = T)


    ## Get the available copy number states at the estimated ploidy:
    proptableOfStates <- prop.table(table(log2(assayData(QDNASeq)$copynumber/clone.ploidy)))
    availableStates <- as.numeric(names(proptableOfStates))
    my.par <- def.par[rep(1,length(availableStates)),]
    my.par$mu <- availableStates

    my.par <- my.par[order(my.par$mu),] 
    my.par[my.par$mu==0,] <- def.par[which.min(abs(def.par$mu)),]

    my.par$eta <- 0
    my.par$m <- my.par$mu
    my.par$e <- 0.999999999999999
    my.par$strength <- 1e+50
    
    if(any(my.par$m==-Inf)){
      my.par$eta[my.par$mu==-Inf] <- 5
      if(min(availableStates[is.finite(availableStates)]) < 0){
        my.par$mu[my.par$mu==-Inf] <-my.par$m[my.par$m==-Inf] <- min(availableStates[is.finite(availableStates)])*3
      } else {
        my.par$mu[my.par$mu==-Inf] <-my.par$m[my.par$m==-Inf] <- min(pseudobulk$copy, na.rm=TRUE)
      }
    }
    
    rownames(my.par) <- sort(round(2^(my.par$mu)*clone.ploidy))
    
    pseudobulk.hmmres <- HMMsegment(pseudobulk, param = my.par, verbose = FALSE)


    if(bootstrap){
      bootRes <- bootstrapPseudobulkRes(HMMCopyCellList, my.par, nBoot)
    } else {
      bootRes <- list("Bootstrap Not Performed")
    }

    return(list("pseudobulk"=pseudobulk, "hmmres"=pseudobulk.hmmres, "boot" = bootRes))
}


bootstrapPseudobulkRes <- function(HMMCopyCellList, my.par, nBoot=49){
    require(R.utils)
    message("Boostrapping...")
    bootResList <- list()
    # pb <- R.utils::ProgressBar(max=nBoot)
    # R.utils::reset(pb)
    bootResList <- foreach(jj = seq_len(nBoot)) %dopar% {
      boot.sample <- sample(length(HMMCopyCellList), replace = TRUE)
      pseudobulk <- HMMCopyCellList[[1]][, 1:6]
      for (i in boot.sample) {
        if (i == 1) next
        pseudobulk$reads <- pseudobulk$reads + HMMCopyCellList[[i]]$reads
      }
      pseudobulk <- pseudobulk[order(factor(chr, levels = mixedsort(unique(chr))), start)]
      pseudobulk <- correctReadcount2(pseudobulk, verbose = FALSE)
      pseudobulk.hmmres <- HMMsegment(pseudobulk, param = my.par, verbose = FALSE)
      pseudobulk.hmmres$rho <-  matrix()
      return(list("sampled_cells" = boot.sample, hmmcopy.res = pseudobulk.hmmres))
    # R.utils::increase(pb)
    }
    return(bootResList)
}


inferClonalCNVfromSubcloneMedian <- function(QDNASeq, HMMCopyCellList){#, minCNVStatePrevelance=1e-3){
  
  stopifnot("All cells don't have the same ploidy in a clone" = 
              length(unique(QDNASeq$ploidy.mod)) == 1)
  
  clone.ploidy = unique(QDNASeq$ploidy.mod)

  
  ## Make the pseudobulk object
  pseudobulk <- HMMCopyCellList[[1]][,1:6]
  
  pseudobulk$reads <- apply(sapply(HMMCopyCellList, `[[`, "reads"),1,median, na.rm=TRUE) ## QUESTION: should I do this before or after correction??
  
  pseudobulk <- pseudobulk[order(factor(chr, levels = mixedsort(unique(chr))), start)]
  pseudobulk <- correctReadcount2(pseudobulk, verbose = FALSE)
  
  def.par <- HMMsegment(pseudobulk, getparam = T)
  
  proptableOfStates <- prop.table(table(log2(assayData(QDNASeq)$copynumber/clone.ploidy)))
  # proptableOfStates <- proptableOfStates[proptableOfStates >= minCNVStatePrevelance]
  availableStates <- as.numeric(names(proptableOfStates))
  # availableStates <- availableStates[is.finite(availableStates)]
  my.par <- def.par[rep(1,length(availableStates)),]
  my.par$mu <- availableStates
  
  my.par <- my.par[order(my.par$mu),] 
  my.par[my.par$mu==0,] <- def.par[which.min(abs(def.par$mu)),]
  
  my.par$eta <- 5
  my.par$m <- my.par$mu
  
  if(any(my.par$m==-Inf)){
    if(min(availableStates[is.finite(availableStates)]) < 0){
      my.par$mu[my.par$mu==-Inf] <-my.par$m[my.par$m==-Inf] <- min(availableStates[is.finite(availableStates)])*3
    } else {
      my.par$mu[my.par$mu==-Inf] <-my.par$m[my.par$m==-Inf] <- min(pseudobulk$copy, na.rm=TRUE)
    }
  }
  
  rownames(my.par) <- sort(round(2^(my.par$mu)*clone.ploidy))
  
  pseudobulk.hmmres <- HMMsegment(pseudobulk, param = my.par, verbose = FALSE)
  
  return(list("pseudobulk"=pseudobulk, "hmmres"=pseudobulk.hmmres))
}


inferClonalCNVfromSubcloneMedian2 <- function(QDNASeq, HMMCopyCellList){#, minCNVStatePrevelance=1e-3){
  
  stopifnot("All cells don't have the same ploidy in a clone" = 
              length(unique(QDNASeq$ploidy.mod)) == 1)
  
  clone.ploidy = unique(QDNASeq$ploidy.mod)
  # 
  # for(i in seq_along(HMMCopyCellList)){
  #   HMMCopyCellList[[i]] <- correctReadcount2(HMMCopyCellList[[i]])
  # }
  
  ## Make the pseudobulk object
  pseudobulk <- HMMCopyCellList[[1]]
  
  pseudobulk$copy <- apply(sapply(HMMCopyCellList, `[[`, "copy"),1,median, na.rm=TRUE) ## QUESTION: should I do this before or after correction??
  
  pseudobulk <- pseudobulk[order(factor(chr, levels = mixedsort(unique(chr))), start)]
  # pseudobulk <- correctReadcount2(pseudobulk, verbose = FALSE)
  
  def.par <- HMMsegment(pseudobulk, getparam = T)
  
  proptableOfStates <- prop.table(table(log2(assayData(QDNASeq)$copynumber/clone.ploidy)))
  # proptableOfStates <- proptableOfStates[proptableOfStates >= minCNVStatePrevelance]
  availableStates <- as.numeric(names(proptableOfStates))
  # availableStates <- availableStates[is.finite(availableStates)]
  my.par <- def.par[rep(1,length(availableStates)),]
  my.par$mu <- availableStates
  
  my.par <- my.par[order(my.par$mu),] 
  my.par[my.par$mu==0,] <- def.par[which.min(abs(def.par$mu)),]
  
  my.par$eta <- 5
  my.par$m <- my.par$mu

  
  
  if(any(my.par$m==-Inf)){
    if(min(availableStates[is.finite(availableStates)]) < 0){
      my.par$mu[my.par$mu==-Inf] <-my.par$m[my.par$m==-Inf] <- min(availableStates[is.finite(availableStates)])*3
    } else {
      my.par$mu[my.par$mu==-Inf] <-my.par$m[my.par$m==-Inf] <- min(pseudobulk$copy, na.rm=TRUE)
    }
  }
  
  rownames(my.par) <- sort(round(2^(my.par$mu)*clone.ploidy))
  
  pseudobulk.hmmres <- HMMsegment(pseudobulk, param = my.par, verbose = FALSE)
  
  return(list("pseudobulk"=pseudobulk, "hmmres"=pseudobulk.hmmres))
}



addChrOffsetToDF <- function(correctOutput){
  correctOutput$chr <- factor(correctOutput$chr , mixedsort(unique(correctOutput$chr)))

  chromosome_lengths <- aggregate(data=correctOutput, end~chr , max)
  
  chromosome_offsets <- c(0,cumsum(as.numeric(chromosome_lengths$end))[-nrow(chromosome_lengths)])
  names(chromosome_offsets) <- mixedsort(unique(correctOutput$chr))
  
  correctOutput$len <- correctOutput$end - correctOutput$start -1
  correctOutput$start.x <- correctOutput$start 
  correctOutput$end.x <- correctOutput$end -1
  correctOutput[,end.x := end.x + chromosome_offsets[chr]]
  correctOutput[,start.x := start.x + chromosome_offsets[chr]]
  
  return(list(correctOutput, chromosome_offsets))

}

require(cowplot)
plotClonalProfilesAtK <- function(clonal.cnv.res, k=1){


  ploidy.plot.list <- list()

  for(ploidy in names(clonal.cnv.res)){

    ploidy.res <- clonal.cnv.res[[ploidy]]
    ploidy.num <- as.numeric(strsplit(ploidy, "_")[[1]][2])
    
    k.res <- ploidy.res[[paste0("k_", k)]]
    
    clusterPlots <- lapply(names(k.res), function(clust) {
        clust.res <- k.res[[clust]]
    
        if(is.character(clust.res[[1]])){
          return(NA_character_)
        }
        if(is.matrix(clust.res$hmmres$mu)){
          stateCorrected <- round(2^clust.res$hmmres$mu[,1]*ploidy.num)
        } else {
          stateCorrected <- round(2^clust.res$hmmres$mu*ploidy.num)
        }
        toPlot <- cbind(clust.res$pseudobulk, "state"=stateCorrected[clust.res$hmmres$state])
        tmp <- addChrOffsetToDF(toPlot)
        toPlot <- tmp[[1]]
        chromosome_offsets <- tmp[[2]]
        toPlot$unlog_copy <- 2^toPlot$copy*ploidy.num
        toLabel <- data.frame(x=chromosome_offsets+toPlot[,max(start)/2, chr][,V1], chr = names(chromosome_offsets))
        p <- ggplot(toPlot, aes(x=start.x)) + geom_point(mapping=aes(y=unlog_copy), size=0.2, alpha=0.2) + geom_step(aes(y=state,group=chr), color="green") +
                  theme_classic() + xlab("Chromosome Position") + ylab("Copy Number") +
                  geom_vline(xintercept=chromosome_offsets[-1]) + scale_y_continuous(limits = c(0,8)) + 
          ggtitle(paste0("Clone ", clust, " for ploidy ", ploidy, " with N=", length(clust.res$cells), " cells.")) +
          geom_text(data=toLabel, aes(x=x, y=7.5, label=chr))
        return(p)
    })    
    ploidy.plot.list[[ploidy]] <- clusterPlots[!sapply(clusterPlots, anyNA)]
  }
  plot_grid(plotlist=unlist(ploidy.plot.list, recursive = FALSE), ncol = 1)
}


## Note: I don't think that mapCellsToCommonBreakpoints(mapCellsToCommonBreakpointsInner(x), mapCellsToCommonBreakpointsInner(y)) 
## will return the same as mapCellsToCommonBreakpointsInner(cbind(x,y)), because "chains" of breakpoints are different across the two...
## Need to think about whether/how to fix this...

quantifyDifferencesForKClusters <- function(clonal.cnv.res, k=1){
  
  ploidyClusterStates <- list()
  ploidyClusterStatesNormalized <- list()
  
  for(ploidy in names(clonal.cnv.res)){
    
    ploidy.res <- clonal.cnv.res[[ploidy]]
    ploidy.num <- as.numeric(strsplit(ploidy, "_")[[1]][2])
    
    k.res <- ploidy.res[[paste0("k_", k)]]
    
    
    k.res <- k.res[sapply(k.res, function(x) return(!is.character(x[[1]])))]
    
    clusterStates <- sapply(k.res, function(res){
      
      if(is.matrix(res$hmmres$mu)){
        stateCorrected <- round(2^res$hmmres$mu[,1]*ploidy.num)
      } else {
        stateCorrected <- round(2^res$hmmres$mu*ploidy.num)
      }
      
      
      out <- stateCorrected[res$hmmres$state]
      
      names(out) <- paste0(res$pseudobulk$chr, ":", res$pseudobulk$start, "-", res$pseudobulk$end)
      out
    })
    colnames(clusterStates) <- paste0(colnames(clusterStates), "_", ploidy)
    ploidyClusterStates[[ploidy]] <- clusterStates
    ploidyClusterStatesNormalized[[ploidy]] <- clusterStates/ploidy.num
    
  }

  
  
  commonClusterSegs <- mapCellsToCommonBreakpointsInner(do.call(cbind,ploidyClusterStates), median(as.numeric(sapply(strsplit(names(k.res), "_"), `[[`,2))))
  commonClusterSegsNormed <- mapCellsToCommonBreakpointsInner(do.call(cbind,ploidyClusterStatesNormalized), 
                                                              median(as.numeric(sapply(strsplit(names(k.res), "_"), `[[`,2))))
  
  return(list("unnormed_dists" = dist(t(commonClusterSegs[,-c(1,2,3), with=FALSE]), "manhattan"),
         "ploidy_normed_dists" =   dist(t(commonClusterSegsNormed[,-c(1,2,3), with=FALSE]), "manhattan")))
  
  
}

#metricsToReturn = c("RMSE", "AIC", "BIC")
computeGOFMetrics <- function(clonal.cnv.res){
  
  
  ploidyClusterOut <- list()
  
  for(ploidy in names(clonal.cnv.res)){
    
    ploidy.res <- clonal.cnv.res[[ploidy]]
    ploidy.num <- as.numeric(strsplit(ploidy, "_")[[1]][2])
    
    k.out <- list()
    
    for(k.char in names(ploidy.res)){
      
      k.res <- ploidy.res[[k.char]]
      
      perClusterStats <- sapply(names(k.res), function(clust){
        clust.res <- k.res[[clust]]
        if(is.character(clust.res[[1]])){
          return(c("rmse" = NA_real_, "rmse_log" = NA_real_, "logLik" = NA_real_, numPar = NA_real_, nBins = NA_real_))
        }
        if(is.matrix(clust.res$hmmres$mu)){
          stateCorrected <- round(2^clust.res$hmmres$mu[,ncol(clust.res$hmmres$mu)]*ploidy.num)
        } else {
          stateCorrected <- round(2^clust.res$hmmres$mu*ploidy.num)
        }
        state <- stateCorrected[clust.res$hmmres$state]
        unlog_copy <- 2^clust.res$pseudobulk$copy*ploidy.num
        
        rmse <- sqrt(mean((state-unlog_copy)^2, na.rm=TRUE))
        if(is.matrix(clust.res$hmmres$mu)){
          rmse_log <- sqrt(mean((clust.res$hmmres$mu[clust.res$hmmres$state,ncol(clust.res$hmmres$mu)]-clust.res$pseudobulk$copy)^2, na.rm=TRUE))
        } else {
          rmse_log <- sqrt(mean((clust.res$hmmres$mu[clust.res$hmmres$state]-clust.res$pseudobulk$copy)^2, na.rm=TRUE))
        }
        
        logLik <- clust.res$hmmres$loglik[length(clust.res$hmmres$loglik)]
        numPar <- nrow(clust.res$hmmres$segs)*2 # parameters estimated are: # of breakpoints (1), location of each breakpoint (#segs-1), and state to assign to each seg (#segs)
        return(c("rmse" = rmse, "rmse_log" = rmse_log, "logLik" = logLik, numPar = numPar, nBins = length(unlog_copy)))
      })
      
      # perClusterStats <- perClusterStats[,!is.na(perClusterStats[1,])]
      totalRMSE <- sqrt(mean(perClusterStats["rmse", ]^2, na.rm=TRUE))
      totalRMSElog <- sqrt(mean(perClusterStats["rmse_log", ]^2, na.rm=TRUE))
      totalLogLik <- sum(perClusterStats["logLik", ], na.rm=TRUE)
      totalNumPar <- sum(perClusterStats["numPar", ], na.rm=TRUE)
      totalNumBins <- sum(perClusterStats["nBins", ], na.rm=TRUE)
      totalAIC <- 2*totalNumPar - 2*totalLogLik
      totalBIC <- totalNumPar*log(totalNumBins) - 2*totalLogLik
      k.out[[k.char]] <- list('totalRes' = c("RMSE" = totalRMSE, "logRMSE" = totalRMSElog, "logLik" = totalLogLik, numPar = totalNumPar, nBins = totalNumBins, AIC = totalAIC, BIC = totalBIC), 
                              'perClustRes' = perClusterStats)
    }
    ploidyClusterOut[[ploidy]] <- k.out
  }
  return(ploidyClusterOut)
}


