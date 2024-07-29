getWithinClusterBootDistances <- function(clust.res, distance="manhattan") {
  boot.res <- clust.res$boot
  boot.states <- sapply(boot.res, function(x) return(x$hmmcopy.res$state))
  return(as.vector(dist(t(boot.states))))
}

getWithinClusterVariance <- function(clust.res) {
  boot.res <- clust.res$boot
  boot.states <- sapply(boot.res, function(x) return(x$hmmcopy.res$state))
  withinClustDistances <- as.vector(dist(t(boot.states)))
  withinClustVar <- 1/(length(boot.res)*(length(boot.res)-1))*sum(withinClustDistances^2) 
  return(withinClustVar)
}


getWithinClusterSumSquareDev <- function(clust.res) {
  boot.res <- clust.res$boot
  boot.states <- sapply(boot.res, function(x) return(x$hmmcopy.res$state))
  mean.state <- rowMeans(boot.states)
  
  withinClustSSD <- sum(colSums((mean.state-boot.states)^2)) 
  return(withinClustSSD)
}

getNumFstat <- function(k.res){
    clust.mean.states <- sapply(k.res, function(x) {
      if(!is.character(x[[1]])) {
        boot.res <- x$boot
        boot.states <- sapply(boot.res, function(x) return(x$hmmcopy.res$state))
        return(rowMeans(boot.states))
      } else{
        return(NULL)
      }
      })

    mean.state <- rowMeans(clust.mean.states)
    
    res <- sum(apply(clust.mean.states, 2, function(x) {
      return(sum((x-mean.state)^2)*length(k.res[[1]]$boot))
    }))/(ncol(clust.mean.states)-1)
    
    return(res)
}

getBetweenClusterVariance <- function(k.res){
  
  pairs.of.clusts <- combn(k.res, 2)
  
  betweenClustDistances <- sapply(seq_len(ncol(pairs.of.clusts)), function(x) {
    getBetweenClusterDistances(pairs.of.clusts[, x])
  })
  
  return(1/(length(k.res)*(length(k.res)-1))*sum(betweenClustDistances^2))
}

getBetweenClusterDistances <- function(k.res){
  clust.mean.states <- sapply(k.res, function(x) {if(!is.character(x[[1]])) return(x$hmmres$state) else return(NULL)})
  
  return(as.vector(dist(t(clust.mean.states))))
}


# 
# 
# returnBootstrapTlikeDistsForK <- function(k.res){
    # k.res <- k.res[sapply(k.res, function(x) {
    #     return(!is.character(x[[1]]))
    # })]
#     if(length(k.res)<2) return(NULL)
#     pairs.of.clusts <- combn(k.res, 2)
# 
#     return(sapply(seq_len(ncol(pairs.of.clusts)), function(col){
# 
#         betweenClustDist <- getBetweenClusterDistances(pairs.of.clusts[, col])
#         varianceBetweenClusts <- 1/(length(betweenClustDist))*sum(betweenClustDist^2)
#         # withinClustDists <- unlist(lapply(pairs.of.clusts[, col], getWithinClusterBootDistances))
#         denom <- sqrt(sum(sapply(lapply(pairs.of.clusts[, col], getWithinClusterBootDistances), mean)^2))
#     
#         betweenClustDist / denom
#    }))
# 
# }


returnBootstrapFStats <- function(k.res){
  k.res <- k.res[sapply(k.res, function(x) {
    return(!is.character(x[[1]]))
  })]
  if(length(k.res)<2) return(NULL)
  
  betweenClustVar <- getNumFstat(k.res)
  withinClustVar <- sum(sapply(k.res, getWithinClusterSumSquareDev))
  
  totalSamples <- sum(sapply(k.res, \(x) return(length(x$boot))))

  f.stat <- betweenClustVar*(totalSamples-length(k.res))/withinClustVar
  
  return(c(f.stat=f.stat, pval = pf(f.stat, length(k.res)-1, totalSamples-length(k.res), lower.tail=FALSE)))
}


