require(gtools)
require(Biobase)






## This function merges breakpoints by looking at chains where bin starts are within mergeWidth*binSize apart, and then picking the median breakpoint 
## Should think about adding a minimum number of cells supporting breakpoint filter here
## should think about checking that I am not merging breakpoints in different directions? Or maybe not - could be mirror patterns??
getConsensusBreakpoints <- function(cn.res, mergeWidth = 3){
  
  binSize <- as.numeric(names(sort(table(sapply(strsplit(rownames(cn.res), ":|-"), function(x) return(as.numeric(x[[3]])-as.numeric(x[[2]])+1))), decreasing = TRUE))[[1]])
  
  breakpoint_locs <- lapply(apply(cn.res[2:nrow(cn.res),, drop=FALSE] != cn.res[1:(nrow(cn.res)-1),, drop=FALSE], 2, which), function(x) return(rownames(cn.res)[x]))
  # 
  # unique_breakpoints <- mixedsort(unique(unlist(breakpoint_locs)))
  # unique_breakpoints_df <- data.frame(do.call(rbind, strsplit(unique_breakpoints, ":|-")))
  # 
  # unique_breakpoints_df[,2] <- as.numeric(unique_breakpoints_df[,2])
  # unique_breakpoints_df[,3] <- as.numeric(unique_breakpoints_df[,3])
  # 
  # colnames(unique_breakpoints_df) <- c("chr", "start", "end")
  # 
  # 
  
  all_breakpoints <- mixedsort((unlist(breakpoint_locs)))
  all_breakpoints_df <- data.frame(do.call(rbind, strsplit(all_breakpoints, ":|-")))
  
  all_breakpoints_df[,2] <- as.numeric(all_breakpoints_df[,2])
  all_breakpoints_df[,3] <- as.numeric(all_breakpoints_df[,3])
  
  colnames(all_breakpoints_df) <- c("chr", "start", "end")
  
  
  breakpoint_groups <- list()
  i = 1
  while(i <= length(all_breakpoints)){
    
    active_breakpoints <- list(all_breakpoints_df[i,])
    j <- i+1
    while(j <= length(all_breakpoints)){
      if(all_breakpoints_df[j,"chr"] == all_breakpoints_df[j-1,"chr"] && all_breakpoints_df[j,"start"] - all_breakpoints_df[j-1,"start"] <= binSize*mergeWidth){
        active_breakpoints <- c(active_breakpoints, list(all_breakpoints_df[j,]))
        j <- j + 1
      } else {
        break
      }
    }
    active_breakpoints <- rbindlist(active_breakpoints)
    breakpoint_groups <- c(breakpoint_groups, list(active_breakpoints))
    i <- j
  }
  median_locs <- sapply(sapply(breakpoint_groups, `[[`, "end"), median) ## end here because we took the left side of the breakpoint
  
  breakpoint_groups <- lapply(seq_along(breakpoint_groups), function(i) {
    return(list("location" = list(chr=unique(breakpoint_groups[[i]][["chr"]]), breakpoint=median_locs[[i]]), all_breakpoints_mapped = unique(breakpoint_groups[[i]])))
  })
}




roundTowardsTarget <- function(x, target){
  ifelse(x<=target, ceiling(x), floor(x))
}

## NOTE THIS FUNCTION ASSUMES INTEGER COPY STATES!!
getBestCopyState <- function(x, cellPloidy=2){
  # This function will first find the mode of the distributions of copy states passed to it. If that is not unique, it will use the 
  # copy state with the longest uninterrupted run out of the modes. If that is still not unique, it will take the median, and round towards the cell ploidy, 
  # making a conservative assumption of no event. 
  
  x.nona <- as.vector(na.omit(x))
  if(length(x.nona)==0){
    return(NA_real_)
  }
  most.common.values <- which(max(tabulate(x.nona))==tabulate(x.nona))
  if(length(most.common.values)>1){ # if we have a tie, then we take the value with the longest run (using rle)
    rlex <- rle(x.nona)
    most.common.values.2 <- rlex$values[rlex$values%in%most.common.values][max(rlex$lengths[rlex$values%in%most.common.values]) == rlex$lengths[rlex$values%in%most.common.values]]
    if(length(most.common.values.2)>1){
      most.common.values <- roundTowardsTarget(median(most.common.values.2), cellPloidy)
    } else {
      most.common.values <- most.common.values.2
    }
  }
  return(most.common.values)
}




mapCellsToCommonBreakpointsInner <- function(cn.res,cell.ploidies, mergeWidth=3){
  
  coord.matrix <- data.frame(do.call(rbind, strsplit(rownames(cn.res), ":|-")))
  
  coord.matrix[,2] <- as.numeric(coord.matrix[,2])
  coord.matrix[,3] <- as.numeric(coord.matrix[,3])
  
  colnames(coord.matrix) <- c("chr", "start", "end")
  coord.matrix <- data.table(coord.matrix)
  
  chr.starts <- rep(1, length(unique(coord.matrix$chr)))
  chr.ends <- coord.matrix[,max(end),chr][,V1]
  names(chr.ends) <- unique(coord.matrix$chr)
  
  common_breakpoints <- rbindlist(lapply(getConsensusBreakpoints(cn.res, mergeWidth), `[[`, "location"))
  chr.seg.list <- list()
  
  for(chrom in unique(coord.matrix$chr)){
    chr.segs <- common_breakpoints[chr==chrom]
    colnames(chr.segs) <- c("chr", "start")
    chr.segs <- rbind(list(chrom, 0),chr.segs)
    if(chr.segs[,max(start)] != chr.ends[[chrom]]){
      
      chr.segs[,end := c(start[-1],chr.ends[[chrom]])]
      chr.segs[,start := start+1]
      
    } else {
      endpoints <- chr.segs$start
      chr.segs <- chr.segs[-.N,]
      chr.segs[,end := endpoints[-1]]
      chr.segs[,start := start+1]
    }
    chr.seg.list[[chrom]] <- chr.segs
  }
  
  chr.segs <- rbindlist(chr.seg.list)
  
  segementedSections <- sapply(seq_len(nrow(chr.segs)),function(ii){
    cur.seg <- chr.segs[ii]

    myx <- coord.matrix[,which(chr==cur.seg$chr & start >= cur.seg$start & end <= cur.seg$end)]
    sapply(seq_len(ncol(cn.res)), function(i){
      getBestCopyState(cn.res[myx,i], cellPloidy = cell.ploidies[i])
    })
  })
  ## Dealing with degenerate case of 1 sample
  if(is.null(dim(segementedSections))){
    out <- cbind(chr.segs, segementedSections)
    colnames(out)[ncol(out)] <- colnames(cn.res)
    return(out)
  } else {
    segementedSections <- t(segementedSections)
    colnames(segementedSections) <- colnames(cn.res)
    out <- cbind(chr.segs, segementedSections)
    return(out)
  }
  
}


mapCellsToCommonBreakpoints <- function(QDNAobject, mergeWidth=3, assay="copynumber"){
  cn.res <- assayData(QDNAobject)[[assay]]
  cell.ploidies <- QDNAobject$ploidy.mod
  return(mapCellsToCommonBreakpointsInner(cn.res, cell.ploidies, mergeWidth = mergeWidth))
  
}
