require(scales)
require(GenomicRanges)

plotComplexHeatmapCommonSegsOld <- function(commonSegmentation, clust.rows=TRUE, clust.dist = "manhattan", clust.method = "ward.D2", annot_rows = NULL, row_split=NULL, column_title = NULL){
  
  chrNames <- commonSegmentation$chr
  
  toPlot <- t(data.matrix(commonSegmentation[,-c(1:3), with=FALSE]))
  colors = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(16)
  
  brks <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4,5,6,7,8,9)
  
  ht = Heatmap(toPlot, cluster_columns  = FALSE, cluster_rows = clust.rows ,col=circlize::colorRamp2(colors, breaks = brks),
               clustering_distance_rows =  clust.dist, clustering_method_rows = clust.method, 
               show_row_names=FALSE, show_column_names = FALSE, column_split = factor(chrNames, levels=unique(chrNames)), 
               column_gap=unit(c(rep(1,times=12), rep(2,times=2), rep(2,times=3), rep(3,4), rep(2, times=2)),"mm"), 
               column_title_gp=gpar(fontsize=10), right_annotation = annot_rows, 
               row_split = row_split, column_title = column_title )
  
  ComplexHeatmap::draw(ht, background="grey70")
  
}

plotComplexHeatmapCommonSegs <- function(commonSegmentation, clust.rows=TRUE, clust.dist = "manhattan", clust.method = "ward.D2", annot_rows = NULL, row_split=NULL, column_title = NULL, ...){
  
  # chrNames <- commonSegmentation$chr
  
  chrLengths <- commonSegmentation[,max(end),chr][,V1]
  
  names(chrLengths) <- commonSegmentation[,max(end),chr][,chr]
  
  plotCoords <- tileGenome(chrLengths, tilewidth=1e5, cut.last.tile.in.chrom = TRUE)
  commonGRanges <- makeGRangesFromDataFrame(commonSegmentation, keep.extra.columns = TRUE)
  
  mappedSegments <- as.data.frame(mergeByOverlaps(plotCoords,commonGRanges)[,-c(1,2)])
  chrNames <- as.character(seqnames(mergeByOverlaps(plotCoords,commonGRanges)[,1]))
  
  toPlot <- t(data.matrix(mappedSegments))
  colors = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(16)
  
  brks <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4,5,6,7,8,9)
  
  ht = Heatmap(toPlot, cluster_columns  = FALSE, cluster_rows = clust.rows ,col=circlize::colorRamp2(colors, breaks = brks),
               clustering_distance_rows =  clust.dist, clustering_method_rows = clust.method, 
               show_row_names=FALSE, show_column_names = FALSE, column_split = factor(chrNames, levels=unique(chrNames)), 
               column_gap=unit(c(rep(1,times=12), rep(2,times=2), rep(2,times=3), rep(3,4), rep(2, times=2)),"mm"), 
               column_title_gp=gpar(fontsize=10), right_annotation = annot_rows, 
               row_split = row_split, column_title = column_title,... )
  
  ComplexHeatmap::draw(ht, background="grey70")
  return(invisible(ht))
}


plotComplexHeatmapCommonSegs3 <- function(commonSegmentation, clust.rows=TRUE, clust.dist = "manhattan", clust.method = "ward.D2", annot_rows = NULL, row_split=NULL, column_title = NULL){
  
  # chrNames <- commonSegmentation$chr
  
  chrStarts <- cumsum(c(0,commonSegmentation[,max(end),chr][,V1][-1]))
  
  names(chrStarts) <- commonSegmentation[,max(end),chr][,chr]
  
  
  toPlot <- data.frame(start.x = chrStarts[commonSegmentation$chr] + commonSegmentation$start,
    end.x = chrStarts[commonSegmentation$chr] + commonSegmentation$end)
  
  toPlot <- t(data.matrix(commonSegmentation[,-c(1:3), with=FALSE]))
  colors = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(16)
  
  brks <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4,5,6,7,8,9)
  
  ht = Heatmap(toPlot, cluster_columns  = FALSE, cluster_rows = clust.rows ,col=circlize::colorRamp2(colors, breaks = brks),
               clustering_distance_rows =  clust.dist, clustering_method_rows = clust.method, 
               show_row_names=FALSE, show_column_names = FALSE, column_split = factor(chrNames, levels=unique(chrNames)), 
               #column_gap=unit(c(rep(1,times=12), rep(2,times=2), rep(2,times=3), rep(3,4), rep(2, times=2)),"mm"), 
               column_title_gp=gpar(fontsize=10), right_annotation = annot_rows, 
               row_split = row_split, column_title = column_title )
  
  ComplexHeatmap::draw(ht, background="grey70")
  
}



plotComplexHeatmapChromosomes <- function(QDNAobject, distance = "manhattan", clustering_method_rows = "complete", 
                                          cluster_rows = TRUE, annot_rows = NULL, 
                                          background = "grey70"){

  toPlot <- t(QDNAobject@assayData$copynumber[complete.cases(QDNAobject@assayData$copynumber),])

  chrNames <- gsub(colnames(toPlot), pattern = "\\:[0-9]+\\-[0-9]+$", rep="")

  colors = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(16)
  
  brks <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4,5,6,7,8,9)
  
  ht = Heatmap(toPlot, cluster_columns  = FALSE, 
               col=circlize::colorRamp2(colors, breaks = brks), 
               clustering_distance_rows =  distance, 
               clustering_method_rows = clustering_method_rows, 
               show_row_names=FALSE, show_column_names = FALSE, 
               right_annotation = annot_rows, 
               column_split = factor(chrNames, levels=unique(chrNames)), 
               column_gap=unit(c(rep(1,times=12), rep(2,times=2), rep(2,times=3), rep(3,4), rep(2, times=3)),"mm"), 
               column_title_gp=gpar(fontsize=10))

  draw(ht, background=background)

}


plotComplexHeatmapChromosomesNoColors <- function(QDNAobject, distance = "manhattan", clustering_method_rows = "complete"){
  
  toPlot <- t(QDNAobject@assayData$copynumber[complete.cases(QDNAobject@assayData$copynumber),])
  
  chrNames <- gsub(colnames(toPlot), pattern = "\\:[0-9]+\\-[0-9]+$", rep="")
  
  # colors = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(16)
  
  # brks <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4,5,6,7,8,9)
  
  ht = Heatmap(toPlot, cluster_columns  = FALSE, clustering_distance_rows =  distance, clustering_method_rows = clustering_method_rows, show_row_names=FALSE, show_column_names = FALSE, column_split = factor(chrNames, levels=unique(chrNames)), column_gap=unit(c(rep(1,times=12), rep(2,times=2), rep(2,times=3), rep(3,4), rep(2, times=3)),"mm"), column_title_gp=gpar(fontsize=10))
  
  draw(ht, background="grey70")
  
}

plotComplexHeatmapPloidyNorm <- function(QDNAobject){

  toPlot <- t(QDNAobject@assayData$copynumber[complete.cases(QDNAobject@assayData$copynumber),]/QDNAobject$ploidy.mod*2)

  chrNames <- gsub(colnames(toPlot), pattern = "\\:[0-9]+\\-[0-9]+$", rep="")

  colors = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(8)
  
  brks <- c(0, 0.5, 1, 1.5, 2, 3, 4,5)
  
  ht = Heatmap(toPlot, cluster_columns  = FALSE, col=circlize::colorRamp2(colors, breaks = brks), clustering_distance_rows =  "manhattan", clustering_method_rows = "complete", show_row_names=FALSE, show_column_names = FALSE, column_split = factor(chrNames, levels=unique(chrNames)), column_gap=unit(c(rep(1,times=12), rep(2,times=2), rep(2,times=3), rep(3,4), rep(2, times=3)),"mm"), column_title_gp=gpar(fontsize=10))

  draw(ht, background="grey70")

}




plotSegments2 <- function (correctOutput, segmentOutput, chromosome=1 , plotMedians=FALSE, range, ... ) {
  if (is.null(segmentOutput$segs)) {
    warning("Processed segments not found, automatically processing")
    segmentOutput$segs <- HMMcopy:::processSegments(segments$segs, 
                                          correctOutput$chr, correctOutput$start, correctOutput$end, 
                                          correctOutput$copy)
  }
  segs <- segmentOutput$segs
  correctOutput$state <- segmentOutput$state
  cols <- rev(brewer_pal('div', palette = 7)(length(unique(segs$state))))
  if(missing(range)){
    range <- quantile(correctOutput$copy, na.rm = TRUE, prob = c(0.01, 0.99))
  }
  a <- subset(correctOutput, chr == chromosome)
  b <- subset(segs, chr == chromosome)
  
  if(is.matrix(segmentOutput$mus)){
    est.mus <- segmentOutput$mus[,ncol(segmentOutput$mus)]
  } else {
    est.mus <- segmentOutput$mus
  }

  
  plot(a$start, a$copy, col = cols[as.numeric(as.character(a$state))], 
       ylim = range, ...)
  
  if(plotMedians){
    for (k in 1:nrow(b)) {
      lines(c(b$start[k], b$end[k]), rep(b$median[k], 2), lwd = 3, 
          col = "green")
    }
  } else {
    for (k in 1:nrow(b)) {
      lines(c(b$start[k], b$end[k]), rep(est.mus[as.numeric(b$state[k])], 2), lwd = 3, 
          col = "green")
    }
  }
}


plotSegmentsAllChr <- function (correctOutput, segmentOutput, plotMedians=FALSE,range, ... )
{
  if (is.null(segmentOutput$segs)) {
    warning("Processed segments now found, automatically processing")
    segmentOutput$segs <- processSegments(segments$segs, 
                                          correctOutput$chr, correctOutput$start, correctOutput$end, 
                                          correctOutput$copy)
  }
  segs <- segmentOutput$segs
  correctOutput$state <- segmentOutput$state

  cols <- rev(brewer_pal('div', palette = 7)(length(unique(segs$state))))
  if(missing(range)){
    range <- quantile(correctOutput$copy, na.rm = TRUE, prob = c(0.01, 0.999))
    
  }
  
  correctOutput$chr <- factor(correctOutput$chr , mixedsort(unique(correctOutput$chr)))

  chromosome_lengths <- aggregate(data=correctOutput, end~chr , max)
  
  chromosome_offsets <- c(0,cumsum(as.numeric(chromosome_lengths$end))[-nrow(chromosome_lengths)])
  names(chromosome_offsets) <- mixedsort(unique(correctOutput$chr))
  
  correctOutput$len <- correctOutput$end - correctOutput$start -1
  correctOutput$start.x <- correctOutput$start 
  correctOutput$end.x <- correctOutput$end -1
  correctOutput[,end.x := end.x + chromosome_offsets[chr]]
  correctOutput[,start.x := start.x + chromosome_offsets[chr]]
  
  segs$start.x <- segs$start + chromosome_offsets[segs$chr]
  segs$end.x <- segs$end + chromosome_offsets[segs$chr]
  
  
  a <- correctOutput
  b <- segs
  
  est.mus <- segmentOutput$mus[,ncol(segmentOutput$mus)]
  
  plot(a$start.x, a$copy, col = cols[as.numeric(as.character(a$state))], 
       ylim = range, ...)
  
  if(plotMedians){
    for (k in 1:nrow(b)) {
      lines(c(b$start.x[k], b$end.x[k]), rep(b$median[k], 2), lwd = 3, 
          col = "green")
    }
  } else {
    for (k in 1:nrow(b)) {
      lines(c(b$start.x[k], b$end.x[k]), rep(est.mus[as.numeric(b$state[k])], 2), lwd = 3, 
          col = "green")
    }
  }

  for (k in seq_along(chromosome_offsets[-1])) {
      abline(v=chromosome_offsets[-1][k])
    }
  
}


plotCellHeatmap <- function(all.segments, cols=c("dodgerblue4", "dodgerblue2", "white", "salmon", "firebrick1", "firebrick"), chr.subset=NA){
  
  all.segments$chr <- factor(all.segments$chr , mixedsort(unique(all.segments$chr)))
  
  if(is.na(chr.subset)){
    # chr.subset <- paste0('chr',unique(all.segments$chr))
    chr.subset <- unique(all.segments$chr)
  }
  
  chromosome_lengths <- aggregate(data=all.segments, end~chr , max)
  
  chromosome_offsets <- c(0,cumsum(as.numeric(chromosome_lengths$end))[-nrow(chromosome_lengths)])
  names(chromosome_offsets) <- mixedsort(unique(all.segments$chr))
  
  all.segments$len <- all.segments$end - all.segments$start -1
  all.segments$start.x <- all.segments$start 
  all.segments$end.x <- all.segments$end -1
  all.segments[,end.x := end.x + chromosome_offsets[chr]]
  all.segments[,start.x := start.x + chromosome_offsets[chr]]
  
  cell_ids <- as.character(unique(all.segments$cell))
  
  #a1 <- all.segments[all.segments$cell=='A1',]
  #plot.all.chr.median.per.cell(a1)
  all.segments$cell <- factor(all.segments$cell , mixedsort(unique(all.segments$cell)))
  all.segments <- all.segments[order(all.segments$cell),]
  cell_pos <- data.frame(cell=levels(all.segments$cell), id=1:length(unique(all.segments$cell)) )
  all.segments <- merge(all.segments, cell_pos)

  chr.starts <- aggregate(data=all.segments, start.x~chr , min)
  chr.end <- aggregate(data=all.segments, end.x~chr , max)
  chr.coord <- merge(chr.starts, chr.end)
  colnames(chr.coord) <- c('chromosome', 'start.chr', 'end.chr')
  chr.coord$chromosome <- chr.coord$chromosome
  chr.coord$chromosome <- factor( chr.coord$chromosome , levels=mixedsort( chr.coord$chromosome ))
  chr.coord <- chr.coord[order(chr.coord$chromosome),]
  chr.coord$mid <- (chr.coord$start.chr + chr.coord$end.chr)/2
  
  ## Make discrete TCN scale
  
  all.segments[,state := factor(c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP")[TCN], levels=c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP"))]

  p <- ggplot(all.segments[chr%in%chr.subset] ) + geom_rect(aes(ymin=id-0.5, ymax=id+0.5, xmin=start.x, xmax=end.x, fill=state), size=1) +
    scale_fill_manual(limits = c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP"), values = cols) + 
    geom_vline(aes(xintercept = start.chr), data=chr.coord, color='grey30') + theme_bw() +
    geom_text(data=chr.coord[chr.coord$chr%in%chr.subset,], aes(label= chromosome, y = length(cell_ids) +min(5, length(cell_ids)*0.6), x=mid ), size=5, angle = 90)+ ylab("")+
    scale_x_continuous(limits = c(min(chr.coord[chr.coord$chromosome%in%chr.subset, 'start.chr']),max(chr.coord[chr.coord$chromosome%in%chr.subset, 'end.chr'])), expand = c(0, 0)) 
  return(p)
  
}

