library(data.table)


res <- list()

for(sample in c("MB243-Nuclei_clones", "ST1R-PDX_clones", "STP-Nuclei_clones", "STP-PDX_clones")){
  
  fl <- (list.files(file.path("infercnv/infercnv_MB/scRNA_scDNA/", sample), pattern="3rd_added.txt", full.names = T))
  fl_unfiltered <- (list.files(file.path("infercnv/infercnv_MB/scRNA_scDNA/", sample), pattern="_scDNA_clones_filtered_cells.txt", full.names = T))

  
  
  tbl_filtered <- fread(fl)
  tbl_filtered$Sample <- sample
  
  tbl_unfiltered <- fread(fl_unfiltered)
  tbl_unfiltered$Sample <- sample
  
  print(sample)
  print("Unfiltered cells")
  tbl_unfiltered[,table(clone_id)]
  print("Aligned cells per clone")
  print(tbl_unfiltered[,table(padj<0.05, clone_id)])
  print("Confident cells per clone")
  print(tbl_filtered[distance > 0.025, table(clone_id)])
  
  print("Final cells per clone")
  print(tbl_filtered[keep.cell == TRUE, table(merged_clone_id)])
  
}

res <- rbindlist(res)
