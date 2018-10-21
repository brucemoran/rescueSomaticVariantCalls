#! /usr/bin/R

##read VCFs into Granges

libs <- c("ensemblVEP", "customProDB", "tidyverse")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})
strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f,sepn)[[1]]})
}
strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

vcfParseGR <- function(vcfIn,germline){

  vcf <- readVcf(vcfIn)
  gr <- suppressWarnings(InputVcf(vcfIn))

  ##parse info
  infor <- info(header(vcf))

  ##somatic
  somName <- names(gr)[names(gr)!=germline]
  print(paste0("Working on: ",somName))
  som <- gr[[somName]]
  seqinfo(som) <- seqinfo(vcf)[seqlevels(som)]

  return(som)
}
