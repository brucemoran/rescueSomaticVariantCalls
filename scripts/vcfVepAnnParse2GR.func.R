#! /usr/bin/R

##read VCFs into GrangesList, take overlap based on 'consensus' files (list)
##use for plotting
libs <- c("ensemblVEP", "org.Hs.eg.db", "customProDB", "tidyverse")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})
strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f,sepn)[[1]]})
}
strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

vcfVepAnnParseGR <- function(vcfIn, germline){

  vcf <- readVcf(vcfIn)
  if(dim(vcf)[1] != 0){
    gr <- suppressWarnings(InputVcf(vcfIn))

    ##parse info
    infor <- info(header(vcf))

    ##VEP annotation naming
    annNames <- unlist(strSplitFun(infor[rownames(infor)=="ANN",]$Description,"\\|"))

    ##somatic
    somName <- names(gr)[names(gr)!=germline]
    print(paste0("Working on: ",somName))
    som <- gr[[somName]]
    seqinfo(som) <- seqinfo(vcf)[seqlevels(som)]

    ##annotation by CANONICAL, and add to mcols
    somAnnDf <- t(as.data.frame(lapply(strSplitFun(som$ANN,"\\|"),function(ff){
      if(ff[annNames=="CANONICAL"]=="YES"){
        if(is.null(ff)){ff<-rep("",length(annNames))}
        if(length(ff)!=length(annNames)){
          lengExtra <- length(annNames)-length(ff)
          ff<-c(ff,rep("",lengExtra))}
        return(ff)}
        else{
          return(rep("",length(annNames)))
        }
      })))
    colnames(somAnnDf) <- annNames

    if(sum(dim(somAnnDf)) != 0){
      values(som) <- cbind(as.data.frame(mcols(som)),somAnnDf)
      som$ANN <- NULL
    }
    som <-unique(som)

    return(som)
  }
  else{
    print("No variants found")
    return(GRanges())
  }
}
