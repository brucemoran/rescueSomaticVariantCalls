#! /usr/bin/R

##functions for SNV_consensus generation and plotting

##load libraries
libs <- c("ensemblVEP", "org.Hs.eg.db", "customProDB", "GenomicRanges", "tidyverse", "bio3d", "plyr", "pheatmap", "data.table")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})

strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f,sepn)[[1]]})
}
strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

##parse VCF into GR, for use on un-annotated VCFs
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

#parse VEP annotated VCF into GR, uses CANONICAL transcript
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

##create single-letter HGVS protein annotation (VEP outputs 3-letter)
##take vector, gsub out aa3 for aa1
subHGVSp <- function(inVec){
  lib <- c("bio3d")
  loadedLib <- lapply(lib,function(l){suppressMessages(library(l, character.only = TRUE))})

  aa1 <- bio3d::aa.table$aa1
  ##amino acid 3 letter to gsub HGVSp
  aa3 <- unlist(lapply(bio3d::aa.table$aa3,function(f){
    sp <- strsplit(f,"")[[1]];
    paste0(sp[1], tolower(sp[2]),tolower(sp[3]))
  }))

  ##include * for Ter
  aa1 <-c(aa1,"*")
  aa3 <- c(aa3, "Ter")

  unlist(lapply(inVec,function(f){
    #check matches (should be none or two)
    a3 <- aa3[!is.na(unlist(stringi::stri_match_all(f,regex=aa3)))]
    a1 <- aa1[!is.na(unlist(stringi::stri_match_all(f,regex=aa3)))]
    ##beauty:
    #https://stackoverflow.com/questions/19424709/r-gsub-pattern-vector-and-replacement-vector
    if(length(a3)>0){
      names(a1) <- a3
      str_replace_all(f,a1)
    }
    else{
      return("")
    }
  }))
}

##create GRsupersets from list of variants
#varList = nested list of [[caller]][[samples1..n]]
##NB hardcode via =NULL
GRsuperSet <- function(varList, impacts=NULL, mcolsWant=NULL, nameCallers=NULL){

  ##set NULL vars
  if(is.null(impacts)){impacts <- c("HIGH","MODERATE")}
  if(is.null(mcolsWant)){mcolsWant <- c("AD", "AD.1", "AF", "Consequence", "IMPACT", "SYMBOL", "HGVSc", "HGVSp", "HGVSp1", "CLIN_SIG")}
  if(is.null(nameCallers)){nameCallers <- c("mutect2", "lancet")}

  if(length(nameCallers) != 2){
    print("Require only 2 callers, no more and no less!")
    exit;
  }

  ##set up output
  GRsuper <- as.list(names(varList[[1]]))
  callerNames <- names(varList)

  ##read first entry
  call1 <- varList[[nameCallers[1]]]
  call2 <- varList[[nameCallers[2]]]

  #exclude MT, GL
  seqwant <- c(seq(from=1,to=22,by=1), "X")

  ##iterate over samples in callerset
  for (x in seq_along(call1)){
    print(paste0("Working on: ",names(call1)[x]))

    calls1 <- call1[[x]]
    calls2 <- call2[[x]]

    calls1$HGVSp1 <- subHGVSp(calls1$HGVSp)
    calls2$HGVSp1 <- subHGVSp(calls2$HGVSp)

    ##sets of call1, 2 and the difference
    gr11 <- calls1[calls1$IMPACT %in% impacts, names(mcols(calls1)) %in% mcolsWant]
    gr22 <- calls2[calls2$IMPACT %in% impacts, names(mcols(calls2)) %in% mcolsWant]
    gr12 <- suppressWarnings(setdiff(gr22, gr11))

    if(length(gr11)!=0 & length(gr12)!=0){
      mcols(gr11) <- mcols(gr11)[,mcolsWant]
      mcols(gr12) <- mcols(gr12)[,mcolsWant]
    }

    ##add 1 and difference (the superset of variants)
    GRsuper[[x]] <- suppressWarnings(c(gr11,gr12))
  }
  return(GRsuper)
}

##find consensus in at least two callers, therefore in GRsuper
##this produces a set of positions to plot
atLeastTwo <- function(varList, GRsuper, tag, tmb=NULL){

  ##run TMB?
  if(is.null(tmb)){tmb <- "snv"}

  ##set vars
  callers <- names(varList)
  samples <- names(varList[[1]])

  ##iterate over list of callers
  GRplots <- lapply(seq_along(samples), function(x){

      sample <- samples[x]
      print(paste0("Working on: ",sample))

      ##all possible combinations of intersects of callers
      ##output to clean GR
      upl <- GRanges()
      upl1 <- apply(t(combn(length(callers), m=2)), 1, function(xx){
        print(paste(callers[xx[1]]," vs. ",callers[xx[2]]))
        gr1 <- varList[[names(varList)[xx[1]]]][[x]]
        gr2 <- varList[[names(varList)[xx[2]]]][[x]]
        gri <- suppressWarnings(GenomicRanges::intersect(gr1,gr2))
        upl <- suppressWarnings(c(upl, gri))
      })
      upl2 <- upl1[[1]]
      for(xx in 2:length(upl1)){
        upl2 <- suppressWarnings(c(upl2, upl1[[xx]]))
      }
      GRplot <- suppressWarnings(GRsuper[[x]][GRsuper[[x]] %in% unique(upl2)])

      #TMB
      fileOut <- paste0(sample,".",tag,".consensus.tab")
      if(tmb == "snv"){
          tmbOut <- exomeTumourMutationBurden(GRplot)
          fileOut <- paste0(sample,".",tag,".TMB_",tmbOut,"_SNV-Mb.consensus.tab")
        }
      write.table(GRplot,file=fileOut,quote=F,row=F,col=T,sep="\t")
      return(GRplot)
    })

    names(GRplots) <- samples
    return(GRplots)
}

##create two plots: all consensus, those in 2+ samples
plotConsensusList <- function(plotList, rawList, tag, includeOrder=NULL){

  ##remove hyphens
  samples1 <- gsub("-","_",names(plotList))
  if(is.null(includeOrder)){
    includeOrder <- samples1
  }
  # ##plots are per mutype
  # mutype <- names(plotList)[x]
  # print(paste0("Working on: ", mutype))

  ##combined set of all samples
  combSet <- suppressWarnings(unique(do.call("c", unname(plotList))))
  seqlevels(combSet) <- sort(seqlevels(combSet))
  combSet <- sort(combSet)
  combDF <- as.data.frame(combSet)

  ##labels for plot
  hgvsp <- unlist(lapply(combSet$HGVSp1,function(f){strsplit(f,"\\.")[[1]][3]}))
  uniqLabels <- paste0(names(combSet),":",combDF$SYMBOL,":",combDF$Consequence, ":", hgvsp)

  ##take those positions, then query raw calls
  ##allows 'rescue' of those falling out from arbitrary filters
  ##enough support previously to allow re-entry AFAICS
  plotDFrawout <- as.data.frame(lapply(rawList,function(ff){
          afs <- rep(0,length(combSet))
          lapply(ff,function(fff){
          seqlevels(fff) <- sort(seqlevels(fff))
          ffs <- sort(fff)
          ffsi <- ffs[ffs %in% combSet]
          afs[combSet %in% ffsi] <- as.numeric(mcols(ffsi)$AF)
          return(afs)
        })}))

  plotDFrawout <- do.call(cbind,lapply(samples1, function(ss){
    apply(plotDFrawout[,grep(ss, colnames(plotDFrawout))],1,max)
  }))

  plotDFrawout <- as.data.frame(plotDFrawout)
  colnames(plotDFrawout) <- samples1
  rownames(plotDFrawout) <- uniqLabels

  ##which samples to include, and order
  ##remove those with all 0 frequency
  plotDForder <- plotDFrawout[rowSums(plotDFrawout)!=0, colnames(plotDFrawout) %in% includeOrder]

  ##reduce any frequency >50% to 50% (somatic should not be >50%)
  ##and plot is purely representative
  plotDForder[plotDForder > 0.5] <- 0.5

  ##find all 0, count to allow separation, order
  notoDF <- do.call(cbind,list(apply(plotDForder,2,function(f){f!=0})))
  notoVec <- apply(notoDF,1,function(f){table(f)[[1]]})
  notoVec <- notoVec[!is.na(notoVec)]
  notoVec1 <- notoVec[notoVec==1]
  notoVec2 <- notoVec[notoVec>1]

  plotDF <- plotDForder[!rownames(plotDForder) %in% names(notoVec1), includeOrder1]
  plotDF2 <- plotDForder[!rownames(plotDForder) %in% names(notoVec2), includeOrder1]

  ##ordering
  plotDF$labels <- rownames(plotDF)
  orderDF <- dplyr::arrange_(plotDF,.dots=includeOrder1)
  rownames(orderDF) <- orderDF$labels
  plotDFordered <- orderDF[,-c(grep("labels",colnames(orderDF)))]
  orderDF2 <- do.call(order, c(data.frame(plotDF2[, 1:(ncol(plotDF2)-1)], plotDF2[, ncol(plotDF2)])))

  ##exit if no variants
  if(is.null(plotDFordered)){
    print("No consensus variants found, exiting")
    break
  }

  else{


    ##no shared variants
    if(!is.null(orderDF2)){
      plotDF2ordered <- plotDF2[orderDF2,]
      plotDF2ordered <- rbind(plotDF2ordered, plotDFordered)
    }
    plotVec <- c()
    plotTag <- c()
    if(!is.null(orderDF2)){
      plotVec <- list(plotDFordered, plotDF2ordered)
      plotTag <- list("only-shared", "private-shared")
      print("Plotting private and shared variants")
    }
    if(is.null(orderDF2)){
      plotVec <- list(plotDFordered)
      plotTag <- list("only-private")
      print("No shared variants, plotting private only")
    }
    if(dim(plotDFordered)[1] == 0){
      print("No variants to plot")
      break
    }
    for (pl in 1:length(plotVec)){
      plotLabels <- rep("",dim(plotVec[[pl]])[1])
      row_fontsize <- 1
      colz <- colorRampPalette(c("lightgrey","dodgerblue","blue"))
      if(dim(plotVec[[pl]])[1] < 120){
        if(dim(plotVec[[pl]])[1]<20){row_fontsize=8}
        if(dim(plotVec[[pl]])[1]<50){row_fontsize=6}
        if(dim(plotVec[[pl]])[1]>50 & dim(plotVec[[pl]])[1]<100){row_fontsize=4}
        if(dim(plotVec[[pl]])[1]>100){row_fontsize=2}
        plotLabels <- rownames(plotVec[[pl]])
      }

      pdf(paste0(tag,".",plotTag[[pl]],".consensus.onlyOverlap.pdf"),onefile=F)
      pheatmap(plotVec[[pl]][,c(1:length(includeOrder1))],
         breaks=seq(from=0,to=0.5,length.out=101),
         color=colz(100),
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         clustering_distance_rows=NA,
         cellwidth=12,
         legend=TRUE,
         fontsize_row=row_fontsize,
         labels_row=plotLabels,
         border_color="lightgrey",
         gaps_col=c(1:length(includeOrder1))
      )
      dev.off()
    }
  }
}

plotConsensusSingle <- function(plotList, rawList, tag, includeOrder=NULL){

  ##remove hyphens
  samples1 <- gsub("-","_",names(plotList))

  ##combined set of all samples
  combSet <- suppressWarnings(unique(do.call("c", unname(plotList))))
  seqlevels(combSet) <- sort(seqlevels(combSet))
  combSet <- sort(combSet)
  combDF <- as.data.frame(combSet)

  ##labels for plot
  hgvsp <- unlist(lapply(combSet$HGVSp1,function(f){strsplit(f,"\\.")[[1]][3]}))
  uniqLabels <- paste0(names(combSet),":",combDF$SYMBOL,":",combDF$Consequence, ":", hgvsp)

  ##take those positions, then query raw calls
  ##allows 'rescue' of those falling out from arbitrary filters
  ##enough support previously to allow re-entry AFAICS
  plotDFrawout <- as.data.frame(lapply(rawList,function(ff){
          afs <- rep(0,length(combSet))
          lapply(ff,function(fff){
          seqlevels(fff) <- sort(seqlevels(fff))
          ffs <- sort(fff)
          ffsi <- ffs[ffs %in% combSet]
          afs[combSet %in% ffsi] <- as.numeric(mcols(ffsi)$AF)
          return(afs)
        })}))

  plotDFrawout <- do.call(cbind,lapply(samples1, function(ss){
    apply(plotDFrawout[,grep(ss, colnames(plotDFrawout))],1,max)
  }))

  plotDFrawout <- as.data.frame(plotDFrawout)
  colnames(plotDFrawout) <- samples1
  rownames(plotDFrawout) <- uniqLabels

  ##which samples to include, and order
  ##remove those with all 0 frequency
  plotDForder <- plotDFrawout[rowSums(plotDFrawout)!=0, colnames(plotDFrawout) %in% includeOrder]

  ##reduce any frequency >50% to 50% (somatic should not be >50%)
  ##and plot is purely representative
  plotDForder[plotDForder > 0.5] <- 0.5

  ##find all 0, count to allow separation, order
  orderPlotDF <- order(plotDForder)
  plotDF <- plotDForder[orderPlotDF]
  orderUniqLabels <- uniqLabels[orderPlotDF]

  ##ordering
  plotDFordered <- data.frame(plotDF, row.names=orderUniqLabels)

  ##exit if no variants
  if(is.null(plotDFordered)){
    print("No consensus variants found, exiting")
    break
  }

  if(!is.null(plotDFordered)){
    plotVec <- plotDFordered
    plotTag <- "variants"

    plotLabels <- rep("",times=dim(plotVec)[1])
    row_fontsize <- 1
    colz <- colorRampPalette(c("lightgrey","dodgerblue","blue"))
    if(dim(plotVec)[1] < 120){
      if(dim(plotVec)[1]<20){row_fontsize=8}
      if(dim(plotVec)[1]<50){row_fontsize=6}
      if(dim(plotVec)[1]>50 & dim(plotVec)[1]<100){row_fontsize=4}
      if(dim(plotVec)[1]>100){row_fontsize=2}
      plotLabels <- rownames(plotVec)
    }

    pdf(paste0(tag,".",plotTag,".consensus.onlyOverlap.pdf"),onefile=F)
    pheatmap(plotVec,
       breaks=seq(from=0,to=0.5,length.out=101),
       color=colz(100),
       cluster_rows=FALSE,
       cluster_cols=FALSE,
       clustering_distance_rows=NA,
       cellwidth=12,
       legend=TRUE,
       fontsize_row=row_fontsize,
       labels_row=plotLabels,
       labels_col=tag,
       border_color="lightgrey"
    )
    dev.off()
  }
}

exomeTumourMutationBurden <- function(GRplot){

  ##get exome for Illumina Nextera Rapid
  exomeBed <- fread("https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed", showProgress=FALSE, data.table=FALSE)

  ##triage
  exomeBed[,1] <- gsub("chr","",exomeBed[,1])
  colnames(exomeBed) <- c("seqname", "start", "end")
  exomeGR <- makeGRangesFromDataFrame(exomeBed, ignore.strand=TRUE)
  exomeSize <- sum(width(exomeGR))/1000000

  ##overlap with input
  hits <- as.data.frame(findOverlaps(GRplot, exomeGR))
  ##output
  GRplotExome <- GRplot[hits$queryHits]
  TMB <- round(length(GRplotExome)/exomeSize,digits=1)
  return(TMB)
}
