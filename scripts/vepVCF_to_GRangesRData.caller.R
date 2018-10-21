#! R

options(stringAsFactors=FALSE)
args <- commandArgs(trailingOnly = TRUE)
LIBPATH <- args[1]
SCRIPTDIR <- args[2]
WORKDIR <- args[3]
GERMLINE <- args[4]
##specific format: caller.mutype.any.thing.you.like.vcf
VEPVCFPATTERN <- args[5]
RAWVCFPATTERN <- args[6]
if(length(args)<6){ RAWVCFPATTERN <- "raw.vcf$" }

.libPaths(LIBPATH)

##these are made into a GRangesList object and saved according to args[5][-1] split on "/"
source(paste0(SCRIPTDIR,"/R/functions/GRanges/vcfVepAnnParse2GR.func.R"))
source(paste0(SCRIPTDIR,"/R/functions/GRanges/vcfParse2GR.func.R"))

setwd(WORKDIR)
##all VCF to be in WORKDIR
vcfList <- dir(pattern=VEPVCFPATTERN)
##ensure all end in .vcf
vcfNames <- unlist(lapply(vcfList, function(f){strSplitFun(f, "\\.")}[[1]][1]))
vcfList <- as.list(vcfList[grep(".vcf$",vcfList)])
outCaller <- strsplit(VEPVCFPATTERN, "\\.")[[1]][1]
outExt <- strsplit(VEPVCFPATTERN,"\\.vcf")[[1]][1]
outName <- gsub(".pass","",VEPVCFPATTERN)
outName <- strsplit(outName, "\\.")[[1]][2]

##run function to make GRangesList
grl <- lapply(vcfList,function(f){
  suppressWarnings(vcfVepAnnParseGR(f, GERMLINE))
  })
names(grl) <- vcfNames

#assign output
assignedName <- paste0(outCaller,".",outName)
assign(assignedName, value=grl)

#save to workdir
saveFile <- paste0(outExt,".RData")
save(list=assignedName, file=saveFile)

##also take the *raw.vcf
rawList <- dir(pattern=paste0(outCaller,".",RAWVCFPATTERN))
rawNames <- unlist(lapply(rawList, function(f){strSplitFun(f, "\\.")}[[1]][1]))
rawExt <- c("raw")

dirtest <- dir(pattern=paste0(outCaller,".raw.RData"))
if(length(dirtest)==0){
  print("Processing raw VCFs")
  grr <- lapply(rawList,function(f){
    vcfParseGR(f, GERMLINE)
    })
  names(grr) <- rawNames

  ##assign output
  assignedNamer <- paste0(outCaller,".",rawExt)
  assign(assignedNamer, value=grr)

  ##save to current dir
  saveFiler <- paste0(outCaller,".",rawExt,".RData")
  save(list=assignedNamer, file=saveFiler)
}
