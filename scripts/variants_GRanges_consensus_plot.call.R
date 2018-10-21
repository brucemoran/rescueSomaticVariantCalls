#! R

options(stringAsFactors=FALSE)
args <- commandArgs(trailingOnly = TRUE)

##where to source functions from
source(args[1])

##string to tag output files
TAG <- args[2]
INCLUDEDORDER <- args[3]

##inputs
tag <- TAG ##tag
includeOrder <- strSplitVec(INCLUDEDORDER,",")[,1] ##comma delim args[3]

##load set of RData GRanges, not raw
print("Loading VEP annotated GRanges")
dirIn <- dir(pattern=".RData")
grIn <- dirIn[grep("raw", dirIn, invert=T)]
loadedGR <- lapply(grIn, function(x){
  load(x,envir=.GlobalEnv)
  })
names(loadedGR) <- loadedGR
mutypes <- unique(strSplitVec(unlist(loadedGR), "\\.")[2,])

#load raw GRanges data
print("Loading raw VCF GRanges")
grawIn <- dir(pattern="raw.RData")
loadedGRaw <- lapply(grawIn, function(x){
  load(x, envir=.GlobalEnv)
  })
names(loadedGRaw) <- loadedGRaw

##get set of RData GRanges into list
mutypeList <- lapply(mutypes, function(f){
  wantGR <- grep(f, names(loadedGR))
  lapply(loadedGR[wantGR], function(g){get(g)})
  })
names(mutypeList) <- mutypes
samples <- names(mutypeList[[1]][[1]])

##get raw GRnages into list
rawList <- lapply(loadedGRaw, function(f){
  get(f)
  })

##get GRanges superset per mutype
GRsuper <- GRsuperSet(mutypeList, impacts=c("HIGH", "MODERATE"))

##get list to plot from
plotList <- atLeastTwo(mutypeList, GRsuper, paste0(tag, ".HM"))

##plot
plotConsensusList(plotList, rawList, paste0(tag, ".HM"), includeOrder)

##run to get all impacts, print but not plot
##get GRanges superset per mutype
GRsuperALL <- GRsuperSet(mutypeList, impacts=c("HIGH", "MODERATE", "MODIFIER", "LOW"))

##get list to plot from
plotListAll <- atLeastTwo(mutypeList, GRsuperALL, paste0(tag, ".ALL"), tmb="snv", mutypeWant="snv")

##plot
plotConsensusList(plotListAll, rawList, paste0(tag, ".ALL"), includeOrder)
