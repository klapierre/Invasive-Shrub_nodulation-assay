library(plyr)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(nlme)
library(lme4)
# library(dist.R)
library(lsmeans)
library(multcomp)
library(multcompView)
library(bipartite)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\nodulation experiment')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20, margin=margin(b=20)), legend.text=element_text(size=16))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  


#################################################
#################################################

#source data management code
source('Invasive-Shrub_nodulation-assay\\La Pierre_invasive shrub_nodulation assay_2015_data management.R')

#bonferroni correction for multiple tests (3 tests, so corrected by 3)
bonferroniCorrection <- 1-0.05/3

#keep only invasives
nodBinaryInv <- subset(nodBinary, subset=(plant_status=='invasive'))

#do conspecifics nodulate more than allospecifics? do analysis for each plant separately (we don't care about differences among the three legumes)
summary(conspecificGEMOmodel <- glm(nod_binary ~ is_conspecific, data=subset(nodBinaryInv, plant=='GEMO'), family=binomial))
summary(conspecificSPJUmodel <- glm(nod_binary ~ is_conspecific, data=subset(nodBinaryInv, plant=='SPJU'), family=binomial))
summary(conspecificULEUmodel <- glm(nod_binary ~ is_conspecific, data=subset(nodBinaryInv, plant=='ULEU'), family=binomial))

  #get confidence intervals around conspecific traits
  print(conspecificCI <- exp(confint(conspecificGEMOmodel, 'is_conspecificTRUE', level=0.95)))
  print(conspecificCI <- exp(confint(conspecificSPJUmodel, 'is_conspecificTRUE', level=0.95)))
  print(conspecificCI <- exp(confint(conspecificULEUmodel, 'is_conspecificTRUE', level=0.95)))

#do native or invasive allospecifics nodulate more the three invasive plants?
#create table of the numbers of strains that are allospecific nodulators for each invader
allospecificTableGEMO <- xtabs(~ host_match_factor + nod_binary, data=droplevels(filter(subset(nodBinaryInv, plant=='GEMO'), host_match %in% c('native', 'invader'))))
#fisher test for native species
print(nativeTest <- fisher.test(allospecificTableGEMO[,], conf.level=bonferroniCorrection))
print(invasiveTest <- fisher.test(allospecificTableGEMO[,][c(2,1),], conf.level=bonferroniCorrection))

allospecificTableSPJU <- xtabs(~ host_match_factor + nod_binary, data=droplevels(filter(subset(nodBinaryInv, plant=='SPJU'), host_match %in% c('native', 'invader'))))
#fisher test for native species
print(nativeTest <- fisher.test(allospecificTableSPJU[,], conf.level=bonferroniCorrection))
print(invasiveTest <- fisher.test(allospecificTableSPJU[,][c(2,1),], conf.level=bonferroniCorrection))

allospecificTableULEU <- xtabs(~ host_match_factor + nod_binary, data=droplevels(filter(subset(nodBinaryInv, plant=='ULEU'), host_match %in% c('native', 'invader'))))
#fisher test for native species
print(nativeTest <- fisher.test(allospecificTableULEU[,], conf.level=bonferroniCorrection))
print(invasiveTest <- fisher.test(allospecificTableULEU[,][c(2,1),], conf.level=bonferroniCorrection))

###heat map depicting binary nodulation outcomes
#order the strains by original host type and species
nodBinaryInv <- nodBinaryInv[order(nodBinaryInv$original_host, nodBinaryInv$original_status),]

#sort out plant status and strain origin
nameFunction <- function(data) {
  result <- data[[2]]
  names(result) <- as.character(data[[1]])
  result
}

plantNames <- (
  nodBinaryInv%>%
    group_by(plant)%>%
    summarize(status={stopifnot(length(unique(plant_status))==1); plant_status[1]})%>%
    nameFunction()
)

bacteriaNames <- (
  nodBinaryInv%>%
    group_by(strain_label)%>%
    summarize(status={stopifnot(length(unique(original_status))==1); original_status[1]})%>%
    nameFunction()
)

#generate heatmap
heatmapNodulation <- nodBinaryInv%>%
  select(plant, strain_label, nod_binary, original_status, original_host)%>%
  spread(plant, nod_binary)
heatmapNodulation <- heatmapNodulation[order(heatmapNodulation$original_status, heatmapNodulation$original_host),]%>%
  select(strain_label, GEMO, SPJU, ULEU)
heatmapNodulation <- heatmapNodulation[c('strain_label', 'GEMO', 'SPJU', 'ULEU')]

heatmapNodulation[is.na(heatmapNodulation)] <- -1

heatmapMatrix <- as.matrix(heatmapNodulation[,-1])
rownames(heatmapMatrix) <- heatmapNodulation[,1]

reverseNAdist <- function(x,method) {
  x[x==-1] <- NA
  dist(x, method=method)
}

#binary heatmap
binaryHeatmapMatrix <- heatmapMatrix
binaryColors <- c('white', 'grey', 'black')
binaryBreaks <- c(-2,-0.5,0.5,2)
heatmap(binaryHeatmapMatrix, Rowv=NA, Colv=NA, margins=c(5,18), breaks=binaryBreaks, col=binaryColors, scale='none',
        RowSideColors=ifelse(bacteriaNames[rownames(heatmapMatrix)]=='native', '#009900', '#FF9900'),
        distfun=function(x) reverseNAdist(x, method='binary'))
legend('topright',
       legend=c('N/A', 'No', 'Yes'),
       title='Nodulation Success',
       fill=binaryColors,
       bty='n')
legend('right',
      legend=c('Native', 'Invasive'),
      title='Rhizobia Origin',
      fill=c('#009900', '#FF9900'),
      bty='n')
#export at 1200x800 and edit strain labels and rearrange components to make more clear