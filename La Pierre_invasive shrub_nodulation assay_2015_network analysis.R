library(plyr)
library(ggplot2)
library(grid)
library(bipartite)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\nodulation experiment')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

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


###make an interaction web for total nodules and functional nodules

#total nodules
matrixTotal <- interactionTotal
rownames(matrixTotal) <- matrixTotal$plant
matrixTotal <-  matrixTotal%>%
  select(-plant_status, -plant)

#network test from both directions
networkIndices <- specieslevel(t(matrixTotal))

#plot the interaction web

plantcolors<-c("#00990099", "#00990099", "#FF990099",  "#00990099", "#FF990099", "#FF990099")
straincolors<-c(rep("black"))


plotweb(t(matrixTotal), bor.col.interaction=plantcolors,  col.interaction=plantcolors, arrow="no", method="normal", labsize=1.5, abuns.type="independent",
        text.rot=90, low.lab.dis=NULL,
        y.width.high=0.005,y.width.low=0.005, col.low="black", col.high="black",
        low.abun.col="black", high.abun.col="black",
        bor.col.high="black", bor.col.low="black",
        sequence=list(seq.high=c('ACGL', 'ACWR', 'LUAR', 'GEMO', 'SPJU', 'ULEU'),
                 seq.low=c('ACGL', 'ACHE', 'ACST', 'ACWR', 'LUAR', 'LUBI', 'LUNA', 'GEMO', 'MEPO', 'SPJU', 'ULEU', 'Vicia')))


#functional nodules
matrixFunc <- interactionFunc
rownames(matrixFunc) <- matrixFunc$plant
matrixFunc <-  matrixFunc%>%
  select(-plant_status, -plant)

#network test from both directions
networkIndices <- specieslevel(t(matrixFunc))

#plot the interaction web

plantcolors<-c("#00990099", "#00990099", "#FF990099",  "#00990099", "#FF990099", "#FF990099")
straincolors<-c(rep("black"))


plotweb(t(matrixFunc), bor.col.interaction=plantcolors,  col.interaction=plantcolors, arrow="no", method="normal", labsize=1.5, abuns.type="independent",
        text.rot=90, low.lab.dis=NULL,
        y.width.high=0.005,y.width.low=0.005, col.low="black", col.high="black",
        low.abun.col="black", high.abun.col="black",
        bor.col.high="black", bor.col.low="black",
        sequence=list(seq.high=c('ACGL', 'ACWR', 'LUAR', 'GEMO', 'SPJU', 'ULEU'),
                      seq.low=c('ACGL', 'ACHE', 'ACST', 'ACWR', 'LUAR', 'LUBI', 'LUNA', 'GEMO', 'MEPO', 'SPJU', 'ULEU', 'Vicia')))























