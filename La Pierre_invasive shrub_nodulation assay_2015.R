library(plyr)
library(ggplot2)
library(reshape2)
library(grid)
library(nlme)
library(lsmeans)
library(multcomp)
library(multcompView)
library(reshape2)
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

#import datafile
nodRaw <- read.csv('La Pierre_invasive shrub_noduation assay_2015.csv')%>%
  #add total nodule number and functional nodule number variables
  mutate(nod_total=clear+white+pink+brown+purple+red+green, nod_func=pink+purple+red)%>%
  #add total biomass and root:shoot ratio variables
  mutate(total_biomass=shoots_g+roots_g, rootshoot=roots_g/shoots_g)%>%
  #remove plants that did not get inoculated or that died before harvest and control plants
  filter(harvest_date!='NA', host_match!='control')%>%
  #remove CYSC plants for not, because not inoculated with own rhizobia, and LUBI because so many died
  filter(plant!='CYSC', plant!='LUBI')

#get proportion of plants that nodulated for each category
nodCore <- nodRaw%>%
  mutate(nodulated=ifelse(nod_total>0, 1, 0))%>%
  select(pot, plant, plant_status, strain_label, host_match, nodulated, nod_total)%>%
  group_by(plant, plant_status, host_match)%>%
  summarise(nod_N=length(nod_total), nod_count=sum(nodulated))%>%
  mutate(nod_proportion=nod_count/nod_N)














