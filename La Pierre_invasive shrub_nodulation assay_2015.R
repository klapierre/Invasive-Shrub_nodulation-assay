library(plyr)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(nlme)
library(lme4)
library(dist.R)
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

#do conspecifics nodulate more than allospecifics? do natives nodulate more than invasives?
summary(conspecificPlantStatusModel <- glm(nod_binary ~ plant_status + is_conspecific, data=nodBinary, family=binomial))

  #get confidence intervals around plant status and conspecific traits
  print(conspecificCI <- exp(confint(conspecificPlantStatusModel, 'is_conspecificTRUE')))
  print(plantStatusCI <- exp(confint(conspecificPlantStatusModel, 'plant_statusnative')))

#is there an interactions between plant status and rhizobial status (con vs allospecific) in nodulation success?
#test for interaction with nested model chi-squared test
summary(conspecificPlantStatusInteraction <- update(conspecificPlantStatusModel, . ~ . + plant_status:is_conspecific))
print(conspecificPlantStatusInteractionPvalue <- anova(conspecificPlantStatusModel, conspecificPlantStatusInteraction, test='Chisq')[2,'Pr(>Chi)'])

#do native allospecifics nodulate more with native plants and invasive allospecifics nodulate more with invasive plants?
#create table of the numbers of strains that are allospecific nodulators for natives and invasives
allospecificTable <- xtabs(~ plant_status + host_match_factor + nod_binary, data=droplevels(filter(nodBinary, host_match %in% c('native', 'invader'))))
#fisher test for native species
print(nativeTest <- fisher.test(allospecificTable['native',,]))
print(invasiveTest <- fisher.test(allospecificTable['invasive',,][c(2,1),]))




###heat map depicting binary nodulation outcomes
#order the strains by original host type and species
nodBinary <- nodBinary[order(nodBinary$original_host, nodBinary$original_status),]

#sort out plant status and strain origin
nameFunction <- function(data) {
  result <- data[[2]]
  names(result) <- as.character(data[[1]])
  result
}

plantNames <- (
  nodBinary%>%
    group_by(plant)%>%
    summarize(status={stopifnot(length(unique(plant_status))==1); plant_status[1]})%>%
    nameFunction()
)

bacteriaNames <- (
  nodBinary%>%
    group_by(strain_label)%>%
    summarize(status={stopifnot(length(unique(original_status))==1); original_status[1]})%>%
    nameFunction()
)

#generate heatmap
heatmapNodulation <- nodBinary%>%
  select(plant, strain_label, nod_binary, original_status, original_host)%>%
  spread(plant, nod_binary)
heatmapNodulation <- heatmapNodulation[order(heatmapNodulation$original_status, heatmapNodulation$original_host),]%>%
  select(strain_label, ACGL, ACWR, GEMO, LUAR, SPJU, ULEU)

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
heatmap(binaryHeatmapMatrix, Rowv=NA, margins=c(5,18), breaks=binaryBreaks, col=binaryColors, scale='none',
        RowSideColors=ifelse(bacteriaNames[rownames(heatmapMatrix)]=='native', '#595959', '#A2A2A2'),
        ColSideColors=ifelse(plantNames[colnames(heatmapMatrix)]=='native', '#595959', '#A2A2A2'),
        distfun=function(x) reverseNAdist(x, method='binary'))
legend('topright',
       legend=c('N/A', 'No', 'Yes'),
       title='Nodulation Success',
       fill=binaryColors,
       bty='n')
legend('right',
      legend=c('Native', 'Invasive'),
      title='Plant Status and\nRhizobia Origin',
      fill=c('#595959', '#A2A2A2'),
      bty='n')
#export at 1200x800 and edit strain labels and rearrange components to make more clear

#boxplot of probability of nodulation by plant status
ggplot(data=barGraphStats(data=predProbInfo, variable='pred_prob', byFactorNames=c('plant_status', 'plant', 'host_match')), aes(x=plant_status, y=mean, fill=host_match, label=plant)) +
  geom_boxplot() +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  # geom_text(hjust='left', vjust='center', nudge_x=0.05, size=6) +
  scale_x_discrete(limits=c('native', 'invasive')) +
  scale_y_continuous(breaks=seq(0.5, 1, 0.1), name="Probability of Nodulation") +
  scale_fill_manual(name='Isolate Origin',
                    breaks=c('invader', 'native', 'original'),
                    labels=c('invasive allospecific', 'native allospecific', 'conspecific'),
                    values=c("#636363", "#bdbdbd", "#FFFFFF")) +
  coord_cartesian(ylim=c(0.45, 1)) +
  xlab("Test-Host Status") +
  theme(legend.justification=c(0,0), legend.position=c(0,0))
#export at 700x700