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

#mixed effects model for nodulation success (binary nodulation variable)
nodBinaryModel <- glmer(nod_binary ~ plant_status*host_match + (1|plant), data=nodBinary, family=binomial)
summary(nodBinaryModel)
lsmeans(nodBinaryModel, cld~plant_status*host_match)
#get table of estimates with 95% CI
se <- sqrt(diag(vcov(nodBinaryModel))) #standard errors
nodBinaryModelResults <- exp(cbind(estimate=fixef(nodBinaryModel), lower=fixef(nodBinaryModel)-1.96*se, upper=fixef(nodBinaryModel)+1.96*se)) #exponentiate coeffecients to get odds ratios

#get predicted probabilities from glmer model
jvalues <- with(nodBinary, list(host_match))
predProb <- lapply(jvalues, function(j){
  nodBinary$host_match <- j
  predict(nodBinaryModel, newdata=nodBinary, type='response')
})
#bind predicted probabilities to main dataset
predProbInfo <- cbind(nodBinary, predProb)
names(predProbInfo)[names(predProbInfo)=='structure(c(0.842005311020883, 0.842005311020883, 0.842005311020883, '] <- 'pred_prob'
#keep only relevant information
predProbInfo <- predProbInfo%>%  
  select(plant, plant_status, host_match, concatenated_OTU, ITS_OTU, nifd_OTU, original_status, pred_prob)

#predicted probabilities of nodulation for genotyped strains
predProbGeno <- predProbInfo%>%
  group_by(plant, concatenated_OTU)%>%
  summarise(mean_pred_prob=mean(pred_prob))


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

#for host_match: can test overall if there is a difference in prob of nodulation within native or invasive for the three host_match categories
#for plant_status: within each host_match category, is there a difference between the plant_status categories

#can test model accuracy by leo bryman approach: pull out 20% of data, run model on the 80% and see how well it predicts outcome in remaining 20%
#can do this first on the plant status question, then do it with the added information of host-match and see does the prediction improve a lot?
#

# #plot probability of nodulation by original host status and plant status
# ggplot(data=barGraphStats(data=predProbInfo, variable='pred_prob', byFactorNames=c('plant_status', 'host_match')), aes(x=plant_status, y=mean, fill=host_match)) +
#   geom_bar(stat='identity', position=position_dodge(), colour='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   scale_x_discrete(limits=c('native', 'invasive')) +
#   scale_fill_manual(name='Isolate Source',
#                     breaks=c('native', 'invader', 'original'),
#                     labels=c('native allospecific', 'invasive allospecific', 'conspecific'),
#                     values=c("#636363", "#bdbdbd", "#FFFFFF")) +
#   xlab('Test-Legume Status') +
#   ylab('Probability of Nodulation')

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




# #subset out only the plant species relevent to field data from proportional nodulation data
# nodPropField <- nodProp%>%
#   # filter(plant=='ACGL' | plant=='GEMO' | plant=='LUAR' | plant=='SPJU' | plant=='ULEU')%>%
#   filter(host_match!='control')


# #native vs invasive, without conspecific strains
# nodPropFieldAllospp <- lme(nod_proportion_total ~ plant_status*host_match, random=~1|plant, data=subset(nodPropField, host_match!='original'))
# summary(nodPropFieldAllospp)
# anova(nodPropFieldAllospp)
# lsmeans(nodPropFieldAllospp, cld~plant_status*host_match)
# 
# allosppFig <- ggplot(data=barGraphStats(data=subset(nodPropField, host_match!='original'), variable='nod_proportion_total', byFactorNames=c('plant_status', 'host_match')), aes(x=plant_status, y=mean, fill=host_match)) +
#   geom_bar(stat='identity', position=position_dodge(), colour='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   scale_y_continuous(limits=c(0,1)) +
#   scale_x_discrete(limits=c('native', 'invasive')) +
#   scale_fill_manual(breaks=c('native', 'invader'),
#                     labels=c('native', 'invasive'),
#                     values=c("#636363", "#bdbdbd")) +
#   xlab('Plant Status') +
#   ylab('Proportion Strains Nodulating') +
#   theme(axis.title.y=element_text(margin=margin(r=15)),
#         axis.title.x=element_text(margin=margin(t=10))) +
#   annotate('text', x=0.45, y=1, label='(a)', size=8, hjust='left')

# #mixed effects model for proportion nodulating (total nodules)
# nodPropFieldModel <- lme(nod_proportion_total ~ plant_status*host_match, random=~1|plant, data=nodPropField)
# summary(nodPropFieldModel)
# anova(nodPropFieldModel)
# lsmeans(nodPropFieldModel, cld~plant_status*host_match)
# 
# #plot proportion strains nodulating by host match and plant status
# ggplot(data=barGraphStats(data=nodPropField, variable='nod_proportion_total', byFactorNames=c('plant_status', 'host_match')), aes(x=plant_status, y=mean, fill=host_match)) +
#   geom_bar(stat='identity', position=position_dodge(), colour='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   scale_x_discrete(limits=c('native', 'invasive')) +
#   scale_fill_manual(name='Isolate Source',
#                     breaks=c('native', 'invader', 'original'),
#                     labels=c('native allospecific', 'invasive allospecific', 'conspecific'),
#                     values=c("#636363", "#bdbdbd", "#FFFFFF")) +
#   xlab('Test-Legume Status') +
#   ylab('Proportion Isolates Nodulating')
# 
# 
# #boxplot of proportion strains nodulating by plant status
# ggplot(data=barGraphStats(data=nodPropField, variable='nod_proportion_total', byFactorNames=c('plant_status', 'plant')), aes(x=plant_status, y=mean, label=plant)) +
#   geom_boxplot() +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
#   geom_text(hjust='left', vjust='center', nudge_x=0.05, size=6) +
#   scale_x_discrete(limits=c('native', 'invasive')) +
#   scale_y_continuous(breaks=seq(0, 1, 0.2), name="Proportion Isolates Nodulating") +
#   coord_cartesian(ylim=c(0, 1)) +
#   xlab("Test-Host Status")


# #exploratory analysis
# 
# ###mixed effects model for height with plant status (invasive/native) and host match (original host/local Bay Area host/away invasive host) and plant species as a random factor
# 
# #mixed effects model
# nodRegHeight <- lme(height_cm ~ plant_status*host_match, random=~1|plant, data=nodRaw)
# summary(nodRegHeight)
# anova(nodRegHeight)
# lsmeans(nodRegHeight, cld~plant_status*host_match)
# 
# # #testing assumptions - mostly normal data
# # plot(ranef(nodRegHeight))
# # resNodRegHeight<-residuals(nodRegHeight)
# # plot(resNodRegHeight)
# # qqnorm(resNodRegHeight)
# # qqline(resNodRegHeight)
# # plot(nodRegHeight)
# 
# ggplot(data=barGraphStats(data=nodRaw, variable="height_cm", byFactorNames=c("plant_status", "host_match")), aes(x=plant_status, y=mean, fill=host_match)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Height (cm)') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# 
# 
# ###mixed effects model for nodule number with plant status (invasive/native) and host match (original host/local Bay Area host/away invasive host) and plant species as a random factor
# 
# #mixed effects model
# nodRegTot <- lme(nod_total ~ plant_status*host_match, random=~1|plant, data=nodRaw)
# summary(nodRegTot)
# anova(nodRegTot)
# lsmeans(nodRegTot, cld~plant_status*host_match)
# 
# # #testing assumptions - not quite normal data
# # plot(ranef(nodRegTot))
# # resNodRegTot<-residuals(nodRegTot)
# # plot(resNodRegTot)
# # qqnorm(resNodRegTot)
# # qqline(resNodRegTot)
# # plot(nodRegTot)
# 
# ggplot(data=barGraphStats(data=nodRaw, variable="nod_total", byFactorNames=c("plant_status", "host_match")), aes(x=plant_status, y=mean, fill=host_match)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Nodule Number') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# 
# ggplot(data=barGraphStats(data=nodRaw, variable="nod_func", byFactorNames=c("plant_status", "host_match")), aes(x=plant_status, y=mean, fill=host_match)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Functional Nodule Number') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# 
# ggplot(data=barGraphStats(data=nodRaw, variable="total_biomass", byFactorNames=c("plant_status", "host_match")), aes(x=plant_status, y=mean, fill=host_match)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Total Biomass (g)') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# 
# ggplot(data=barGraphStats(data=nodRaw%>%filter(host_match!='control'), variable="total_biomass", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Total Biomass (g)') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# 
# #for only those strains that were compatible
# ggplot(data=barGraphStats(data=subset(nodRaw, subset=(nod_total>0)), variable="total_biomass", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Total Biomass (g)')
# 
# 
# ggplot(data=barGraphStats(data=nodProp%>%filter(host_match!='control'), variable="nod_proportion", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Proportion of Strains Nodulating') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# 
# ggplot(data=barGraphStats(data=nodProp%>%filter(host_match!='control'), variable="nod_proportion_total", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Proportion of Strains Nodulating') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# ggplot(data=barGraphStats(data=nodProp%>%filter(host_match!='control'), variable="nod_proportion_func", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Proportion of Strains\nMaking Functional Nodules') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# 
# ggplot(data=barGraphStats(data=nodRaw%>%filter(host_match!='control'), variable="nod_total", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Status') + ylab('Number of Nodules') +
#   scale_shape_discrete(name='Strain Origin') +
#   theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!
# 
# 
# 
# 
# 
# 
# ###strain by species interaction figures
# 
# ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   xlab('Strain') + ylab('Total Biomass (g)') +
#   theme(axis.text.x=element_text(angle=90)) +
#   ggtitle('ACGL')
# 
# ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   xlab('Strain') + ylab('Total Biomass (g)') +
#   theme(axis.text.x=element_text(angle=90)) +
#   ggtitle('ACWR')
# 
# ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   xlab('Strain') + ylab('Total Biomass (g)') +
#   theme(axis.text.x=element_text(angle=90)) +
#   ggtitle('GEMO')
# 
# ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   xlab('Strain') + ylab('Total Biomass (g)') +
#   theme(axis.text.x=element_text(angle=90)) +
#   ggtitle('LUAR')
# 
# ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   xlab('Strain') + ylab('Total Biomass (g)') +
#   theme(axis.text.x=element_text(angle=90)) +
#   ggtitle('SPJU')
# 
# ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   xlab('Strain') + ylab('Total Biomass (g)') +
#   theme(axis.text.x=element_text(angle=90)) +
#   ggtitle('ULEU')
# 
# 
# ggplot(data=barGraphStats(data=nodSize, variable="bio", byFactorNames=c("plant", "ino_type")), aes(x=plant, y=mean, fill=as.factor(ino_type))) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Plant Species') + ylab('Total Biomass (g)')
# 
# 
# ggplot(data=barGraphStats(data=nodSize%>%filter(plant=='ACWR'), variable="bio", byFactorNames=c("original_host", "noded")), aes(x=original_host, y=mean, fill=as.factor(noded))) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2)
# 
# #from rhizobial side of things
# #total nodules formed
# ggplot(data=barGraphStats(data=rhizProp%>%filter(host_status!='control'), variable="rhiz_proportion_total", byFactorNames=c("original_host", "host_status")), aes(x=host_status, y=mean, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Rhizobial Strain Origin') + ylab('Proportions of Plants Nodulated')
# 
# ggplot(data=barGraphStats(data=rhizProp%>%filter(host_status!='control'), variable="rhiz_proportion_total", byFactorNames=c("host_status")), aes(x=host_status, y=mean)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Rhizobial Strain Origin') + ylab('Proportions of Plants Nodulated')
# 
# 
# #functional nodules formed
# ggplot(data=barGraphStats(data=rhizProp%>%filter(host_status!='control'), variable="rhiz_proportion_func", byFactorNames=c("original_host", "host_status")), aes(x=host_status, y=mean, fill=original_host)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Rhizobial Strain Origin') + ylab('Proportions of Plants Nodulated')
# 
# ggplot(data=barGraphStats(data=rhizProp%>%filter(host_status!='control'), variable="rhiz_proportion_func", byFactorNames=c("host_status")), aes(x=host_status, y=mean)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
#   xlab('Rhizobial Strain Origin') + ylab('Proportions of Plants Nodulated')



data <- nodBinary

# Sort out plant status and bacteria origin for later use
named.vector.from.data.frame <- function(data) {
  result <- data[[2]]
  names(result) <- as.character(data[[1]])
  result
}
plant.to.status <- (
  data
  %>% group_by(plant)
  %>% summarize(status={stopifnot(length(unique(plant_status)) == 1); plant_status[1]})
  %>% named.vector.from.data.frame()
)
bacteria.to.origin <- (
  data
  %>% group_by(strain_label)
  %>% summarize(status={stopifnot(length(unique(original_status)) == 1); original_status[1]})
  %>% named.vector.from.data.frame()
)

# Heatmap preliminaries
heatmap.data <- data%>%
  select(plant, strain_label, nod_binary, original_status, original_host)%>%
  spread(plant, nod_binary)
heatmap.data <- heatmap.data[order(heatmap.data$original_status, heatmap.data$original_host),]%>%
  select(strain_label, ACGL, ACWR, GEMO, LUAR, SPJU, ULEU)

heatmap.matrix <- as.matrix(heatmap.data[, -1])
rownames(heatmap.matrix) <- heatmap.data[, 1]
heatmap.matrix[is.na(heatmap.matrix)] <- -1

reverse.na.dist <- function(x, method) {
  x[x == -1] <- NA
  dist(x, method=method)
}
draw.heatmap <- function(heatmap.matrix, distance.method, colors, breaks) {
  row.colors <- ifelse(bacteria.to.origin[rownames(heatmap.matrix)] == 'native', 'green', 'blue')
  col.colors <- ifelse(plant.to.status[colnames(heatmap.matrix)] == 'native', 'green', 'blue')
  
  heatmap(
    heatmap.matrix,
    distfun=function(x) reverse.na.dist(x, method=distance.method),
    scale='none',
    col=colors,
    breaks=breaks,
    RowSideColors=row.colors,
    ColSideColors=col.colors,
    margins=c(5,18),
    Rowv=NA
  )
  legend(
    "right",
    legend=c('Native', 'Invasive'),
    title='Plant status and\nbacteria origin',
    fill=c('green', 'blue'),
    bty='n'
  )
}

# Count heatmap
colors <- c('black', brewer.pal(9, 'YlOrRd'))
node.count.quantiles <- quantile(
  data$nod_total[data$nod_total > 0],
  seq(0, 1, length.out=length(colors) - 2),
  na.rm=TRUE
)
breaks <- c(-2, -0.5, 0, node.count.quantiles)
draw.heatmap(heatmap.matrix, 'manhattan', colors, breaks)
legend(
  'topright',
  legend=c('Missing', '0', sprintf('(%.1f, %.1f]', breaks[-c(1, 2, length(breaks))], breaks[-(1:3)])),
  title='Nodule count',
  fill=colors,
  bty='n'
)

# Binary heatmap
binary.heatmap.matrix <- ifelse(heatmap.matrix > 0, 1, heatmap.matrix)
binary.heatmap.matrix <- binary.heatmap.matrix[order(bacteria.to.origin[rownames(binary.heatmap.matrix)]),]
binary.colors <- c('white', 'grey', 'black')
binary.breaks <- c(-2, -0.5, 0.5, 2)
draw.heatmap(binary.heatmap.matrix, 'binary', binary.colors, binary.breaks)
legend(
  'topright',
  legend=c('Missing', 'No', 'Yes'),
  title='Nodulated?',
  fill=binary.colors,
  bty='n'
)
