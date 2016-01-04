library(plyr)
library(ggplot2)
library(grid)
library(nlme)
library(lsmeans)
library(multcomp)
library(multcompView)
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

###mixed effects model for height with plant status (invasive/native) and host match (original host/local Bay Area host/away invasive host) and plant species as a random factor

#mixed effects model
nodRegHeight <- lme(height_cm ~ plant_status*host_match, random=~1|plant, data=nodRaw)
summary(nodRegHeight)
anova(nodRegHeight)
lsmeans(nodRegHeight, cld~plant_status*host_match)

# #testing assumptions - mostly normal data
# plot(ranef(nodRegHeight))
# resNodRegHeight<-residuals(nodRegHeight)
# plot(resNodRegHeight)
# qqnorm(resNodRegHeight)
# qqline(resNodRegHeight)
# plot(nodRegHeight)

ggplot(data=barGraphStats(data=nodRaw, variable="height_cm", byFactorNames=c("plant_status", "host_match")), aes(x=plant_status, y=mean, fill=host_match)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Height (cm)') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!



###mixed effects model for nodule number with plant status (invasive/native) and host match (original host/local Bay Area host/away invasive host) and plant species as a random factor

#mixed effects model
nodRegTot <- lme(nod_total ~ plant_status*host_match, random=~1|plant, data=nodRaw)
summary(nodRegTot)
anova(nodRegTot)
lsmeans(nodRegTot, cld~plant_status*host_match)

# #testing assumptions - not quite normal data
# plot(ranef(nodRegTot))
# resNodRegTot<-residuals(nodRegTot)
# plot(resNodRegTot)
# qqnorm(resNodRegTot)
# qqline(resNodRegTot)
# plot(nodRegTot)

ggplot(data=barGraphStats(data=nodRaw, variable="nod_total", byFactorNames=c("plant_status", "host_match")), aes(x=plant_status, y=mean, fill=host_match)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Nodule Number') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!


ggplot(data=barGraphStats(data=nodRaw, variable="nod_func", byFactorNames=c("plant_status", "host_match")), aes(x=plant_status, y=mean, fill=host_match)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Functional Nodule Number') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!


ggplot(data=barGraphStats(data=nodRaw, variable="total_biomass", byFactorNames=c("plant_status", "host_match")), aes(x=plant_status, y=mean, fill=host_match)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Total Biomass (g)') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!


ggplot(data=barGraphStats(data=nodRaw%>%filter(host_match!='control'), variable="total_biomass", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Total Biomass (g)') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!


#for only those strains that were compatible
ggplot(data=barGraphStats(data=subset(nodRaw, subset=(nod_total>0)), variable="total_biomass", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Total Biomass (g)')


ggplot(data=barGraphStats(data=nodProp%>%filter(host_match!='control'), variable="nod_proportion", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Proportion of Strains Nodulating') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!


ggplot(data=barGraphStats(data=nodProp%>%filter(host_match!='control'), variable="nod_proportion_total", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Proportion of Strains Nodulating') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!

ggplot(data=barGraphStats(data=nodProp%>%filter(host_match!='control'), variable="nod_proportion_func", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Proportion of Strains\nMaking Functional Nodules') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!


ggplot(data=barGraphStats(data=nodRaw%>%filter(host_match!='control'), variable="nod_total", byFactorNames=c("plant", "original_host")), aes(x=plant, y=mean, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Status') + ylab('Number of Nodules') +
  scale_shape_discrete(name='Strain Origin') +
  theme(legend.title=element_text('Strain Origin') ) #can't get legend title to work!






###strain by species interaction figures

ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  xlab('Strain') + ylab('Total Biomass (g)') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('ACGL')

ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  xlab('Strain') + ylab('Total Biomass (g)') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('ACWR')

ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  xlab('Strain') + ylab('Total Biomass (g)') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('GEMO')

ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  xlab('Strain') + ylab('Total Biomass (g)') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('LUAR')

ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  xlab('Strain') + ylab('Total Biomass (g)') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('SPJU')

ggplot(data=nodRaw%>%filter(plant=='ACGL'), aes(x=strain_label, y=total_biomass, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  xlab('Strain') + ylab('Total Biomass (g)') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('ULEU')


ggplot(data=barGraphStats(data=nodSize, variable="bio", byFactorNames=c("plant", "ino_type")), aes(x=plant, y=mean, fill=as.factor(ino_type))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Plant Species') + ylab('Total Biomass (g)')


ggplot(data=barGraphStats(data=nodSize%>%filter(plant=='ACWR'), variable="bio", byFactorNames=c("original_host", "noded")), aes(x=original_host, y=mean, fill=as.factor(noded))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2)

#from rhizobial side of things
#total nodules formed
ggplot(data=barGraphStats(data=rhizProp%>%filter(host_status!='control'), variable="rhiz_proportion_total", byFactorNames=c("original_host", "host_status")), aes(x=host_status, y=mean, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Rhizobial Strain Origin') + ylab('Proportions of Plants Nodulated')

ggplot(data=barGraphStats(data=rhizProp%>%filter(host_status!='control'), variable="rhiz_proportion_total", byFactorNames=c("host_status")), aes(x=host_status, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Rhizobial Strain Origin') + ylab('Proportions of Plants Nodulated')


#functional nodules formed
ggplot(data=barGraphStats(data=rhizProp%>%filter(host_status!='control'), variable="rhiz_proportion_func", byFactorNames=c("original_host", "host_status")), aes(x=host_status, y=mean, fill=original_host)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Rhizobial Strain Origin') + ylab('Proportions of Plants Nodulated')

ggplot(data=barGraphStats(data=rhizProp%>%filter(host_status!='control'), variable="rhiz_proportion_func", byFactorNames=c("host_status")), aes(x=host_status, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab('Rhizobial Strain Origin') + ylab('Proportions of Plants Nodulated')














  

