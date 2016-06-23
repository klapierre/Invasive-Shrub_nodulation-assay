library(plyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(nlme)
library(lsmeans)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\nodulation experiment\\isotope data')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
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

isotope <- read.csv('La Pierre KJL_nod_pilots 0116.csv')%>%
  filter(d15N!='NA')

Ndfa <- isotope%>%
  select(per15N, species, inoculation_status)%>%
  group_by(species, inoculation_status)%>%
  summarise(per15N=mean(per15N))%>%
  spread(key=inoculation_status, value=per15N)%>%
  mutate(perNdfa=(100*(ctl-ino))/(ctl))

#plot percent 15N for each species*inoculation combination
ggplot(data=barGraphStats(data=isotope, variable='per15N', byFactorNames=c('species', 'inoculation_status')), aes(x=species, y=mean, fill=inoculation_status)) +
  geom_bar(stat='identity', position=position_dodge(), colour='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  scale_fill_manual(breaks=c('ctl', 'ino'),
                    labels=c('uninoculated', 'inoculated'),
                    values=c("#0072B2", "#D55E00")) +
  xlab('Plant Species') +
  ylab('Percent 15N')

#plot percent 15N for each species*inoculation combination as point graph, with notes about nodules annotated
ggplot(data=isotope, aes(x=species, y=per15N, colour=inoculation_status, shape=inoculation_status, label=nodule_notes)) +
  geom_point(size=5) +
  geom_text_repel(segment.color='white', size=5) +
  scale_colour_manual(breaks=c('ctl', 'ino'),
                    labels=c('uninoculated', 'inoculated'),
                    values=c("#0072B2", "#D55E00")) +
  xlab('Plant Species') +
  ylab('Percent 15N')

#plot percent N derived from atmosphere for each species
ggplot(data=Ndfa, aes(x=species, y=perNdfa)) +
  geom_bar(stat='identity', colour='black', fill='white') +
  xlab('Plant Species') +
  ylab('Percent N Derived from Atmosphere')




























