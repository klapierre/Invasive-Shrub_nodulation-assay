library(plyr)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\nodulation experiment')

#################################################
#################################################

#import datafile
nodRaw <- read.csv('La Pierre_invasive shrub_noduation assay_2015.csv')%>%
  #add total nodule number and functional nodule number variables
  mutate(nod_total=clear+white+pink+brown+purple+red+green, nod_func=pink+purple+red)%>%
  #add total biomass and root:shoot ratio variables
  mutate(total_biomass=shoots_g+roots_g, rootshoot=roots_g/shoots_g)%>%
  #remove plants that did not get inoculated or that died before harvest and control plants
  filter(harvest_date!='NA')%>%
  #remove CYSC plants for not, because not inoculated with own rhizobia, and LUBI because so many died
  filter(plant!='CYSC', plant!='LUBI')

#get proportion of plants that nodulated for each category
nodProp <- nodRaw%>%
  mutate(nodulated_total=ifelse(nod_total>0, 1, 0), nodulated_func=ifelse(nod_func>0, 1, 0))%>%
  select(pot, plant, plant_status, strain_label, original_host, host_match, nodulated_total, nodulated_func, nod_total, nod_func)%>%
  group_by(plant, plant_status, host_match, original_host)%>%
  summarise(nod_N=length(nod_total), nod_total_count=sum(nodulated_total), nod_func_count=sum(nodulated_func))%>%
  mutate(nod_proportion_total=nod_total_count/nod_N, nod_proportion_func=nod_func_count/nod_N)

#get proportion of plants each strain can affiliate with
rhizProp <- nodRaw%>%
  mutate(total_nod_formed=ifelse(nod_total>0, 1, 0))%>%
  mutate(func_nod_formed=ifelse(nod_func>0, 1, 0))%>%
  select(pot, plant, plant_status, strain_label, original_host, host_match, total_nod_formed, func_nod_formed, nod_total, nod_func)%>%
  group_by(strain_label, original_host)%>%
  summarise(rhiz_count_total=sum(total_nod_formed), rhiz_count_func=sum(func_nod_formed))%>%
  mutate(rhiz_proportion_total=rhiz_count_total/6, rhiz_proportion_func=rhiz_count_func/6)%>%
  mutate(host_status=ifelse(original_host=='GEMO', 'invasive', ifelse(original_host=='MEPO', 'invasive', ifelse(original_host=='SPJU', 'invasive', ifelse(original_host=='ULEU', 'invasive', ifelse(original_host=='Vicia', 'invasive', ifelse(original_host=='ACGL', 'native', ifelse(original_host=='ACHE', 'native', ifelse(original_host=='ACST', 'native', ifelse(original_host=='ACWR', 'native', ifelse(original_host=='LUAR', 'native', ifelse(original_host=='LUBI', 'native', ifelse(original_host=='LUNA', 'native', 'control')))))))))))))

###make interaction matrix for plant by strain origin interactions (with proportion of the strains nodulating as the data)

#total nodules
interactionTotal <- nodProp%>%
  filter(host_match!='control')%>%
  ungroup()%>%
  select(plant, plant_status, original_host, nod_proportion_total)%>%
  spread(key=original_host, value=nod_proportion_total, fill=NA)

#functional nodules
interactionFunc <- nodProp%>%
  filter(host_match!='control')%>%
  ungroup()%>%
  select(plant, plant_status, original_host, nod_proportion_func)%>%
  spread(key=original_host, value=nod_proportion_func, fill=NA)


#figure out which strains didn't nodulate anything
nodNone <- nodRaw%>%
  group_by(strain_label)%>%
  summarise(tot_nodules=sum(nod_total))

#make a column with control vs not control
nodSize <- nodRaw%>%
  mutate(total_nod=ifelse(nod_total==0, 0, 1), func_nod=ifelse(nod_func==0, 0, 1), ino_status=ifelse(host_match=='control', 'ctl', 'ino'), ino_type=paste(total_nod, func_nod, ino_status, sep='::'))%>%
  group_by(plant, original_host, func_nod, total_nod, ino_type)%>%
  summarise(bio=mean(total_biomass))