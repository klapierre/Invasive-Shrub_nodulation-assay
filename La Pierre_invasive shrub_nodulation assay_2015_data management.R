library(plyr)
library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\bigcb\\invasive shrubs project\\nodulation experiment')

#################################################
#################################################

#import strain information
strain <- read.csv('La Pierre_invasive shrub_noduation assay_2015_strains used.csv')%>%
  select(-original_host)

#import datafile
nodRaw <- read.csv('La Pierre_invasive shrub_noduation assay_2015.csv')%>%
  #merge with strain information
  left_join(strain)%>%
  #add total nodule number and functional nodule number variables
  mutate(nod_total=clear+white+pink+brown+purple+red+green, nod_func=pink+purple+red)%>%
  #add total biomass and root:shoot ratio variables
  mutate(total_biomass=shoots_g+roots_g, rootshoot=roots_g/shoots_g)%>%
  #remove plants that did not get inoculated or that died before harvest and control plants
  filter(harvest_date!='NA')%>%
  #remove CYSC plants for not, because not inoculated with own rhizobia, and LUBI because so many died
  filter(plant!='CYSC', plant!='LUBI')%>%
  #remove two pots that had been inoculated with the same strains as other pots
  filter(pot!=545, pot!=469)

#create binary table of what plants were nodulated
nodBinary <- nodRaw%>%
  select(pot, block, plant, plant_status, strain_label, original_host, host_match, nod_total, concatenated_OTU, ITS_OTU, nifd_OTU, strain_taxonomy)%>%
  #remove controls
  filter(host_match!='control')%>%
  #disregard whether the original host plant is the same as the test plant
  mutate(original_status=ifelse(original_host %in% c('ACGL', 'ACWR', 'LUAR', 'LUBI', 'LUNA', 'ACHE', 'ACST'),  'native', 'invasive'))%>%
  #create binary nodulation variable
  mutate(nod_binary=ifelse(nod_total>0, 1, 0))%>%
  mutate(is_conspecific=(host_match == 'original'))%>%
  mutate(host_match_factor=plyr::mapvalues(host_match, c('original', 'native', 'invader'), c('Conspecific', 'Native allospecific', 'Invasive allospecific')))

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
  mutate(host_status=ifelse(original_host %in% c('GEMO', 'MEPO', 'SPJU', 'ULEU', 'Vicia'), 'invasive', ifelse(original_host %in% c('ACGL', 'ACHE', 'ACST', 'ACWR', 'LUAR', 'LUBI', 'LUNA'), 'native', 'control')))

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

#subset out only plants inoculated with strains that have molecular data (i.e., assigned OTUs)
nodOTU <- nodRaw%>%
  filter(concatenated_OTU!='NA', concatenated_OTU!='')

#write OTU file
nodOTUtable <- nodOTU%>%
  select(plant, strain_label, concatenated_OTU, ITS_OTU, nifd_OTU, nod_total)
write.csv(nodOTUtable, 'La Pierre_invasive shrub_nodulation assay_OTU strains and field spp.csv')
