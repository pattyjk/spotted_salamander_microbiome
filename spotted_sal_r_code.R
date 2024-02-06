library(vegan)
library(ggplot2)
library(tidyverse)
library(stringi)
library(reshape2)
library(ggpubr)
library(dplyr)
library(plyr)

#read in metadata for each experiment
juv_meta<-read.delim('spotted_salamander_microbiome/juvenile_exp_meta.txt', header=T)
larv_meta<-read.delim('spotted_salamander_microbiome/larval_exp_meta.txt', header=T)
life_meta<-read.delim('spotted_salamander_microbiome/life_stages_metadata.txt', header=T)

#read in asv table
asv.tbl<-read.delim("spotted_salamander_microbiome/spot_sal_asv_table.txt", row.names=1, header=T)
dim(asv.tbl)
#2106  241

#remove low sequencing depth samples
asv.tbl<-asv.tbl[,-which(colSums(asv.tbl)<400)]
dim(asv.tbl)
#2106  225

#rarefy table
set.seed(555)
asv.rare<-rrarefy(t(asv.tbl), sample = 400)
dim(asv.tbl)

#create asv table for life stages experiment
life_stg_tbl<-asv.rare[row.names(asv.rare) %in% life_meta$SampleID, ]

#create asv table for larval experiment
larv_tbl<-asv.rare[row.names(asv.rare) %in% larv_meta$SampleID, ]

#create asv table for juvenile experiment
juv_tbl<-asv.rare[row.names(asv.rare) %in% juv_meta$SampleID, ]

#calculate beta diversity
life.pcoa<-capscale(life_stg_tbl~1, method='bray')
larv.pcoa<-capscale(larv_tbl~1, method='bray')

#calculate percent variation explained
100*round(life.pcoa$CA$eig[1]/sum(life.pcoa$CA$eig), 3)
#23.9
100*round(life.pcoa$CA$eig[2]/sum(life.pcoa$CA$eig), 3)
#13.6

100*round(larv.pcoa$CA$eig[1]/sum(larv.pcoa$CA$eig), 3)
#42.4
100*round(larv.pcoa$CA$eig[2]/sum(larv.pcoa$CA$eig), 3)
#22.9


#extract the coordinates
larv.scores<-scores(larv.pcoa)
larv.coords<-as.data.frame(larv.scores$sites)
larv.coords$SampleID<-rownames(larv.coords)
larv.coords<-merge(larv.coords, larv_meta, by=c('SampleID'))

life.scores<-scores(life.pcoa)
life.coords<-as.data.frame(life.scores$sites)
life.coords$SampleID<-rownames(life.coords)
life.coords<-merge(life.coords, life_meta, by=c('SampleID'))
life.coords$LifeStage<-gsub("PondWaterControl", "Pond Water", life.coords$LifeStage)

larv.pcoa2<- 
  
  ggplot(larv.coords, aes(MDS1, MDS2, color=Treatment2))+  
  geom_point(size=4)+
  theme_bw()+
    #geom_mark_circle(aes(color = Treatment2))+
    scale_color_manual(values=c("#ABDAF4", "#F4CF7F", "#7F7F7F"))+
  xlab("PC1-23.9%")+
      ylab("PC2-11.8%")+
  theme(legend.position = c(0.8, 0.2))+
  theme(text = element_text(size=14),
        axis.text = element_text(size=14), legend.text=element_text(size=14))

life.pcoa2<-
 
  ggplot(life.coords, aes(MDS1, MDS2,color=LifeStage))+  
  geom_point(size=4)+
  theme_bw()+
  #geom_mark_circle(aes(color = LifeStage))+
  xlab("PC1-22.6%")+
  theme(legend.position = c(0.8, 0.2))+
  ylab("PC2-13.6%")+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  scale_color_manual(values=c("#FBF1A1", "#80CFBA", "#B0D9F1", "#7F7F7F"))+
  theme(text = element_text(size=14),
        axis.text = element_text(size=14), legend.text=element_text(size=14))

#PERMANOVA
larv.dis<-vegdist(larv_tbl, method='bray')
larv.dis<-as.data.frame(as.matrix(larv.dis))
larv.dis$SampleID<-row.names(larv.dis)
larv.dis<-merge(larv.dis, larv_meta[,c(1,6)], by= 'SampleID', all.y=F)
larv.dis<-larv.dis[,-1]
larv.permanova<-adonis(larv.dis[,-16] ~ larv.dis$Treatment2, permutations = 10000)
larv.permanova$aov.tab
#treatment is significant, R2=.33

life.dis<-vegdist(life_stg_tbl, method='bray')
life.dis<-as.data.frame(as.matrix(life.dis))
life.dis$SampleID<-row.names(life.dis)
life.dis<-merge(life.dis, life_meta[,c(1,4)], by='SampleID', all.y=F)
life.dis<-life.dis[,-1]
life.permanova<-adonis(life.dis[,-44] ~ life.dis$LifeStage, permutations=10000)
life.permanova$aov.tab
#life stages are signficant, R2=0.61

#####################################################
#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(larv_tbl))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, larv_meta, by='SampleID')
larv.alph$Richness<-as.numeric(larv.alph$`specnumber(larv_tbl)`)
pairwise.t.test(larv.alph$Richness, larv.alph$Treatment2, p.adjust.method = 'hochberg')
#        Bacillus Bsal  
#Bsal    0.3938   -     
#Control 0.0122   0.0037

life.alpha<-as.data.frame(specnumber(life_stg_tbl))
life.alpha$SampleID<-row.names(life.alpha)
life.alpha<-merge(life.alpha, life_meta, by='SampleID')
pairwise.t.test(life.alpha$`specnumber(life_stg_tbl)`, life.alpha$LifeStage, p.adjust.method = 'hochberg')
#no significant differences in richness

juv.alpha<-as.data.frame(specnumber(juv_tbl))
juv.alpha$SampleID<-row.names(juv.alpha)
juv.alpha<-merge(juv.alpha, juv_meta, by='SampleID')

#plot it
life.alpha2<-ggplot(life.alpha, aes(LifeStage, `specnumber(life_stg_tbl)`, fill=LifeStage))+
  geom_boxplot()+
  scale_fill_manual(values=c("#FBF1A1", "#80CFBA", "#B0D9F1", "#7F7F7F"))+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("")+
  ylab("ITS sOTU Richness")

larv.alpha2<-ggplot(larv.alph, aes(Treatment2, `specnumber(larv_tbl)`, fill=Treatment2))+
  geom_boxplot()+
  scale_fill_manual(values=c("#ABDAF4", "#F4CF7F", "#7F7F7F"))+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("")+
  ylab("ITS OTU Richness")

juv.alpha.fig<-ggplot(juv.alpha, aes(as.factor(SampleDay), `specnumber(juv_tbl)`, fill=Treatment))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c("#9991C3", "#D1CFA2", "#88BB9A", "#7F7F7F", "#B3CDE6", "#D7A2CE", "#E6B3BA", "#B3CDE6", "#B18980"))+
  facet_wrap(~Treatment, ncol=9)+
  ylab("ITS OTU Richness")+
  xlab("")+
  coord_cartesian(ylim=c(0,50))+
  theme(legend.position = "none")


#######################Analyze taxonomy

#read in taxonomy 
tax<-read.delim("spotted_salamander_microbiome/taxonomy.tsv")

#read in metadata for each experiment
juv_meta<-read.delim('spotted_salamander_microbiome/juvenile_exp_meta.txt', header=T)
larv_meta<-read.delim('spotted_salamander_microbiome/larval_exp_meta.txt', header=T)
life_meta<-read.delim('spotted_salamander_microbiome/life_stages_metadata.txt', header=T)

#read in asv table
asv.tbl<-read.delim("spotted_salamander_microbiome/spot_sal_asv_table.txt", row.names=1, header=T)
dim(asv.tbl)
#2106  234

#remove low sequencing depth samples
asv.tbl<-asv.tbl[,-which(colSums(asv.tbl)<400)]
dim(asv.tbl)
#2106  230

#convert from abundance to relative abundance
asv.tbl<-sweep(asv.tbl, 2, colSums(asv.tbl), '/')

#create asv table for life stages experiment
life_stg_tbl<-asv.tbl[,names(asv.tbl) %in% life_meta$SampleID]
colSums(life_stg_tbl)

#create asv table for larval experiment
larv_tbl<-asv.tbl[,names(asv.tbl) %in% larv_meta$SampleID]
colSums(larv_tbl)

#create asv table for juvenile experiment
juv_tbl<-asv.tbl[,names(asv.tbl) %in% juv_meta$SampleID]
colSums(juv_tbl)

#add taxonomy to asv table
larv_tbl$Feature.ID<-row.names(larv_tbl)
larv_tbl<-merge(larv_tbl, tax[,-3], by='Feature.ID', all.x=T)

juv_tbl$Feature.ID<-row.names(juv_tbl)
juv_tbl<-merge(juv_tbl, tax[,-3], by='Feature.ID')

life_stg_tbl$Feature.ID<-row.names(life_stg_tbl)
life_stg_tbl<-merge(life_stg_tbl, tax[,-3], by='Feature.ID')

#split taxonomy into different groups
juv_split<-separate(juv_tbl, Taxon, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
larv_split<-separate(larv_tbl, Taxon, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
life_split<-separate(life_stg_tbl, Taxon, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#shape the dataframe
juv_melt<-melt(juv_split)
larv_melt<-melt(larv_split)
life_melt<-melt(life_split)

#calculate the relative abundace for each class for each sample
juv_sum<-ddply(juv_melt, c("variable", 'Class'), summarize, rel_abun=sum(value))
larv_sum<-ddply(larv_melt, c("variable", 'Class'), summarize, rel_abun=sum(value))
life_sum<-ddply(life_melt, c("variable", 'Class'), summarize, rel_abun=sum(value))

#add metadata
life_sum<-merge(life_sum, life_meta, by.x = 'variable', by.y='SampleID')
larv_sum<-merge(larv_sum, larv_meta, by.x = 'variable', by.y='SampleID')
juv_sum<-merge(juv_sum, juv_meta, by.x = 'variable', by.y='SampleID')

#summarize only for life stage data
life_stage_sum<-ddply(life_sum, c("LifeStage", "Class"), summarize, rel_abun=mean(rel_abun))
life_cast<-dcast(life_stage_sum, Class ~ LifeStage)
life_cast$sum<-rowSums(life_cast[,-1])
life_cast<-life_cast[order(life_cast$sum, decreasing = T),]
life_cast<-life_cast[1:15,]
life_cast<-as.data.frame(life_cast[,-grep('sum', names(life_cast))])

#get 'others' category (things not in top20)
others<-as.data.frame(t(1-colSums(life_cast[,-1])))
others$tax<-"Others"
others<-others %>% select(tax, everything())
names(others)<-names(life_cast)
life_cast<-rbind(life_cast, others)
write.table(life_cast, 'life_cast.txt', row.names = F, quote=F, sep='\t')

#reread in fix taonomy
lifestages_class<-read.delim("spotted_salamander_microbiome/lifestages_class.txt", header=T)
lifestages_melt<-melt(lifestages_class)

#read in the best pallette
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

#plot it
life_taxonomy<-
  
  ggplot(lifestages_melt, aes(variable, value, fill=Class))+
  geom_bar(stat='identity')+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=pal)+
  guides(fill=guide_legend(ncol=1))+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Relative Abundance")+
  theme_bw()+
  theme(text = element_text(size=14))

###Analyze juvenile experimental data for beta diversity (alpha diversity above)
#read in metadata for each experiment
juv_meta<-read.delim('spotted_salamander_microbiome/juvenile_exp_meta.txt', header=T)

#read in asv table
asv.tbl<-read.delim("spotted_salamander_microbiome/spot_sal_asv_table.txt", row.names=1, header=T)
dim(asv.tbl)
#2106  241

#remove low sequencing depth samples
asv.tbl<-asv.tbl[,-which(colSums(asv.tbl)<400)]
dim(asv.tbl)
#2106  225

#rarefy table
set.seed(555)
asv.rare<-rrarefy(t(asv.tbl), sample = 400)

#subset mapping file by treatment
juv_meta_split<-split(juv_meta, juv_meta$Treatment)

#use subset maps to subset asv table into each experiment
agitated_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$Agitated$SampleID,]
bac_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$Bacillus$SampleID,]
chry_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$Chryseobacterium$SampleID,]
cont_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$Control$SampleID,]
highbs_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$HighBsal$SampleID,]
highrep_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$Highrepeated$SampleID,]
highthy_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$HighThyroid$SampleID,]
lowbs_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$LowBsal$SampleID,]
pene_asv<-asv.rare[row.names(asv.rare) %in% juv_meta_split$Penicillium$SampleID,]

#calculate beta diversity
agi.pcoa<-capscale(agitated_asv~1, method='bray')
bac.pcoa<-capscale(bac_asv~1, method='bray')
chry.pcoa<-capscale(chry_asv~1, method='bray')
cont.pcoa<-capscale(cont_asv~1, method='bray')
highbs.pcoa<-capscale(highbs_asv~1, method='bray')
highrep.pcoa<-capscale(highrep_asv~1, method='bray')
highthy.pcoa<-capscale(highthy_asv~1, method='bray')
lowbs.pcoa<-capscale(lowbs_asv~1, method='bray')
pene.pcoa<-capscale(pene_asv~1, method='bray')


#calculate percent variation explained
100*round(agi.pcoa$CA$eig[1]/sum(agi.pcoa$CA$eig), 3)
#33
100*round(agi.pcoa$CA$eig[2]/sum(agi.pcoa$CA$eig), 3)
#20

100*round(chry.pcoa$CA$eig[1]/sum(chry.pcoa$CA$eig), 3)
#34.6
100*round(chry.pcoa$CA$eig[2]/sum(chry.pcoa$CA$eig), 3)
#20.2

100*round(bac.pcoa$CA$eig[1]/sum(bac.pcoa$CA$eig), 3)
#24.9
100*round(bac.pcoa$CA$eig[2]/sum(bac.pcoa$CA$eig), 3)
#21

100*round(cont.pcoa$CA$eig[1]/sum(cont.pcoa$CA$eig), 3)
#30.9
100*round(cont.pcoa$CA$eig[2]/sum(cont.pcoa$CA$eig), 3)
#17.5

100*round(highbs.pcoa$CA$eig[1]/sum(highbs.pcoa$CA$eig), 3)
#27.7
100*round(highbs.pcoa$CA$eig[2]/sum(highbs.pcoa$CA$eig), 3)
#19.5

100*round(highrep.pcoa$CA$eig[1]/sum(highrep.pcoa$CA$eig), 3)
#32.4
100*round(highrep.pcoa$CA$eig[2]/sum(highrep.pcoa$CA$eig), 3)
#17.5

100*round(highthy.pcoa$CA$eig[1]/sum(highthy.pcoa$CA$eig), 3)
#25.4
100*round(highthy.pcoa$CA$eig[2]/sum(highthy.pcoa$CA$eig), 3)
#20.9

100*round(lowbs.pcoa$CA$eig[1]/sum(lowbs.pcoa$CA$eig), 3)
#47.6
100*round(lowbs.pcoa$CA$eig[2]/sum(lowbs.pcoa$CA$eig), 3)
#16.8

100*round(pene.pcoa$CA$eig[1]/sum(pene.pcoa$CA$eig), 3)
#42.6
100*round(pene.pcoa$CA$eig[2]/sum(pene.pcoa$CA$eig), 3)
#15.4


#extract the coordinates
agi.scores<-scores(agi.pcoa)
agi.coords<-as.data.frame(agi.scores$sites)
agi.coords$SampleID<-rownames(agi.coords)
agi.coords<-merge(agi.coords, juv_meta, by=c('SampleID'), all.y=F)

bac.scores<-scores(bac.pcoa)
bac.coords<-as.data.frame(bac.scores$sites)
bac.coords$SampleID<-rownames(bac.coords)
bac.coords<-merge(bac.coords, juv_meta, by=c('SampleID'), all.y=F)

chry.scores<-scores(chry.pcoa)
chry.coords<-as.data.frame(chry.scores$sites)
chry.coords$SampleID<-rownames(chry.coords)
chry.coords<-merge(chry.coords, juv_meta, by=c('SampleID'), all.y=F)

cont.scores<-scores(cont.pcoa)
cont.coords<-as.data.frame(cont.scores$sites)
cont.coords$SampleID<-rownames(cont.coords)
cont.coords<-merge(cont.coords, juv_meta, by=c('SampleID'), all.y=F)

highbs.scores<-scores(highbs.pcoa)
highbs.coords<-as.data.frame(highbs.scores$sites)
highbs.coords$SampleID<-rownames(highbs.coords)
highbs.coords<-merge(highbs.coords, juv_meta, by=c('SampleID'), all.y=F)

highrep.scores<-scores(highrep.pcoa)
highrep.coords<-as.data.frame(highrep.scores$sites)
highrep.coords$SampleID<-rownames(highrep.coords)
highrep.coords<-merge(highrep.coords, juv_meta, by=c('SampleID'), all.y=F)

highthy.scores<-scores(highthy.pcoa)
highthy.coords<-as.data.frame(highthy.scores$sites)
highthy.coords$SampleID<-rownames(highthy.coords)
highthy.coords<-merge(highthy.coords, juv_meta, by=c('SampleID'), all.y=F)

lowbs.scores<-scores(lowbs.pcoa)
lowbs.coords<-as.data.frame(lowbs.scores$sites)
lowbs.coords$SampleID<-rownames(lowbs.coords)
lowbs.coords<-merge(lowbs.coords, juv_meta, by=c('SampleID'), all.y=F)

pene.scores<-scores(pene.pcoa)
pene.coords<-as.data.frame(pene.scores$sites)
pene.coords$SampleID<-rownames(pene.coords)
pene.coords<-merge(pene.coords, juv_meta, by=c('SampleID'), all.y=F)


#plot it
all_coords<-rbind(agi.coords, bac.coords, chry.coords, cont.coords, highbs.coords, highrep.coords, lowbs.coords, highthy.coords, pene.coords)
all_coords$SampleDay<-as.factor(all_coords$SampleDay)

juv.coords<-
  
  ggplot(all_coords, aes(MDS1, MDS2, color=Treatment, shape=SampleDay))+  
  geom_point(aes(size=2))+
  theme_bw()+
  scale_color_manual(values=c("#9991C3", "#D1CFA2", "#88BB9A", "#7F7F7F", "#B3CDE6", "#D7A2CE", "#E6B3BA", "#B3CDE6", "#B18980"))+
  xlab("PC1")+
  guides(fill='none', color='none', size='none')+
  #theme(legend.position = "none")+
 #geom_mark_circle(aes(color = as.factor(SampleDay)))+
  ylab("PC2")+
  theme(legend.position = c(0.05, 0.15))+
  #theme(text = element_text(size=14),axis.text = element_text(size=14), legend.text=element_text(size=14))+
  facet_wrap(~Treatment, nrow = 1)
  
  
#make manuscript plots for larval/life stages data
ggarrange(juv.alpha.fig, juv.coords, ncol=1, labels=c("A", "B"), widths=c(1,10))
ggarrange(life.alpha2, life.pcoa2, life_taxonomy, labels=c("A", "B", "C"), ncol=3, nrow=1, widths = c(1.5, 1.5, 2))
ggarrange(larv.alpha2, larv.pcoa2, ncol=2, labels=c("A", "B"))


