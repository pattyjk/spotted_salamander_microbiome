library(vegan)
library(ggplot2)
library(reshape2)

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
juv.pcoa<-capscale(juv_tbl~1, method='jaccard')
life.pcoa<-capscale(life_stg_tbl~1, method='jaccard')
larv.pcoa<-capscale(larv_tbl~1, method='jaccard')

#calculate percent variation explained
100*round(juv.pcoa$CA$eig[1]/sum(juv.pcoa$CA$eig), 3)
#22.6
100*round(juv.pcoa$CA$eig[2]/sum(juv.pcoa$CA$eig), 3)
#12.8

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

juv.scores<-scores(juv.pcoa)
juv.coords<-as.data.frame(juv.scores$sites)
juv.coords$SampleID<-rownames(juv.coords)
juv.coords<-merge(juv.coords, juv_meta, by=c('SampleID'))


ggplot(juv.coords, aes(MDS1, MDS2, color=Treatment))+  
  geom_point(aes(size=2))+
  theme_bw()+
  xlab("PC1-%")+
  ylab("PC2-%")+
  theme(text = element_text(size=14),
        axis.text = element_text(size=14), legend.text=element_text(size=14))

ggplot(larv.coords, aes(MDS1, MDS2, color=Treatment2))+  
  geom_point(aes(size=2))+
  theme_bw()+
  xlab("PC1-%")+
  ylab("PC2-%")+
  theme(text = element_text(size=14),
        axis.text = element_text(size=14), legend.text=element_text(size=14))

ggplot(life.coords, aes(MDS1, MDS2, color=LifeStage))+  
  geom_point(aes(size=2))+
  theme_bw()+
  xlab("PC1-23.9%")+
  ylab("PC2-13.6%")+
  theme(text = element_text(size=14),
        axis.text = element_text(size=14), legend.text=element_text(size=14))

#######################analyze taxonomy

#read in taxonomy 
tax<-read.delim("spotted_salamander_microbiome/taxonomy.tsv")

#read in metadata
meta<-read.delim("spotted_salamander_microbiome/ITS_FINAL_RUN_Mapping_File.txt", header=T)

#read in asv table
asv.tbl<-read.delim("spotted_salamander_microbiome/spot_sal_asv_table.txt", row.names=1, header=T)
dim(asv.tbl)
#2106  234

#remove low sequencing depth samples
asv.tbl<-asv.tbl[,-which(colSums(asv.tbl)<400)]
dim(asv.tbl)
#2106  223