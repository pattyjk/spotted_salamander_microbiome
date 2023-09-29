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
asv.rare<-rrarefy(t(asv.tbl), sample = 400)
dim(asv.tbl)

#create asv table for life stages experiment
life_stg_tbl<-asv.rare[row.names(asv.rare) %in% life_meta$SampleID, ]

#create asv table for larval experiment
larv_tbl<-asv.rare[row.names(asv.rare) %in% larv_meta$SampleID, ]

#create asv table for juvenile experiment
juv_tbl<-asv.rare[row.names(asv.rare) %in% juv_meta$SampleID, ]
View(juv_tbl)

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