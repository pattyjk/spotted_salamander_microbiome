library(vegan)
library(ggplot2)
library(reshape2)

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

#rarefy table
asv.rare<-rrarefy(t(asv.tbl), sample = 400)
dim(asv.tbl)



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