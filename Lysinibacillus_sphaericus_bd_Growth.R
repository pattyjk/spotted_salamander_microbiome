####Bd/Bsal growth when exposed to  Lysinibacillus sphaericus metabolites

#read in data
bd_growth<-read.delim("./spotted_salamander_microbiome/lysin_bd_assay.txt", header=T)

#filter data
bd_growth<-bd_growth[-which(bd_growth$Keep == "N"),]

#get averages
library(dplyr)
bd_av<-ddply(bd_growth, c("Chytrid", "Conc"), summarize, mean=mean(Growth), sd=sd(Growth), n=length(Growth), se=sd/n)

#plot data
library(ggplot2)

ggplot(bd_av, aes(Conc, mean, fill=Conc))+
  facet_wrap(~Chytrid)+
  geom_bar(stat='identity')+
  ylab("Average growth %")+
  xlab("Lysinibacillus sphaericus Metabolite Concentration %")+
  guides(fill=F)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  scale_fill_manual(values=c("#55B4E9", "#55B4E9","#55B4E9","#55B4E9","#55B4E9","#55B4E9","#55B4E9","#55B4E9", "#E69F00", "#E69F00"))
  

Pos: #E69F00
  Others: #55B4E9