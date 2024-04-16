cort_data<-read.delim("spotted_salamander_microbiome/cort_data.txt", header=T)
library(ggplot2)


cort_1<-ggplot(cort_data, aes(Treatment, cort_data$CORTpg.g.h))+
  geom_boxplot()+
  theme_bw()+
  ggtitle('(A)- 1 hour post exposure')+
  ylab("Corticosterone release rate (pg g-1 hr-1)")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

cort_data2<-read.delim("spotted_salamander_microbiome/cort_data2.txt", header=T)
cort_2<-
  ggplot(cort_data2, aes(Treatment, cort_data2$CORT.pg.g.h))+
  geom_boxplot()+
  theme_bw()+
  ggtitle('(B)- 30 days post exposure')+
  ylab("Corticosterone release rate (pg g-1 hr-1)")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

library(ggpubr)
ggarrange(cort_1, cort_2)

##Proportional change in mass
mass<-read.delim("spotted_salamander_microbiome/juv_mass.txt")

ggplot(mass, aes(Treatment, Prop_change_mass))+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("Proportional change in mass")

bartlett.test(mass$Prop_change_mass ~ mass$Treatment)
#Bartlett's K-squared = 22.301, df = 8, p-value = 0.004388
aov(mass$Prop_change_mass ~ mass$Treatment)
#                 mass$Treatment Residuals
#Sum of Squares        0.252981  4.653210
#Deg. of Freedom              8        86

TukeyHSD(aov(mass$Prop_change_mass ~ mass$Treatment))
#NO MATH DIFF