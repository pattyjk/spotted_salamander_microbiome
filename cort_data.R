cort_data<-read.delim("spotted_salamander_microbiome/cort_data.txt", header=T)
library(ggplot2)

ggplot(cort_data, aes(Treatment, cort_data$CORTpg.g.h))+
  geom_boxplot()+
  theme_bw()+
  ylab("Corticosterone release rate (pg g-1 hr-1)")+
  xlab("")

ggplot(cort_data, aes(Treatment, cort_data$Cort2))+
  geom_boxplot()

ggplot(cort_data, aes(Treatment, cort_data$Proportional.Change))+
  geom_boxplot()


##Proportional change in mass
mass<-read.delim("spotted_salamander_microbiome/juv_mass.txt")

ggplot(mass, aes(Treatment, Prop_change_mass))+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  ylab("Proportional change in mass")

bartlett.test(mass$Prop_change_mass ~ mass$Treatment)
aov(mass$Prop_change_mass ~ mass$Treatment)
#                 mass$Treatment Residuals
#Sum of Squares        0.252981  4.653210
#Deg. of Freedom              8        86

TukeyHSD(aov(mass$Prop_change_mass ~ mass$Treatment))
#NO MATH DIFF