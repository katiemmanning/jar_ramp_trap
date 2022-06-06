##jar ramp and pitfall trap data analysis

#bring in order data sets from github

pitfall_order <- read.csv("https://raw.githubusercontent.com/katiemmanning/jar_ramp_trap/main/Insect%20ID%202020_pitfall_order.csv",na.strings = NULL)

jar_order <- read.csv("https://raw.githubusercontent.com/katiemmanning/jar_ramp_trap/main/Insect%20ID%202020_jarramp_order.csv",na.strings = NULL)

#taxa_order <- read.csv("")

#add trap type as a column on each data file
pitfall_order$Trap="pitfall"
jar_order$Trap="jar"

#combine order data tables 
library (plyr)
total_order <- rbind.fill (pitfall_order, jar_order)

#############
#NMDS of insect community by order between trap types
library (vegan)

#Create matrix of environmental variables
env.matrix_order<-total_order[c(1:3,16)]
#create matrix of community variables
com.matrix_order<-total_order[c(4:15)]

#ordination by NMDS
NMDS_order<-metaMDS(com.matrix_order, distance="bray", k=2, autotransform=FALSE, trymax=100)
NMDS_order
stressplot(NMDS_order)
#stress=0.19

#order NMDS visualization 

#what taxa to display using "taxa"
flying_order<-as.vector(t(taxa_order[1,]))
flying_order<-flying_order[-1]
crawling_order<-as.vector(t(taxa_order[2,]))
crawling_order<-crawling_order[-1]
include_order<-as.vector(t(taxa_order[3,]))
include_order<-include_order[-1]

#plot order NMDS
plot(NMDS_order, disp='sites', type="n")
title(main="Order", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "pitfall")
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "jar") 
#display ground trap data as solid shapes - pitfall=circle, ramp trap=square, jar=triangle, flying trap as triangle outline
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="pitfall"),pch=19, col="#E69F00")
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="jar"), pch=17, col="#009E73")
#add legend
legend(1.13,1.43, title=NULL, pch=c(19,17), col=c("#E69F00","#009E73"), cex=.7, legend=c("Pitfall", "Jar ramp"))
#add insect taxa as text
ordilabel(NMDS_order, display="species", select =which (include_order==TRUE & crawling_order == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS_order, display="species", select =which (include_order==TRUE & flying_order == TRUE), cex=0.6, col="white", fill="black")

#bootstrapping and testing for differences between the groups (traps)
fit<-adonis(com.matrix_order ~ Trap, data = env.matrix_order, permutations = 999, method="bray")
fit
#P-value = 0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_order)
anova(betadisper(distances_data, env.matrix_order$Trap))
#P-value = 0.1971 -- assume homogeneity of multivariate dispersion

#####
#calculate order Abundance
insects.abun_order <- rowSums(total_order[,4:15])
total_order$abundance <- insects.abun_order

#calculate order Richness
insects.rowsums_order <- rowSums(total_order[,4:15]>0)
total_order$richness <- insects.rowsums_order

#calculate order Shannon diversity
diversity_order <-diversity(total_order[,4:15])
total_order$diversity <-diversity_order

#calculate order Evenness
evenness_order <-diversity_order/log(specnumber(total_order[,4:15]))
total_order$evenness <- evenness_order


#######
#Mixed effects models
library(lme4)
library(lmerTest) #to obtain p values
library (emmeans) #for pairwise comparisons
library (multcompView) #to view letters

#order richness
##AIC 294
richness.model_order<-lmer(richness ~ Trap + Date + (1 | Site), data=total_order)
summary(richness.model_order)
anova(richness.model_order)
AIC(richness.model_order)
#pairwise comparison 
rich.emm_order<-emmeans(richness.model_order,pairwise~Trap)
rich.emm_order
#results: p=0.02 
rich.cld_order<-multcomp::cld(rich.emm_order, alpha = 0.05, Letters = LETTERS)
rich.cld_order

#order abundance
##AIC 613
abundance.model_order<-lmer(abundance ~ Trap + Date + (1 | Site), data=total_order)
summary(abundance.model_order)
anova(abundance.model_order)
AIC(abundance.model_order)
#pairwise comparison 
abun.emm_order<-emmeans(abundance.model_order,pairwise~Trap)
abun.emm_order
#results: p=0.0001
abun.cld_order<-multcomp::cld(abun.emm_order, alpha = 0.05, Letters = LETTERS)
abun.cld_order

#order diversity
##AIC 94
diversity.model_order<-lmer(diversity ~ Trap + Date + (1 | Site), data=total_order)
summary(diversity.model_order)
anova(diversity.model_order)
AIC(diversity.model_order)
#pairwise comparison 
div.emm_order<-emmeans(diversity.model_order,pairwise~Trap)
div.emm_order
#results: p=0.16
div.cld_order<-multcomp::cld(div.emm_order, alpha = 0.05, Letters = LETTERS)
div.cld_order

#order evenness
##AIC -58
evenness.model_order<-lmer(evenness ~ Trap + Date + (1 | Site), data=total_order)
summary(evenness.model_order)
anova(evenness.model_order)
AIC(evenness.model_order)
#pairwise comparison 
even.emm_order<-emmeans(evenness.model_order,pairwise~Trap)
even.emm_order
#results: p=0.056
even.cld_order<-multcomp::cld(even.emm_order, alpha = 0.05, Letters = LETTERS)
even.cld_order

#looking at linear models

richness.model_order<-lm(richness ~ Trap, data=total_order)
summary(richness.model_order)
anova(richness.model_order)
AIC(richness.model_order)

evenness.model_order<-lm(evenness ~ Trap, data=total_order)
summary(evenness.model_order)
anova(evenness.model_order)
AIC(evenness.model_order)


###
library(ggplot2)
#order abundance plot
abundance.plot_order<-ggplot(total_order, aes(x =Trap, y = abundance, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Abundance (log10)")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  geom_text(data=abun.cld_order, aes(y = 600, label = .group))
abundance.plot_order

#order richness plot
richness.plot_order<-ggplot(total_order, aes(x =Trap, y = richness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Richness")+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  geom_text(data=rich.cld_order, aes(y = 25, label = .group))
richness.plot_order

#order diversity plot
diversity.plot_order<-ggplot(total_order, aes(x =Trap, y = diversity, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Shannon Diversity")+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  geom_text(data=div.cld_order, aes(y = 2.5, label = .group))
diversity.plot_order

#order evenness plot
evenness.plot_order<-ggplot(total_order, aes(x =Trap, y = evenness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Evenness")+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  geom_text(data=even.cld_order, aes(y = 1.2, label = .group))
evenness.plot_order

#Mush order plots together
library(ggpubr) 
orderfigure <- ggarrange(richness.plot_order, abundance.plot_order, diversity.plot_order, evenness.plot_order,
                         labels = c("A", "B", "C", "D"),
                         ncol = 2, nrow = 2,
                         common.legend = TRUE, legend = "bottom")
pdf("order.pdf", height=6, width=8) #height and width in inches
orderfigure
dev.off()



#############################

#bring in functional data sets from github

pitfall <- read.csv("https://raw.githubusercontent.com/katiemmanning/jar_ramp_trap/main/Insect%20ID%202020_pitfall_functional.csv",na.strings = NULL)

jar <- read.csv("https://raw.githubusercontent.com/katiemmanning/jar_ramp_trap/main/Insect%20ID%202020_jarramp_functional.csv",na.strings = NULL)

#taxa <- read.csv("")

#add trap type as a column on each data file
pitfall$Trap="pitfall"
jar$Trap="jar"

#combine data tables 
library (plyr)
total <- rbind.fill (pitfall, jar)

#############
#NMDS of insect community by functional classification between trap types
library (vegan)

#Create matrix of environmental variables
env.matrix<-total[c(1:3,43)]
#create matrix of community variables
com.matrix<-total[c(4:42)]

#ordination by NMDS
NMDS<-metaMDS(com.matrix, distance="bray", k=2, autotransform=FALSE, trymax=100)
NMDS
stressplot(NMDS)
#stress=0.21

#functional classification NMDS visualization 

#what taxa to display using "taxa"
flying_func<-as.vector(t(taxa[1,]))
flying_func<-flying_func[-1]
crawling_func<-as.vector(t(taxa[2,]))
crawling_func<-crawling_func[-1]
include_func<-as.vector(t(taxa[3,]))
include_func<-include_func[-1]

#plot functional NMDS
plot(NMDS, disp='sites', type="n")
title(main="Functional", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "pitfall")
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "jar") 
#display ground trap data as solid shapes - pitfall=circle, ramp trap=square, jar=triangle, flying trap as triangle outline
points(NMDS, display="sites", select=which(env.matrix$Trap=="pitfall"),pch=19, col="#E69F00")
points(NMDS, display="sites", select=which(env.matrix$Trap=="jar"), pch=17, col="#009E73")
#add legend
legend(1.0,1.51, title=NULL, pch=c(19,17), col=c("#E69F00","#009E73"), cex=.7, legend=c("Pitfall", "Jar ramp"))

#add insect taxa as text
ordilabel(NMDS, display="species", select =which (include_func==TRUE & crawling_func == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS, display="species", select =which (include_func==TRUE & flying_func == TRUE), cex=0.6, col="white", fill="black")

#bootstrapping and testing for differences between the groups (traps)
fit<-adonis(com.matrix ~ Trap, data = env.matrix, permutations = 999, method="bray")
fit
#P-value = 0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$Trap))
#P-value = .0001 -- cannot assume homogeneity of multivariate dispersion


################
#calculate Abundance
insects.abun <- rowSums(insects[,4:42])
insects$abundance <- insects.abun

#calculate Richness
insects.rowsums <- rowSums(insects[,4:42]>0)
insects$richness <- insects.rowsums

#calculate Shannon diversity
diversity <-diversity(insects[,4:42])
insects$diversity <-diversity

#calculate Evenness
evenness <-diversity/log(specnumber(insects[,4:42]))
insects$evenness <- evenness

#######
#Mixed effects models
library(lme4)
library(lmerTest) #to obtain p values
library (emmeans) #for pairwise comparisons
library (multcompView) #to view letters

#richness
##AIC 718
richness.model<-lmer(richness ~ Trap + Date + (1 | Site), data=insects)
summary(richness.model)
anova(richness.model)
AIC(richness.model)
#pairwise comparison 
rich.emm<-emmeans(richness.model,pairwise~Trap)
rich.emm
#results: jar-pitfall no sig diff (0.0610), sig dif btw all others
rich.cld<-multcomp::cld(rich.emm, alpha = 0.05, Letters = LETTERS)
rich.cld

#abundance
##AIC 1795
abundance.model<-lmer(abundance ~ Trap + Date + (1 | Site), data=insects)
summary(abundance.model)
anova(abundance.model)
AIC(abundance.model)
#pairwise comparison 
abun.emm<-emmeans(abundance.model,pairwise~Trap)
abun.emm
#results: jar-pitfall no sig diff (0.8089), sig dif btw all others
abun.cld<-multcomp::cld(abun.emm, alpha = 0.05, Letters = LETTERS)
abun.cld

#diversity
##AIC 175 (152 w/o Date)
diversity.model<-lmer(diversity ~ Trap + Date + (1 | Site), data=insects)
summary(diversity.model)
anova(diversity.model)
AIC(diversity.model)
#pairwise comparison 
div.emm<-emmeans(diversity.model,pairwise~Trap)
div.emm
#results: no sig diff btw jar-pitfall (0.2016), jar-sticky (0.9540), pitfall-sticky (0.0661); sig diff btw all others 
div.cld<-multcomp::cld(div.emm, alpha = 0.05, Letters = LETTERS)
div.cld

#evenness
##AIC -172 (-193 w/o Date)
evenness.model<-lmer(evenness ~ Trap + Date + (1 | Site), data=insects)
summary(evenness.model)
anova(evenness.model)
AIC(evenness.model)
#pairwise comparison 
even.emm<-emmeans(evenness.model,pairwise~Trap)
even.emm
#results: no sig diff btw jar-pitfall (0.2851) or ramp-sticky (0.0974); sig diff btw all others
even.cld<-multcomp::cld(even.emm, alpha = 0.05, Letters = LETTERS)
even.cld

###########
library(ggplot2)
#abundance plot
abundance.plot<-ggplot(insects, aes(x =Trap, y = abundance, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=abun.cld, aes(y = 600, label = .group))
abundance.plot

#richness plot
richness.plot<-ggplot(insects, aes(x =Trap, y = richness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=rich.cld, aes(y = 25, label = .group))
richness.plot

#diversity plot
diversity.plot<-ggplot(insects, aes(x =Trap, y = diversity, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=div.cld, aes(y = 2.5, label = .group))
diversity.plot

#evenness plot
evenness.plot<-ggplot(insects, aes(x =Trap, y = evenness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=even.cld, aes(y = 1.2, label = .group))
evenness.plot

#Mush order plots together
library(ggpubr) 
functionalfigure <- ggarrange(richness.plot, abundance.plot, diversity.plot, evenness.plot,
                              labels = c("E", "F", "G", "H"),
                              ncol = 2, nrow = 2,
                              common.legend = TRUE, legend = "bottom")
pdf("functional.pdf", height=6, width=8) #height and width in inches
functionalfigure
dev.off()

