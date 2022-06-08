##jar ramp and pitfall trap data analysis

#bring in data sets from github

pitfall <- read.csv("https://raw.githubusercontent.com/katiemmanning/jar_ramp_trap/main/Insect%20ID%202020_pitfall_functional.csv",na.strings = NULL)

jar <- read.csv("https://raw.githubusercontent.com/katiemmanning/jar_ramp_trap/main/Insect%20ID%202020_jarramp_functional.csv",na.strings = NULL)

#taxa <- read.csv("")

#add trap type as a column on each data file
pitfall$Trap="pitfall"
jar$Trap="jar"

#calculate mean richness and abundance of each trap type
#insects.abun <- rowSums(pitfall[,4:42])
#pitfall$abundance <- insects.abun
#insects.rowsums <- rowSums(pitfall[,4:42]>0)
#pitfall$richness <- insects.rowsums

#insects.abun <- rowSums(jar[,4:42])
#jar$abundance <- insects.abun
#insects.rowsums <- rowSums(jar[,4:42]>0)
#jar$richness <- insects.rowsums

#mean(pitfall$abundance) #14.82
#mean(pitfall$richness) #5.28
#mean(jar$abundance) #26.17
#mean(jar$richness) #6.36

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
NMDS<-metaMDS(com.matrix, distance="bray", k=2, autotransform=TRUE, trymax=100)
NMDS
stressplot(NMDS)
#stress=0.26

#what taxa to display using "taxa"
#flying_func<-as.vector(t(taxa[1,]))
#flying_func<-flying_func[-1]
#crawling_func<-as.vector(t(taxa[2,]))
#crawling_func<-crawling_func[-1]
#include_func<-as.vector(t(taxa[3,]))
#include_func<-include_func[-1]

#plot functional NMDS
#8x11
plot(NMDS, disp='sites', type="n")
#title(main="Functional", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "pitfall")
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "jar") 
#display ground trap data as solid shapes - pitfall=circle, ramp trap=square, jar=triangle, flying trap as triangle outline
points(NMDS, display="sites", select=which(env.matrix$Trap=="pitfall"),pch=19, col="#E69F00")
points(NMDS, display="sites", select=which(env.matrix$Trap=="jar"), pch=17, col="#009E73")
#add legend
legend(1.123,1.06, title=NULL, pch=c(19,17), col=c("#E69F00","#009E73"), cex=.7, legend=c("Pitfall", "Jar ramp"))

#add insect taxa as text
#ordilabel(NMDS, display="species", select =which (include_func==TRUE & crawling_func == TRUE), cex=0.6, col="black", fill="white")
#ordilabel(NMDS, display="species", select =which (include_func==TRUE & flying_func == TRUE), cex=0.6, col="white", fill="black")

#bootstrapping and testing for differences between the groups (traps)
fit<-adonis(com.matrix ~ Trap, data = env.matrix, permutations = 999, method="bray")
fit
#P-value = 0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$Trap))
#P-value = .5233 -- assumes homogeneity of multivariate dispersion


################
#calculate Abundance
insects.abun <- rowSums(total[,4:42])
total$abundance <- insects.abun

#calculate Richness
insects.rowsums <- rowSums(total[,4:42]>0)
total$richness <- insects.rowsums

#calculate Shannon diversity
diversity <-diversity(total[,4:42])
total$diversity <-diversity

#calculate Evenness
evenness <-diversity/log(specnumber(total[,4:42]))
total$evenness <- evenness

#######
#Mixed effects models
library(lme4)
library(lmerTest) #to obtain p values
library (emmeans) #for pairwise comparisons
library (multcompView) #to view letters
library (wiqid) #for AICc

#richness
##AICc 340
richness.model<-lmer(richness ~ Trap + Date + (1 | Site), data=total)
##AICc 348
richness.model<-lmer(richness ~ Trap + (1 | Site), data=total)
##AICc 354
richness.model<-lm(richness ~ Trap, data=total)
##AICc 353
richness.model<-lm(richness ~ Trap + Date, data=total)
##AICc 343
richness.model<-lm(richness ~ Trap + Site, data=total)
summary(richness.model)
anova(richness.model)
AICc(richness.model)
#pairwise comparison 
rich.emm<-emmeans(richness.model,pairwise~Trap)
rich.emm
#results: 
rich.cld<-multcomp::cld(rich.emm, alpha = 0.05, Letters = LETTERS)
rich.cld

#abundance
##AICc 617
abundance.model<-lmer(abundance ~ Trap + Date + (1 | Site), data=total)
##AICc 657
abundance.model<-lmer(abundance ~ Trap + (1 | Site), data=total)
##AICc 661
abundance.model<-lm(abundance ~ Trap, data=total)
##AICc 649
abundance.model<-lm(abundance ~ Trap + Date, data=total)
##AICc 665
abundance.model<-lm(abundance ~ Trap + Site, data=total)
summary(abundance.model)
anova(abundance.model)
AICc(abundance.model)
#pairwise comparison 
abun.emm<-emmeans(abundance.model,pairwise~Trap)
abun.emm
#results: 
abun.cld<-multcomp::cld(abun.emm, alpha = 0.05, Letters = LETTERS)
abun.cld

#diversity
##AICc 118
diversity.model<-lmer(diversity ~ Trap + Date + (1 | Site), data=total)
##AICc 102
diversity.model<-lmer(diversity ~ Trap + (1 | Site), data=total)
##AICc 105
diversity.model<-lm(diversity ~ Trap, data=total)
##AICc 111
diversity.model<-lm(diversity ~ Trap + Date, data=total)
##AICc 90
diversity.model<-lm(diversity ~ Trap + Site, data=total)
summary(diversity.model)
anova(diversity.model)
AICc(diversity.model)
#pairwise comparison 
div.emm<-emmeans(diversity.model,pairwise~Trap)
div.emm
#results: 
div.cld<-multcomp::cld(div.emm, alpha = 0.05, Letters = LETTERS)
div.cld

#evenness
##AICc -46
evenness.model<-lmer(evenness ~ Trap + Date + (1 | Site), data=total)
##AICc -76
evenness.model<-lmer(evenness ~ Trap + (1 | Site), data=total)
##AICc -87
evenness.model<-lm(evenness ~ Trap, data=total)
##AICc -83
evenness.model<-lm(evenness ~ Trap + Date, data=total)
##AICc -91
evenness.model<-lm(evenness ~ Trap + Site, data=total)
summary(evenness.model)
anova(evenness.model)
AICc(evenness.model)
#pairwise comparison 
even.emm<-emmeans(evenness.model,pairwise~Trap)
even.emm
#results: 
even.cld<-multcomp::cld(even.emm, alpha = 0.05, Letters = LETTERS)
even.cld

#

##AICc 354
richness.model<-lm(richness ~ Trap, data=total)
summary(richness.model)
anova(richness.model)
AICc(richness.model)
#pairwise comparison 
rich.emm<-emmeans(richness.model,pairwise~Trap)
rich.emm
#results: p = 0.0238 -- difference in richness
rich.cld<-multcomp::cld(rich.emm, alpha = 0.05, Letters = LETTERS)
rich.cld

##AICc 661
abundance.model<-lm(abundance ~ Trap, data=total)
summary(abundance.model)
anova(abundance.model)
AICc(abundance.model)
#pairwise comparison 
abun.emm<-emmeans(abundance.model,pairwise~Trap)
abun.emm
#results: p = 0.0005 -- difference in abundance
abun.cld<-multcomp::cld(abun.emm, alpha = 0.05, Letters = LETTERS)
abun.cld

##AICc 105
diversity.model<-lm(diversity ~ Trap, data=total)
summary(diversity.model)
anova(diversity.model)
AICc(diversity.model)
#pairwise comparison 
div.emm<-emmeans(diversity.model,pairwise~Trap)
div.emm
#results: 0.1590 -- no difference in diversity
div.cld<-multcomp::cld(div.emm, alpha = 0.05, Letters = LETTERS)
div.cld

##AICc -87
evenness.model<-lm(evenness ~ Trap, data=total)
summary(evenness.model)
anova(evenness.model)
AICc(evenness.model)
#pairwise comparison 
even.emm<-emmeans(evenness.model,pairwise~Trap)
even.emm
#results: 0.1514 -- no difference in evenness
even.cld<-multcomp::cld(even.emm, alpha = 0.05, Letters = LETTERS)
even.cld


###########
library(ggplot2)
#abundance plot
abundance.plot<-ggplot(total, aes(x =Trap, y = abundance, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Abundance")+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  geom_text(data=abun.cld, aes(y = 80, label = .group))
abundance.plot

#richness plot
richness.plot<-ggplot(total, aes(x =Trap, y = richness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Richness")+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  geom_text(data=rich.cld, aes(y = 15, label = .group))
richness.plot

#diversity plot
diversity.plot<-ggplot(total, aes(x =Trap, y = diversity, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Shannon diversity")+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  geom_text(data=div.cld, aes(y = 2.5, label = .group))
diversity.plot

#evenness plot
evenness.plot<-ggplot(total, aes(x =Trap, y = evenness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Evenness")+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  geom_text(data=even.cld, aes(y = 1.1, label = .group))
evenness.plot

#Mush order plots together
library(ggpubr) 
metrics<- ggarrange(richness.plot, abundance.plot, diversity.plot, evenness.plot,
                              #labels = c("A", "B", "C", "D"),
                              ncol = 2, nrow = 2,
                              common.legend = TRUE, legend = "bottom")
pdf("boxplots.pdf", height=6, width=8) #height and width in inches
metrics
dev.off()
metrics

