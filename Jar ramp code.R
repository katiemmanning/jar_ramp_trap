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
