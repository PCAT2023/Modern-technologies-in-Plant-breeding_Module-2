###############################################################################################
###                 Scientific Training Center for Plant Biotechnologies                    ###
###                                   REGULAR LEVEL                                         ###
###                                       ----------                                        ###
###                               Hands-on: QTL detection using Rqtl                        ###
###                                          ----                                           ###
###             Script edited by Pr Laurent Gentzbittel & Pr Cécile Ben, Skoltech           ###
###############################################################################################

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
#On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
#Choose the directory containing the datafile and the associated R script.

### II. Installation R packages needed for the analysis on RStudio:
#Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
#Comment #1: R package installation requires a connection to internet
#Comment #2: Once packages have been installed, no need to re-install them again when you close-open again RStudio.
library(qtl)#For QTL analysis
library(lattice) #For graphs


### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# use of the constraint 'set-to-zero' for ANOVAs ## will see later in this script
options(contrasts=c('contr.treatment','contr.poly'))
#can also use 'contr.sum' for a 'sum-to-zero' constraint

# postscript figures with colours
trellis.device(device = postscript,color = T)

###################
# Data import in R
###################

QTL_LR5RalstoDataRaw<-read.cross("csvs", genfile="LR5_Geno.csvs", phefile="RalstoLR5_Scoredpi.csvs", sep=";", dec=",")
#genfile: file containing genotypes, phefile:file containing phenotypes.

QTL_LR5RalstoDataRaw<- convert2riself(QTL_LR5RalstoDataRaw)
#The function convert2riself converts the type of crossing indicating that the analysis was made on RILs obtained by self-fertilization.
summary(QTL_LR5RalstoDataRaw)
str(QTL_LR5RalstoDataRaw)


#################
########	A. ANALYSIS OF THE QUALITY OF PHENOTYPIC AND GENOTYPIC DATA
#################

###
#A.1. Graphic visualization
###

x11()
plot(QTL_LR5RalstoDataRaw)
#Graph 1: Evaluation of the quality of genotypic data: visualization of missing data.
#Graph 2: Visualization of the genetic map
#Graph 3 to 15 (phe 1 to phe 13): Distribution of phenotypic data for the different traits evaluated within the population of RILs
#Warning: the graph 'phe1' does not represent a real phenotype but only the number of RILs. It therefore has no interest!

# Detail of the genetic map
x11()
plot.map(QTL_LR5RalstoDataRaw, show.marker.names=TRUE)

# Details of missing genotypic data
x11()
plotMissing(QTL_LR5RalstoDataRaw, reorder=TRUE)

###
# A.2. More detailed analysis of the quality of genotypes
###

# A.2.a. Calculation of recombination fractions
QTL_LR5RalstoDataRaw<-est.rf(QTL_LR5RalstoDataRaw)
x11()
plotRF(QTL_LR5RalstoDataRaw)
# Estimated recombination fractions for all pairs of markers, along with LOD scores for the test of r = 1/2. 
# The recombination fractions are in the upper left triangle; LOD scores are in the lower right triangle.
#Red indicates a large LOD score or a small recombination fraction, while blue is the opposite.

#Conclusion: This graph does not reveal any problem with the positioning of the markers on the chromosomes.

# A.2.b. Calculation of probabilities of genotyping errors
QTL_LR5RalstoDataRaw<-calc.errorlod(QTL_LR5RalstoDataRaw, error.prob=0.1)
#This function calculates, for each individual at each marker, a LOD score measuring the strength of evidence for a genotyping error, as described by Lincoln and Lander (1992)
#Evidence for a genotype being in error is considered assuming that all other genotypes (for that individual, on that chromosome) are correct
#If markers are extremely tightly linked,recombination events can give large error LOD scores. The error LOD scores should not be trusted blindly, but should be viewed as a tool for identifying genotypes deserving further study.

top.errorlod(QTL_LR5RalstoDataRaw)

x11()
plotGeno(QTL_LR5RalstoDataRaw)
#Observed marker genotype data for RIL population. White and black circles correspond to the genotypes AA, and BB (AB are ignored due to convert2riself), respectively. 
#Genotypes flagged by red squares would be indicated to be likely genotyping errors. Here there is none.
#Blue crosses indicate recombination events.

###
# A.2. More detailed analysis of phenotypes
###
# Prior to the QTL analysis, it is necessary to show by ANOVA that there is an effect of the "line" factor on the phenotypic trait analyzed,
# and to compute adjusted means corrected from the block or environmental effect.
# This part of the analysis of the raw phenotypic data is not processed in this hands-on tutorial
# It uses the concepts & methods that were discussed in the 'Advanced Experimental Plans' sessions.

Pheno_Score7dpi<-data.frame(indice=QTL_LR5RalstoDataRaw$pheno$T7)
Pheno_Score7dpi
str(Pheno_Score7dpi)
shapiro.test(Pheno_Score7dpi$indice)

Pheno_cfu<-data.frame(indice=QTL_LR5RalstoDataRaw$pheno$Moyenne.logcfu)
Pheno_cfu
str(Pheno_cfu)
shapiro.test(Pheno_cfu$indice)

#################
######## B. DETECTION OF QTL BY ANOVA: 'MANUAL' METHOD
#################

#For example, this analysis is carried out on 2 markers (MTE117 located on chr3 and MTE33 located on chr5)
#and a phenotypic trait: the disease score at 7 days after inoculation (T7)
#Genotype Code: 1=A, 2=B

ExManuel<-data.frame(marker1=QTL_LR5RalstoDataRaw$geno$LG3$data[,"MTE117"], marker2=QTL_LR5RalstoDataRaw$geno$LG5$data[,"MTE33"]
,indice=QTL_LR5RalstoDataRaw$pheno$T7)
ExManuel

x11()
par(mfrow=c(1,2))
boxplot(ExManuel$indice~ExManuel$marker1,main='MTE117', ylab='Disease score at 7dpi')
boxplot(ExManuel$indice~ExManuel$marker2,main='MTE33', ylab='Disease score at 7dpi')

#Analysis for the MTE117 marker located on chr3
summary(aov(indice ~ factor(marker1),data=ExManuel))
#Conclusion?

#Analysis for the MTE33 marker located on chr5
summary(aov(indice ~ factor(marker2),data=ExManuel))
#Conclusion?

#Graph of the effect of the genotype with the MTE33 marker on the disease index at 7 days after inoculation
x11()
plotPXG(QTL_LR5RalstoDataRaw, pheno.col=5,"MTE33")
# Empty circles represent raw data. The genotypes of individuals that were not typed is inferred based on the genotype at linked markers via a single imputation from sim.geno; these points are plotted in red.
# Also shown are the disease score means for each genotype as well as the confidence intervals of the mean.

#################
######## D. DETECTION DE QTL PAR SCANONE: mr method ie Marker regression or ANOVA at each marker (Single QTL genome scan)
#################
# the previous method is programmed for all markers

###
## Example 1: Analysis on the disease score at 7dpi
###

QTLT7_LR5Ralsto.mr<-scanone(QTL_LR5RalstoDataRaw, pheno.col=5,method="mr")
# pheno.col = 5 corresponds to the disease index at 7dpi (T7)
# method = "mr" indicates the method by marker regression

QTLT7_LR5Ralsto.mr
summary(QTLT7_LR5Ralsto.mr) # only the lod max on each chromosome

###
## Example 2: Analysis of the quantity of bacteria (in log (cfu)) present in the aerial part of the plant at the end of the experiment (cfu)
###

QTLcfu_LR5Ralsto.mr<-scanone(QTL_LR5RalstoDataRaw, pheno.col=2,method="mr") #cfu
QTLcfu_LR5Ralsto.mr
summary(QTLcfu_LR5Ralsto.mr) # only the lod max on each chromosome

#################
######## E. DETECTION OF QTL with SCANONE: em method ie Simple interval mapping by maximum likelihood (Single QTL genome scan)
#################

## E.1. For this method, it is essential to impute the missing genotypes.
QTL_LR5RalstoData<- calc.genoprob(QTL_LR5RalstoDataRaw, step=1, error.prob=0.01)
#The function calc.genoprob calculates conditional genotype probabilities given the multipoint marker data.Uses the hidden Markov model technology 
#to calculate the probabilities of the true underlying genotypes given the observed multipoint marker data, with possible allowance for genotyping errors.
#i.e. for an individual whose genotype is unknown to a marker, gives the probabilities of the genotypes to this marker according to the haplotype
# error.prob: probability of imputing a genotype A to the marker when it is B
#step: Maximum distance (in cM) between positions at which the genotype probabilities are calculated, 
#though for step = 0, probabilities are calculated only at the marker locations.

###
## E.2. Example 1: Analysis on the disease score at 7dpi
###

## QTL detection
QTLT7_LR5Ralsto.em<-scanone(QTL_LR5RalstoData, pheno.col=5,method="em") #T7
QTLT7_LR5Ralsto.em
summary(QTLT7_LR5Ralsto.em) # only the lod max on each chromosome

## Graphic visualization
x11()
plot(QTLT7_LR5Ralsto.em)

## Graphic visualization of the results of the genotype effect at the marker closest to the major QTL peak
max(QTLT7_LR5Ralsto.em)
find.marker(QTL_LR5RalstoData, chr="LG5", pos=29.0) #Goal: find the marker closest to the QTL's lod peak on ch5
x11()
effectplot(QTL_LR5RalstoData, mname1="LG5@29", pheno.col=5, main="LR5 Ralstonia - QTL LG5 effect")
x11()
plotPXG(QTL_LR5RalstoData, marker="MTE33", pheno.col=5, main="LR5 Ralstonia - QTL LG5 effect")

## Calculation of QTL confidence intervals (CI) and identification of flanking markers
# By default, decrease by 1.5 of the lod score otherwise use the argument drop = 1 for example.

lodint(QTLT7_LR5Ralsto.em, chr="LG5", expandtomarkers=TRUE)
lodint(QTLT7_LR5Ralsto.em, chr="LG5", expandtomarkers=TRUE, drop=1) #20cM between MTE32 and MTE34
# the CI is given in positions, not in markers:
ICfinT7<-lodint(QTLT7_LR5Ralsto.em, chr="LG5", drop=1)#12cM

## Graphic visualization listing all the information concerning the major QTL
x11()
plot(QTLT7_LR5Ralsto.em,chr="LG5",
     show.marker.names=T,
	xlab="Map Position on LG5 (cM)",
     ylab="LOD Score LR5Ralsto_ScoreT7")
abline(h=3,col="darkgreen",lty=2)
abline(v=c(ICfinT7$pos), col=c("blue", "red", "blue"), lty=c(2,1,2))

###
## Example 2: Analysis of the quantity of bacteria (in log (cfu)) present in the aerial part of the plant at the end of the experiment (cfu)
###

### QTL detection
QTLcfu_LR5Ralsto.em<-scanone(QTL_LR5RalstoData, pheno.col=2,method="em")#cfu
QTLcfu_LR5Ralsto.em
summary(QTLcfu_LR5Ralsto.em) # only the lod max on each chromosome

## Graphic visualization
x11()
plot(QTLcfu_LR5Ralsto.em)

## Graphic visualization of the results of the genotype effect at the marker closest to the major QTL peak
max(QTLcfu_LR5Ralsto.em)
find.marker(QTL_LR5RalstoData, chr="LG5", pos=23.0) #Goal: find the marker closest to the qtl lod peak on ch.5
x11()
effectplot(QTL_LR5RalstoData, mname1="LG5@23", pheno.col=2, main="LR5 Ralstonia - QTL LG5 effect")
x11()
plotPXG(QTL_LR5RalstoData, marker="MTE32", pheno.col=2, main="LR5 Ralstonia - QTL LG5 effect")

## Calculation of QTL confidence intervals and identification of flanking markers
# By default, decrease by 1.5 of the lod score otherwise use the argument drop = 1 for example

lodint(QTLcfu_LR5Ralsto.em, chr="LG5", expandtomarkers=TRUE)
lodint(QTLcfu_LR5Ralsto.em, chr="LG5", expandtomarkers=TRUE, drop=1) #42.2cM between MTE30 and MTE34
# the CI is given in positions, not in markers:
ICfincfu<-lodint(QTLcfu_LR5Ralsto.em, chr="LG5", drop=1) #19cM
ICfincfu

## Graphic visualization listing all the information concerning the major QTL
x11()
plot(QTLcfu_LR5Ralsto.em,chr="LG5",
     show.marker.names=T,
	xlab="Map Position on LG5 (cM)",
     ylab="LOD Score LR5Ralsto_cfu")
abline(h=3,col="darkgreen",lty=2)
abline(v=c(ICfincfu$pos), col=c("blue", "red", "blue"), lty=c(2,1,2))

###
### Major QTL IC comparisons on LG5 for cfu and at T7
###

x11()
plot(QTLcfu_LR5Ralsto.em,QTLT7_LR5Ralsto.em, chr="LG5",
     show.marker.names=T,
	xlab="Map Position on LG5 (cM)",
     ylab="LOD Score LR5Ralsto")
abline(h=3,col="darkgreen",lty=2)
abline(v=c(ICfincfu$pos), col=c("black"), lty=c(2,1,2))
abline(v=c(ICfinT7$pos), col=c("red"), lty=c(2,1,2))


#################
######## F. DETECTION OF QTL with SCANONE: imp method ie imputation method or anova after imputation of missing genotypes (Single QTL genome scan)
#################

## For this method, it is essential to impute the missing genotypes.
QTL_LR5RalstoDataImp<- sim.geno(QTL_LR5RalstoData, step=2, n.draws=16, error.prob=0.01)    
#sim.geno simulates sequences of genotypes from their joint distribution, given the observed marker data.
#Uses the hidden Markov model technology to simulate from the joint distribution Pr(g | O) 
#where g is the underlying genotype vector and O is the observed multipoint marker data, with possible allowance for genotyping errors.
#n.draws=16 = Number of imputations
# error.prob: probability to impute a genotype A to the marker while it is B
# step = 2, every 2cM
#For each individual with a missing genotype, this method creates 16 individuals with the same phenotype as the parent individual and with the probable genotypes for this marker according to the distance between the 2 markers
#This method virtually increases the size of the population therefore increases the power of detection


###
## Example 1: Analysis on the disease score at 7dpi
###

## QTL detection
QTLT7_LR5Ralsto.imp<-scanone(QTL_LR5RalstoDataImp, pheno.col=5,method="imp")#T7
QTLT7_LR5Ralsto.imp
summary(QTLT7_LR5Ralsto.imp) # only the lod max on each chromosome

## Graphic visualization comparing the major QTL on LG5 detected by the imp and em method.
x11()
plot(QTLT7_LR5Ralsto.imp, QTLT7_LR5Ralsto.em, chr="LG5", col=c("green", "pink"), ylab="LodScore Ralsto LR5 T7")

###
## Example 2: Analysis of the quantity of bacteria (in log (cfu)) present in the aerial part of the plant at the end of the experiment (cfu)
###

#It's your turn...

#################
######## G. DETECTION OF QTL with SCANONE (Single QTL genome scan): hk method Haley-Knott regression (Haley-Knott, 1992)
#################

###
## Example: Analysis on all the phenotypic traits measured.
###

#By default, Rqtl only works on the 1st phenotype (here, name of the lines !!)
#It must be indicated which phenotype(s) to consider here pheno.col = 2: 13

## Detection of QTL with the Haley-Knott regression method (Haley-Knott, 1992)
QTLall_LR5Ralsto.hk <- scanone(QTL_LR5RalstoData, pheno.col=2:13,method="hk")      # all cfu + scoredpi phenotypes
summary(QTLall_LR5Ralsto.hk)
#Warning: this summary of results table is misleading because it only displays for each chromosome the marker with the highest lod score for the first phenotype analyzed (here Moyenne.logcfu).

## Graphic visualization
x11()
plot(QTLall_LR5Ralsto.hk, lodcolumn=1:12, col=c("blue", "red", "green"))
# Problem: Only the first 3 phenotypes are displayed (Moyenne.logcfu, T3, T5) and not all of the 10 phenotypes

# We display the graphs by subset of 3 phenotypes
## Graphic visualization
x11()
plot(QTLall_LR5Ralsto.hk, lodcolumn=1:3, col=c("blue", "red", "green"))
x11()
plot(QTLall_LR5Ralsto.hk, lodcolumn=2:4, col=c("blue", "red", "green"))

x11()
plot(QTLall_LR5Ralsto.hk, lodcolumn=4:6, col=c("pink", "orange", "grey"))
x11()
plot(QTLall_LR5Ralsto.hk, lodcolumn=7:9, col=c("violet", "cyan", "yellow"))
x11()
plot(QTLall_LR5Ralsto.hk, lodcolumn=10:12, col=c("black", "darkgreen", "purple"))

## the listings of QTL positions, by phenotype ....
# rq: it also edits the lod values for the other phenotypes at the peaks considered - 1 single peak per chr
summary(QTLall_LR5Ralsto.hk,format="onepheno",lodcolumn=1)  # cfu
summary(QTLall_LR5Ralsto.hk,format="onepheno",lodcolumn=2)  # T3
summary(QTLall_LR5Ralsto.hk,format="onepheno",lodcolumn=3)  # T5
summary(QTLall_LR5Ralsto.hk,format="onepheno",lodcolumn=4)  # T7
summary(QTLall_LR5Ralsto.hk,format="onepheno",lodcolumn=11)  # T24

###
## Calculation of lod thresholds by permutation
###
#WARNING!!!!!!! VERY LONG CALCULATION PER 1000 PERMUTATIONS. DO FOR 200 AT THE FIRST TIME as part of this tutorial.QTL_LR5RalstoData.operm.hk<-scanone(QTL_LR5RalstoData, pheno.col=2:13, method="hk", n.perm=200)

#save(QTL_LR5RalstoData.operm.hk,file="QTL_LR5RalstoData.operm200hk.RData")    to save data and calculations
#load("QTL_LR5RalstoData.operm200hk.RData")             # to retrieve the data of this calculation: "load" function 

summary(QTL_LR5RalstoData.operm.hk, alpha=0.05)
#LOD thresholds (200 permutations)
#   Moyenne.logcfu  T3  T5   T7  T10  T12  T14  T17  T19  T21  T24  T26
#5%           2.64 2.5 2.4 2.38 2.36 2.45 2.59 2.56 2.49 2.42 2.53 2.52
summary(QTL_LR5RalstoData.operm.hk, alpha=0.1)
summary(QTL_LR5RalstoData.operm.hk, alpha=0.01)


summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=1, alpha=0.05, pvalues=TRUE)#cfu
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=2, alpha=0.05, pvalues=TRUE)#T3 pas significatif
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=3, alpha=0.05, pvalues=TRUE)#T5
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=4, alpha=0.05, pvalues=TRUE)#T7
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=5, alpha=0.05, pvalues=TRUE)#T10
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=6, alpha=0.05, pvalues=TRUE)#T12
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=7, alpha=0.05, pvalues=TRUE)#T14
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=8, alpha=0.05, pvalues=TRUE)#T17
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=9, alpha=0.05, pvalues=TRUE)#T19
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=10, alpha=0.05, pvalues=TRUE)#T21
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=11, alpha=0.05, pvalues=TRUE)#T24
summary(QTLall_LR5Ralsto.hk, perms=QTL_LR5RalstoData.operm.hk, lodcolumn=12, alpha=0.05, pvalues=TRUE)#T26

## For the rest of the analysis, we continue on the cfu and T5 dataset (most significant dpi score)

#### Calculation of QTL confidence intervals and identification of flanking markers
# By default, decrease by 1.5 of the lod score otherwise use the argument drop = 1 for example.

lodint(QTLall_LR5Ralsto.hk, lodcolumn=1, chr="LG5", expandtomarkers=TRUE, drop=1)
lodint(QTLall_LR5Ralsto.hk, lodcolumn=1, chr="LG5", drop=1)

lodint(QTLall_LR5Ralsto.hk, lodcolumn=3, chr="LG5", expandtomarkers=TRUE, drop=1)
lodint(QTLall_LR5Ralsto.hk, lodcolumn=3, chr="LG5", drop=1)

#################
######## H. DETECTION OF QTL with SCANONE (Single QTL genome scan): DETECTION OF QTL with covariates, Haley & Knott regression
#################

###
## Example 1: Analysis on the disease score at 5dpi (T5)
###


max(QTLall_LR5Ralsto.hk, lodcolumn=3)#T5
QTLT5_LR5Ralsto.hk<-scanone(QTL_LR5RalstoData, pheno.col=4, method="hk")

x11()
plot(QTLT5_LR5Ralsto.hk)

max(QTLT5_LR5Ralsto.hk)
find.marker(QTL_LR5RalstoData, chr="LG5", pos=30.2) #Goal: find the marker closest to the qtl lod peak on ch.5
#Marker flanking the nearest qtl MTE33

## Use MTE33 as a cofactor
g<-pull.geno(QTL_LR5RalstoData) [, "MTE33"]
mean(is.na(g)) # ~ 3.4% of individuals are not genotyped on this marker
gnum<-cbind(as.numeric(g==1), as.numeric(g==2))
#g<-pull.geno(fill.geno((QTL_LR5RalstoData)) [, "MTE33"] # We complete the genotypes by imputation. pb: 1 single imputation. 
QTLT5_LR5Ralsto.hk.MTE33 <- scanone(QTL_LR5RalstoData, pheno.col=4, method="hk", addcovar=gnum)
summary(QTLT5_LR5Ralsto.hk.MTE33)

x11()
plot(QTLT5_LR5Ralsto.hk, QTLT5_LR5Ralsto.hk.MTE33, col=c("blue", "red"))
# conclusion ?
lodint(QTLT5_LR5Ralsto.hk, chr="LG5", expandtomarkers=TRUE, drop=1)

### can be a threshold problem? Let's calculate permutations
QTLT5_LR5Ralsto.hk.MTE33operm <- scanone(QTL_LR5RalstoData, pheno.col=4, method="hk", addcovar=gnum, n.perm=200)
summary(QTLT5_LR5Ralsto.hk.MTE33operm)

## Use MTE34 as a cofactor
g4<-pull.geno(QTL_LR5RalstoData) [, "MTE34"]
mean(is.na(g4)) # ~ 6.9% of individuals are not genotyped on this marker

QTLT5_LR5Ralsto.hk.MTE34<-scanone(QTL_LR5RalstoData, pheno.col=4, method="hk", addcovar=g4)
summary(QTLT5_LR5Ralsto.hk.MTE34)

x11()
plot(QTLT5_LR5Ralsto.hk, QTLT5_LR5Ralsto.hk.MTE34, col=c("blue", "red"))

x11()
plot(QTLT5_LR5Ralsto.hk, QTLT5_LR5Ralsto.hk.MTE34, chr="LG5",col=c("blue", "red"))

x11()
plot(QTLT5_LR5Ralsto.hk, QTLT5_LR5Ralsto.hk.MTE33, QTLT5_LR5Ralsto.hk.MTE34, chr="LG5",col=c("blue", "darkgreen", "red"))

lodint(QTLT5_LR5Ralsto.hk, chr="LG5", drop=1)
lodint(QTLT5_LR5Ralsto.hk.MTE33, chr="LG5", drop=1)
lodint(QTLT5_LR5Ralsto.hk.MTE34, chr="LG5", drop=1)
## We don't even reduce the confidence interval!
## To go further, it is necessary to work with a denser LG5 genetic map, because here the choice of cofactors interferes too much on the detection of the QTL. 
## The possible cofactors are 'too' close and there is not enough recombination in the population between these cofactors and the potential QTL.


### Graphic representation of QTL on genetic map
qtlGraphMap<-makeqtl(QTL_LR5RalstoData, chr=c("LG5"), pos=30.2, what="prob")
x11()
plot(qtlGraphMap, main="LR5Ralsto T5")
x11()
effectplot(QTL_LR5RalstoData, mname1="LG5@30.2", pheno.col=4, main="LR5 ralsto T5 QTL LG5 effect")
x11()
plotPXG(QTL_LR5RalstoData, marker="MTE33", pheno.col=4, main="LR5 ralsto T5 QTL LG5 effect")



###
## Example 2: Analysis on the cfu
###

QTLcfu_LR5Ralsto.hk<-scanone(QTL_LR5RalstoData, pheno.col=2, method="hk")
max(QTLcfu_LR5Ralsto.hk)
find.marker(QTL_LR5RalstoData, chr="LG5", pos=22.7) #Goal: find the marker closest to the qtl lod peak on ch.5

g2<-pull.geno(QTL_LR5RalstoData) [, "MTE32"]
mean(is.na(g2)) # ~ 2.3% of individuals are not genotyped on this marker
#g2<-pull.geno(fill.geno((QTL_LR5RalstoData)) [, "MTE32"] # We complete the genotypes by imputation. pb: 1 only imputation. 
QTLcfu_LR5Ralsto.hk.MTE32<-scanone(QTL_LR5RalstoData, pheno.col=2, method="hk", addcovar=g2)
summary(QTLcfu_LR5Ralsto.hk.MTE32)

x11()
plot(QTLcfu_LR5Ralsto.hk, QTLcfu_LR5Ralsto.hk.MTE32, col=c("blue", "red"))

lodint(QTLcfu_LR5Ralsto.hk, chr="LG5", expandtomarkers=TRUE, drop=1)
#      chr  pos      lod
#MTE30 LG5  0.0 0.474508
#MTE32 LG5 22.7 2.891881
#MTE33 LG5 30.2 1.737840
lodint(QTLcfu_LR5Ralsto.hk.MTE32, chr="LG5", expandtomarkers=TRUE, drop=1)

# Same case as for the previous phenotype T5.

#################
######## I. DETECTION OF QTL with Multiple QTL Mapping
#################
#MQM is an automatic three-stage procedure in which, 
#in the first stage, missing data is 'augmented. In other words, rather than guessing one likely genotype, multiple genotypes are modeled with their estimated probabilities. 
#In the second stage important markers are selected by multiple regression and backward elimination. 
#In the third stage a QTL is moved along the chromosomes using these pre-selected markers as cofactors, except for the markers in the window around the interval under study. 
#QTL are (interval) mapped using the most 'informative' model through maximum likelihood. 

## first stage: 'Filling in' empty genotypes
#one may perform rough genome scans by marker regression without having to drop individuals with missing genotype data.
x11()
geno.image(QTL_LR5RalstoDataRaw)

### Several possible methods
#m.imp <- mqmaugment(QTL_LR5RalstoDataRaw)  ## version 'by default', often impossible to use
#Pb: if too much missing data doesn't work
# Impossible to do this calculation for the set of possible combinations with the classic mqmaugment command as below:
# minprob = The default value of 0.1 will drop all genotypes that are less likely of occurring than 1 in 10, compared against the most likely genotype.
#Calculation much too long !!!

# Impose only one possible imputation (the most probable), as recommended in the tutorial mqm.
m.imp <- mqmaugment(QTL_LR5RalstoData, minprob=1.0)

## Other methods : 
#m.imp <- fill.geno(QTL_LR5RalstoDataRaw,method='argmax')
#argmax: for each individual the most probable sequence of genotypes, given the observed data

#With method="no_dbl_XO", non-recombinant intervals are filled in; recombinant intervals are left missing. For example, a sequence of genotypes like A---A---H---H---A (with A and H corresponding to genotypes AA and AB, respectively, and with - being a missing value) will be filled in as AAAAA---HHHHH---A.

#With method="maxmarginal", the conditional genotype probabilities are calculated with calc.genoprob, and then at each marker, the most probable genotype is determined. This is taken as the imputed genotype if it has probability greater than min.prob; otherwise it is made missing.

#With method="no_dbl_XO" and method="maxmarginal", some missing genotypes likely remain. With method="maxmarginal", some observed genotypes may be made missing.

x11()
geno.image(m.imp)

###
## Example 1: Analysis on the disease score at 5dpi (T5)
###

abrev <- 'Disease index T5'        
phenotype <- 4
colonne <- phenotype-1

## sans cofacteurs
Phenotype.MQMnocofactor <- mqmscan(m.imp, pheno.col=phenotype, window.size=10, step.size=1)

# comparison mqm without cofactors and ehk: must be similar
summary(Phenotype.MQMnocofactor, threshold=1.5)
### comparison with hk
summary(QTLall_LR5Ralsto.hk,threshold=1.5, lodcolumn=colonne)

## both QTLs at the same time:
x11()
plot(Phenotype.MQMnocofactor, 
         chr='LG5',
     show.marker.names=T, col=c('red'),
#	xlab="Map Position on LG5 (cM)",
     ylab="LOD Score")

plot(QTLall_LR5Ralsto.hk,chr="LG5",
     col=c('darkgreen'), add=T)

abline(h=2.5,lty=2)

legend('topleft',col=c('darkgreen','red'),lwd=2,legend=c(paste(abrev,'; ehk'),paste(abrev,'; mqm')))
#dev.off()

########
######## we start a 'forward' search
########

#######
# indexes of useful cofactors, 1 by chromosome
summary(Phenotype.MQMnocofactor)
find.marker(m.imp, chr= summary(Phenotype.MQMnocofactor)$chr, pos= summary(Phenotype.MQMnocofactor)$pos)

index1 <- find.markerindex(m.imp, find.marker(m.imp, chr= summary(Phenotype.MQMnocofactor)$chr, pos= summary(Phenotype.MQMnocofactor)$pos))
index1

# the set of cofactors, defined on the set of m.imp genotypes thanks to the indexes of index1:
set1 <- mqmsetcofactors(m.imp,cofactors=index1)
set1

### detection with cofactors:
Phenotype.MQMset1 <- mqmscan(m.imp , set1, cofactor.significance=0.99, pheno.col=phenotype, window.size=10,step.size=1)


## cofactor.significance = 0.99 is a 'joke': this means keeping the cofactors from the moment the pvalue (of Ho: the marker has no effect) is less than 0.99
summary(Phenotype.MQMset1)

x11()
plot(mqmgetmodel(Phenotype.MQMset1))


### we compare the three methods: ehk, mqm without cofactors and mqm with forward cofactors
x11()
plot(Phenotype.MQMnocofactor,Phenotype.MQMset1,
     show.marker.names=T, col=c('blue','red'),
#	xlab="Map Position on LG5 (cM)",
     ylab="LOD Score"
     )
plot(QTLall_LR5Ralsto.hk, lodcolumn=colonne,
     show.marker.names=T, col=c('darkgreen'),
#	xlab="Map Position on LG5 (cM)",
#     ylab="LOD Score",
     add=T
     )
abline(h=3,lty=2)
legend('topleft',col=c('blue','red','darkgreen'),lwd=2,legend=c(paste(abrev,'; mqm'),paste(abrev,'; mqm with cofactors forward'),paste(abrev,'; ehk')))



##########
# the indexes of 'useful' cofactors: based on guesstimate with a lod of 2.5 (defined in threshold)
seuil <- 2.5

summary(Phenotype.MQMnocofactor,threshold=seuil)

index2 <- find.markerindex(m.imp, find.marker(m.imp, chr= summary(Phenotype.MQMnocofactor, threshold=seuil)$chr, pos= summary(Phenotype.MQMnocofactor, threshold=seuil)$pos))
index2

# the set of cofactors, defined on the set of m.imp genotypes thanks to the indexes of index2:
set2 <- mqmsetcofactors(m.imp,cofactors=index2)


### detection with cofactors, the 'strongest' in LOD:
Phenotype.MQMset2 <- mqmscan(m.imp , set2, cofactor.significance=0.02,pheno.col=phenotype, window.size=10, step.size=1)

summary(Phenotype.MQMset2)


x11()
plot(mqmgetmodel(Phenotype.MQMset2))


### we compare the four methods: ehk, mqm without cofactors and mqm with forward cofactors according to two methods
x11()
plot(Phenotype.MQMnocofactor,Phenotype.MQMset1,Phenotype.MQMset2,
     show.marker.names=T, col=c('darkgreen','blue','red'),
#	xlab="Map Position on LG5 (cM)",
     ylab="LOD Score"
     )
plot(QTLall_LR5Ralsto.hk, lodcolumn=colonne,
     show.marker.names=T, col=c('black'),
#	xlab="Map Position on LG5 (cM)",
#     ylab="LOD Score",
     add=T
     )
abline(h=3,lty=2)
legend('topleft',col=c('darkgreen','blue','red','black'),lwd=2,legend=c(paste(abrev,'; mqm'),paste(abrev,'; mqm avec cofacteurs forward larges'),paste(abrev,'; mqm avec cofacteurs forward restreints'),paste(abrev,'; ehk')))





####################
#
# Unsupervised cofactor selection through backward elimination: 10 cofactors
  autocofactorsPhenotype <- mqmautocofactors(m.imp, 10)

########
## another function that can be used to choose the cofactors (see instructions for use)
## autocofactorsPhenotype<-mqmsetcofactors(m.imp, 20)
########


# by running mqmautocofactors several times, we see that different cofactors are chosen each time
x11()
mqmplot.cofactors(m.imp, autocofactorsPhenotype, justdots=TRUE)

# detection of QTL by mqm, with cofactors determined by unsupervised, but which will be tested for their significance
Phenotype.MQMsetauto <- mqmscan(m.imp, autocofactorsPhenotype, pheno.col=phenotype, window.size=8, step.size=1, cofactor.significance=0.05)
summary(Phenotype.MQMsetauto)


### what are the selected cofactors?
x11()
plot(mqmgetmodel(Phenotype.MQMsetauto))


## comparison of the three mqm: base, forward cofactors and backward cofactors
x11()
#png('MSa+P_MQM_forbackward_5.png')
plot(Phenotype.MQMnocofactor,Phenotype.MQMset2,Phenotype.MQMsetauto,
     show.marker.names=T, col=c('darkgreen','red','blue'),
#	xlab="Map Position on LG5 (cM)",
     ylab="LOD Score"
     )
abline(h=3.5,lty=2)
legend('topleft',col=c('darkgreen','red','blue'),lwd=2,legend=c(paste(abrev,'; mqm'),paste(abrev,'; mqm with restricted forward cofactors'), paste(abrev,'; mqm with backward cofactors')))
#dev.off()





#########
######### data permutations, to determine the threshold LODs in 'data-driven' mode
#########
#require(snow)         # to use multiprocessors on multiprocessor computers.
nbpermut <- 150       # ridiculously low, but pedagogically acceptable


# We take the indexes of the cofactors retained in backward
index2 <- find.markerindex(m.imp, mqmgetmodel(Phenotype.MQMsetauto)$name )
set2 <- mqmsetcofactors(m.imp, cofactors=index2)

# calculations of permutations
resultspermut <- mqmpermutation(m.imp, pheno.col=phenotype, scanfunction=mqmscan , window.size=10, step.size=1, cofactors= set2, n.perm=nbpermut, multicore=TRUE,batchsize=5 )

resultsummary <- mqmprocesspermutation(resultspermut)

summary(resultsummary)


## Shape of permutations.
x11()
mqmplot.permutations(resultspermut)

seuil <- round(summary(resultsummary)[1,1],1)
seuil

QTLs <- summary(Phenotype.MQMsetauto,threshold=seuil)
QTLs

## how many groups and which ones?
Pics <- QTLs$chr
Pics


##
## a close-up on the groups of interest - forward and backward comparison
##

x11()
  icfin<-lodint(Phenotype.MQMsetauto, chr=Pics[1], drop=1)

  plot(Phenotype.MQMset2, Phenotype.MQMsetauto ,chr=Pics[1] ,
       show.marker.names=T,
       col=c('blue','red'),lwd=3,
       xlab=paste("Map Position on", Pics[1] ,"(cM)"),
       ylab=paste("LOD Score", abrev),ylim=c(0,max(QTLs[,3]+0.5))
       )
  legend('topleft',col=c('blue','red'),lwd=2,legend=c(paste(abrev,'; mqm forward restreint'), paste(abrev,'; mqm backward')))
  abline(h= seuil ,col="darkgreen",lty=2)
  abline(v=c(icfin$pos), col=c("blue", "red", "blue"), lty=c(2,1,2))



### QTL figure at whole genome level
x11()

par(mar=c(6,6,1,2)) # margins modification

plot(Phenotype.MQMsetauto,
#     show.marker.names=T,
     col=c('black'),lwd=4,     
#	xlab="Map Position on LG5 (cM)",
     ylab="LOD Score",
     cex.lab=3
     )
abline(h=seuil,lty=2)
legend('topleft',col=c('black'),lwd=2,legend=c(paste(abrev,'; mqm with backward cofactors')))


