###############################################################################################
###                 Scientific Training Center for Plant Biotechnologies                    ###
###                                   REGULAR LEVEL                                         ###
###                                       ----------                                        ###
###                               Hands-on: Building Genetic Map using Rqtl                 ###
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
library(qtl)#For genetic map construction & QTL analysis
library(broman) # includes myround() and holmans triangle stuff


### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())

# postscript figures with colours
trellis.device(device = postscript,color = T)

##################################################################
### R code from vignette source 'geneticmaps.Rnw'

################################################################
################################################################
###PRELIMINARY STEPS TO PREPARE A SIMPLE 'TOY' DATASET
################################################################
################################################################

###################################################
### code chunk number 1: geneticmaps.Rnw:38-40
###################################################
options(width=87, digits=3, scipen=4)
set.seed(61777369)


###################################################
### code chunk number 2: myround
###################################################
#library(broman) # includes myround() and holmans triangle stuff


###################################################
### code chunk number 3: loaddata
###################################################
#library(qtl)
data(mapthis)

#In practice, one would use the function read.cross() to load a data set into R. 
#See the help file (type ?read.cross), look at the sample data sets at https://rqtl.org/sampledata,
#and consider Chapter 2 of Broman and Sen (2009) A Guide to QTL Mapping with R/qtl (https://rqtl.org/book).

str(mapthis)

###################################################
#To remove allele-miscalling problems in initial dataset and create a simple toy dataset
toswitch <- c("C3M16", "C3M3",  "C3M4",  "C3M2",  "C3M10", "C3M6",  "C5M8",  "C5M5",  "C5M1",  "C5M6", "C2M16", "C2M24", "C2M13", "C1M23", "C1M16", "C1M1",  "C1M36")
mapthis <- switchAlleles(mapthis, toswitch)

################################################################
################################################################
###DATA EXPLORATION & QUALITY TRIMMING
################################################################
################################################################


###################################################
### code chunk number 5: summarycross
###################################################
summary(mapthis)

###################################################
### code chunk number 6: plotmissing (eval = FALSE)
###################################################
x11()
plotMissing(mapthis)
#The result indicates several individuals with a great deal of missing data (horizontal lines), 
#as well as several markers with a great deal of missing data (vertical lines)

###################################################
### code chunk number 9: plotntypedplot
###################################################
x11()
par(mfrow=c(1,2), las=1, cex=0.8)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",
     main="No. genotypes by marker")
#The function ntyped() provides the numbers of genotyped markers for each individual (or the number of genotyped individuals for each marker).

#six individuals missing almost all genotypes, and there are four markers that are missing genotypes at about half of the individuals. 
#Such appreciable missing data often indicates a problem (either bad DNAs or difficult-to-call markers) and can lead to difficulties in the genetic map construction
#so it is best, at this stage, to omit them,
#though one might consider adding them back in later

###################################################
### code chunk number 10: dropind
###################################################
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>50))
summary(mapthis)

#For further analysis, we only keep individual genotyped on at least 50 markers.

###################################################
### code chunk number 11: dropmarkers
###################################################
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 200])
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)

#To omit the markers with lots of missing data (herein, less than 200 genotyping datapoints available over 300 individuals),
#we first need to identify the names of the markers. 
#We then use drop.markers().

###################################################
### code chunk number 13: comparegenoplot
###################################################

#Identify duplicate individuals
#useful to compare the genotypes between all pairs of individuals, in order to reveal pairs with unusually similar genotypes, 
#which may indicate sample duplications or monozygotic twins. 
#In either case, we will want to remove one individual from each pair.

cg <- comparegeno(mapthis)
#output is a matrix, whose contents are the proportions of markers at which the pairs have matching genotypes

x11()
par(mar=c(4.1,4.1,0.1,0.6),las=1)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes",
     main="")
rug(cg[lower.tri(cg)])
#a pair of F2 siblings typically share genotypes at 40% of the markers. 
#But there are some pairs with well over 90% matching genotypes.

###################################################
### code chunk number 14: matchingpairs
###################################################
#identify these pairs that share genotypes at >90% of the markers
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

###################################################
### code chunk number 15: matchinggenotypes
###################################################
#inspect the genotype matches for these pairs
g <- pull.geno(mapthis) #1:A, 2:H, 3:B
table(g[144,], g[292,])
table(g[214,], g[216,])
table(g[238,], g[288,])
#in each case, the pairs have matching genotypes at all but 2 or 3 markers.
#simply omit one individual from each pair. See below.. 

###################################################
### code chunk number 16: dropmismatches
###################################################
#First we will omit the genotypes that mismatch, as these are indicated to be errors in one or the other individual (or both).
for(i in 1:nrow(wh)) {
  tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
  mapthis$geno[[1]]$data[wh[i,1],tozero] <- NA
}

###################################################
### code chunk number 17: omitdup
###################################################
#Now, we omit one individual from each pair
mapthis <- subset(mapthis, ind=-wh[,2])

###################################################
### code chunk number 18: finddupmar
###################################################
#useful to look for duplicate markers (that is, markers with identical genotypes).
#This is particularly true for the case of very large sets of markers
#multiple markers with identical genotypes will invariably map to the same location, 
#and so one might as well thin out the markers so that there are no such duplicates, 
#as the extra markers simply slow down all of our analyses.

print(dup <- findDupMarkers(mapthis, exact.only=FALSE))
#no markers with matching genotypes.

###################################################
### code chunk number 19: genotable: Look for markers with distorted segregation patterns
###################################################
#In F2, We expect the genotypes to appear with the frequencies 1:2:1. 
#Moderate departures from these frequencies are not unusual and may indicate the presence of partially lethal alleles. 
#Gross departures from these frequencies often indicate problematic markers that should be omitted, at least initially: 
#monomorphic markers (that is, where the two parental lines actually have the same allele), 
#or markers that are especially difficult to call (e.g., AA are often called as AB). 

gt <- geno.table(mapthis)
#geno.table fct calculates the genotype frequencies and also a P-value for a test of departure from the 1:2:1 expected ratios. 

#We will focus on those markers that show significant distortion at the 5% level, after a Bonferroni correction for the multiple tests.
gt[gt$P.value < 0.05/totmar(mapthis),]

#The first of these markers, C4M2, is not terrible. 
#The others appear to be monomorphic with a few errors (C2M9 and C2M15) or have one genotype that is quite rare, 
#which likely indicates difficulties in genotyping. 
#It would be best to omit the worst of these markers (i.e. with p-value<1.e-10).See below...

###################################################
### code chunk number 20: dropbadmarkers
###################################################
todrop <- rownames(gt[gt$P.value < 1e-10,])
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)

###################################################
### code chunk number 22: plotgenofreqbyind
###################################################
#Just as we expect the markers to segregate 1:2:1, we expect the individuals to have genotype frequencies that are in ~ those proportions. 
#Studying such genotype frequencies may help to identify individuals with high genotyping error rates 
#or some other labeling or breeding mistake

g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
x11()
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
       ylim=c(0,1))
#No particular problems, though the small number of short chromosomes result in considerable variability, 
#including one individual with no BB genotypes, and the frequencies of AB genotypes varies from 10-89%. 
#However, if there was an individual with ??? 100% AA or BB genotypes (like one of the parental strains), we would see it here.

###################################################
### code chunk number 23: triangleplot
###################################################
#more clear view of these genotype frequencies is obtained by representing them as points in an equilateral triangle.
#For any point within an equilateral triangle, the sum of the distances to the three sides is constant, and so one may represent a trinomial distribution as a point within the triangle.

x11()
par(mar=rep(0.1,4), pty="s")
triplot(labels=c("AA","AB","BB"))
tripoints(gfreq, cex=0.8)
tripoints(c(0.25, 0.5, 0.25), col="red", lwd=2, cex=1, pch=4)

#red X in the center of the figure corresponds to the expected genotype frequencies (1/4, 1/2, 1/4).
#Each black dot corresponds to an individual's genotype frequencies.
#one point along the right edge corresponds to the individual with no BB genotypes.

################################################################
################################################################
###GENETIC MAP CONSTRUCTION
################################################################
################################################################

###################################################
###################################################
###1. CREATING LINKAGE GROUPS
###################################################
###################################################

###################################################
### code chunk number 24: pairwiselinkage
###################################################
mapthis <- est.rf(mapthis)

#function to estimate the recombination fraction between each pair
#Calculation of a LOD score for a test of rf = 1/2

#NB: In the presence of segregation distortion, unlinked markers can appear to be linked.
#If segregation distortion is rampant in a dataset (and such things do happen), the usual tests of pairwise linkage will give distorted results, 
#and so it is best to instead use a simple chi-square or likelihood ratio test, to assess the association between markers. 
#The function markerlrt() behaves just like est.rf(), but uses a general likelihood ratio test in place of the usual test of pairwise linkage

###################################################
### code chunk number 25: checkAlleles
###################################################
checkAlleles(mapthis, threshold=5)
#In our case, "no apparent problem"

#In some cases: markers may be tightly associated with other markers, but with recombination fractions well above 1/2, 
#which indicates that some markers likely have their alleles switched
#You will be able to identify those problematic markers analysing diff.in.max.LOD 
#i.e. difference between the maximum LOD score for the cases where the estimated recombination fraction is > 1/2
#and the maximum LOD score for the cases where the estimated recombination fraction is <1/2. 

###################################################
### code chunk number 27: lodvrfplot
###################################################
#plot of the LOD scores against the estimated recombination fractions for all marker pairs

rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
x11()
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex=0.8)
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#As expected, all marker pairs with large LOD scores have estimated recombination fractions <1/2.So no pb!

###################################################
### code chunk number 28: forminitialgroups
###################################################
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6)

#Two markers will be placed in the same linkage groups if they have estimated recombination fraction ??? max.rf and LOD score ??? min.lod. 
#if markers a and b are linked and markers b and c are linked, then all three are placed in the same linkage group. 
#generally start with min.lod relatively large (say 6 or 8 or even 10). 
#The appropriate value depends on the number of markers and chromosomes (and individuals)
#The aim is to get as many truly linked markers together, but avoid completely putting unlinked markers in the same linkage group. 
#It is usually easier to combine linkage groups after the fact rather than to have to split linkage groups apart.

table(lg[,2])

#The inferred linkage groups are numbered in decreasing order of size (so that linkage group 1 has the largest number of markers).

###################################################
### code chunk number 40b: formgroupsagain
###################################################
lg <- formLinkageGroups(mapthis, max.rf=0.20, min.lod=2)
table(lg[,2])

###################################################
### code chunk number 40c: formgroupsagain
###################################################
lg <- formLinkageGroups(mapthis, max.rf=0.25, min.lod=8)
table(lg[,2])

###################################################
### code chunk number 29: reorganizemarkers
###################################################
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)

###################################################
### code chunk number 31: plotrfplot
###################################################
x11()
par(mar=c(4.1,4.1,2.1,2.1), las=1)
plotRF(mapthis, main="", alternate.chrid=TRUE)
# Estimated recombination fractions for all pairs of markers, along with LOD scores for the test of r = 1/2. 
# The recombination fractions are in the upper left triangle; LOD scores are in the lower right triangle.
#Red indicates a large LOD score or a small recombination fraction, while blue is the opposite.

#plot of the pairwise recombination fractions and LOD scores looks about as expected. 
#The markers within linkage group 1 are linked to each other. The pattern within the group is a bit random, 
#but that is because we have not yet ordered the markers in any way

#we have five clear linkage groups, with markers within a group linked to one another, 
#and markers in distinct groups showing no evidence of linkage. The results couldn't possibly be cleaner! ;-)

###################################################
### code chunk number 33: plotrfonemarkerplot
###################################################
#picking out a marker from one group and studying the recombination fractions and LOD scores for that marker against all others. 
#Let us pick the third marker in linkage group 4.

x11()
par(mar=c(4.1,4.1,1.1,0.6), las=1)
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
mn4 <- markernames(mapthis, chr=4)
par(mfrow=c(2,1))
plot(rf, mn4[3], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, mn4[3], bandcol="gray70", alternate.chrid=TRUE)

#No problem: High lod associated to low rf with markers allocated to the same lg.

###################################################
### code chunk number 34: genocrosstab
###################################################
geno.crosstab(mapthis, mn4[3], mn4[1])
mn5 <- markernames(mapthis, chr=5)
geno.crosstab(mapthis, mn4[3], mn5[1])

###################################################
###################################################
###2. LOCAL REORDERING OF MOLECULAR MARKERS WITHIN EACH LINKAGE GROUP
###################################################
###################################################

#Starting with the smallest linkage groups gathering the fewest molecular markers.
#fewer possible orders and so the process is quicker.

###################################################
### code chunk number 46: chrfivemap
###################################################
mapthis <- orderMarkers(mapthis, chr=5)
#orderMarkers() may be used to establish an initial order. 
#It picks an arbitrary pair of markers, and then adds an additional marker (chosen at random), one at a time, in the best supported position. 

# pull.map() to inspect the result.
pull.map(mapthis, chr=5)

#Since the marker names were chosen to indicate the true marker order, we can see that we got the order exactly right. 
#Though note that the second marker, C5M2 was omitted at some point during the course of our analysis.
##The whole chromosome is flipped, but we have no information, from the genotype data, to orient the chromosome.
#Of course, we wouldn't know this, and so we should make an attempt to explore alternate orders, in case another order might be seen to be an improvement.

###################################################
### code chunk number 48: ripplechr5 (eval = FALSE)
###################################################
rip5 <- ripple(mapthis, chr=5, window=7)

#function ripple() is used after the addition of each marker, to consider all possible orders in a sliding window of markers, to attempt to improve marker order. 
#The argument window defines the length of the window. Larger values will explore more possible orders but will require considerably more computation time. 
#The value window=7 is usually about the largest one would ever consider.
#If window is equal to the number of markers on the chromosome, than all possible orders are considered.
# function ripple()chooses among orders in order to minimize the number of obligate crossovers.

###################################################
### code chunk number 49: summaryripple5
###################################################
summary(rip5)
#The good orders generally are those that result in the smallest numbers of crossovers. 
#comparing orders to minimize the number of obligate crossovers (with method="countxo", the default, in ripple()), 
#this was already done when we called orderMarkers(), and it is not necessary to run it again.
#Nevertheless, the results will indicate how close, in terms of number of crossovers, the next-best marker order is to the inferred one.

#Switching markers 2 and 1 results in one additional obligate crossover. 
#If, instead, one switches markers 5 and 4, there are two additional obligate crossovers.

###################################################
### code chunk number 51: ripplechr5lik (eval = FALSE)
###################################################
#The more refined, but considerably slower approach, is to compare the likelihoods of the different orders. 
#(The likelihood for a given marker order is the probability of the data assuming that order is correct and plugging in estimates of the inter-marker distances.) 
#We may also indicate a genotyping error probability (through error.prob).

#To study the likelihood of different orders, though we will want to greatly reduce the window argument, so that it can be accomplished in a reasonable amount of time.

rip5lik <- ripple(mapthis, chr=5, window=4, method="likelihood",
                   error.prob=0.005)

###################################################
### code chunk number 52: summaryripple5lik
###################################################
summary(rip5lik)
#switching markers 1 and 2 has a LOD score (that is, log10 likelihood, relative to the initial order) of 0.1, which indicates that it is slightly preferred. 
#However, the estimated chromosome length is slightly longer (38.5 versus 38.2 cM).

#We know, from the marker names, that the initial order was the true order, but in practice we would not have such information, and we would probably want to switch markers 1 and 2. 
#Usually we are looking to have the estimated chromosome length be as short as possible,
#but if we trust the likelihood calculation, we should go with the alternate order. 
#Of course,the two orders are not really distinguishable; we can't really say, on the basis of these data, whether the correct order is the first or the second.

###################################################
### code chunk number 81: comparexo2likplot
###################################################
#interesting to compare the numbers of obligate crossovers for different orders to the LOD scores. 
#We have LOD scores for a much smaller number of orders (114 versus 13680 total orders), since we had used a smaller value for window, 
#but we can pull out the obligate crossover counts for those orders that we evaluated in terms of likelihood.
x11()
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
pat5 <- apply(rip5[,1:9], 1, paste, collapse=":")
pat5lik <- apply(rip5lik[,1:9], 1, paste, collapse=":")
rip5 <- rip5[match(pat5lik, pat5),]
plot(rip5[,"obligXO"], rip5lik[,"LOD"], xlab="obligate crossover count",
     ylab="LOD score")
#clear negative relationship between the crossover counts and the likelihoods, 
#though the relationship is not perfect, particularly for the orders that are not well supported.

###################################################
### code chunk number 53: compareorder
###################################################
#results can be sensitive to the assumed genotyping error rate. how sensitive the results are to the assumed error rate?
compareorder(mapthis, chr=5, c(2, 1,3:9), error.prob=0.01)
compareorder(mapthis, chr=5, c(2, 1,3:9), error.prob=0.001)
compareorder(mapthis, chr=5, c(2, 1,3:9), error.prob=0)

#with smaller assumed genotyping error rates, the evidence in favor of switching markers 1 and 2 increases somewhat. 
#Note also that the map length increases quite a bit.
#If we were looking at these data blindly, I would likely go with switching markers 1 and 2, so let's go ahead and do that here.

###################################################
### code chunk number 54: switchorder
###################################################
mapthis <- switch.order(mapthis, chr=5, c(2,1,3:9), error.prob=0.005)
pull.map(mapthis, chr=5)
#Note that the map is slightly smaller than what we had seen above, after running orderMarkers(), as we had used the default value for error.prob in that function, and now we are using error.prob=0.005. 
#Also note that markers C5M10 and C5M9 are quite close together and are separated from the next marker by 24.6 cM. This explains why it is difficult to assess the appropriate order for these two markers.

###################################################
### code chunk number 55: orderchrfour (eval = FALSE)
###################################################
mapthis <- orderMarkers(mapthis, chr=4)
pull.map(mapthis, chr=4)
#The marker names tell us the true order, and so we see that the inferred order is correct.

###################################################
### code chunk number 58: ripplechr4 (eval = FALSE)
###################################################
rip4 <- ripple(mapthis, chr=4, window=7)

###################################################
### code chunk number 59: summaryripple4
###################################################
summary(rip4)
#The order with markers 9 and 10 switched gives the same number of obligate crossovers,
#and so we can not distinguish between these two marker orders. 

###################################################
### code chunk number 61: ripplechr4lik (eval = FALSE)
###################################################
#We turn to the likelihood comparison.
rip4lik <- ripple(mapthis, chr=4, window=4, method="likelihood",
                  error.prob=0.005)

###################################################
### code chunk number 62: summaryripple4lik
###################################################
summary(rip4lik)
#There is reasonably no good evidence to switch markers 9 and 10
#Let's keep the initial order which is nice, as this results in the true marker order

###################################################
### code chunk number 64: orderchrthree (eval = FALSE)
###################################################
mapthis <- orderMarkers(mapthis, chr=3)
pull.map(mapthis, chr=3)
#Note, from the marker names, that the marker order is the true order. 
#The whole chromosome is flipped, but we have no information, from the genotype data, to orient the chromosome.

###################################################
### code chunk number 67: ripplechr3 (eval = FALSE)
###################################################
rip3 <- ripple(mapthis, chr=3, window=7)

###################################################
### code chunk number 68: summaryripple3
###################################################
summary(rip3)
#The next-best order, with markers 2 and 3 switched, results in an additional 6 obligate crossovers.

###################################################
### code chunk number 70: ripplechr3lik (eval = FALSE)
###################################################
#We turn to the likelihood comparison.
#rip3lik <- ripple(mapthis, chr=3, window=4, method="likelihood",
#                   error.prob=0.005)


###################################################
### code chunk number 71: summaryripple3lik
###################################################
#summary(rip3lik)
#The next-best order (with markers 2 and 3 switched) is considerably worse than our initial order.

###################################################
### code chunk number 72: orderchrtwo (eval = FALSE)
###################################################

mapthis <- orderMarkers(mapthis, chr=2)
pull.map(mapthis, chr=2)

#as seen from the marker names, the inferred order is the true one.

#NB: It may happen that the function starts badly, drifts and finally does not give the expected result. 
#In our case, it is easy to realize this because we know the expected order. 
#Alternatively, you may identify strangely large genetic distances between 2 markers (i.e. >100cM). 
#I advise you to run the function several times and keep the order that minimizes the distance between markers and therefore the total size of the linkage group.

mapthis1 <- orderMarkers(mapthis, chr=2)
pull.map(mapthis, chr=2)

mapthis2 <- orderMarkers(mapthis, chr=2)
pull.map(mapthis, chr=2)

###################################################
### code chunk number 75: ripplechr2 (eval = FALSE)
###################################################
rip2 <- ripple(mapthis, chr=2, window=7)

###################################################
### code chunk number 76: summaryripple2
###################################################
summary(rip2)
#The next-best order, with markers 3 and 4 switched, has 6 additional obligate crossovers

###################################################
### code chunk number 78: ripplechr2lik (eval = FALSE)
###################################################
#We turn to the likelihood comparison
#rip2lik <- ripple(mapthis, chr=2, window=4, method="likelihood",
#                   error.prob=0.005)

###################################################
### code chunk number 79: summaryripple2lik
###################################################
#summary(rip2lik)
#The next-best order, in terms of likelihood, has markers 21 and 22 switched, and is considerably worse than our initial order.

###################################################
### code chunk number 82: orderchrone (eval = FALSE)
###################################################
mapthis <- orderMarkers(mapthis, chr=1)
pull.map(mapthis, chr=1)
#Again, the inferred order is precisely the true order.

###################################################
### code chunk number 85: ripplechr1 (eval = FALSE)
###################################################
rip1 <- ripple(mapthis, chr=1, window=7)
#The ripple() analysis is excellent for comparing nearby marker orders. 
#And for chromosomes with a modest number of markers (like chromosome 5, here), one can consider all possible orders exhaustively. 
#But with many markers (as with chromosome 1), we can investigate only a small proportion of the possible orders. 
#For example, with 33 markers (as on chromosome 1), there are 33!/2 ??? 1036 possible marker orders, and we investigated only 117,360 of them using ripple with window=7.

###################################################
### code chunk number 86: summaryripple1
###################################################
summary(rip1)
#Two alternate orders (switching markers 28 and 29, or switching markers 31 and 32) have 2 additional obligate crossovers, compared to our initial order. 

###################################################
### code chunk number 88: ripplechr1lik (eval = FALSE)
###################################################
#We turn to the likelihood comparison.
#rip1lik <- ripple(mapthis, chr=1, window=4, method="likelihood",
#                   error.prob=0.005)

###################################################
### code chunk number 89: summaryripple1lik
###################################################
#summary(rip1lik)
#The next best order (switching markers 28 and 29) results in a moderately large decrease in likelihood, 
#and so we stick with the initial order produced by orderMarkers().

###################################################
#NB: Unexpected good performance of fct orderMarkers
#quite surprised to see how well orderMarkers() performed with these data. 
#Of course, the data are quite clean (by design) and comprise quite a large number of individuals. 
#In practice, one can't expect orderMarkers() to perform so well, and it is worthwhile to study the estimated map and
#the pairwise linkage information closely. 

###################################################
### code chunk number 90: summarymap
###################################################
#look at the actual map: are there large gaps between markers, indicating adjacent markers that are only weakly linked?
summaryMap(mapthis)
#summary of the average inter-marker distance and the largest gap on each chromosome.
#By default, the map is estimated assuming no crossover interference, but a map function is used to derive the genetic distances (though, by default, the Haldane map function is used).
#chromosome 1 has a 45.6 cM gap; the other chromosomes have gaps no larger than 25 cM. 
#The large gap on chromosome 1 is suspicious, but it is not terrible.
#With gross mistakes in marker order, one will often see much larger gaps, say > 100 cM.

###################################################
### code chunk number 91: savesummarymap
###################################################
firstsummary <- summaryMap(mapthis)

###################################################
### code chunk number 93: plotmapplot
###################################################
x11()
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
plotMap(mapthis, main="", show.marker.names=TRUE)
#Most of marker names are unreadable, because the markers are densely spaced, 
#but at least we learn the identity of markers that are surrounded by large gaps.

###################################################
### code chunk number 95: plotrfonemoretimeplot
###################################################
x11()
par(mar=c(4.1,4.1,1.6,1.6), las=1)
plotRF(mapthis, main="")
# Estimated recombination fractions for all pairs of markers, along with LOD scores for the test of r = 1/2. 
# The recombination fractions are in the upper left triangle; LOD scores are in the lower right triangle.
#Red indicates a large LOD score or a small recombination fraction, while blue is the opposite.

#close to what we want to see: nearby markers show clear association and no distant markers show any association: 
#red along the diagonal, dissipating to blue away from the diagonal.

###################################################
### code chunk number 97: plotrfafterreorderplot
###################################################
#If there were gross problems with marker order, we might see groups of distantly placed markers that are more highly associated than more closely placed markers. 
#Just as an example, suppose we split chromosome 1 into three pieces and moved the latter piece into the middle.

x11()
par(mar=c(4.1,4.1,1.6,1.6), las=1, pty="s", cex=0.8)
messedup <- switch.order(mapthis, chr=1, c(1:11,23:33,12:22),
                         error.prob=0.005)
plotRF(messedup, chr=1, main="")

#sort of pattern one should expect with gross problems in marker order: 
#the markers in the first and last segments are more tightly associated to each other than they are to the markers in the middle segment

###################################################
### code chunk number 99: plotmapmessedupplot
###################################################
x11()
par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.8)
plotMap(messedup, main="", show.marker.names=TRUE)

#the genetic map shows a 50 cM gap between the first and middle segments and an 80 cM gap between the middle and last segments.
#The map is not as telling as the pairwise linkage information, but these are the sorts of things to look at, in trying to decide whether there are gross problems that need to be corrected.
#If features such as those in the 2 previous figures were seen, 
#1. one should identify the segments of markers that need to be moved around, and then use switch.order() to reorganize the markers. 
#2. One might first use compareorder() to compare the current order to the reorganized one, to see that the new order gave a clear improvement in likelihood. 
#3. In addition, one will often need to alternate between ripple() and switch.order() until the final marker order is established.

###################################################
###################################################
###2. REFINING THE GENETIC MAP
###################################################
###################################################
#Problem statement:
#The large gaps in the genetic map on chromosome 1 remain a concern. 
#While such gaps may indicate problems with the order of markers, they might also indicate a high genotyping error rate at certain markers. 
#If an individual marker is more prone to genotyping errors than others, it would often (in the sort of analyses performed above) be placed at one end of the chromosome or the other, 
#but in some cases (particularly if the genotyping error rate is not high and there are a large number of individuals in the cross) 
#it may be placed at approximately the correct position but result in reasonably large gaps surrounding the marker.

###################################################
### code chunk number 100: droponemarker (eval = FALSE)
###################################################
#One approach for identifying such problematic markers is to drop one marker at a time and
#investigate the change in chromosome length and the change in log likelihood.

dropone <- droponemarker(mapthis, error.prob=0.005)

###################################################
### code chunk number 103: plotdroponeplot
###################################################
x11()
par(mar=c(4.1,4.1,1.6,0.1), mfrow=c(2,1), cex=0.8)
plot(dropone, lod=1, ylim=c(-100,0))
plot(dropone, lod=2, ylab="Change in chr length (cM)")

#The top panel contains LOD scores; positive values would indicate that dropping a marker improves the likelihood. 
#The bottom panel indicates the decrease in estimated chromosome length (in cM) following dropping a marker.

#there is no one marker for which its omission results in an increase in likelihood, 
#but there are a number of markers that give an appreciable decrease in chromosome length, of 15-25 cM, when omitted. 
#Markers at the ends of chromosomes will often result in a smaller estimated chromosome length, 
#but that is just because such terminal markers hang off some distance from the rest of the markers, and so these changes can often be discounted. 
#Interior markers that result in a big change, (and there appears to be one on each of chromosomes 1, 2 and 3), 
#might indicate error-prone markers that are best omitted.

###################################################
### code chunk number 104: worstmarkers
###################################################
summary(dropone, lod.column=2)

#One may identify the marker on each chromosome whose omission results in the largest decreases in chromosome length 
#For chromosomes 4 and 5, these are terminal markers. 
#The markers on chromosomes 1, 2 and 3 are all interior markers. 
#One should probably study the pairwise linkage between these markers and surrounding markers before proceeding, 
#but we will go ahead and remove these markers without any further investigations. See below...

###################################################
### code chunk number 105: dropbadmarkers
###################################################
badmar <- rownames(summary(dropone, lod.column=2))[1:3]
mapthis <- drop.markers(mapthis, badmar)

###################################################
### code chunk number 106: reestimatemap
###################################################
#One should re-estimate the genetic map. We use replace.map() to insert it into the cross object.

newmap <- est.map(mapthis, error.prob=0.005)
mapthis <- replace.map(mapthis, newmap)
summaryMap(mapthis)
#Removing those three markers resulted in a decrease in the overall map length from 650 cM to 524 cM.

###################################################
### code chunk number 107: savenewsummary
###################################################
secondsummary <- summaryMap(mapthis)

###################################################
### code chunk number 109: countxoplot  - Look for problem individuals
###################################################
#Pb statement:
#Now that we have the markers in their appropriate order, it is a good idea to return to the question of whether there are particular individuals whose data are problematic. 
#If a particular individual showed considerable genotyping errors or did not actually belong to the cross under investigation (e.g. labeling or breeding error), 
#what sort of aberrations might be seen in the data? 
#amount of missing genotype data for each individual and the individuals' genotype frequencies previously studied.
#Another feature to investigate is the observed number of crossovers in each individual. 

x11()
par(mar=c(4.1,4.1,0.6,0.6), cex=0.8)
plot(countXO(mapthis), ylab="Number of crossovers")
thecounts <- countXO(mapthis)
countCObyInd <- rev(sort(thecounts, decreasing=TRUE))
countCObyInd
worst <- rev(sort(thecounts, decreasing=TRUE)[1:2])
worst
#The crossover counts clearly indicate two problematic individuals, with 73 and 86 crossovers; 
#the other individuals have 3-20 crossovers. We should remove these individuals. See below.

###################################################
### code chunk number 110: drophighxoind
###################################################
mapthis <- subset(mapthis, ind=(countXO(mapthis) < 50))
#We only keep individuals with less than 50 CO
#Ideally, we would now revisit the entire process again

###################################################
### code chunk number 111: rip5again
###################################################
#after removing these problematic individuals, is there evidence that marker order needs to be changed? 
#Let us at least look at chromosome 5.
summary(rip <- ripple(mapthis, chr=5, window=7))
summary(rip <- ripple(mapthis, chr=5, window=2, method="likelihood",
                      error.prob=0.005))
#There is good evidence for switching markers 1 and 2

###################################################
### code chunk number 112: switchchr5again
###################################################
mapthis <- switch.order(mapthis, chr=5, c(2, 1, 3:9), error.prob=0.005)
pull.map(mapthis, chr=5)

#One should investigate the other four chromosomes similarly,
#Such analysis provided no evidence for further changes in marker order (data not shown). 

###################################################
### code chunk number 113: reestmapagain
###################################################
#We should, finally, re-estimate the genetic map.
newmap <- est.map(mapthis, error.prob=0.005)
mapthis <- replace.map(mapthis, newmap)
summaryMap(mapthis)
#The overall map length has decreased further, from 524 cM to 510 cM

###################################################
### code chunk number 114: savethirdsummary
###################################################
thirdsummary <- summaryMap(mapthis)


###################################################
### code chunk number 115: studyerrorrate (eval = FALSE)- Estimate genotyping error rate 
###################################################
#Above, we had generally assumed a genotyping error rate of 5/1000 (in estimating map distances and in likelihood calculations comparing different marker orders); 
#(error rate value used to simulate the genotype data).

#estimation of the genotyping error rate from the data:
#the function est.map() not only estimates the inter-marker distances, but also calculates the log likelihood for each chromosome. 
#Thus, if we run est.map() with different assumed values for the genotyping error rate (specified with the error.prob argument), 
#one can identify the maximum likelihood estimate of the error rate.

loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
for(i in seq(along=err)) {
   cat(i, "of", length(err), "\n")
   tempmap <- est.map(mapthis, error.prob=err[i])
   loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)


###################################################
### code chunk number 118: ploterrorratelikplot
###################################################
x11()
par(mar=c(4.1,4.1,0.6,0.6), las=1)
plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02),
     ylab=expression(paste(log[10], " likelihood")))

#The log10 likelihood indicates that the MLE is approximately 0.005. 
#Error rates of 0.0025 and 0.0075 have log10 likelihoods that are 3 less than that of the MLE. 
#We might investigate additional assumed error rates, or even use the R function optimize() to refine our estimate, but we won't pursue that here.

###################################################
### code chunk number 119: errorlod (eval = FALSE)-Look for genotyping errors
###################################################
#While the methods for estimating inter-marker distances allow for the presence of genotyping errors at a fixed rate, 
#it is nevertheless worthwhile to look for, and ideally correct, potential genotyping errors in the data. 
#Such errors may be identified through apparent tight double-crossovers, with a single marker being out of phase with its adjacent markers.

#most convenient approach for identifying such double-crossovers is to calculate genotyping error LOD scores,
#The LOD score compares the likelihood for a genotype being in error versus it not being in error.

#One must assume a genotyping error rate (error.prob), but the results are almost identical for a wide range of values.

mapthis <- calc.errorlod(mapthis, error.prob=0.005)

###################################################
### code chunk number 121: toperrorlod
###################################################
print(toperr <- top.errorlod(mapthis, cutoff=6))
#List of the genotypes with the largest error LOD scores. 
#One may generally focus on those with quite large values, say at least 4-5. 
#Here will we look at just those genotypes with error LOD ??? 6 (see cutoff)

###################################################
### code chunk number 123: plotgenoplot
###################################################
x11()
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex.axis=0.9)
plotGeno(mapthis, chr=1, ind=toperr$id[toperr$chr==1], main="", cex=0.8,
          include.xo=FALSE, cutoff=6)
#Genotypes on chromosome 1 for individuals with some potential errors flagged by red squares. White, gray and black circles correspond to AA, AB and BB genotypes, respectively.
#All of the flagged genotypes are cases with an exchange from one homozygote to the other and then back again 
#(thus, two crossovers in each interval flanking a single marker)

###################################################
### code chunk number 124: dropgenotypes
###################################################
#One might zero out these suspicious genotypes (that is, make them missing). 
#Even better would be to revisit the raw genotyping information, or even re-genotype these instances. 
#But we are talking about just 12 genotypes out of 25,147, and with our allowance for genotyping errors in the map estimation, they have little influence on the results.
#If we did wish to delete these genotypes, we could do so as follows.

mapthis.clean <- mapthis
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  mapthis.clean$geno[[chr]]$data[mapthis$pheno$id==id, mar] <- NA
}

###################################################
### code chunk number 126: plotsegdis-Revisit segregation distortion
###################################################
gt <- geno.table(mapthis, scanone.output=TRUE)
x11()
par(mar=c(4.1,4.1,0.6,0.6), las=1, mfrow=c(2,1), cex=0.8)
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")

#The top panel contains ???log10 P-values from tests of 1:2:1 segregation at each marker. 
#The bottom panel contains the observed genotype frequencies at each marker (with black, blue and red corresponding to AA, AB and BB genotypes, respectively)

#The greatest departure from 1:2:1 segregation is on chromosome 4, 
#with somewhat more AA genotypes and somewhat fewer BB genotypes. 
#If we apply a Bonferroni correction for the 88 tests (88 is the total number of markers we have retained in the data), 
#we would look for P ??? 0.05/88 which corresponds to ???log10 P ??? 3.25, and there is one marker on chromosome 4 that exceeds this.
#There are also some departures from 1:2:1 segregation on chromosomes 2 and 3, but these appear to be within the range of what would be expected by chance

#The aberrant segregation pattern on chromosome 4 is not too worrisome. 
#Multipoint estimates of genetic map distances are little affected by segregation distortion, 
#and the pattern of distortion on the chromosome indicates rather smooth changes in genotype frequency. 
#More worrisome would be a single distorted marker in the midst of other markers with normal segregation, which would indicate genotyping errors rather than, 
#for example, the presence of partially lethal alleles.

###################################################
### code chunk number 128: plotfinalmapplot
###################################################
#Finally, we're done. Let us plot the final map
x11()
par(las=1, mar=c(4.6,4.6,0.6,0.6), cex=0.8)
plotMap(mapthis, main="", show.marker.names=TRUE)

#The map is not pretty, as most of the marker names are obscured. 
#R/qtl does not produce production-quality figures automatically; 
#I always go through quite a few extra contortions within R to produce a figure suitable for a paper.

###############
###LinkageMapView: Plot Linkage Group Maps with Quantitative Trait Loci
#Produces high resolution, publication ready linkage maps and quantitative trait loci maps. 
#Input can be output from 'R/qtl', simple text or comma delimited files. Output is currently a portable document file.
###############

#install.packages("LinkageMapView")
library(LinkageMapView)
x11()
lmv.linkage.plot(mapthis, "ForPublication.pdf")
# lmv.linkage.plot: main function to produce linkage group maps and has many parameters to customize the pdf output.

###############
### RgeneticMap
# A R function to draw genetic maps (linkage map), and can be used with R/qtl directly.
# You can export the maps as svg format and edit them with any svg editor (Inkscape is a good FREE svg editor).
###############
