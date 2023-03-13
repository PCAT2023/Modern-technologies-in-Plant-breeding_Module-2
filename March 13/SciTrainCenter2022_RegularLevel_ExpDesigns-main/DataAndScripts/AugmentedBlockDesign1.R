#' ---
#' title:  |
#'   | Case Study  : Augmented Block Design.
#'   | Analysis with Fixed Model
#' author:
#'   - Prof L. Gentzbittel Skoltech, Project Center for Agro Technologies ^[l.gentzbittel@skoltech.ru]
#'   - Prof C. Ben, Skoltech, Project Center for Agro Technologies ^[c.ben@skoltech.ru]
#' date: "Oct. 20th ,  2021 - Skoltech"
#' output: 
#'   pdf_document:
#'     keep_tex: true
#' use_bookdown: TRUE
#' latex_engine: xelatex
#' header-includes:
#'   - \usepackage{bbold}
#'   - \def\+#1{\mathbf{#1}}
#' geometry: left = 2cm, right = 1.5cm, top = 1.5cm, bottom = 1.5cm
#' ---
#' 
#+ echo = FALSE, message = FALSE, warning = FALSE
# just forgot these lines. They are used to generate the printed version
knitr::opts_chunk$set(fig.width = 7, fig.height=6, warning=FALSE, message=FALSE,
                      tidy.opts=list(width.cutoff=80), tidy=TRUE)
# Now we resume to nominal situation
#' 
#' # CASE STUDY PRESENTATION 
#' 
#'  A toy dataset : three blocks, 4 checks and 8 new varieties
#' 
#' 
#' # PREPARATION OF THE WORKING INTERFACE IN R 
#' 
### I. Set working directory  ####
# On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
# Choose the directory containing the Alpha_latticefile and the associated R script.

### II. Possibly, installation of new R packages needed for the analysis on RStudio:
# Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
# Comment #1: R package installation requires a connection to internet
# Comment #2: Once packages have been installed, 
# no need to re-install them again when you close-open again RStudio.

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())


# this is a trick to detect which folder contains the R script and the Alpha_lattice
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)

#'
#' # LOADING REQUIRED METHODS FOR ANALYSIS 
#'

library(agricolae)
library(emmeans)
library(ggplot2) ; library(gridExtra)

# constraints on effects for ANOVA
options(contrasts=c("contr.sum", "contr.poly"))

#' 
#' # ANALYSIS OF THE CASE STUDY 
#' 
#' The case study will ba analysed using a classical method based on fixed model (GLM) 


## ### the data, with map of the field 
Produc <- read.table("AugmentedBlockDesign1b.csv", sep = ";", dec = ".", header = TRUE)
str(Produc)

### the trial ?
with(Produc,
     table(Accession,Bloc))

## visualize the field:
library(desplot)
x11(width = 7, height = 5)
desplot(Bloc ~ x + y, Produc,
        col =  Type, text = Accession, cex = 1,
        out1 = Bloc, out2 = Accession,
        out2.gpar=list(col = "gray50", lwd = 1, lty = 1),
        main = "Augmented Block Design - Example 1")


############################################  
############################################  Analysis 'by hand' using a fixed model  ############################################
############################################  

## type I SS
Ana1a <- aov(weight ~ Bloc + Accession , data = Produc)  # accessions 'corrected for the bloc effect'
summary(Ana1a)


Ana1b <- aov(weight ~ Accession + Bloc , data = Produc)  # blocks 'corrected' for genotype effect 
summary(Ana1b)
## the two ANOVA are different because some combinations are lacking



###################### Analysis of the checks only ###################### 
# this is a simple RCBD

## SS of accessions corrected for block effect
AnaChecks <- aov(weight ~ Bloc + Accession , data = subset(Produc, subset = Type == "check"))
summary(AnaChecks)



###################### Phenotypic values of new entries, "corrected for the block effect" 
#' the *block effect* for correction is estimated in the ANOVA of the checks


## What are the block effects ? 
model.tables(AnaChecks)
model.tables(AnaChecks)$tables$Bloc

### Syntax to correct the phenotypic values of the 'new entries" 
Produc$EffetBloc <-sapply(Produc$Bloc, switch,                         ## use of the swith() fnction
                          B1=model.tables(AnaChecks)$tables$Bloc[1],   ## the Bi block effect is taken from
                          B2=model.tables(AnaChecks)$tables$Bloc[2],   ## the results of ANOVA table
                          B3=model.tables(AnaChecks)$tables$Bloc[3]
                          )

## create the vector of all data "corrected for the block effect" 
Produc$weightCorr <- Produc$weight - Produc$EffetBloc

Produc

x11()
(gr1 <- ggplot(Produc) + aes(x = weight, y = weightCorr) +
                geom_point(color = 'firebrick', cex = 4) +
                geom_abline(slope = 1, intercept = 0) +
                theme_classic()
        )


########################   to verify / validate :

## CAUTION : This verification works ONLY gor 'checks' beacuse the design is balanced for checks
## so the SS of accessions is independnat from the SS of block: we and do the ANOVA in both direction
## and it shall produce the SAME results

Valid1 <- aov(weightCorr ~ Bloc + Accession, data = subset(Produc,subset = Type == "check"))
summary(Valid1)

Valid2 <- aov(weightCorr ~ Accession + Bloc, data = subset(Produc,subset = Type == "check"))
summary(Valid2)

#' Are the results expected ?
#' 

###################### the above procedure can be simply done by computing "marginal means" 
########### also called 'lsmeans'

emmeans(Ana1a, ~ Accession)



###################### whaht accessions differ from checks ? 
################### débobinage : quelles sont les accessions différentes des checks ?

## just use the formulas of the slides: 
c <- 4
r <- 3
ddl <- (r-1)*(c-1)
MSE <- summary(Ana1a)[[1]]$'Mean Sq'[3]  # 3rd term of 'Mean Sq' column of the summary table of ANOVA 

# SE of the difference between 2 adj means of selections in different blocks (cf dias)
Sv <- sqrt(2*(c+1)*MSE/c)

#vSE of the difference between adjusted selection mean and check. c'est la quantité utile à connaitre !
Svc <- sqrt(((r+1)*(c+1)*MSE)/(r*c)) 

# We are , for example, interested to those corrected vaues that outperform the checks.
## the CI is given following a t distribution 

# critical value of t:
qt(0.95,ddl)

## numerical value of Least Significant Interval:
LSI <- qt(0.95,ddl) * Svc
LSI
# so any new entry that outperforms the 'best check' by this quantity is significantly different from it.
Produc
#' So, what is the list of the new entreis that outperfom the best check ?




############################################  
############################################  Analysis using agricolae  ############################################
############################################  

## another elegant way to see the design:
with(Produc,
by(Accession, Bloc, as.character)
)

## the agricolae::DAU() fits augmented design 

modelDAU1 <- DAU.test(Produc$Bloc, Produc$Accession, Produc$weight,  
                      method = "lsd", group = TRUE,
                      console = TRUE)
## The table "ANOVA, Treatment Adjusted" is the same as Ana1a ( Bloc + Accessions)
## the table "ANOVA, Block Adjusted" is the same as Ana1b (Accessions + Bloc)


## il y a une petite diff numérique sur Svc entre agricolae et calculs manuels

options(digits = 2)
modelDAU1$means
### see lsmeans and weightCorr to compare 


## If one wants pairwise comapsrisons détail : si on veut les comparaisons deux à deux, et non pas les groupes.
modelDAU2 <- DAU.test(Produc$Bloc, Produc$Accession, Produc$weight, 
                      method="lsd", group = FALSE, 
                      console = TRUE)

head(modelDAU2$comparison, 12)

## en conclusion : agricolae est utile et juste pour l'analyse des plans en blocs augmentés.
## Il utilise un modèle fixe d'ANOVA.
## on peut raffiner un peu mieux que agricolae.



###################
######################################  Analysis using plantbreeding
###################

install.packages("plantbreeding", repos="http://R-Forge.R-project.org")  ## unclear if maintained
## requires reshape

library(plantbreeding)
Ana3 <- aug.rcb(dataframe = Produc, genotypes = "Accession", block = "Bloc", yvar = "weight")

# analysis of variance
Ana3$anova
## compare with Ana2 and modelDAU1


Ana3$adjusted_values # yield observed and expected value table  
## same as agricolae. Does not provide for checks. 

str(Ana3)
## Ana3$se_geno (Difference between two varieties /entries and a check mean) equates Svc , as expected
##  approximation in agricolae ?


