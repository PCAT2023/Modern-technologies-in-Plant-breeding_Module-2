#' ---
#' title:  |
#'   | Case Study  : Alpha Lattice design
#' author:
#'   - Prof L. Gentzbittel Skoltech, Project Center for Agro Technologies ^[l.gentzbittel@skoltech.ru]
#'   - Prof C. Ben, Skoltech, Project Center for Agro Technologies ^[c.ben@skoltech.ru]
#' date: "Feb 2022 - Skoltech"
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
#' A spring oats trial. There were 24 varieties in 3 replicates, each consisting of 6 incomplete blocks of 4 plots. Planted in a  alpha design.
#'  72 observations on the following 5 variables: PLOT Plot number; REP Replicate code; BLOCK Incomplete block code; GEN Genotype code; YIELD Observed dry matter yield (tonnes/ha)


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


# this is a trick to detect which folder contains the R script and the datafile
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)

#'
#' # LOADING REQUIRED METHODS FOR ANALYSIS 
#'
library(agricolae)
library(tidyverse)
library(ggplot2)
library(openxlsx)   ## to import/export Excel files

#'
#' # ANALYSIS OF THE CASE STUDY 
#' 

#################   
#################   generate an alpha-lattice using agricolae
################# 

# Treatments. Here: num. of genotypes
Genotype <- c(paste("gen0",1:9,sep=""),paste("gen",10:24,sep=""))

g <- length(Genotype) ; g
r <- 3  # two replicates
k <- 4  # block size , num. of plots per block
s <- 6 # num of blocks per replicate

## we have the condition for a alpha-lattice :  g = s * k

( obs <- s * k * r )  # total num. of plots
( b <- s * r )        # total num. of blocks

book <- design.alpha(Genotype, k, r, seed = 3264)
## the Efficiency factor is 75% of a RCBD. 

#' The different quantities associated with this design:
ls(book)
book$parameters
# map of field plots
book$sketch


#' After implementation in the field and the harvest, the yield data are gathered and typed in:
dbook<- read.table('AlphaLatticeOat.csv', sep = ' ', header = TRUE)
str(dbook)

## Function for lattice analysis in agricolae: 
str(PBIB.test)
##function (block, trt, replication, y, k, method = c("REML", "ML", "VC"), test = c("lsd", "tukey"), alpha = 0.05, console = FALSE, group = TRUE)  

#' 
#' the details of this analysis rely on Mixed Linear Models (MLM). We do not go into the details here. This will be covered in 'Advanced Level'
#' 

# require(nlme) if method = REML or LM in PBIB.test; and require(MASS) if method=VC
model <- PBIB.test(dbook$BLOCK, dbook$GEN, dbook$REP, dbook$YIELD, k = 4, method = c("REML"), 
                        test = "tukey", alpha = 0.05, group = TRUE)
#' ``replication`` and ``block in replications`` are considered as random factors. ``Genotypes`` are considered as fixed factors. This is the *intrablock analysis*.

ls(model)

#' Test for differences among genotypes:
model$ANOVA
#' **Adjusted** means:
model$means
#' Groups of adj. means based on Tukey test: 
model$groups

plot(model,las=2)

## How to plot the Alpha_lattice
x11()
bar.group(model$groups, ylim = c(0,6),
          xlab = "variety", ylab = "Yield", las = 2)
x11()
bar.err(model$means, ylim = c(0,6), xlab = "variety", ylab = "Yield", las = 2)

# Save the means in a excel file:  
write.xlsx( as.data.frame(model$means %>% mutate(Genotype = row.names(.))), "AdjustedMeans.xlsx", overwrite = TRUE)


x11()
model$means %>% select('dbook$YIELD.adj', 'dbook$YIELD') %>% rename( Adj = 'dbook$YIELD.adj' , Obs = 'dbook$YIELD' ) %>% 
ggplot() + aes(y = Adj, x = Obs) +
  geom_point(size = 2) + geom_abline(slope = 1, intercept = 0) +
  xlim(3.4, 5.3) + ylim(3.4, 5.3) +
  labs( x='Obs. yield (tons/ha)', y='Adjust. yield (tons/ha)') +
  theme_minimal()


#'  **What is your opinion when comparing observed means/raw data and adjusted means** ?. 
