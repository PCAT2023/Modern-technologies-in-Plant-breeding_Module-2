#' ---
#' title:  |
#'   | Case Study  : Balanced Incomplete Block Design - example 2
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
#'  The design involves 6 varieties of wheat crop in a BIBD with
#'  10 blocks of 3 plots each.
#' 
#' # PREPARATION OF THE WORKING INTERFACE IN R 
#' 
### I. Set working directory  ####
# On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
# Choose the directory containing the datafile and the associated R script.

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


# this is a trick to detect which folder contains the R script and the data
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)

#'
#' # LOADING REQUIRED METHODS FOR ANALYSIS 
#'

library(agricolae)
library(emmeans)


#' # ANALYSIS OF THE CASE STUDY 
#' 
## Load the data and examine 
BIBD2 <- read.table('BIBD2.txt', sep = ' ', header = TRUE)
str(BIBD2)
BIBD2$block <- factor(BIBD2$block)
BIBD2$variety <- factor(BIBD2$variety)
str(BIBD2)
## Draw at the white board. What are the characteristics of this design ?  

BIBD2
BIBD2[order(BIBD2$block), ]


#################   
#################   generate the same BIBD using agricolae
################# 

## function to create BIBD
str(design.bib)  
# function (trt, k, r = NULL, serie = 2, seed = 0, kinds = "Super-Duper", 
#           maxRep = 20, randomization = TRUE)  
  
trt <- c("V1", "V2", "V3", "V4", 'V5', 'V6' )  # treatments
k <- 3   # size of blocks

## generate the exp. design. Remember the properties of a BIBD
outdesign <- design.bib(trt, k, seed = 3264, serie = 2)


print(outdesign$parameters)
book <- outdesign$book
plots <- as.numeric( book[,1] )
matrix(plots, byrow = TRUE, ncol = k)
print(outdesign$sketch)
# write in hard disk
# write.csv(book,"book.csv", row.names=FALSE)
# file.show("book.csv")
#################   


#################   
#################   Analysis of the BIBD using agricolae
################# 

Analysis <- BIB.test(block = BIBD2$block, 
         trt = BIBD2$variety, 
         y = BIBD2$y,
         test = c("tukey"),
         alpha = 0.05, group = TRUE,  # compute groups of means
         console = TRUE)

## use of agricolae::bar.group() function 
x11()
par(mfrow=c(1,2), cex = 1)  ## two graphics within the same figure

bar.group(Analysis$groups, 
          col = "blue", border = "red", density = 6,
         ylim = c(0, 90),
         main = 'Adjusted means'
          )

bar.err(Analysis$means, variation = "SD",
        col = 'grey10', density = 8,
        ylim = c(0,90),
        main = "Uncorrected means \n and standard deviations" )


#################   
#################   Analysis of the BIBD by hand' 
################# 

## Because not all treatmnents are tested in all block, the SS of block and 
## treatment are NOT independant

## ANOVA for treatment adjusted for block effect, block unadjusted
Ana1a <- aov(y ~ block + variety, data = BIBD2)
summary(Ana1a)
## marginal means ie adjusted means
emmeans(Ana1a, ~ variety)

## ANOVA for blocks adjusted for treatment effect, treatment unadjusted
## ** NOT really useful **
Ana1b <- aov(y ~ variety + block, data = BIBD2)
summary(Ana1b)
## marginal means ie adjusted means
emmeans(Ana1b, ~ variety)
