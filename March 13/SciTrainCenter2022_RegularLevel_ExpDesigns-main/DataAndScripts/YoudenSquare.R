#' ---
#' title:  |
#'   | Case Study  : Youden Square or Incomplete latin Square
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
#'  In one of the experiments, the experimenter is interested in making comparisons among 7 
#' treatments and there are 28 experimental units available.  These 28 experimental units are 
#' arranged in a Youden Square design with 4 rows and 7 columns with one observation per 
#' cell. 
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
YoudenSq <- read.table('YoudenSquare1.csv', sep = ',', header = TRUE)
str(YoudenSq)
YoudenSq$row <- factor(YoudenSq$row)
YoudenSq$col <- factor(YoudenSq$col)
YoudenSq$trt <- factor(YoudenSq$trt)
str(YoudenSq)
## What are the characteristics of this design ? 
#'The parameters of the design are v (number of treatments)  = 7, 
#'p (number of rows)= 4, q (number of columns) =7, 
#'r (replication of treatments)  = 4.

YoudenSq



#################   
#################   generate the same Youden square using agricolae
################# 

## function to create Youden Square
str(design.youden)  
# function (trt, r, serie = 2, seed = 0, kinds = "Super-Duper", first = TRUE, 
# randomization = TRUE)  
  
trt <- c("T1", "T2", "T3", "T4", 'T5', 'T6', 'T7' )  # treatments
r  <- 4   # size of blocks

## generate the exp. design. Remember the properties of a BIBD
outdesign <- design.youden(trt, r, seed = 3264, serie = 2)


print(outdesign$parameters)
book <- outdesign$book
plots <- as.numeric( book[,1] )

print(outdesign$sketch)
# write in hard disk
# write.csv(book,"book.csv", row.names=FALSE)
# file.show("book.csv")
#################   

#################   
#################   Analysis of the YS by hand' 
################# 

## Because not all treatmnents are tested in all rows, the SS of rows, cols and 
## treatment are NOT independant.

## ANOVA for treatment adjusted for row and col effect, 
Ana1a <- aov(y ~ row + col + trt, data = YoudenSq)
summary(Ana1a)

## marginal means ie adjusted means
Fitted.Model <- emmeans(Ana1a, ~ trt)
## LSD for adjusted means
pairs(Fitted.Model, adjust = 'Tukey')

## be careful : the below is based on **observed** means so may be incorrect
## in some cases (not here because differences are clear)
outHSD <- HSD.test(Ana1a, "trt",console = TRUE)


## imagine we just give up with rows
Ana2 <- aov(y ~ col + trt, data = YoudenSq)
summary(Ana2)
#' Comparison of the two MSE show we greatly improve the model by accounting for rows
