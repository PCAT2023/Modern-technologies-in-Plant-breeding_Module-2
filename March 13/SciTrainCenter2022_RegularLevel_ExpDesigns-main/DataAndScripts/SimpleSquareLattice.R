#' ---
#' title:  |
#'   | Case Study  : Simple Square Lattice design
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
#' 25 soybean varieties are tested in a Simple Square Lattice
#' with block size = 5 and 2 replicates.
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

#' 
#' # ANALYSIS OF THE CASE STUDY 
#' 


#################   
#################   Analysis of simple square lattice using agricolae
################# 

lattice5x5 <- read.table("soy.txt", header=T)
str(lattice5x5)
lattice5x5$Group <- factor(lattice5x5$Group)
lattice5x5$Block <- factor(lattice5x5$Block)
lattice5x5$Treatmnt <- factor(lattice5x5$Treatmnt)
str(lattice5x5)

## see the field experiment:
lattice5x5
## what are the parameters of this design ?

attach(lattice5x5)
model1 <-PBIB.test(Block, Treatmnt, Group, Yield, k = 5,
                  method = "VC", group = TRUE, console = TRUE)

model1$means
#model1$comparison
model1$groups


## How to plot results of the simple square lattice
x11()
bar.group(model1$groups, ylim = c(0, 30),
          xlab = "treatment", ylab = "Yield")
x11()
bar.err(model1$means, ylim = c(0, 30), xlab = "treatment", ylab = "Yield")


## alternative method to estimate variances
model2 <-PBIB.test(Block, Treatmnt, Group, Yield, k = 5,
                   method="REML", group = TRUE, console = TRUE)


#' What are your conclusions ?

