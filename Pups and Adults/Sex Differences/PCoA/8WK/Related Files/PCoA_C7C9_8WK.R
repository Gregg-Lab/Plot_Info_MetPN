#TITLE
###PCoA Sex Difference Plot

#DESCRIPTION:
###Plots a PCoA plot separating experimental groups further by sex.


#AUTHOR:
###Date Created: 06-04-2019
###Date Last Modified:  06-10-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

library(RColorBrewer) #gives color schemes
library(ggplot2) #data plotter
library(gplots) #data plotter
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(stringr) #loaded to use str_replace_all() which removes all special characters in a string
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)
library(tidyverse) #ggplot2, tibble, purrr, tidyr, readr, dplyr, stringr, forcats
library(ggpubr)

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

# VARIABLES #####################################################################
WORKING.DIRECTORY <- "~/Downloads/Microbiome-June2019/Cohort_7_9/8wk_SexDifference_PCoA"
TAXONOMY.FILE <- "combined_pn.final.0.03.cons.taxonomy"
CORR.FILE <- "combined_pn.final.0.03.pick.0.03.filter.spearman.corr.axes"
PCOA.FILE <- "combined_pn.final.0.03.pick.thetayc.0.03.lt.pcoa.axes"
DESIGN.FILE <- "combined_pn_8wk.design.txt"
SHARED.FILE <- "combined_pn.final.8wk.shared"
TITLE <- "C7C9 8WK Sex Difference PCoA"

#################################################################################


setwd(WORKING.DIRECTORY)

##subset to get significant OTUs
##import data
tax <- read.table(file=TAXONOMY.FILE, row.names = 1,
                  header=TRUE,
                  check.names=FALSE,
                  comment.char="")
shared <- read.table(file=SHARED.FILE, row.names = 2, header = TRUE)
design <- read.table(file=DESIGN.FILE, header = TRUE)

#clean data
tax_clean <- cleanTaxonomyKeepUnclassified(tax)
tax_clean$OTU <- rownames(tax_clean)
tax_clean_unclassified <- cleanUnclassified(tax_clean, column = "Genus")
shared <- cleanShared(shared, design, tax_clean_unclassified)
design <- cleanDesign(design, shared)

corr.axes <- read.table(file=CORR.FILE, header = TRUE)

corr.tax <- merge(corr.axes, tax_clean_unclassified, by = "OTU", all.y = TRUE)
corr.tax.1 <- subset(corr.tax, p.value < 0.001)
dim(corr.tax.1) #22 OTUs
corr.tax.2 <- subset(corr.tax.1, p.value.1 < 0.001)
dim(corr.tax.2) #1 OTUs

##make plot
pcoa.axes <- read.table(file=PCOA.FILE, header = TRUE)
pcoa.pch <- c(2, 1, 17, 16) #2 is open triange, 1 is open circle, 17 is filled triangle, 16 is filled circle
pcoa.pre.final <- merge(design, pcoa.axes, by.x = "Sample", by.y = "group")
pcoa.final <- pcoa.pre.final[order(pcoa.pre.final$Group),]
shape <- data.frame(c(rep(2, 4), rep(1, 6), rep(17, 4), rep(16, 7)))  #labels males as circles, females as triangles, controls as open, and mets as filled
colnames(shape) <- "shape" #maintains column name of "shape"
pcoa.final.shape <- cbind(shape, pcoa.final) #combines shape labels with data table
#write.csv(pcoa.axes.final, file = "pcoa.axes.final.csv") #exports csv file to run in mothur!


#p19 plot (ALL DATA KEPT)
plot(pcoa.final.shape$axis1, 
     pcoa.final.shape$axis2, 
     pch = pcoa.final.shape$shape, #selects shape from shape column
     col = "black",  #selects color from group column
     main = "C7C9 8WK Sex Difference PCoA", 
     xlab ="Axis 1 (34.44%)", 
     ylab = "Axis 2 (29.13%)", 
     ylim = c(-.8, .8),  #sets limits from -60% to 60%
     xlim = c(-.8, .8), #sets limits from -60% to 60%
     cex = 1.3)

#add legend
legend(x = -0.3, #places legend and (0.33, .4)
       y = .7, 
       legend = c("Female Ctrl PN", "Male Ctrl PN", "Female Met PN", "Male Met PN"), 
       pch = pcoa.pch,
       col = c("black"),
       cex = 0.8)

#add ellipse
dataEllipse(pcoa.pre.final$axis1, 
            pcoa.pre.final$axis2,
            center.pch = FALSE, #doesn't plot center ellipse point
            add = TRUE, #adds to existing plot
            ellipse.label = "yes", #doesn't have ellipse label
            plot.points = FALSE, #doesn't plot points again
            groups = as.factor(pcoa.pre.final$shape), #separates points by group
            group.labels = FALSE, #removes ellipse group labels
            col = c("blue", "red"),
            lwd = 0.1, #change ellipse line thickness
            fill = TRUE, #fills ellipse with col color
            fill.alpha = 0.3, #makes ellipse filled color transculent
            levels=c(0.975), #ellipse 97.5% confidence level
            cex = 0.6,
            labels = FALSE)
