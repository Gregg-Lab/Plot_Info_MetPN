#TITLE
###Sex Difference PCoA Plot
###Adapted PCoA plot that labels males and females separately.

#DESCRIPTION:
###Plots PCoA for each group and each sex.
###Removes samples without 1000 reads, OTUs under 200 reads, and all "_unclassified" OTUs at the "OTU" level.


#AUTHOR:
###Date Created: 5-13-2019
###Date Last Modified:  06-10-2019
###Modifications by Zach Carlson, BS
###Original Code by Anne Seekatz, Rcode_erinsubset/erinsubset_Fig4.pcoa.R, viewed on her GitHub account (Dated 2017)
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT
library(car)
library(shape)

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

# VARIABLES #####################################################################
WORKING.DIRECTORY <- "~/Downloads/Microbiome-May2019/Cohort_7_9/P19_SexDifference_PCoA"
TAXONOMY.FILE <- "combined_pn.final.0.03.cons.taxonomy"
CORR.FILE <- "combined_pn.final.0.03.pick.0.03.filter.spearman.corr.axes"
PCOA.FILE <- "combined_pn.final.0.03.pick.thetayc.0.03.lt.pcoa.axes"
DESIGN.FILE <- "combined_pn_P19.design.txt"
SHARED.FILE <- "combined_pn.final.p19.shared"
TITLE <- "C7C9 P19 Sex Difference PCoA"

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
dim(corr.tax.1) #20 OTUs
corr.tax.2 <- subset(corr.tax.1, p.value.1 < 0.001)
dim(corr.tax.2) #2 OTUs

##make plot
pcoa.axes <- read.table(file=PCOA.FILE, header = TRUE)
pcoa.pch <- c(2, 1, 17, 16) #2 is open triange, 1 is open circle, 17 is filled triangle, 16 is filled circle
pcoa.pre.final <- merge(design, pcoa.axes, by.x = "Sample", by.y = "group")
pcoa.final <- pcoa.pre.final[order(pcoa.pre.final$Group),]
shape <- data.frame(c(rep(2, 2), rep(1, 8), rep(17, 4), rep(16, 4)))  #labels males as circles, females as triangles, controls as open, and mets as filled
colnames(shape) <- "shape" #maintains column name of "shape"
pcoa.final.shape <- cbind(shape, pcoa.final) #combines shape labels with data table
#write.csv(pcoa.axes.final, file = "pcoa.axes.final.csv") #exports csv file to run in mothur!


#p19 plot (ALL DATA KEPT)
plot(pcoa.final.shape$axis1, 
     pcoa.final.shape$axis2, 
     pch = pcoa.final.shape$shape, #selects shape from shape column
     col = "black",  #selects color from group column
     main = "C7C9 P19 Sex Difference PCoA", 
     xlab ="Axis 1 (33.00%)", 
     ylab = "Axis 2 (22.96%)", 
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
