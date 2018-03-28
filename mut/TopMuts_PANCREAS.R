### HISTORY ##############################################################################
# Version			Date					Coder						Comments
# 1.0 				16/06/2017				Ema & Jacek & Stefano
#
#   Contact: Emanuela, Jacek, Stefano (e.gadaleta@qmul.ac.uk; j.marzec@qmul.ac.uk, s.pirro@qmul.ac.uk)
#
#   Barts Cancer Institute
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
### DESCRIPTION ##########################################################################
#  Read in the mutation and clinical data downloaded previously from TCGA and plot top/defined mutations
## https://cancergenome.nih.gov/; https://ocg.cancer.gov/programs/target/research
## http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html


#========================================================================================#
# Clear workspace
rm(list=ls());
# Close any open graphics devices
graphics.off();

#========================================================================================#
#    									FUNCTIONS										 #
#========================================================================================#


#========================================================================================#
#    						INSTALL AND LOAD LIBRARIES									 #
#========================================================================================#
# install libraries if not present
TCGA.libs <- c("TCGAbiolinks","SummarizedExperiment");

if (length(setdiff(TCGA.libs, rownames(installed.packages()))) > 0) {
   source("http://bioconductor.org/biocLite.R")
   biocLite(setdiff(TCGA.libs, rownames(installed.packages())))
}

# load libraries
library("TCGAbiolinks")
library("SummarizedExperiment")


#========================================================================================#
#    										MAIN						  				 #
#========================================================================================#
# define arguments
args <- commandArgs();

mutFile <- args[4];
annFile <- args[5];
outFolder <- args[6];
genes <- args[7]; 
 
# VERIFY THE DATA INPUTS
# verify the input parameters
if ( is.na(args[4]) || is.na(args[5]) || is.na(args[6]) || is.na(args[7]) ) {
    cat("\nError: Some arguments are missing! Please provide all required parameters.\n\n");
    q();
}


# READ IN MUTATION AND CLINICAL DATA
annData <- read.table(file = annFile, row.names = 1, header = T, sep = "\t");
mutData <- read.table(file = mutFile, row.names = NULL, header = T, sep = "\t");


# SET THE WORK ENVIRONMENT
# create the directory
if (file.exists(outFolder)){
	cat( paste("\nThis output folder already exists. Files will be written to", outFolder, "\n", sep=" ") );
} else {
	dir.create(outFolder, recursive = TRUE);
    cat( paste("Writing files to ", outFolder, "\n", sep = " ") );
}


# KEEP TRACK OF PARAMETERS
write(args, file = paste(outFolder, "R_parameters.txt", sep=""), append = FALSE, sep="\t");


#========================================================================================#
#    							PREPARE DATA INPUTS										 #
#========================================================================================#
# make life easy
uniqGenes <- unique(mutData$Hugo_Symbol);

# create input data matrix
allMut <- matrix(data=0, nrow=length(uniqGenes), ncol=1, dimnames=list(rownames=uniqGenes, colnames="no_reports"));

# populate the matrix
for(gene in uniqGenes) {

	# identify position of each gene
	genePos <- which(mutData$Hugo_Symbol == gene);
	# populate the matrix with the number of times gene is reported
	allMut[gene, "no_reports"] <- length(genePos);

}

# Order by most reported gene
allMut <- allMut[ order(allMut[, "no_reports"], decreasing=TRUE),];

# identify top mutated genes
topGenes <- names(allMut[1:genes]);

allPos <- NULL;
for( gene in topGenes) {
	genePos <- grep(paste("\\b", gene, "\\b", sep=""), mutData$Hugo_Symbol);
	allPos <- c(allPos, genePos);
}


#========================================================================================#
#    							PLOTTING R OBJECTS										 #
#========================================================================================#
# MUTATION DATA FOR PLOT
#subset mutations to top genes (see args[7])
mutData.top <- mutData[ allPos , ];
topMut.genes <- unique(mutData.top$Hugo_Symbol);

# CLINICAL DATA FOR PLOT
# Subset clinical data to covariates of interest
annCovar <- c("bcr_patient_barcode","stage","gender","race", "surv.stat", "alcohol_history");
annData.slimmed <- annData[ , annCovar ];

# GENERATE ONCOPRINT
TCGAvisualize_oncoprint(
	mut = mutData.top,
	genes = topMut.genes,
	filename = paste(outFolder, "oncoprint_BRCA.pdf", sep=""),
	annotation = annData.slimmed,
	color=c("background"="#CCCCCC","DEL"="purple", "INS"="yellow","SNP"="brown"),
	rows.font.size= 8,
	width = 5,
	heatmap.legend.side = "right",
	dist.col = 0,
	label.font.size = 6
);
