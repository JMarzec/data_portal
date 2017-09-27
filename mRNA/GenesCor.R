################################################################################
#
#   File name: GenesCor.R
#
#   Authors: Jacek Marzec ( j.marzec@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

################################################################################
#
#   Description: Script for calculating coexpression between user-defined genes/probes. Person correlation coefficient is used to measure correaltion between genes' expression. The pairwise correlation between expression of user-defined genes is depicted in a form of a correlation matrix heat map. NOTE: the script allowes to process gene matrix with duplicated gene IDs. It allows to process up to 50 genes
#
#   Command line use example: R --file=./GenesCor.R --args "TCGA_PAAD_normalized.txt" "TCGA_PAAD_target.txt" "Genes_of_interest.txt" "Example_results/PC_GenesCor" "PDAC"
#
#   First arg:      Full path with name of the normalised expression matrix
#   Second arg:     Full path with name of the text file with samples annotation. The file is expected to include the following columns: sample name (1st column) and annotation (3rd column)
#   Third arg:      ID of gene/probe of interest
#   Fourth arg:     Full path with name of the output folder
#   Fifth arg (OPTIONAL):  Samples group to use for the analysis
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

##### Prepare object to write into a file
prepare2write <- function (x) {

	x2write <- cbind(rownames(x), x)
    colnames(x2write) <- c("",colnames(x))
	return(x2write)
}

##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0


##### Deal with the duplicated genes
duplGenes <- function(expData) {

    genesList <- NULL
    genesRepl <- NULL

    for ( i in 1:nrow(expData) ) {

        geneName <- expData[i,1]

        ##### Distingish duplicated genes by adding duplicate number
        if ( geneName %in% genesList ) {

            ##### Report genes with more than one duplicates
            if ( geneName %in% names(genesRepl) ) {

                genesRepl[[ geneName ]] = genesRepl[[ geneName ]]+1

                geneName <- paste(geneName, ".", genesRepl[[ geneName ]], sep="")

            } else {
                genesRepl[[ geneName ]] <- 2

                geneName <- paste(geneName, ".2", sep="")
            }
        }
        genesList <- c(genesList,geneName)
    }

    rownames(expData) <- genesList

    ##### Remove the first column with gene names, which now are used as row names
    expData <- expData[, -1]

    return(expData)
}


#===============================================================================
#    Load libraries
#===============================================================================

library(gplots)
library(plotly)
library(heatmaply)

#===============================================================================
#    Main
#===============================================================================

##### Catch the arguments from command line
args <- commandArgs()

expFile <- args[4]
annFile <- args[5]
genesFile <- args[6]
outFolder <- args[7]


##### Check if all required parameters were provided
if ( is.na(args[4]) || is.na(args[5]) || is.na(args[6]) || is.na(args[7]) ) {

    cat("\nError: Some arguments are missing! Please provide all required parameters.\n\n")
    q()
}


##### Read file with expression data
expData <- read.table(expFile,sep="\t", as.is=TRUE, header=TRUE, row.names=NULL)

##### Deal with the duplicated genes
expData <- duplGenes(expData)


##### Retieve the expression data file name
coreName <- strsplit(expFile, "/")
coreName <- coreName[[1]][length(coreName[[1]])]


##### Read sample annotation file
annData <- read.table(annFile,sep="\t",as.is=TRUE,header=TRUE,row.names=1)
rownames(annData) <- gsub("-", ".", rownames(annData))


##### Read file with the gene list
genes <- unlist(read.table(genesFile,sep="\t",as.is=TRUE,header=FALSE))


##### Allow to process up to 50 genes
if ( length(genes) > 50 ) {

    cat(paste("More than 50 genes were selected! Please make sure that you select up to 50 genes.\n"))
    q()
}
##### Make sure that there are more than 1 user-defined gene
if ( length(genes) < 2 ) {
	cat(paste("Please provide at least 2 genes.\n"))
	q()
}


##### Keep only samples with annotation info
expData <- expData[,colnames(expData) %in% rownames(annData)]
annData <- subset(annData, rownames(annData) %in% colnames(expData))

##### Make sure that the samples order in the data matrix and annotation file is the same
annData <- annData[colnames(expData),]


##### If provided, get the data for the samples of defined group
if ( !is.na(args[8]) ) {

    group <- args[8]

    ##### Check if expresison data is available for that group
    if ( group %in% annData[,2] ) {

        ##### Keep only data for the samples from the group of interest
        expData <- expData[, annData[,2] %in% group ]
        annData <- annData[ annData[,2] %in% group, ]

    } else {
        cat("There is no data available for the",group , "group!\n")
        q()
    }
}


##### Identify genes of interest not present in the expression matrix
absentGenes <- genes[genes %!in% rownames(expData)]


##### Set/create a directory for the output files
if (file.exists(outFolder)){
	cat( paste("The output folder already exists. Writing files to the", outFolder, "\n", sep = " ") )
} else {
	dir.create(outFolder, recursive = TRUE);
    cat( paste("Writing files to the", outFolder, "\n", sep = " ") )
}

##### Change working directory to the project workspace
setwd(outFolder)

##### Report genes of interest not present in the expression matrix
if ( length(absentGenes) > 0 ) {

    write(absentGenes, file = paste(coreName, gene, "absent_genes.txt", sep = "_"), append = FALSE, sep="\t")
}

##### Report used parameters to a file
write(args, file = paste(coreName, "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#     Calculate Pearson correlation coefficients
#===============================================================================

##### Keep only genes of interest
expData <- expData[ rownames(expData) %in% genes, ]

##### Remove gense with standard deviation = 0 (otherwise the cor.test complains)
#####  Keep only genes/probes with variance > 0 across all samples (otherwise the cor.test complains)
rsd<-apply(expData,1,sd)
expData <- expData[rsd>0,]


##### Check if used-defined genes are present in the data
for ( i in 1:length(genes) ) {
	if ( genes[i] %!in% rownames(expData) ) {
    	cat("The genes/probes", genes[i], "are not present in the data!", sep=" ")

			genes <- genes[-i]
		}
}


##### Calculate Pearson correlation coefficient for user-defined genes
corr.res <- matrix(data = 0, nrow = nrow(expData), ncol = (2*length(genes))+1, dimnames = list( rownames(expData), c( "Gene", paste(rep(genes, each = 2), c("Correlation", "P-value")) ) ))

for ( i in 1:length(genes) ) {
	for ( j in 1:nrow(expData) ) {

    #### Pearson correlation coefficient and test P
		corr.res[j, paste0(genes[i], " Correlation")] <- cor.test( as.numeric(expData[ genes[i] ,]), as.numeric(expData[ j, ]), method = "pearson" )$estimate;
		corr.res[j, paste0(genes[i], " P-value")] <- cor.test( as.numeric(expData[ genes[i] ,]), as.numeric(expData[ j, ]), method = "pearson" )$p.value;
	}
}

corr.res <- data.frame(corr.res)
corr.res[,"Gene"] <- rownames(expData)
colnames(corr.res) <- c( "Gene", paste(rep(genes, each = 2), c("Correlation", "P-value")) )

##### Write the results into a file
write.table(corr.res, file=paste(coreName, "corr.txt", sep="_"), sep="\t", row.names=FALSE)


#===============================================================================
#     Pairwise correlation heat map for defined genes
#===============================================================================

##### Extract the the correlation results for the user-defined genes
corr.res.genes <- corr.res[genes, paste(genes, "Correlation")]
colnames(corr.res.genes) <- rownames(corr.res.genes)
corr.res.genes[upper.tri(corr.res.genes)] <- NA


##### Generate heatmap including the top correlated genes (PLOTLY)
p <- heatmaply(data.frame(corr.res.genes), dendrogram="none", colors = colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), scale="none", trace="none", hide_colorbar = FALSE, fontsize_row = 8, fontsize_col = 8) %>%
layout(autosize = TRUE, width = 800, margin = list(l=150, r=50, b=150, t=50, pad=4), showlegend = FALSE)

##### Save the heatmap as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, "corr_heatmap.html", sep = "_"))


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
