################################################################################
#
#   File name: GeneCoExpr.R
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
#   Description: Script for calculating coexpression of user-defined genes across all samples or samples in user-defined group. Person correlation coefficients among samples are depicted in a form of a pair-wise heat map. NOTE: the script allowes to process gene matrix with duplicated gene IDs.
#
#   Command line use example: R --file=./GeneCoExpr.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "Genes_of_interest.txt" "Example_results/PC_GeneCoExpr" "Primary Tumour"
#
#   First arg:      Full path with name of the normalised expression matrix
#   Second arg:     Full path with name of the text file with samples annotation. The file is expected to include the following columns: sample name (1st column) and annotation (3rd column)
#   Third arg:      Full path with name of the text file with listed of genes/probes of interest
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

##### Assign colours to analysed groups
getTargetsColours <- function(targets) {
    
    ##### Predefined selection of colours for groups
    targets.colours <- c("red","blue","green","darkgoldenrod","darkred","deepskyblue", "coral", "cornflowerblue", "chartreuse4", "bisque4", "chocolate3", "cadetblue3", "darkslategrey", "lightgoldenrod4", "mediumpurple4", "orangered3")
    
    f.targets <- factor(targets)
    vec.targets <- targets.colours[1:length(levels(f.targets))]
    targets.colour <- rep(0,length(f.targets))
    for(i in 1:length(f.targets))
    targets.colour[i] <- vec.targets[ f.targets[i]==levels(f.targets)]
    
    return( list(vec.targets, targets.colour) )
}


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


##### Perform pair-wise correlation (report p-values)
cor.prob <- function(X, dfr = nrow(X) - 2) {
    R <- cor(X)
    above <- row(R) < col(R)
    r2 <- R[above]^2
    Fstat <- r2 * dfr / (1 - r2)
    R[above] <- 1 - pf(Fstat, 1, dfr)
    
    cor.mat <- t(R)
    cor.mat[upper.tri(cor.mat)] <- NA
    
    return(cor.mat)
}


##### Perform pair-wise correlation (report correlation coefficients)
cor.R <- function(X, dfr = nrow(X) - 2) {
    R <- cor(X)
    above <- row(R) > col(R)
    r2 <- R[above]^2
    Fstat <- r2 * dfr / (1 - r2)
    R[above] <- 1 - pf(Fstat, 1, dfr)
    
    cor.mat <- t(R)
    cor.mat[upper.tri(cor.mat)] <- NA
    
    return(cor.mat)
}

#===============================================================================
#    Install and load libraries
#===============================================================================
packages <- c("gplots", "plotly", "heatmaply")

# check if packages are present and install
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
   source("http://bioconductor.org/biocLite.R")
   biocLite(setdiff(packages, rownames(installed.packages())))
}
	
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


##### Make sure that at least 2 genes were selected
if ( length(genes) < 2 ) {
    
    cat("Please select at least two genes.\n")
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
    
    write(absentGenes, file = paste(coreName, "absent_genes.txt", sep = "_"), append = FALSE, sep="\t")
}

##### Report used parameters to a file
write(args, file = paste(coreName, "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#     Calculate Pearson correlation coefficients
#===============================================================================

##### Keep only genes/probes of interest
expData <- expData[ genes, ]

##### Remove gene with standard deviation = 0 (otherwise the cor.test complains)
##### Keep only genes/probes with variance > 0 across all samples (otherwise the cor.test complains)
rsd<-apply(expData,1,sd)
expData <- expData[rsd>0,]

##### Give a warning if there are any genes with no expression difference across samples)
if ( nrow(expData) != length(genes) ) {
    
    cat( "Gene(s)", names(rsd[rsd == 0]), " have no expression difference across samples!\n" )
}

    
##### Perform pair-wise correlation using expression values for user-defined genes (report correlation coefficients)
expData.cor.R <- cor.R(expData)

##### Generate heatmap (PLOTLY)
p <- heatmaply(data.frame(expData.cor.R), Rowv=NULL, Colv=NULL, colors = colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), scale="none", trace="none", hide_colorbar = FALSE, fontsize_row = 8, fontsize_col = 8) %>%
layout(autosize = TRUE, width = 800, margin = list(l=100, r=50, b=100, t=50, pad=4))

##### Save the heatmap as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, "heatmap_corr_R.html", sep = "_"))


##### Perform pair-wise correlation using expression values for user-defined genes (report p-values)
expData.cor.prob <- cor.prob(expData)
expData.cor.prob[ expData.cor.prob == 1 ] <- NA

##### Generate heatmap (PLOTLY)
p <- heatmaply(data.frame(expData.cor.prob), Rowv=NULL, Colv=NULL, colors = colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), scale="none", trace="none", hide_colorbar = FALSE, fontsize_row = 8, fontsize_col = 8) %>%
layout(autosize = TRUE, width = 800, margin = list(l=100, r=50, b=100, t=50, pad=4))

##### Save the heatmap as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, "heatmap_corr_P.html", sep = "_"))


##### Write the results into a file
write.table(prepare2write(expData.cor.R), file=paste(coreName, "corr_R.txt", sep="_"), sep="\t", row.names=FALSE)
write.table(prepare2write(expData.cor.prob), file=paste(coreName, "corr_P.txt", sep="_"), sep="\t", row.names=FALSE)


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
