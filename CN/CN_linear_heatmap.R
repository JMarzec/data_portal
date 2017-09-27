################################################################################
#
#   File name: CN_linear_heatmap.R
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
#   Description: Script performing clustering and generating heatmap for user-defined genes using relative linear copy-number values. The script imputes missing data for genes with non-missing data in >50% of samples. Genes with >50% missing data are removed prior clustering. Missing values should be denoted as 'NA'. The heatmap presents only the 500 genes with the highest copy-number variance across samples. NOTE: the script allowes to process gene matrix with duplicated gene IDs.
#
#   Command line use example: R --file=./CN_linear_heatmap.R --args "CCLE_PC_CN_processed.txt" "CCLE_PC_target.txt" "Example_results/PC_CN_linear_heatmap" "Gene_list.txt" "3" "10"
#
#   First arg:      Full path with name of the relative linear copy-number matrix
#   Second arg:     Full path with name of the text file with samples annotation. The file is expected to include the following columns: sample name (1st column) and annotation (3rd column)
#   Third arg:      Full path with name of the output folder
#   Fourth arg (OPTIONAL):   Full path with name of the text file listing genes to be considered for the heatmap. Individual genes are expected to be listed per row
#   Fifth arg (OPTIONAL):    The desired number of groups by which to color the dendrogram’s branches in the columns (samples). If this argument is provided, one is also expected to provide the sixth argument as well
#   Sixth arg (OPTIONAL):    The desired number of groups by which to color the dendrogram’s branches in the rows (genes/probes). If this argument is provided, one is also expected to provide the fifth argument as well
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


#===============================================================================
#    Load libraries
#===============================================================================

library(Amelia)
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
outFolder <- args[6]

##### Check if all required parameters were provided
if ( is.na(args[4]) || is.na(args[5]) || is.na(args[6]) ) {
    
    cat("\nError: Some arguments are missing! Please provide all required parameters.\n\n")
    q()
}

##### Read file with CN data
expData <- read.table(expFile,sep="\t", as.is=TRUE, header=TRUE, row.names=NULL)

##### Deal with the duplicated genes
expData <- duplGenes(expData)


##### Retieve the CN data file name
coreName <- strsplit(expFile, "/")
coreName <- coreName[[1]][length(coreName[[1]])]


##### Read sample annotation file
annData <- read.table(annFile,sep="\t",as.is=TRUE,header=TRUE,row.names=1)
rownames(annData) <- gsub("-", ".", rownames(annData))


##### Keep only samples with annotation info
expData <- expData[,colnames(expData) %in% rownames(annData)]
annData <- subset(annData, rownames(annData) %in% colnames(expData))

##### Make sure that the samples order in the data matrix and annotation file is the same
annData <- annData[colnames(expData),]


##### If provided, get the list of genes of interest
if ( !is.na(args[7]) ) {
    
    genesFile <- args[7]
    
    ##### Read file with the gene list
    genes <- unlist(read.table(genesFile,sep="\t",as.is=TRUE,header=FALSE))
    
} else {
    genes <- rownames(expData)
}

##### If provided, get the desired number of groups by which to color the dendrogram’s branches in the columns (samples) and rows (genes/probes)
if ( !is.na(args[8]) && !is.na(args[9]) ) {
    
    k_col <- as.numeric(args[8])
    k_row <- as.numeric(args[9])
    
} else {
    k_col <- 1
    k_row <- 1
}


##### Identify genes of interest not present in the CN matrix
absentGenes <- genes[genes %!in% rownames(expData)]


##### Keep only genes of interest
expData <- expData[ rownames(expData) %in% genes, ]


##### Set/create a directory for the output files
if (file.exists(outFolder)){
	cat( paste("The output folder already exists. Writing files to the", outFolder, "\n", sep = " ") )
} else {
	dir.create(outFolder, recursive = TRUE);
    cat( paste("Writing files to the", outFolder, "\n", sep = " ") )
}

##### Change working directory to the project workspace
setwd(outFolder)

##### Report genes of interest not present in the CN matrix
if ( length(absentGenes) > 0 ) {
    
    write(absentGenes, file = paste(coreName, "absent_genes.txt", sep = "_"), append = FALSE, sep="\t")
}

##### Report used parameters to a file
write(args, file = paste(coreName, "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#    Imputation of missing data using Ameila package
#===============================================================================

##### Remove genes with missing data in >50% of samples
genes2remove <- NULL
allowedNAs <- round(0.5*ncol(expData), digits = 0)

for ( i in 1:nrow(expData) ) {
    
    if ( sum(is.na(expData[i,])) > allowedNAs ) {
        
        genes2remove <- c(genes2remove, i)
    }
}

if ( !is.null(genes2remove) ) {
    expData <- expData[-genes2remove,]
}


##### Use amelia function to impute missing date. NOTE: the parameter 'empri' is set to the number of rows of the imput data. If error is produces by 'amelia' funciton try to increase this number
##### First check if there is any missing data
if ( any(is.na(expData)) ) {
    
    expDataImput <- amelia(x = expData, empri = nrow(expData))
    expData <- expDataImput$imputations$imp1
}

##### Keep 500 genes with the highest CN values variance across samples
rsd<-apply(expData,1,sd)

if ( nrow(expData) < 500 ) {
    
    sel<-order(rsd, decreasing=TRUE)[1:nrow(expData)]
} else {
    sel<-order(rsd, decreasing=TRUE)[1:500]
}

expData <- expData[sel,]

#===============================================================================
#     DNA copy-number heat map
#===============================================================================

##### Prepare samples annotation info
annot <- as.matrix(annData[,c(2:ncol(annData))])
names(annot) <- names(annData)[c(2:ncol(annData))]

##### Prepare colours for sample groups
targets.colour <- getTargetsColours(annData[,2])


##### Transpose matrix
expData.t <- data.frame(t(expData))

##### Cluster genes
hc <- hclust(as.dist(1-cor(expData.t, method="pearson")), method="ward.D")

##### Cluster samples
hr <- hclust(as.dist(dist(expData.t, method="euclidean")), method="ward.D")

##### Generate heatmap including the top correlated genes
#pdf(paste(coreName, "heatmap.pdf", sep = "_"), width = 10, height = ncol(expData)/5, pointsize = 12)
#heatmap.2(as.matrix(expData.t), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), RowSideColors=targets.colour[[2]], scale="col", margins = c(8, 14), labCol=FALSE, trace="none", key = TRUE, keysize = 1.5, density.info="none")
#####  Add the legend
#legend("topright", legend=levels(factor(annData[,2])),fill=targets.colour[[1]], box.col = "transparent", cex=0.8)
#dev.off()


##### Generate heatmap (PLOTLY)
p <- heatmaply(data.frame(cbind( expData.t, annot )), k_col=k_col, k_row=k_row, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), colors = colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), scale="col", trace="none", hide_colorbar = TRUE, fontsize_row = 8, fontsize_col = 8) %>%
layout(autosize = TRUE, width = 800, margin = list(l=100, r=50, b=150, t=50, pad=4), showlegend = FALSE)

##### Save the heatmap as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, "linear_heatmap.html", sep = "_"))


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
