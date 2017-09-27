################################################################################
#
#   File name: GeneCoExprHeatmap.R
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
#   Description: Script for calculating coexpression between user-defined gene/probe and all other genes/probe. Person correlation coefficient is used to measure correaltion between genes' expression. The top correlated genes' expression is depicted in a form of a heat map. NOTE: the script allowes to process gene matrix with duplicated gene IDs. It allows to process up to 40 genes
#
#   Command line use example: R --file=./GeneCoExprHeatmap.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "KRAS" "20" "Example_results/PC_GeneCoExprHeatmap" "Gene_list.txt"
#
#   First arg:      Full path with name of the normalised expression matrix
#   Second arg:     Full path with name of the text file with samples annotation. The file is expected to include the following columns: sample name (1st column) and annotation (3rd column)
#   Third arg:      ID of gene/probe of interest
#   Fourth arg:     Number of most correlated genes to be displayed in the form of heat map
#   Fifth arg:      Full path with name of the output folder
#   Sixth arg (OPTIONAL):  Full path with name of the text file listing genes to be considered for the analysis. Individual genes are expected to be listed per row
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
gene <- args[6]
topGenesNo <- as.numeric(args[7])

##### Allow to process up to 40 genes
if ( topGenesNo > 40 ) {
    cat(paste("More than 40 genes were selected! Please make sure that you select up to 40 genes."))
    q()
}


outFolder <- args[8]


##### Check if all required parameters were provided
if ( is.na(args[4]) || is.na(args[5]) || is.na(args[6]) || is.na(args[7])  || is.na(args[8])  ) {

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


##### Keep only samples with annotation info
expData <- expData[,colnames(expData) %in% rownames(annData)]
annData <- subset(annData, rownames(annData) %in% colnames(expData))

##### Make sure that the samples order in the data matrix and annotation file is the same
annData <- annData[colnames(expData),]


##### If provided, get the list of genes of interest
if ( !is.na(args[9]) ) {

    genesFile <- args[9]

    ##### Read file with the gene list
    genes <- unlist(read.table(genesFile,sep="\t",as.is=TRUE,header=FALSE))

} else {
    genes <- rownames(expData)
}


##### Check if the queried genes is present in the expression data
if ( gene %!in% rownames(expData) ) {
    cat("The gene/probe", gene, "is not present in the data!", sep=" ")
    q()
}


##### Extract the expression of the gene of inteterest
gene.expr <- expData[gene, ]

##### ... and remove it from the expression matrix
expData <- expData[rownames(expData) %!in% gene, ]
genes <- genes[genes %!in% gene]

##### Identify genes of interest not present in the expression matrix
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

##### Report genes of interest not present in the expression matrix
if ( length(absentGenes) > 0 ) {

    write(absentGenes, file = paste(coreName, gene, "absent_genes.txt", sep = "_"), append = FALSE, sep="\t")
}

##### Report used parameters to a file
write(args, file = paste(coreName, gene, "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#     Calculate Pearson correlation coefficients
#===============================================================================

##### Remove gene with standard deviation = 0 (otherwise the cor.test complains)
#####  Keep only genes/probes with variance > 0 across all samples (otherwise the cor.test complains)
rsd<-apply(expData,1,sd)
expData <- expData[rsd>0,]


##### Get Pearson correlation coefficient for each gene
corr.res <- matrix(data = 0, nrow = nrow(expData), ncol = 3, dimnames = list( rownames(expData), c("Gene","Correlation", "P-value")))

for ( i in 1:nrow(expData) ) {

    #### Pearson correlation coefficient and test P
	corr.res[i, "Correlation"] <- cor.test( as.numeric(gene.expr), as.numeric(expData[ i, ]), method = "pearson" )$estimate;
	corr.res[i, "P-value"] <- cor.test( as.numeric(gene.expr), as.numeric(expData[ i, ]), method = "pearson" )$p.value;
}

corr.res <- data.frame(corr.res)

corr.res[,"Gene"] <- rownames(expData)
corr.res[,"Correlation"] <- as.numeric(corr.res[,"Correlation"])

##### Order results by decreasing correlation coefficient (start with the most positively correlated genes)
corr.sorted.dec <- order(corr.res[,"Correlation"], decreasing = TRUE)
corr.res.pos <- corr.res[corr.sorted.dec, ]

##### Order results by increasing correlation coefficient (start with the most negatively correlated genes)
corr.sorted.inc <- order(corr.res[,"Correlation"], decreasing = FALSE)
corr.res.neg <- corr.res[corr.sorted.inc, ]

##### Write the results into a file
write.table(corr.res.pos, file=paste(coreName, gene, "corr.txt", sep="_"), sep="\t", row.names=FALSE)


##### Draw histogram of correlation coefficients
pdf(paste(coreName, gene, "corr_hist.pdf", sep = "_"), width = 8, height = 6, pointsize = 12)
histInfo <- hist(as.numeric(corr.res.pos[, "Correlation"]), breaks=30, plot=FALSE)
hist(as.numeric(corr.res.pos[, "Correlation"]), col="grey", xlab="Pearson correlation coeffcient", main="", breaks=30, xlim=c(-1,1), labels = TRUE, ylim=c(0,(max(histInfo$counts)+2)))
dev.off()


##### Draw histogram of correlation coefficients (PLOTLY)
p <- plot_ly(x = ~ corr.res.pos[, "Correlation"], type = 'histogram', nbinsx = 50, width = 800, height = 500) %>%
layout(xaxis = list( title = "Pearson correlation coeffcient"), yaxis = list( title = "Frequency"), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F)

##### Save the histogram as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, gene, "corr_hist.html", sep = "_"))


#===============================================================================
#     Expression heat map
#===============================================================================

##### Get expression data for user-defined top correlated genes
topGenesCor.pos <- expData[rownames(corr.res.pos)[1:topGenesNo],]

##### Sort samples based on increasing expression of the gene of interest
gene.expr.order <- order(gene.expr, decreasing = FALSE)
topGenesCor.pos <- topGenesCor.pos[,gene.expr.order]

##### Merge the expression of gene of interest with the expression matrix
topGenesCor.pos <- rbind(gene.expr[gene.expr.order], topGenesCor.pos)

##### Get expression data for user-defined top correlated genes
topGenesCor.neg <- expData[rownames(corr.res.neg)[1:topGenesNo],]

##### Sort samples based on increasing expression of the gene of interest
topGenesCor.neg <- topGenesCor.neg[,gene.expr.order]

##### Merge the expression of gene of interest with the expression matrix
topGenesCor.neg <- rbind(gene.expr[gene.expr.order], topGenesCor.neg)


##### Prepare colours for sample groups
##### Prepare samples annotation info
annData <- annData[gene.expr.order,]
annot <- as.matrix(annData[,c(2:ncol(annData))])
colnames(annot) <- names(annData)[c(2:ncol(annData))]

targets.colour <- getTargetsColours(annData[,2])


##### Transpose matrix
topGenesCor.pos <- data.frame(t(topGenesCor.pos))
topGenesCor.neg <- data.frame(t(topGenesCor.neg))


##### Generate heatmap including the top correlated genes
#pdf(paste(coreName, gene, "heatmap.pdf", sep = "_"), width = topGenesNo/6, height = ncol(expData)/5, pointsize = 12)
##### Report top positively correlated genes
#heatmap.2(as.matrix(topGenesCor.pos), Rowv=FALSE, Colv=FALSE, dendrogram="none", col=colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), RowSideColors=targets.colour[[2]], scale="col",margins = c(14, 14), labCol=paste(colnames(topGenesCor.pos), "  (", round(as.numeric(c(1,corr.res.pos[ c(1:topGenesNo), "Correlation"])), digits=2), ")", sep=""), labRow=rownames(topGenesCor.pos), trace="none", key = TRUE, keysize = 1.5, density.info="none")
#####  Add the legend
#legend("topright", legend=levels(factor(annData[,2])),fill=targets.colour[[1]], box.col = "transparent", cex=0.8)

##### Report top negatively correlated genes
#heatmap.2(as.matrix(topGenesCor.neg), Rowv=FALSE, Colv=FALSE, dendrogram="none", col=colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), RowSideColors=targets.colour[[2]], scale="col",margins = c(14, 14), labCol=paste(colnames(topGenesCor.neg), "  (", round(as.numeric(c(1,corr.res.neg[ c(1:topGenesNo), "Correlation"])), digits=2), ")", sep=""), labRow=rownames(topGenesCor.neg), trace="none", key = TRUE, keysize = 1.5, density.info="none")
#####  Add the legend
#legend("topright", legend=levels(factor(annData[,2])),fill=targets.colour[[1]], box.col = "transparent", cex=0.8)
#dev.off()


##### Generate heatmap including the top correlated genes (PLOTLY)

##### Report top positively correlated genes
p <- heatmaply(data.frame(cbind( topGenesCor.pos, annot )), dendrogram="none", colors = colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), scale="col", labCol=paste(colnames(topGenesCor.pos), "  (", round(as.numeric(c(1,corr.res.pos[ c(1:topGenesNo), "Correlation"])), digits=2), ")", sep=""), trace="none", hide_colorbar = TRUE, fontsize_row = 8, fontsize_col = 8) %>%
layout(autosize = TRUE, width = 800, margin = list(l=150, r=50, b=150, t=50, pad=4), showlegend = FALSE)

##### Save the heatmap as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, gene, "heatmap_pos.html", sep = "_"))


##### Report top negatively correlated genes
p <- heatmaply(data.frame(cbind( topGenesCor.neg, annot )), dendrogram="none", colors = colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), scale="col", labCol=paste(colnames(topGenesCor.neg), "  (", round(as.numeric(c(1,corr.res.neg[ c(1:topGenesNo), "Correlation"])), digits=2), ")", sep=""), trace="none", hide_colorbar = TRUE, fontsize_row = 8, fontsize_col = 8) %>%
layout(autosize = TRUE, width = 800, margin = list(l=150, r=50, b=150, t=50, pad=4), showlegend = FALSE)

##### Save the heatmap as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, gene, "heatmap_neg.html", sep = "_"))


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
