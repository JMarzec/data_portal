################################################################################
#
#   File name: Expr_CN_Profile.R
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
#   Description: Script generating box-plots and bar-plots to visualise expression measurments across samples and groups (as indicated in target file) from normalised expression data for user-defined gene. NOTE: the script allowes to process gene matrix with duplicated gene IDs.
#
#   Command line use example: R --file=./Expr_CN_Profile.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_processed_CN.txt" "CCLE_target.txt" "Target" "KRAS" "Example_results/PC_Expr_CN_Profile"
#
#   First arg:      Full path with name of the normalised expression matrix
#   Second arg:     Full path with name of the relative linear copy-number matrix
#   Third arg:      Full path with name of the text file with samples annotation. The file is expected to include the following columns: sample name (1st column) and annotation (3rd column)
#   Forth arg:      Variable from the samples annotation file to be used for samples colouring
#   Fifth arg:      ID of gene/probe of interest
#   Six arg:        Full path with name of the output folder

#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

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

library(plotly)

#===============================================================================
#    Main
#===============================================================================

##### Catch the arguments from command line
args <- commandArgs()

expFile <- args[4]
cnFile <- args[5]
annFile <- args[6]
target <- args[7]
gene <- args[8]
outFolder <- args[9]

##### Check if all required parameters were provided
if ( is.na(args[4]) || is.na(args[5]) || is.na(args[6]) || is.na(args[7]) || is.na(args[8]) || is.na(args[9]) ) {
    
    cat("\nError: Some arguments are missing! Please provide all required parameters.\n\n")
    q()
}


##### Read file with expression data
expData <- read.table(expFile,sep="\t",as.is=TRUE,header=TRUE,row.names=NULL)

##### Deal with the duplicated genes
expData <- duplGenes(expData)


##### Read file with CN data
cnData <- read.table(cnFile,sep="\t",as.is=TRUE,header=TRUE,row.names=NULL)

##### Deal with the duplicated genes
cnData <- duplGenes(cnData)


##### Keep only samples present in both the expression and CN datasets
absentSamples.cnData <- colnames(expData)[colnames(expData) %!in% colnames(cnData)]
absentSamples.expData <- colnames(cnData)[colnames(cnData) %!in% colnames(expData)]

expData <- expData[,colnames(expData) %in% colnames(cnData)]
cnData <- cnData[,colnames(cnData) %in% colnames(expData)]


##### Make sure that the samples order in the expression and CN matrices are the same
cnData <- cnData[, colnames(expData)]


##### Retieve the expression data file name
coreName <- strsplit(expFile, "/")
coreName <- coreName[[1]][length(coreName[[1]])]


##### Read sample annotation file
annData <- read.table(annFile,sep="\t",as.is=TRUE,header=TRUE,row.names=1)
rownames(annData) <- gsub("-", ".", rownames(annData))


##### Keep only samples with annotation info
expData <- expData[,colnames(expData) %in% rownames(annData)]
cnData <- cnData[,colnames(cnData) %in% rownames(annData)]
annData <- subset(annData, rownames(annData) %in% colnames(expData))

##### Make sure that the samples order in the data matrix and annotation file is the same
annData <- annData[colnames(expData),]


##### Check if the queried genes is present in the expression data
genes <- rownames(expData)

if ( gene %!in% rownames(expData) ) {
    cat("The gene/probe", gene, "is not present in the data!", sep=" ")
    q()
    
##### ... and extract the expression of the gene of inteterest
} else {
    gene.expr <- data.matrix(expData[gene, ])
    gene.cn <- data.matrix(cnData[gene, ])
}


##### Set/create a directory for the output files
if (file.exists(outFolder)){
    cat( paste("The output folder already exists. Writing files to the", outFolder, "\n", sep = " ") )
} else {
    dir.create(outFolder, recursive = TRUE);
    cat( paste("Writing files to the", outFolder, "\n", sep = " ") )
}

##### Change working directory to the project workspace
setwd(outFolder)


##### Report samples not present in the the expression or CN matrices
if ( length(absentSamples.expData) > 0 ) {
    
    write(absentSamples.expData, file = paste(coreName, gene, "absent_in_mRNA_data.txt", sep = "_"), append = FALSE, sep="\t")
}

if ( length(absentSamples.cnData) > 0 ) {
    
    write(absentSamples.cnData, file = paste(coreName, gene, "absent_in_CN_data.txt", sep = "_"), append = FALSE, sep="\t")
}


##### Report used parameters to a file
write(args, file = paste(coreName, gene, "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#     Generate mRNA expression vs DNA copy-number scatterplot
#===============================================================================

targets <- annData[,target]
targets.colour <- getTargetsColours(targets)

##### Calculate Pearson correlation coefficient
expr_cn.corr <- round(cor.test( as.numeric(gene.expr), as.numeric(gene.cn), method = "pearson" )$estimate, digits=2)


##### Generate scatter plot
pdf(paste(coreName, gene, "mRNA_vs_CN_plot.pdf", sep = "_"), pointsize = 12)
plot(gene.cn, gene.expr, main=gene, cex=1, ylim=c(min(gene.expr) ,max(gene.expr)+1), col=targets.colour[[2]], pch=16, xlab=paste0(gene, " relative linear copy-number values"), ylab=paste0(gene, " mRNA expression") )
points(gene.cn, gene.expr, main=gene, cex=1, pch=1 )

#####  Adding Pearson correlation coefficient
legend("topright", legend=paste0("Pearson's r = ", expr_cn.corr), box.col = "transparent")

#####  Adding the legend
legend("topleft", legend=levels(factor(targets)),pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()


##### Generate scatter plot (PLOTLY)
##### Prepare data frame
gene.df <- data.frame(targets, as.numeric(gene.cn), as.numeric(gene.expr))
colnames(gene.df) <- c("Target", "CN", "mRNA")

p <- plot_ly(gene.df, x = ~CN, y = ~mRNA, color = ~Target, text=colnames(gene.expr), colors = targets.colour[[1]], type='scatter', mode = "markers", marker = list(size=10, symbol="circle"), width = 800, height = 600) %>%
layout(title = paste0("Pearson's r = ", expr_cn.corr), xaxis = list(title = paste0(gene, " relative linear copy-number values")), yaxis = list(title = paste0(gene, " mRNA expression")), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'v', y = 0.5), showlegend=TRUE)

##### Save the scatterplot as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, gene, "mRNA_vs_CN_plot.html", sep = "_"))


#===============================================================================
#     Calculate putative copy-number alterations
#===============================================================================

##### Draw histogram of the relative linear CN values
pdf(paste(coreName, gene, "hist.pdf", sep = "_"), width = 8, height = 6, pointsize = 12)
histInfo <- hist(gene.cn, breaks=5, plot=FALSE)

hist(gene.cn, col="grey", xlab="Relative linear copy-number values", main="", breaks=5, xlim=c(-3,3), labels = TRUE, ylim=c(0,(max(histInfo$counts)+1)))
dev.off()


##### Draw histogram of correlation coefficients (PLOTLY)
p <- plot_ly(x = ~as.numeric(gene.cn), type = 'histogram', width = 800, height = 500) %>%
layout(xaxis = list( title = paste0(gene, " relative linear copy-number values")), yaxis = list( title = "Frequency"), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F)

##### Save the histogram as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, gene, "hist.html", sep = "_"))


##### Assign gain for linear CN values above 0.5 and loss for linear CN values below -0.5
gene.cn[ gene.cn > 0.5 ] <- 1
gene.cn[ gene.cn < -0.5 ] <- -1
gene.cn[ gene.cn <= 0.5 & gene.cn >= -0.5 ] <- 0


#===============================================================================
#     Generate mRNA expression vs putative DNA copy-number alterations box-plot
#===============================================================================


##### Preprare dataframe
gene.df <- data.frame(targets, rep(unique(targets)[1],length(targets)), as.numeric(gene.cn), as.numeric(gene.expr))
colnames(gene.df) <- c("Target", "Box", "CN", "mRNA")


##### Generate box-plot
pdf(paste(coreName, gene, "mRNA_vs_CN_boxplot.pdf", sep = "_"), pointsize = 12 ,width = 2*length(unique(as.numeric(gene.cn))), height = 5)
#par(mar=c(3, 3, 3, 2))
boxplot(data=gene.df, mRNA~CN,col="grey", cex=0.7, ylim=c(min(gene.expr)-0.5 ,max(gene.expr)+1), xlab=paste0(gene, "  putative copy-number alterations"), ylab=paste0(gene, "  mRNA expression"), xaxt="n" )

##### Adding x-axis labels
axis(1, at=c(1,2,3), labels=c("Loss", "Diploid", "Gain"))

##### Add data points
points(factor(gene.df$CN), gene.df$mRNA, cex=1, col=targets.colour[[2]], pch=16)
points(factor(gene.df$CN), gene.df$mRNA, cex=1.1, pch=1)

#####  Adding the legend
legend("topleft", legend=levels(factor(targets)),pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()


gene.cn[ gene.cn == 1 ] <- "(1) Gain"
gene.cn[ gene.cn == -1 ] <- "(-1) Loss"
gene.cn[ gene.cn == 0 ] <- "(0) Diploid"

gene.df <- data.frame(targets, rep(unique(targets)[1],length(targets)), data.frame(t(gene.cn)), as.numeric(gene.expr))
colnames(gene.df) <- c("Target", "Box", "CN", "mRNA")

##### Generate box-plot  (PLOTLY)
p <- plot_ly(gene.df, x = ~CN, y = ~mRNA, color = ~Target, colors = targets.colour[[1]], type='scatter', mode = "markers", marker = list(size=10, symbol="circle"), width = 800, height = 600, text=colnames(gene.expr) ) %>%

add_boxplot(gene.df, x= ~CN, y= ~mRNA, color = ~Box, key=FALSE, line = list(color = "grey"), showlegend=FALSE ) %>%


layout(title = "", xaxis = list(title = paste0(gene, " relative linear copy-number values")), yaxis = list(title = paste0(gene, " mRNA expression")), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'v', y = 0.5), showlegend=TRUE)


##### Save the box-plot as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, gene, "mRNA_vs_CN_boxplot.html", sep = "_"))



##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
