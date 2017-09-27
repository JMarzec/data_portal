################################################################################
#
#   File name: ExprProfile.R
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
#   Command line use example: R --file=./ExprProfile.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "Derived_From" "KRAS" "Example_results/BC_ExprProfile"
#
#   First arg:      Full path with name of the normalised expression matrix
#   Second arg:     Full path with name of the text file with samples annotation. The file is expected to include the following columns: sample name (1st column) and annotation (3rd column)
#   Third arg:      Variable from the samples annotation file to be used for samples grouping
#   Forth arg:      ID of gene/probe of interest
#   Fifth arg:      Full path with name of the output folder

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


##### Assign colours to analysed groups (for plots filling)
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
annFile <- args[5]
target <- args[6]
gene <- args[7]
outFolder <- args[8]

##### Check if all required parameters were provided
if ( is.na(args[4]) || is.na(args[5]) || is.na(args[6]) || is.na(args[7]) || is.na(args[8]) ) {
    
    cat("\nError: Some arguments are missing! Please provide all required parameters.\n\n")
    q()
}


##### Read file with expression data
expData <- read.table(expFile,sep="\t",as.is=TRUE,header=TRUE,row.names=NULL)

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


##### Check if the queried genes is present in the expression data
genes <- rownames(expData)

if ( gene %!in% rownames(expData) ) {
    cat("The gene/probe", gene, "is not present in the data!", sep=" ")
    q()
    
##### ... and extract the expression of the gene of inteterest
} else {
    gene.expr <- expData[gene, ]
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


##### Report used parameters to a file
write(args, file = paste(coreName, gene, "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#     Generate box-plot and bar-plot
#===============================================================================

targets <- annData[,target]
targets.colour <- getTargetsColours(targets)


##### Order samples accordingly to investigated groups
dataByGroups <- NULL
targetByGroups <- NULL
colourByGroups <- NULL
expr <- list()
averegeExpr <- NULL

for (i in 1:length(unique(targets))) {
    
    ##### Separate samples accordingly to investigated groups
    expr.tmp <- gene.expr[ targets %in% unique(sort(targets))[i] ]
    averegeExpr <- c(averegeExpr, rep(mean(as.numeric(expr.tmp)), length(expr.tmp)))
    colour <- targets.colour[[2]][ targets %in% unique(sort(targets))[i] ]
    
    ##### Order samples within each group based on the expression level
    expr.tmp <- expr.tmp[order(expr.tmp)]
    colour <- colour[order(expr.tmp)]
        
    expr[[i]] <- as.numeric(expr.tmp)
    names(expr)[[i]] <- unique(sort(targets))[i]
    dataByGroups <- c(dataByGroups, expr.tmp)
    targetByGroups <- c(targetByGroups, targets[ targets %in% unique(sort(targets))[i] ])
    colourByGroups <- c(colourByGroups, colour)
}

dataByGroups <- unlist(dataByGroups)
    
##### Generate box-plot
pdf(paste(coreName, gene, "expression_profile.pdf", sep = "_"), pointsize = 12 ,width = 2*length(unique(targets)), height = 4)
boxplot(expr,col="grey", cex=0.7, ylim=c(min(gene.expr)-1 ,max(gene.expr)+1), ylab=paste0(gene, "  mRNA expression") )

##### Add data points
for (i in 1:length(expr)) {
        
    points(rep(i, length(expr[[i]])), expr[[i]], cex=0.7)
}

##### Generate bar-plot
df.bar <- barplot(dataByGroups,col=colourByGroups, las = 2, ylim=c(min(gene.expr)-0.2 ,max(gene.expr)+2), ylab=paste0(gene, "  mRNA expression"), xpd = FALSE, axisnames = FALSE)

##### Add mean expression for each group
lines(x = df.bar, y = averegeExpr, lwd=5, col="white")
lines(x = df.bar, y = averegeExpr, lwd=2)

legend("topleft", legend = unique(sort(targets)), fill = unique(colourByGroups), bty = "n", horiz=TRUE, cex=1)
dev.off()


##### Generate box-plot  (PLOTLY)
##### Prepare data frame
gene.expr.df <- data.frame(targets, as.numeric(gene.expr))
colnames(gene.expr.df) <- c("Group", "Expression")


p <- plot_ly(gene.expr.df, y= ~Expression, color = ~Group, type = 'box', jitter = 0.3, pointpos = 0, boxpoints = 'all', marker = list(color = colourByGroups), line = list(color = unique(colourByGroups)), width = 800, height = 600) %>%
layout(yaxis = list( title = paste0(gene, "  mRNA expression")), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'v', y = 0.5), showlegend=TRUE)

##### Save the box-plot as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, gene, "expression_profile_boxplot.html", sep = "_"))


##### Generate bar-plot (PLOTLY)
##### Prepare data frame
dataByGroups.df <- data.frame(targetByGroups, names(dataByGroups), as.numeric(dataByGroups))
colnames(dataByGroups.df) <- c("Group","Sample", "Expression")

##### The default order will be alphabetized unless specified as below
dataByGroups.df$Sample <- factor(dataByGroups.df$Sample, levels = dataByGroups.df[["Sample"]])

p <- plot_ly(dataByGroups.df, x = ~Sample, y = ~Expression, color = ~Group, type = 'bar',  marker = list(color = c(colourByGroups)), width = 800, height = 400) %>%
layout(title = "", xaxis = list(title = ""), yaxis = list(title = paste0(gene, "  mRNA expression")), margin = list(l=50, r=50, b=100, t=50, pad=4), autosize = F, legend = list(orientation = 'v', y = 0.5), showlegend=TRUE)

##### Save the bar-plot as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, gene, "expression_profile_barplot.html", sep = "_"))


##### Write expression data into a file
cat(paste("Writing expression data to", paste(coreName, gene, "expression_profile.exp", sep = "_")))
write.table(t(dataByGroups), file=paste(coreName, gene, "expression_profile.exp", sep = "_"),sep="\t", row.names=FALSE)
    


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
