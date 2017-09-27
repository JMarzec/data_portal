################################################################################
#
#   File name: PCA.R
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
#   Description: Script performing principal component analysis using normalised expression data. It generates scree plot, as well as 2- and 3-dimensional plots of user-defined principal component. NOTE: the script allowes to process gene matrix with duplicated gene IDs.
#
#   Command line use example: R --file=./PCA.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "Derived_From" "1" "Example_results/PC_PCA"
#
#   First arg:      Full path with name of the normalised expression matrix
#   Second arg:     Full path with name of the text file with samples annotation. The file is expected to include the following columns: sample name (1st column) and annotation (3rd column)
#   Third arg:      Variable from the samples annotation file to be used for samples colouring
#   Forth arg:      The principal component to be plotted together with the two subsequent most prevalent principal components
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
annFile <- args[5]
target <- args[6]
PC1 <- as.numeric(args[7])
PC2 <- PC1 + 1
PC3 <- PC1 + 2
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
write(args, file = paste(coreName, "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#     Principal components analysis
#===============================================================================

##### Assign colours according to defined sample annotation
targets <- annData[,target]
targets.colour <- getTargetsColours(targets)


##### Keep only probes with variance > 0 across all samples
rsd<-apply(expData,1,sd)
expData <- expData[rsd>0,]

##### Perform principal components analysis
expData_pca <- prcomp(t(expData), scale=FALSE)


##### Generate scree plot
pdf(paste(coreName, "PCA_screeplot.pdf", sep = "_"), width = 14, height = 6)
screeplot(expData_pca, npcs = length(expData_pca$sdev), main="The variances captured by principal components")
dev.off()


##### Get variance importance for all principal components
importance_pca <- summary(expData_pca)$importance[2,]
importance_pca <- paste(round(100*importance_pca, 2), "%", sep="")


##### Generate PCA plot
pdf(paste(coreName, "PCA.pdf", sep = "_"), pointsize = 12)
plot(expData_pca$x[,PC1], expData_pca$x[,PC2], type="p", xlab = paste("PC ", PC1, " (",importance_pca[PC1],")",sep=""), ylab = paste("PC ", PC2, " (",importance_pca[PC2],")",sep=""), ylim=c( min(expData_pca$x[,PC2]),max(expData_pca$x[,PC2]) + (abs(min(expData_pca$x[,PC2]))+abs(max(expData_pca$x[,PC2])))/4 ), pch=16, col=targets.colour[[2]])

##### Adding colours corresponding to defined sample annotation
points(expData_pca$x[,PC1], expData_pca$x[,PC2], type="p", pch=1)

#####  Adding the legend
legend("topleft", legend=levels(factor(targets)),pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()



##### Generate bar-plot (PLOTLY)
##### Prepare data frame
expData_pca.df <- data.frame(paste0("PC ", c(1:length(expData_pca$sdev))), expData_pca$sdev)
colnames(expData_pca.df) <- c("PC", "Variances")

##### The default order will be alphabetized unless specified as below
expData_pca.df$PC <- factor(expData_pca.df$PC, levels = expData_pca.df[["PC"]])

p <- plot_ly(expData_pca.df, x = ~PC, y = ~Variances, type = 'bar', width = 800, height = 600) %>%
layout(title = "The variances captured by principal components", xaxis = list(title = ""), margin = list(l=50, r=50, b=100, t=100, pad=4), autosize = F)

##### Save the box-plot as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, "PCA_screeplot.html", sep = "_"))


##### Generate PCA plot (PLOTLY)
##### Prepare data frame
expData_pca.df <- data.frame(targets, expData_pca$x[,PC1], expData_pca$x[,PC2], expData_pca$x[,PC3])
colnames(expData_pca.df) <- c("Target", "PC1", "PC2", "PC3")
rownames(expData_pca.df) <- annData$File_name

p <- plot_ly(expData_pca.df, x = ~PC1, y = ~PC2, color = ~Target, text=rownames(expData_pca.df), colors = targets.colour[[1]], type='scatter', mode = "markers", marker = list(size=10, symbol="circle"), width = 800, height = 600) %>%
layout(title = "", xaxis = list(title = paste("PC ", PC1, " (",importance_pca[PC1],")",sep="")), yaxis = list(title = paste("PC ", PC2, " (",importance_pca[PC2],")",sep="")), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'v', y = 0.5), showlegend=TRUE)

##### Save the box-plot as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, "PCA.html", sep = "_"))



##### Generate PCA 3-D plot (PLOTLY)
p <- plot_ly(expData_pca.df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Target, text=rownames(expData_pca.df), colors = targets.colour[[1]], type='scatter3d', mode = "markers", marker = list(size=8, symbol="circle"), width = 800, height = 800) %>%
layout(scene = list(xaxis = list(title = paste("PC ", PC1, " (",importance_pca[PC1],")",sep="")), yaxis = list(title = paste("PC ", PC2, " (",importance_pca[PC2],")",sep="")), zaxis = list(title = paste("PC ", PC3, " (",importance_pca[PC3],")",sep="")) ), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'v', y = 0.5), showlegend=TRUE)

##### Save the box-plot as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste(coreName, "PCA_3D.html", sep = "_"))


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
