## Visualising pre-computed gene expression data

This repository contains script for visualising pre-computed normalised gene expression data. All script allow to process gene matrices with duplicated gene IDs.

Script | Description
------------ | ------------
*[ExprProfile.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/ExprProfile.R)* | Generates box-plots and bar-plots to visualise expression measurments across samples and groups for user-defined gene
*[GenesHeatmap.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GenesHeatmap.R)* | Performs clustering and generates heatmap for all or user-defined genes
*[PCA.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/PCA.R)* | Performs principal component analysis
*[GeneCoExpr.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GeneCoExpr.R)* | Calculates co-expression of user-defined genes across all samples or samples in user-defined group and presents correlation coefficients between samples as well as associated p-values in a form of correlation matrix heatmap
*[GenesCor.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GenesCor.R)* | Performs pairwise comparisons of expression profiles between multiple user-defined genes across all samples or samples in user-defined group and presents correlation coefficients between genes in a form of correlation matrix heatmap
*[GeneCoExprHeatmap.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GeneCoExprHeatmap.R)* | Calculates co-expression between user-defined gene and all other genes and generates heatmaps of top positively and negatively correlated genes
<br />


###   Command line usage with example data

*[ExprProfile.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/ExprProfile.R)*
```
      R --file=./ExprProfile.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "Derived_From" "KRAS" "Example_results/PC_ExprProfile"
```
<br />

*[GenesHeatmap.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GenesHeatmap.R)*
```
      R --file=./GenesHeatmap.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "Example_results/PC_GenesHeatmap" "Gene_list.txt" "3" "3"
```
<br />

*[PCA.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/PCA.R)*
```
      R --file=./PCA.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "Derived_From" "1" "Example_results/PC_PCA"
```
<br />

*[GeneCoExpr.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GeneCoExpr.R)*
```
      R --file=./GeneCoExpr.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "Genes_of_interest.txt" "Example_results/PC_GeneCoExpr" "Primary Tumour"
```
<br />

*[GenesCor.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GenesCor.R)*
```
      R --file=./GenesCor.R --args "TCGA_PAAD_normalized.txt" "TCGA_PAAD_target.txt" "Genes_of_interest.txt" "Example_results/PC_GenesCor" "PDAC"
```
<br />

*[GeneCoExprHeatmap.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GeneCoExprHeatmap.R)*
```
      R --file=./GeneCoExprHeatmap.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_target.txt" "KRAS" "20" "Example_results/PC_GeneCoExprHeatmap" "Gene_list.txt"
```
<br />


###   Output files for the data portal

Script | Output files
------------ | ------------
*[ExprProfile.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/ExprProfile.R)* | <li> [*input_file_name*]\_[*gene_name*]\_expression_profile_barplot.html</li> <li>[*input_file_name*]\_[*gene_name*]\_expression_profile_boxplot.html</li>
*[GenesHeatmap.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GenesHeatmap.R)* | <li> [*input_file_name*]\_heatmap.html</li>
*[PCA.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/PCA.R)* | <li> [*input_file_name*]\_PCA_screeplot.html</li> <li> [*input_file_name*]\_PCA_3D.html</li>
*[GeneCoExpr.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GeneCoExpr.R)* | <li> [*input_file_name*]\_heatmap_corr_R.html</li> <li> [*input_file_name*]\_heatmap_corr_P.html</li>
*[GenesCor.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GenesCor.R)* | <li> [*input_file_name*]\_corr_heatmap.html</li> <li> [*input_file_name*]\_corr.txt</li>
*[GeneCoExprHeatmap.R](https://github.com/JMarzec/data_portal/tree/master/mRNA/GeneCoExprHeatmap.R)* | <li> [*input_file_name*]\_[*gene_name*]\_corr_hist.html</li> <li> [*input_file_name*]\_[*gene_name*]\_heatmap_pos.html</li> <li> [*input_file_name*]\_[*gene_name*]\_heatmap_neg.html</li>
<br />
