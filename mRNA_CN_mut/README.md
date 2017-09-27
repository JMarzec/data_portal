## Visualising pre-computed mRNA expression, DNA copy-number and mutation data

This repository contains script for visualising pre-computed per-gene mRNA expression, DNA copy-number and mutation data. All script allow to process gene matrices with duplicated gene IDs.

Script | Description
------------ | ------------
*[Expr_CN_mut_Profile.R](https://github.com/JMarzec/data_portal/tree/master/mRNA_CN_mut/Expr_CN_mut_Profile.R)* | Generates scatter-plot and box-plot to visualise expression measurments against relative linear copy-number values and putative copy-number alterations, respectively, across samples with indicated mutation status for the gene of interest
<br />


###   Command line usage with example data

*[Expr_CN_mut_Profile.R](https://github.com/JMarzec/data_portal/tree/master/mRNA_CN_mut/Expr_CN_mut_Profile.R)*
```
      R --file=./Expr_CN_mut_Profile.R --args "CCLE_PC_processed_mRNA.txt" "CCLE_PC_processed_CN.txt" "CCLE_PC_processed_mut.txt" "CCLE_PC_target.txt" "KRAS" "Example_results/PC_Expr_CN_mut_Profile"
```
<br />


###   Output files for the data portal

Script | Output files
------------ | ------------
*[Expr_CN_mut_Profile.R](https://github.com/JMarzec/data_portal/tree/master/mRNA_CN_mut/Expr_CN_mut_Profile.R)* | <li> [*input_file_name*]\_[*gene_name*]\_hist.html</li> <li> [*input_file_name*]\_[*gene_name*]\_mRNA_vs_CN_mut_plot.html</li> <li> [*input_file_name*]\_[*gene_name*]\_mRNA_vs_CN_mut_boxplot.html</li>
<br />