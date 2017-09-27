## Visualising pre-computed DNA copy-number data

This repository contains script for visualising pre-computed per-gene DNA copy-number data. All script allow to process gene matrices with duplicated gene IDs.

Script | Description
------------ | ------------
*[CN_linear_heatmap.R](https://github.com/JMarzec/data_portal/tree/master/CN/CN_linear_heatmap.R)* | Performs clustering and generates heatmap for user-defined genes using relative linear copy-number values
*[CN_binary_heatmap.R](https://github.com/JMarzec/data_portal/tree/master/CN/CN_binary_heatmap.R)* | Calls gained and lost genes by estimating the putative copy-number alterations from the linear copy-number values, performs clustering and generates heatmap for user-defined genes
<br />


###   Command line usage with example data

*[CN_linear_heatmap.R](https://github.com/JMarzec/data_portal/tree/master/CN/CN_linear_heatmap.R)*
```
      R --file=./CN_linear_heatmap.R --args "CCLE_PC_CN_processed.txt" "CCLE_PC_target.txt" "Example_results/PC_CN_linear_heatmap" "Gene_list.txt" "3" "10"
```
<br />

*[CN_binary_heatmap.R](https://github.com/JMarzec/data_portal/tree/master/CN/CN_binary_heatmap.R)*
```
      R --file=./CN_binary_heatmap.R --args "CCLE_PC_CN_processed.txt" "CCLE_PC_target.txt" "Example_results/PC_CN_binary_heatmap" "Gene_list.txt" "3" "10"
```
<br />


###   Output files for the data portal

Script | Output files
------------ | ------------
*[CN_linear_heatmap.R](https://github.com/JMarzec/data_portal/tree/master/CN/CN_linear_heatmap.R)* | <li> [*input_file_name*]\_linear_heatmap.html</li>
*[CN_binary_heatmap.R](https://github.com/JMarzec/data_portal/tree/master/CN/CN_binary_heatmap.R)* | <li> [*input_file_name*]\_hist.html</li> <li> [*input_file_name*]\_binary_heatmap.html</li>
<br />
