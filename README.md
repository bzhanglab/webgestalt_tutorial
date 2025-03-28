# webgestalt_tutorial

Tutorial of how to start from abundance data to using WebGestalt for Enrichment Analysis. You can download all the code from this repository [here](https://github.com/bzhanglab/webgestalt_tutorial/archive/refs/heads/main.zip).


This project uses RNA, Protein, and Metabolite data from [Quirós  et. al](https://pubmed.ncbi.nlm.nih.gov/28566324/).

> Quirós PM, Prado MA, Zamboni N, D'Amico D, Williams RW, Finley D, Gygi SP, Auwerx J. Multi-omics analysis identifies ATF4 as a key regulator of the mitochondrial stress response in mammals. J Cell Biol. 2017 Jul 3;216(7):2027-2045. doi: 10.1083/jcb.201702058. Epub 2017 May 31. PMID: 28566324; PMCID: PMC5496626.


## Structure

```
input_data/
├─ abundance data and output of rank files and volcano plots
output_webgestalt/
├─ output from run_webgestalt.R
in.txt
out.txt
limma_processing.R
run_webgestalt.R
```

## 1. Installation

This requires R, and the `limma`, `readr`, `dplyr`, and `enchancedvolcano` packages.

You can install this by running the following R script

```R
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

  BiocManager::install(c('EnhancedVolcano', 'limma', 'readr', 'dplyr'))
```

If you intend to use the WebGestaltR package for enrichment analysis, follow the [installation structures](https://bzhanglab.github.io/WebGestaltR/articles/Installation.html).

## 2. Differential Expression Analysis

This is done through `limma` and can be ran by executing the `limma_processing.R` script (`Rscript limma_processing.R`).

This script uses the samples identified in the `in.txt` and `out.txt` files to determine the positive (in) and negative (out) associated samples. This is currently configured to compare the FCCP condition vs. the Control condition from the example data.

This will create the `.rnk` files for the three omic types (protein, rna, and metabolites). You can copy and paste these files into [webgestalt.org](https://www.webgestalt.org/), or install the `WebGestaltR` package and run the analysis locally.

## 3. WebGestaltR (optional)

If you would rather use the R package to run the analysis locally, you can run the `run_webgestalt.R` script. This will run a multi-omic GSEA analysis using all three omic-types. You can modify the different parameters in the R file to optimize for your analysis.