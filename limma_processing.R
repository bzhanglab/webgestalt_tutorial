#!/usr/bin/env Rscript

# Load required libraries
library(limma)
library(readr)
library(dplyr)
library(EnhancedVolcano)

# Function to perform limma analysis for a single omics dataset
perform_limma_analysis <- function(data_file, phenotype_file_in, phenotype_file_out) {
    # Read the data file
    data <- read.delim(data_file, check.names = FALSE)

    # Identify the first column name
    first_col_name <- names(data)[1]

    # Remove duplicate rows, keeping the last occurrence
    data <- data[!duplicated(data[[first_col_name]], fromLast = TRUE), ]

    # Read phenotype files
    samples_in <- read_lines(phenotype_file_in)
    samples_out <- read_lines(phenotype_file_out)

    # Find matching samples
    matching_samples <- intersect(samples_in, names(data))
    matching_samples_out <- intersect(samples_out, names(data))

    # Check if we have samples
    if (length(matching_samples) == 0 || length(matching_samples_out) == 0) {
        stop(paste(
            "No matching samples found in", data_file,
            "\nSamples in:", paste(samples_in, collapse = ", "),
            "\nSamples out:", paste(samples_out, collapse = ", "),
            "\nActual data columns:", paste(names(data), collapse = ", ")
        ))
    }

    # Create design matrix
    all_samples <- c(matching_samples, matching_samples_out)
    phenotype <- rep(c(1, -1), c(length(matching_samples), length(matching_samples_out)))
    design <- model.matrix(~phenotype)
    rownames(design) <- all_samples

    # Prepare expression matrix
    expr_matrix <- data[, all_samples]
    rownames(expr_matrix) <- data[[first_col_name]]

    # Perform limma analysis
    fit <- lmFit(expr_matrix, design)
    fit <- eBayes(fit)

    # Extract results
    results <- topTable(fit, coef = 2, number = Inf)

    # Create signed -log(p) column
    results$signed_log_p <- -log10(results$P.Value) * sign(results$logFC)


    # Retain first column name and signed log p
    results$analyte <- rownames(results)
    results$signed_log_p <- results$signed_log_p

    # get base file name and use as output file name plus title

    base_file_name <- tools::file_path_sans_ext(data_file)
    # remove the folder
    hard_base <- basename(base_file_name)
    volcano_plot <- EnhancedVolcano(
        results,
        lab = results$analyte,
        subtitle = "",
        x = "logFC",
        y = "P.Value",
        title = hard_base,
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 3.0,
        labSize = 3.0
    )
    # Save the plot
    ggsave(paste0(base_file_name, "_volcano.pdf"), plot = volcano_plot, width = 8, height = 6)
    return(results[, c("analyte", "signed_log_p")])
}

# Function to save rnk file
save_rnk_file <- function(results, output_file) {
    write.table(results, output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
}

# Omics datasets to process
omics_datasets <- c("input_data/rna.tsv", "input_data/protein.tsv", "input_data/metabolites.tsv")

# Process each dataset
for (dataset in omics_datasets) {
    # Derive output filename
    base_name <- tools::file_path_sans_ext(dataset)
    rnk_file <- paste0(base_name, ".rnk")

    # Perform analysis
    results <- perform_limma_analysis(
        data_file = dataset,
        phenotype_file_in = "in.txt",
        phenotype_file_out = "out.txt"
    )

    # Save rnk file
    save_rnk_file(results, rnk_file)

    # Print confirmation
    cat("Processed", dataset, "- Generated", rnk_file, "\n")
}

cat("Analysis complete. RNK files generated for each omics dataset.\n")
