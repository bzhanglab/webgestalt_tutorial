library(WebGestaltR)

files <- c("input_data/protein.rnk", "input_data/rna.rnk", "input_data/metab.rnk")
file_names <- c("protein", "rna", "metab")
id_types <- c("genesymbol", "genesymbol", "hmdb")
output_dir <- "output_webgestalt"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
# Perform WebGestalt analysis for each file
WebGestaltRMultiOmics(analyteListFiles = files,
                      listNames = file_names,
                      analyteTypes = id_types,
                      outputDirectory = output_dir,
                      organism = "hsapiens",
                      enrichMethod = "GSEA",
                      enrichDatabase = "pathway_WikiPathways",
                      minNum = 5,
                      maxNum = 5000,
                      fdrMethod = "BH",
                      fdrThr = 0.05,
                      topThr = 50,)
