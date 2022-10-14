#!/usr/bin/env Rscript

# ----- LOG FILE -----
log_file <- file("log.txt", open = "wt")
sink(log_file)
sink(log_file, type = "output")

# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep = "/"))
}
#Import the nmr package
source_local("utils.r")
pkgs <- c("devtools", "MassSpecWavelet")
loadAndDisplayPackages(pkgs)
cat("\n\n");
devtools::install_github("EBI-Metabolights/icl_nmr_R")

# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = F, quote = F, sep = "\t")
cat("\n\n")
print(args)
cat("\n\n")
cat("hit final line")