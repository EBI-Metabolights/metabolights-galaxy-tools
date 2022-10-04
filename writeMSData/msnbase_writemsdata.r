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
source_local("lib.r")

pkgs <- c("MSnbase", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");

# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = F, quote = F, sep = "\t")

cat("\n\n")

# ----- LOAD RDATA -----
if (!exists("singlefile")) singlefile <- NULL
filePath <- retrieveSingleFileInWorkingDir(singlefile, args)
singlefile <- filePath$singlefile
files <- filePath$files

md5sumList <- list("origin" = getMd5sum(files))

cat("\n\n")

# this assumes that the RData has the same structure as at the point
# it is saved in write_ms_data
rdata <- load(singlefile)
raw_data <- rdata@raw_data
# ----- MAIN PROCESSING INFO / EXPORT -----


cat("\tMSnExp OBJECT INFO\n")
print(raw_data)
cat("\t\tphenoData\n")
print(raw_data@phenoData@data)
cat("\n\n")

writeMSData(rdata, paste0("write_ms", singlefile), verbose <- TRUE)

cat("\tDONE\n")
