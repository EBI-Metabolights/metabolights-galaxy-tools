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
args <- parseCommandArgs(FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names = F, quote = F, sep = "\t")

cat("\n\n")

# ----- LOAD RDATA -----
if (!exists("singlefile")) singlefile <- NULL
filePath <- retrieveSingleFileInWorkingDir(singlefile, args, "")
singlefile <- filePath$singlefile
files <- filePath$files

md5sumList <- list("origin" = getMd5sum(files))

cat("\n\n")

# this assumes that the RData has the same structure as at the point
# it is saved in write_ms_data
cat(paste0(unlist(singlefile)))
cat("\n\n")
cat(files[1])
cat("\n\n")
cat(paste0(unlist(files)))
cat("\n\n")
rdata <- get(load(files[1]))
cat("\n\n")
cat(str(rdata))
raw_data <- rdata$raw_data
# ----- MAIN PROCESSING INFO / EXPORT -----

cat("\n\n")
cat("\tMSnExp OBJECT INFO\n")
print(raw_data)
cat("\t\tphenoData\n")
print(raw_data["phenoData"]["data"])
cat("\n\n")
print(str(rdata@processingData@files))
cat("\n\n")
cat("length of samples")
cat((length(rdata@experimentData@samples)))
cat("\n\n")
ls()

# remove the record of the filename.mzML to prevent msnbase from trying to locate it.
rdata@processingData@files <- rdata@processingData@removedPeaks

newobj <- new("MSnExp")
newobj@assayData = rdata@assayData
newobj@phenoData = rdata@phenoData
newobj@featureData = rdata@featureData
newobj@experimentData = rdata@experimentData
newobj@protocolData = rdata@protocolData
newobj@processingData = rdata@processingData


#zzfil <- tempfile(pattern = str(rdata@processingData@files), fileext = ".mzml")
#zz <- file(zzfil, "w")  # open an output file connection

#file.create(rdata@processingData@files)
writeMSData(newobj, file = paste0(tempfile(), ".mzML"), copy = FALSE)

cat("\tDONE\n")
