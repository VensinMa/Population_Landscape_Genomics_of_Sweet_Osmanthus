#!/usr/bin/env Rscript

# Ensure necessary libraries are installed and loaded
if (!requireNamespace("argparser", quietly = TRUE)) install.packages("argparser")
if (!requireNamespace("pophelper", quietly = TRUE)) install.packages("pophelper")

library(argparser)
library(pophelper)

# Create argument parser with a description
p <- arg_parser("Draw structure figure for admixture result",
                description = "This script plots structure figures from admixture result files. It requires a directory containing Q matrix files, a sample order file (e.g., .nosex), and an output prefix for saving the figures.")

# Add command line arguments
p <- add_argument(p, "dir", help="Input: Directory containing Q matrix files", type="character")
p <- add_argument(p, "sample", help="Input: Sample order file (e.g., .nosex)", type="character")
p <- add_argument(p, "output", help="Output prefix for saved figures", type="character")

# Parse the command line arguments
argv <- parse_args(p)

# Validate input
if (!dir.exists(argv$dir)) {
  cat(sprintf("Error: The specified directory '%s' does not exist.\n", argv$dir))
  print_help(p)
  quit(status = 1)
}

if (!file.exists(argv$sample)) {
  cat(sprintf("Error: The specified sample file '%s' does not exist.\n", argv$sample))
  print_help(p)
  quit(status = 1)
}

# Proceed with the script logic
dir <- argv$dir
sample_file <- argv$sample
output_prefix <- argv$output

Qfiles <- list.files(dir, pattern="Q", full.names=TRUE)
qlist <- readQ(Qfiles, indlabfromfile=FALSE)

sample_labels <- read.table(sample_file, header=FALSE)
for (i in 1:length(qlist)) {
  rownames(qlist[[i]]) <- sample_labels$V1
}

# Plot structure figure using pophelper
plotQ(qlist,
      sortind="all",
      imgtype="pdf",
      ordergrp=FALSE,
      imgoutput="join",
      width=20,
      height=4,
      exportpath=getwd(),
      outputfilename=output_prefix,
      useindlab=TRUE,
      sharedindlab=FALSE,
      showindlab=TRUE)

cat("succcedusr/bin/env Rscript

# Ensure necessary libraries are installed and loaded
if (!requireNamespace("argparser", quietly = TRUE)) install.packages("argparser")
if (!requireNamespace("pophelper", quietly = TRUE)) install.packages("pophelper")

library(argparser)
library(pophelper)

# Create argument parser with a description
p <- arg_parser("Draw structure figure for admixture result",
                description = "This script plots structure figures from admixture result files. It requires a directory containing Q matrix files, a sample order file (e.g., .nosex), and an output prefix for saving the figures.")

# Add command line arguments
p <- add_argument(p, "dir", help="Input: Directory containing Q matrix files", type="character")
p <- add_argument(p, "sample", help="Input: Sample order file (e.g., .nosex)", type="character")
p <- add_argument(p, "output", help="Output prefix for saved figures", type="character")

# Parse the command line arguments
argv <- parse_args(p)

# Validate input
if (!dir.exists(argv$dir)) {
  cat(sprintf("Error: The specified directory '%s' does not exist.\n", argv$dir))
  print_help(p)
  quit(status = 1)
}

if (!file.exists(argv$sample)) {
  cat(sprintf("Error: The specified sample file '%s' does not exist.\n", argv$sample))
  print_help(p)
  quit(status = 1)
}

# Proceed with the script logic
dir <- argv$dir
sample_file <- argv$sample
output_prefix <- argv$output

Qfiles <- list.files(dir, pattern="Q", full.names=TRUE)
qlist <- readQ(Qfiles, indlabfromfile=FALSE)

sample_labels <- read.table(sample_file, header=FALSE)
for (i in 1:length(qlist)) {
  rownames(qlist[[i]]) <- sample_labels$V1
}

# Plot structure figure using pophelper
plotQ(qlist,
      sortind="all",
      imgtype="pdf",
      ordergrp=FALSE,
      imgoutput="join",
      width=20,
      height=4,
      exportpath=getwd(),
      outputfilename=output_prefix,
      useindlab=TRUE,
      sharedindlab=FALSE,
      showindlab=TRUE)

cat("Structure plot has been saved with the prefix:", output_prefix, "\n")

