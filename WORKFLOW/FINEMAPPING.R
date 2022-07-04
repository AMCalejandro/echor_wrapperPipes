#!/bin/Rscript

## DESCRIPTION ##
# This is the workflow to automate the finemapping from echocolocatoR

# This script takes a txt file with the following informations
# RS_PATH
# SS_path


## Example usage
#Rscript FINEMAPPING.R [path_to_metadata_file.txt]


# Activate echoR conda env
#system("conda activate echoR")


# Load libraries
library(data.table) # Efficient multiprocessor data reading/writing/processing
library(echolocatoR) # To fo the finemapping
library(tidyverse) # For data wrangling / munging
library(tools) # To perform format changes in the SS file
library(rlang) # Data masking and more
library(vautils) # Use util to find nearest gene from an input SNP
library(here) # Efficient path management within the R project

# Read file with metadata to perform finemapping
args <- commandArgs(trailingOnly = TRUE)
metadata = fread(args[1])
source(here::here("R", "utils_v2.R"))


# Create / Localice RS dir
make_results_dir(metadata[2])


# Perform minor processing in the input df
SS_check(metadata[1])

# Get the input data, take the lead variant on each (nominal) significant 
# locus, then we find the nearest gene to it, and write out the top_SNPs file to be used in echolocatoR
lead_variants = gwas_lead_snps(metadata[1])

# Then we make the top_SNPs file
get_topSNPs(lead_variants, build = "hg19")





## At this point

# We have the SS GWAS updates if needed, RS updated if needed, top_SNPs file created

# We are ready to run finemap_loci

















# Get SS data
# Get SNP_GENE dat ato analyse
# Run echolocatoR
