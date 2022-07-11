

library(ggplot2)
library(data.table) # Efficient multiprocessor data reading/writing/processing
library(echolocatoR) # To fo the finemapping
library(tidyverse) # For data wrangling / munging

#metadata2 = readLines("/mnt/rreal/RDS/acarrasco/TRIALS/ECHOLOCATOR_WORKFLOW/metadata_2.txt")

args <- commandArgs(trailingOnly = TRUE)
metadata2 = readLines(args[1])
# Making sure I remove any problematic white spaces
metadata2 = metadata2[which(metadata2!="")]
metadata2 = trimws(gsub("\\s+", "", metadata2))


# Find the number of
loci_to_plot = list.files(metadata2[1])
ld_ref <- "UKB"
plot.zoom = c("All", "4x", "10x", "20x")
show_plot=T

source("/mnt/rreal/RDS/acarrasco/TOOLS/echor_wrapperPipes/R/finemap_plotting.R")
map(loci_to_plot, ~get_regional_plots(.x, 
                                      plot.zoom = plot.zoom,
                                      ld_ref = ld_ref))

