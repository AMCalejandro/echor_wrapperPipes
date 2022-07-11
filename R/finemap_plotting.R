


# DESCRIPTION
# This is the script whwere I keep the function I will use to run finemape
finemapping_wrapper = function(top_SNPs = top_SNPs, 
                               study_name = "PD_GWAS",
                               study_type = "motor_progression",
                               build = "hg19",
                               finemap_tools = c("ABF", "FINEMAP", "SUSIE", "POLYFUN_SUSIE"),
                               mean_SS = NULL,
                               ld_ref = "UKB",
                               .metadata = metadata) {
  
  
  # In a future, I should create two independent functions.
  # Finemapping_cc
  # Finemappinc_quant
  echolocatoR::finemap_loci(top_SNPs = top_SNPs, loci = top_SNPs$Locus,
                dataset_name = study_name, dataset_type = study_type,
                force_new_subset = F, force_new_LD = T,
                force_new_finemap = T, remove_tmps = F,
               
               # SUMMARY STATS ARGUMENTS
               fullSS_genome_build = build,
               fullSS_path = .metadata[2], results_dir = .metadata[1], query_by = "tabix",
               chrom_col = "CHR", position_col = "POS", snp_col = "SNP",
               pval_col = "Pval", effect_col = "Effect", stderr_col = "StdErr",
               freq_col = "Freq1", MAF_col = "calculate",
               A1_col = "Allele1", A2_col = "Allele2",
               #N_cases_col = "TotalSampleSize",
               #N_controls = 0,

               # FILTERING ARGUMENTS
               ## It's often desirable to use a larger window size
               ## (e.g. 2Mb which is bp_distance=500000*2),
               ## but we use a small window here to speed up the process.
               bp_distance = 500000*2,
               min_MAF = 0.001,
               trim_gene_limits = F,
               # FINE-MAPPING ARGUMENTS
               ## General
               finemap_methods = finemap_tools,
               n_causal = 5,
               PP_threshold = .95,
               consensus_threshold = 2,
               # LD ARGUMENTS
               LD_genome_build = build,
               LD_reference = 
               superpopulation = "EUR",
               download_method = "axel",
                             
              # Additional arguments - My arguments
              case_control = F,
              nThread = 15,
              sample_size = mean_SS,
                   
              # PLOT ARGUMENTS
              ## general
              plot.types = c("fancy"),
              server = T,
              ## Generate multiple plots of different window sizes;
              ### all SNPs, 4x zoomed-in, and a 50000bp window
              plot.zoom = c("All", "1x","4x","10x"),
              verbose = TRUE,
              # ENVIRONMENT ARGS
              conda_env= "echoR" )
  
  # Then I write the full path to be used for the plotting workflow
  writeLines(paste(.metadata[1],study_type,study_name, sep = "/"),
         paste(.metadata[1], "metadata_2.txt", sep = "/") )
  


}





get_regional_plots = function(locus, 
                              ld_ref = "UKB",
                              plot.zoom = c("All", "4X", "10x", "20x"),
                              show_plot = T,
                              .metadata = metadata2) {
  
  finemap_file = list.files(paste(.metadata, locus,"Multi-finemap", sep = "/"))
  multifinemap = fread(paste(.metadata,
                             locus,
                             "Multi-finemap", 
                             finemap_file, sep = "/")) %>%
    echolocatoR::find_consensus_SNPs(.)
  
  ldmatrix_index = list.files(paste(.metadata, locus,"LD", sep = "/")) %>% 
    grep(pattern="*RDS", x = .)
  
  ldmatrix_file = list.files(paste(.metadata, locus,"LD", sep = "/"))[ldmatrix_index]
  ldmatrix = readRDS(paste(.metadata, 
                           locus,
                           "LD", 
                           ldmatrix_file, sep = "/"))
  
  locus_dir = paste(.metadata, locus, sep = "/")
  
  # To get the finemapping plot
  trk_plot <- PLOT.locus(finemap_dat=multifinemap, 
                         LD_matrix=as.matrix(ldmatrix), 
                         LD_reference=ld_ref,
                         locus_dir=locus_dir,  
                         save_plot=T,
                         show_plot=show_plot,
                         plot.zoom=plot.zoom, 
                         conda_env = "echoR")
  
  names_old = list.files(locus_dir, pattern="^multiview")
  names_new = paste0("finemapplot_",list.files(locus_dir, pattern="^multiview"))
  file.rename(paste(locus_dir, names_old, sep = "/"),
              paste(locus_dir, names_new, sep = "/"))
  
  # To get the locus plot with FANTOM_ENHANCER data
  trk_plot <- PLOT.locus(finemap_dat=multifinemap, 
                         LD_matrix=ldmatrix, 
                         LD_reference=ld_ref,
                         locus_dir=locus_dir,  
                         save_plot=T,
                         show_plot=show_plot,
                         plot.zoom=plot.zoom,
                         nThread = parallel::detectCores()-1,
                         conda_env = "echoR",
                         
                         XGR_libnames=c("FANTOM5_Enhancer"))
  
  names_old = list.files(locus_dir, pattern="^multiview")
  names_new = paste0("XGRplot_",list.files(locus_dir, pattern="^multiview"))
  file.rename(paste(locus_dir, names_old, sep = "/"),
              paste(locus_dir, names_new, sep = "/"))
  
  
  # To get the locus plot with Roadmap brain tissues data
  trk_plot <- PLOT.locus(finemap_dat=multifinemap, 
                         LD_matrix=ldmatrix, 
                         LD_reference=ld_ref,
                         locus_dir=locus_dir,  
                         save_plot=T,
                         show_plot=show_plot,
                         plot.zoom=plot.zoom,
                         nThread = parallel::detectCores()-1,
                         conda_env = "echoR",
                         
                         Roadmap=T, 
                         Roadmap_query=c("brain"))
  
  names_old = list.files(locus_dir, pattern="^multiview")
  names_new = paste0("RoadmapBrainplot_",list.files(locus_dir, pattern="^multiview"))
  file.rename(paste(locus_dir, names_old, sep = "/"),
              paste(locus_dir, names_new, sep = "/"))
  
  
  # To get the locus plot with Nott et al data
  trk_plot <- PLOT.locus(finemap_dat=multifinemap, 
                         LD_matrix=ldmatrix, 
                         LD_reference=ld_ref,
                         locus_dir=locus_dir,  
                         save_plot=T,
                         show_plot=show_plot,
                         plot.zoom=plot.zoom,
                         nThread = parallel::detectCores()-1,
                         conda_env = "echoR",
                         
                         Nott_epigenome=T,  
                         Nott_binwidth = 200,
                         Nott_regulatory_rects = T, 
                         Nott_show_placseq = T)
  
  
  names_old = list.files(locus_dir, pattern="^multiview")
  names_new = paste0("NOTTplot_",list.files(locus_dir, pattern="^multiview"))
  file.rename(paste(locus_dir, names_old, sep = "/"),
              paste(locus_dir, names_new, sep = "/"))
}
