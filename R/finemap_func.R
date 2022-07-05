


# DESCRIPTION
# This is the script whwere I keep the function I will use to run finemape
finemapping_wrapper = function(top_SNPs = top_SNPs, 
                               study_name = "PD_GWAS",
                               study_type = "motor_progression",
                               build = "hg19",
                               finemap_tools = c("ABF", "FINEMAP", "SUSIE", "POLYFUN_SUSIE"),
                               mean_SS = NULL,
                               .metadata = metadata) {
  
  
  # In a future, I should create two independent functions.
  # Finemapping_cc
  # Finemappinc_quant
  echolocatoR::finemap_loci(top_SNPs = top_SNPs, loci = top_SNPs$Locus,
                dataset_name = study_name, dataset_type = study_type,
                force_new_subset = F, force_new_LD = F,
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
               LD_reference = "UKB",
               superpopulation = "EUR",
               download_method = "axel",
                             
              # Additional arguments - My arguments
              case_control = F,
              nThread = 15,
              sample_size = mean_SS,
                   
              # PLOT ARGUMENTS
              ## general
              plot.types = c("simple"),
              server = T,
              ## Generate multiple plots of different window sizes;
              ### all SNPs, 4x zoomed-in, and a 50000bp window
              plot.zoom = c("1x","4x","10x"),

              plot.Nott_epigenome=F,
              plot.Nott_show_placseq = F,
              verbose = TRUE,
              # ENVIRONMENT ARGS
              conda_env= "echoR" )
  
  
  # Then I write the full path to be used for the plotting workflow
  fwrite(list(paste(.metadata[1],study_type,study_name, sep = "/")),
         paste(.metadata[1], "metadata_2.txt", sep = "/"),
         quote= F, sep ="", col.names=F, row.names=F)
}