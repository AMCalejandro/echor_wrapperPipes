#TODO
# Add function descriptions
SS_check = function(data, metadata) {
  mydata = data
  firstchar = substr(colnames(mydata)[1],1,1)
  
  
  # Sort input file based on CHR POS - Tabix crashes otherwise
  # Also, use scipen option to pevent scientific notation
  options(scipen=999)
  mydata = mydata %>% arrange(CHR, POS)
  

  if (firstchar != "#") {
    colnames(mydata)[1] <- paste0("#", colnames(mydata)[1])
  }
  
  
  message("Saving data in vcf format and tsv extension in the working dir")
  new_path = paste(metadata[2], metadata[3], sep="/")
  fwrite(mydata, new_path, col.names = T, row.names = F, sep ="\t", quote = F)
  
  # Apending new path to metadata file
  metadata = c(metadata, new_path)
  return(metadata)
}

gwas_lead_snps = function(data,
                       pval_thres = 5e-8,
                       pval_col = "Pval",
                       chr_col = "CHR",
                       pos_col = "POS") {
  
  chrvar <- rlang::sym(chr_col)
  posvar <-  rlang::sym(pos_col)
  pvvar = rlang::sym(pval_col)

  data = data %>%
    dplyr::filter(!!pvvar <= pval_thres) %>%
    group_by(!!chrvar) %>%
    dplyr::arrange(!!pvvar) %>%
    dplyr::filter(row_number() == 1) %>%
    ungroup()

  return(data)
}



# get_topSNPs(prueba, pval_col= "P-value", 
#             pval_thres = pval_threshold,
#             chr_col = chr_col,
#             bp_col = bp_col)


make_results_dir <- function(results_path){
  if(!dir.exists(results_path)){
    dir.create(results_path)
  } else {
    message(stringr::str_c(results_path, " directory already exists.. \n\n"))
  }
  return(results_path)
}





# Nominate the closest gene


# Check this two sources to come up with the closest gene
#https://jef.works/blog/2016/12/06/mapping-snps-and-peaks-to-genes-in-R/
#install_github('git@github.com:oyhel/vautils.git')



make_topSNPs = function(data,
                        build = "hg19", 
                        write.out = T,
                        custom_gene = NULL,
                        .metadata_file = metadata) {
  # I need to do this in advance to enter vautils function
  data = data %>%
    dplyr::rename(rsid = `#SNP`, 
                  chromosome = CHR,
                  position = POS)

  snps_mapped = vautils::find_nearest_gene(as.data.frame(data),
                                           build=build,
                                           collapse=F,
                                           snp = "rsid",
                                           flanking=1000)
  

  if(any(snps_mapped$distance == "intergenic")) {
    snps_mapped = snps_mapped %>%
      mutate(distance = recode(distance, "intergenic" = "0"))
  }
  
  snps_mapped_filt = snps_mapped %>%
    mutate(distance = abs(as.numeric(distance))) %>%
    arrange(distance) %>%
    group_by(rsid) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  if (!is.null(custom_gene)) {
     custom_gene_filt = snps_mapped %>%
       dplyr::filter(GENE %in% custom_gene)
     snps_mapped_filt = rbind(snps_mapped_filt, custom_gene_filt)
  }
  
  snps_mapped_filt = data.frame(Locus = snps_mapped_filt$GENE,
                     Gene =snps_mapped_filt$GENE,
                     CHR = snps_mapped_filt$chromosome,
                     POS = snps_mapped_filt$position,
                     SNP = snps_mapped_filt$rsid)
  
  
  tmp_2 = data %>% 
    dplyr::filter(rsid %in% snps_mapped_filt$SNP) %>%
    dplyr::select(SNP = rsid, Effect, P = Pval)
  mydf = snps_mapped_filt %>% inner_join(tmp_2)
  
  
  if (write.out) {
    fwrite(mydf, paste(.metadata_file[1],"top_SNPs.txt", sep = "/"),
           col.names = T,
           row.names = F,
           sep ="\t",
           quote = F)
  }
  return(mydf)
}












  






