

#TODO
# Add function descriptions




SS_check = function(file) {
  myfile = fread(file, header=TRUE)
  firstchar = substr(colnames(myfile)[1],1,1)
  extension = tools::file_ext(file)
  
  if (firstchar != "#") {
    colnames(myfile)[1] <- paste0("#", colnames(myfile)[1])
    vcf_reformat = T
  }
  if (extension != "tsv") {
    change_format = T
  }
  if (vcf_reformat | change_format) {
    print("\n\n REFORMATTING SS \n\n ")
    new = gsub(extension, "tsv", file)
    fwrite(myfile, new, colnames = T, row.names = F, sep ="\t", quote = F)
    fullRS_path = new
  } else {
    print("\n\n ALL GOOD IN SS. NO CHANGES MADE \n\n")
  }
}


gwas_lead_snps = function(data,
                       pval_thres = 5e-8,
                       pval_col = Pval,
                       chr_col = CHR,
                       bp_col = BP) {
  
  chrvar <- rlang::sym(chr_col)
  bpvar <-  rlang::sym(bp_col)
  pvvar = rlang::sym(pval_col)
  
  data = data %>%
    dplyr::filter(!!pvvar <= pval_thres) %>%
    group_by(!!chrvar, !!bpvar) %>%
    dplyr::arrange(!!pval_col) %>%
    dplyr::filter(row_number() == 1)
  
  return(data)
}



# get_topSNPs(prueba, pval_col= "P-value", 
#             pval_thres = pval_threshold,
#             chr_col = chr_col,
#             bp_col = bp_col)


make_results_dir <- function(results_path){
  if(!dir.exists(results_dir_path)){
    dir.create(results_dir_path)
  } else {
    print(stringr::str_c(results_dir_path, " \n\n directory already exists.. \n\n"))
  }
  return(results_dir_path)
}





# Nominate the closest gene


# Check this two sources to come up with the closest gene
#https://jef.works/blog/2016/12/06/mapping-snps-and-peaks-to-genes-in-R/
#install_github('git@github.com:oyhel/vautils.git')

mydf = data.frame(matrix(ncol = 7, nrow = 0))
colnames(mydf) <- c("Locus","Gene","CHR","POS","SNP","P","Effect")


data = fread("EARLY_PD/POST_GWAS/ECHOLOCATOR/for_echolocatoR_axialOutcome_3.tsv")
data = data %>% arrange(Pval)
data = head(data)


make_topSNPs = function(data, build = "hg19") {
  

  snps_mapped = vautils::find_nearest_gene(data,
                                           build=build,
                                           collapse=F,
                                           snp ="#MarkerName",
                                           flanking=1000,
                                           chr = "CHR",
                                           bp = "POS" )
  
  if(any(snps_mapped$distance == "intergenic")) {
    snps_mapped_intergenic = snps_mapped %>%
      filter(distance == "intergenic")
  }
  
  snps_mapped_nearest = snps_mapped %>%
    filter(!rsid %in% c(snps_mapped_intergenic %>% pull(rsid))) %>%
    mutate(distance = abs(as.numeric(distance))) %>%
    arrange(distance) %>%
    group_by(rsid) %>%
    filter(row_number() == 1)
  
  tmp = data.frame(Locus = c(snps_mapped_intergenic$GENE, snps_mapped_nearest$GENE),
                   Gene =c(snps_mapped_intergenic$GENE, snps_mapped_nearest$GENE),
                   CHR = c(snps_mapped_intergenic$chromosome, snps_mapped_nearest$chromosome),
                   POS = c(snps_mapped_intergenic$position, snps_mapped_nearest$position),
                   SNP = c(snps_mapped_intergenic$rsid, snps_mapped_nearest$rsid))
  
  tmp_2 = data %>% filter(`#MarkerName` %in% tmp$SNP) %>%
    select(SNP = `#MarkerName`, Effect, P = Pval)
  
  mydf = tmp %>% inner_join(tmp_2)
  
  message("\n\n Writing top_SNPs... \n\n")
  fwrite(mydf, paste0(metadata[1],"/topSNPs.txt"),
         quote = F, col.names=T, row.names=F, sep = "\t")
  #print(mydf)
  
}











  






