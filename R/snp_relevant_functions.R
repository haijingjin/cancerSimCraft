# snp_relevant_functions_v3.R

vcf_to_snp_list <- function(phased_vcf_table, sample_name){

  tmp <- phased_vcf_table[,c('CHROM', 'POS', 'REF', 'ALT', sample_name)]
  tmp <- tmp %>% tidyr::separate(sample_name, into = c("GT", "DS"), sep = ":") %>% dplyr::select(CHROM, POS, REF, ALT, GT) %>% dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1)

  snp_list <- list()
  snp_list[['all']] <- tmp
  snp_list[['1|1']] <- tmp %>% dplyr::filter(GT == '1|1')
  snp_list[['1|0']] <- tmp %>% dplyr::filter(GT == '1|0')
  snp_list[['0|1']] <- tmp %>% dplyr::filter(GT == '0|1')

  return(snp_list)
}


check_ref_snp_match <- function(seg_names, snp_list, ref_genome){

  for(i in 1:length(seg_names)){
    ref_seq <- ref_genome[[ seg_names[i] ]]
    snp_table <- snp_list[[ seg_names[i] ]]$all

    # construct IRanges instance
    ir <- IRanges(start = snp_table$POS, end = snp_table$POS)
    # extract seq from ref genome
    seq_ref <- extractAt(x = ref_seq, at = ir)
    # extract seq from snp table
    snp_ref <- snp_table$REF
    # check if two refs are identifcal
    if(!all(seq_ref == snp_ref)){
      warning("seq_ref and snp_ref don't match!")
    }else{
      print(paste0("seq_ref and snp_ref matches in ", seg_names[i], "!"))
    }
  }
}


insert_snps_to_genome <- function(seg_names, snp_list, ref_genome){

  sim_genome <- list(maternal = list(), paternal = list())
  snp_info <- list(maternal = list(), paternal = list())

  for(i in 1:length(seg_names)){
    ref_seq <- ref_genome[[ seg_names[i] ]]
    homo_snp <- snp_list[[ seg_names[i] ]]$`1|1`
    hetero_snp_maternal <- snp_list[[ seg_names[i] ]]$`1|0`
    hetero_snp_paternal <- snp_list[[ seg_names[i] ]]$`0|1`

    # identify duplicated snp positions
    all_pos <- c(homo_snp$POS, hetero_snp_maternal$POS, hetero_snp_paternal$POS)
    count_pos <- table(all_pos)
    dup_pos <- names(count_pos[count_pos > 1])


    # Function to insert SNPs into a genome sequence
    insert_snps_to_hap_segment <- function(ref_seq, homo_snp, hap_hetero_snp, dup_pos) {

      # combine homo and snp
      combined_snp <- rbind(homo_snp, hap_hetero_snp)

      # retrieve unique snp index (remove snps with duplicated positions)
      unique_snp_idx <- which(! combined_snp$POS %in% dup_pos)
      combined_snp <- combined_snp[unique_snp_idx,]

      hap_genome <- replaceLetterAt(x = ref_seq,
                                    at = combined_snp$POS,
                                    letter = combined_snp$ALT)

      hap_snp <- combined_snp

      return(list(genome = hap_genome, snp = hap_snp))
    }

    # insert SNPs into the genome
    for (haplotype in c("maternal", "paternal")) {

      if (haplotype == "maternal") {
        hap_hetero_snp <- hetero_snp_maternal
      } else if(haplotype == "paternal") {
        hap_hetero_snp <- hetero_snp_paternal
      }

      print(paste0("Insert snp to the ", seg_names[i], " of the ", haplotype, " genome."))
      result <- insert_snps_to_hap_segment(ref_seq = ref_seq,
                                           homo_snp = homo_snp,
                                           hap_hetero_snp = hap_hetero_snp,
                                           dup_pos = dup_pos)
      sim_genome[[haplotype]][[seg_names[i]]] <- result$genome
      snp_info[[haplotype]][[seg_names[i]]] <- result$snp
    }

  }

  return(list(sim_genome = sim_genome, snp_info = snp_info))
}


check_alt_snp_match <- function(seg_names, snp_list, sim_genome){

  for(i in 1:length(seg_names)){

    homo_snp <- snp_list[[ seg_names[i] ]]$`1|1`
    hetero_snp_maternal <- snp_list[[ seg_names[i] ]]$`1|0`
    hetero_snp_paternal <- snp_list[[ seg_names[i] ]]$`0|1`

    all_pos <- c(homo_snp$POS, hetero_snp_maternal$POS, hetero_snp_paternal$POS)

    # remove duplicated genotypes
    count_pos <- table(all_pos)
    dup_pos <- names(count_pos[count_pos > 1])

    for(haplotype in c("maternal", "paternal")){

      if (haplotype == "maternal") {
        hap_hetero_snp <- hetero_snp_maternal
      } else if(haplotype == "paternal") {
        hap_hetero_snp <- hetero_snp_paternal
      }


      combined_snp <- rbind(homo_snp, hap_hetero_snp)
      unique_snp_idx <- which(! combined_snp$POS %in% dup_pos)
      combined_snp <- combined_snp[unique_snp_idx, ]

      combined_snp_alt <- combined_snp$ALT

      # construct IRanges instance
      ir <- IRanges(start = combined_snp$POS, end = combined_snp$POS)

      # extract seq from sim genome
      sim_seq <- sim_genome[[haplotype]][[ seg_names[i] ]]
      sim_seq_alt <- extractAt(x = sim_seq, at = ir)

      # check if two refs are identifcal
      if(!all(sim_seq_alt == combined_snp$ALT)){
        warning("seq_alt and snp_alt don't match in ", haplotype, " ", seg_names[i], "!")
      }else{
        print(paste0("seq_alt and snp_alt matches in ", haplotype, " ", seg_names[i], "!"))
      }

    }
  }
}



