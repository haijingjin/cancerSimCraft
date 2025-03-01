

#' Convert Phased VCF Table to SNP Lists by Genotype
#'
#' @description
#' Takes a phased VCF table and converts it into a list of filtered SNP data frames,
#' separated by genotype categories (homozygous alternative and heterozygous variants).
#'
#' @param phased_vcf_table A data frame containing phased VCF data with at least the following columns:
#'   CHROM, POS, REF, ALT, and a sample-specific genotype column
#' @param sample_name Character string specifying the name of the sample column in the VCF table
#'
#' @return A list containing four data frames:
#'   \itemize{
#'     \item all: All SNPs (single-nucleotide variants only)
#'     \item 1|1: Homozygous alternative variants
#'     \item 1|0: Heterozygous variants with alternative allele on first haplotype
#'     \item 0|1: Heterozygous variants with alternative allele on second haplotype
#'   }
#'
#' @details
#' The function filters for single-nucleotide variants only (where REF and ALT are
#' single characters) and separates the genotype field (GT) from the dosage field (DS)
#' in the sample column.
#'
#' @importFrom tidyr separate
#' @importFrom dplyr select filter
#'
#' @examples
#' vcf_df <- data.frame(
#'   CHROM = c("chr1", "chr1"),
#'   POS = c(1000, 2000),
#'   REF = c("A", "C"),
#'   ALT = c("G", "T"),
#'   sample1 = c("1|0:0.5", "0|1:0.5")
#' )
#' snp_lists <- vcf_to_snp_list(vcf_df, "sample1")
#' @export
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


#' Check Reference Genome and SNP List Compatibility
#'
#' This function verifies that the reference nucleotides in the SNP list
#' match the corresponding positions in the reference genome. It is used to
#' ensure that the reference genome and SNP list are compatible before
#' introducing SNPs into the genome.
#'
#'
#' @param seg_names Character vector of segment/chromosome names to check
#' @param snp_list A nested list where first level contains segment names and second level
#'   contains a data frame named 'all', which includes all SNPs having at least POS and REF columns
#' @param ref_genome A DNAStringSet or similar object containing reference sequences,
#'   with names matching seg_names
#'
#' @return No return value, called for side effects:
#'   \itemize{
#'     \item Prints confirmation message for each matching segment
#'     \item Issues warning if mismatches are found
#'   }
#'
#' @details
#' The function performs the following steps for each segment:
#' 1. Creates an IRanges object from SNP positions
#' 2. Extracts reference sequences at those positions
#' 3. Compares extracted sequences with SNP reference alleles
#' 4. Reports matches/mismatches via messages and warnings
#'
#' @importFrom IRanges IRanges
#' @importFrom Biostrings extractAt
#'
#' @examples
#' \dontrun{
#' # Example reference genome
#' ref_genome <- DNAStringSet(c(
#'   chr1 = "ACTGACTGACTG",
#'   chr2 = "GTCAGTCAGTCA"
#' ))
#'
#' # Example SNP list
#' snp_list <- list(
#'   chr1 = list(all = data.frame(POS = c(1,5), REF = c("A","C"))),
#'   chr2 = list(all = data.frame(POS = c(2,6), REF = c("T","T")))
#' )
#'
#' check_ref_snp_match(c("chr1", "chr2"), snp_list, ref_genome)
#' }
#' @export
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

#' Insert SNPs into Reference Genome
#'
#'
#' #' This function introduces single nucleotide polymorphisms (SNPs) into a reference genome
#' based on a list of SNPs. The SNPs are categorized as homozygous (`1|1`) or heterozygous (`1|0` or `0|1`),
#' and the function ensures that SNPs at duplicated positions are handled appropriately.
#' The result is a synthetic genome with SNPs inserted for both maternal and paternal haplotypes.
#'
#'
#' @param seg_names Character vector of segment/chromosome names to process
#' @param snp_list A nested list containing SNP information for each segment, with
#'   sublists for different genotypes ('1|1', '1|0', '0|1'), each containing
#'   data frames with at least POS and ALT columns
#' @param ref_genome A DNAStringSet or similar object containing reference sequences,
#'   with names matching seg_names
#'
#' @return A list with two components:
#'   \itemize{
#'     \item sim_genome: List containing maternal and paternal haplotype sequences
#'     \item snp_info: List containing SNP information for each haplotype
#'   }
#'
#' @details
#' The function processes each segment as follows:
#' 1. Combines homozygous SNPs with appropriate heterozygous SNPs for each haplotype
#' 2. Identifies and removes SNPs with duplicated positions
#' 3. Inserts alternative alleles into the reference sequence
#' 4. Stores both modified sequences and SNP information
#'
#' Maternal haplotypes receive '1|0' heterozygous variants, while paternal
#' haplotypes receive '0|1' variants. Both haplotypes receive all homozygous
#' alternative ('1|1') variants.
#'
#' @importFrom Biostrings replaceLetterAt
#'
#' @examples
#' \dontrun{
#' # Example reference genome
#' ref_genome <- DNAStringSet(c(
#'   chr1 = "ACTGACTGACTG",
#'   chr2 = "GTCAGTCAGTCA"
#' ))
#'
#' # Example SNP list
#' snp_list <- list(
#'   chr1 = list(
#'     "1|1" = data.frame(POS = 1, ALT = "G"),
#'     "1|0" = data.frame(POS = 5, ALT = "T"),
#'     "0|1" = data.frame(POS = 8, ALT = "C")
#'   ),
#'   chr2 = list(
#'     "1|1" = data.frame(POS = 2, ALT = "A"),
#'     "1|0" = data.frame(POS = 6, ALT = "G"),
#'     "0|1" = data.frame(POS = 9, ALT = "T")
#'   )
#' )
#'
#' result <- insert_snps_to_genome(c("chr1", "chr2"), snp_list, ref_genome)
#' }
#' @export
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

#' Check Alternative Allele Matches Between SNP List and Simulated Genome
#'
#' @description
#' This function verifies that the alternate alleles (ALT) in the SNP list match
#' the corresponding positions in the synthetic genome. It ensures that the SNPs
#' have been correctly inserted into the synthetic genome for both maternal and
#' paternal haplotypes. The function checks each segment (e.g., chromosome) and
#' issues warnings if mismatches are found.
#'
#' @param seg_names Character vector of segment/chromosome names to check
#' @param snp_list A nested list where first level contains segment names and second level
#'   contains data frames for different genotypes ('1|1', '1|0', '0|1') with POS and ALT columns
#' @param sim_genome A nested list containing maternal and paternal sequences for each segment,
#'   structured as sim_genome$maternal[[seg_name]] and sim_genome$paternal[[seg_name]]
#'
#' @return No return value, called for side effects:
#'   \itemize{
#'     \item Prints confirmation message for each matching segment and haplotype
#'     \item Issues warning if mismatches are found
#'   }
#'
#' @details
#' For each segment and haplotype, the function:
#' 1. Combines appropriate homozygous and heterozygous SNPs
#' 2. Removes SNPs with duplicated positions
#' 3. Extracts sequences at SNP positions from simulated genome
#' 4. Compares extracted sequences with expected alternative alleles
#' 5. Reports matches/mismatches via messages and warnings
#'
#' Maternal haplotype checks include '1|1' and '1|0' variants, while paternal
#' haplotype checks include '1|1' and '0|1' variants.
#'
#' @importFrom IRanges IRanges
#' @importFrom Biostrings extractAt
#'
#' @examples
#' \dontrun{
#' # Example simulated genome
#' sim_genome <- list(
#'   maternal = list(chr1 = DNAStringSet("GCTGACTGACTG")),
#'   paternal = list(chr1 = DNAStringSet("ACTGACTCACTG"))
#' )
#'
#' # Example SNP list
#' snp_list <- list(
#'   chr1 = list(
#'     "1|1" = data.frame(POS = 1, ALT = "G"),
#'     "1|0" = data.frame(POS = 5, ALT = "T"),
#'     "0|1" = data.frame(POS = 8, ALT = "C")
#'   )
#' )
#'
#' check_alt_snp_match("chr1", snp_list, sim_genome)
#' }
#' @export
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



