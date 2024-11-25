# ground_truth_check_v1.R

process_vcf <- function(vcf) {
  # Convert the vcf object to a data frame with genotype information
  new_vcf <- as.data.frame(cbind(vcf@fix, vcf@gt))

  # Split the genotype column into separate columns and calculate REF_AD, ALT_AD
  new_vcf <- new_vcf %>%
    tidyr::separate(col = colnames(.)[ncol(.)], into = c("GT", "PL", "DP", "AD"), sep = ":") %>%
    mutate(REF_AD = as.numeric(sub(",.*", "", AD)),
           ALT_AD = as.numeric(sub(".*,", "", AD)))

  # Calculate B-allele frequency (BAF)
  new_vcf <- new_vcf %>%
    mutate(BAF = ALT_AD / (REF_AD + ALT_AD))

  new_vcf$DP <- as.integer(new_vcf$DP)
  new_vcf$POS <- as.integer(new_vcf$POS)

  return(new_vcf)
}
