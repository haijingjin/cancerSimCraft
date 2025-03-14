% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp_relevant_functions.R
\name{check_alt_snp_match}
\alias{check_alt_snp_match}
\title{Check Alternative Allele Matches Between SNP List and Simulated Genome}
\usage{
check_alt_snp_match(seg_names, snp_list, sim_genome)
}
\arguments{
\item{seg_names}{Character vector of segment/chromosome names to check}

\item{snp_list}{A nested list where first level contains segment names and second level
contains data frames for different genotypes ('1|1', '1|0', '0|1') with POS and ALT columns}

\item{sim_genome}{A nested list containing maternal and paternal sequences for each segment,
structured as sim_genome$maternal[\link{seg_name}] and sim_genome$paternal[\link{seg_name}]}
}
\value{
No return value, called for side effects:
\itemize{
\item Prints confirmation message for each matching segment and haplotype
\item Issues warning if mismatches are found
}
}
\description{
This function verifies that the alternate alleles (ALT) in the SNP list match
the corresponding positions in the synthetic genome. It ensures that the SNPs
have been correctly inserted into the synthetic genome for both maternal and
paternal haplotypes. The function checks each segment (e.g., chromosome) and
issues warnings if mismatches are found.
}
\details{
For each segment and haplotype, the function:
\enumerate{
\item Combines appropriate homozygous and heterozygous SNPs
\item Removes SNPs with duplicated positions
\item Extracts sequences at SNP positions from simulated genome
\item Compares extracted sequences with expected alternative alleles
\item Reports matches/mismatches via messages and warnings
}

Maternal haplotype checks include '1|1' and '1|0' variants, while paternal
haplotype checks include '1|1' and '0|1' variants.
}
\examples{
\dontrun{
# Example simulated genome
sim_genome <- list(
  maternal = list(chr1 = DNAStringSet("GCTGACTGACTG")),
  paternal = list(chr1 = DNAStringSet("ACTGACTCACTG"))
)

# Example SNP list
snp_list <- list(
  chr1 = list(
    "1|1" = data.frame(POS = 1, ALT = "G"),
    "1|0" = data.frame(POS = 5, ALT = "T"),
    "0|1" = data.frame(POS = 8, ALT = "C")
  )
)

check_alt_snp_match("chr1", snp_list, sim_genome)
}
}
