% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp_relevant_functions.R
\name{check_ref_snp_match}
\alias{check_ref_snp_match}
\title{Check Reference Genome and SNP List Compatibility}
\usage{
check_ref_snp_match(seg_names, snp_list, ref_genome)
}
\arguments{
\item{seg_names}{Character vector of segment/chromosome names to check}

\item{snp_list}{A nested list where first level contains segment names and second level
contains a data frame named 'all', which includes all SNPs having at least POS and REF columns}

\item{ref_genome}{A DNAStringSet or similar object containing reference sequences,
with names matching seg_names}
}
\value{
No return value, called for side effects:
\itemize{
\item Prints confirmation message for each matching segment
\item Issues warning if mismatches are found
}
}
\description{
This function verifies that the reference nucleotides in the SNP list
match the corresponding positions in the reference genome. It is used to
ensure that the reference genome and SNP list are compatible before
introducing SNPs into the genome.
}
\details{
The function performs the following steps for each segment:
\enumerate{
\item Creates an IRanges object from SNP positions
\item Extracts reference sequences at those positions
\item Compares extracted sequences with SNP reference alleles
\item Reports matches/mismatches via messages and warnings
}
}
\examples{
\dontrun{
# Example reference genome
ref_genome <- DNAStringSet(c(
  chr1 = "ACTGACTGACTG",
  chr2 = "GTCAGTCAGTCA"
))

# Example SNP list
snp_list <- list(
  chr1 = list(all = data.frame(POS = c(1,5), REF = c("A","C"))),
  chr2 = list(all = data.frame(POS = c(2,6), REF = c("T","T")))
)

check_ref_snp_match(c("chr1", "chr2"), snp_list, ref_genome)
}
}
