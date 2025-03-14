% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sanity_check_functions.R
\name{segments_sanity_check}
\alias{segments_sanity_check}
\title{Perform Sanity Checks on Genomic Segments}
\usage{
segments_sanity_check(segment_list, chr_lengths)
}
\arguments{
\item{segment_list}{A list containing two data frames (maternal and paternal),
each with columns:
\itemize{
\item chrom - Chromosome name
\item ref_start - Segment start position
\item ref_end - Segment end position
\item CN_change - Copy number change (-1 indicates segments to exclude)
}}

\item{chr_lengths}{A list containing two named numeric vectors (maternal and paternal)
where names are chromosome identifiers and values are chromosome lengths}
}
\value{
A list where:
\itemize{
\item Names are concatenated strings of haplotype and chromosome (e.g., "maternal chr1")
\item Values are error messages for failed checks
\item Empty list indicates all checks passed
}
}
\description{
This function performs a sanity check on chromosome segments to ensure that the total length of the
segments matches the expected chromosome lengths for both maternal and paternal haplotypes. It checks each
chromosome and haplotype combination and reports any discrepancies.
}
\details{
The function performs the following checks for each chromosome in both haplotypes:
\enumerate{
\item Excludes segments with CN_change == -1
\item Sums the lengths of remaining segments (ref_end - ref_start + 1)
\item Compares total segment length with expected chromosome length
\item Records any mismatches in the returned list
}
}
\examples{
\dontrun{
# Create sample segment list
segments <- list(
  maternal = data.frame(
    chrom = c("chr1", "chr1"),
    ref_start = c(1, 101),
    ref_end = c(100, 200),
    CN_change = c(0, 1)
  ),
  paternal = data.frame(
    chrom = c("chr1", "chr1"),
    ref_start = c(1, 101),
    ref_end = c(100, 200),
    CN_change = c(0, 0)
  )
)

# Create sample chromosome lengths
chr_lengths <- list(
  maternal = c(chr1 = 200),
  paternal = c(chr1 = 200)
)

# Run sanity checks
results <- segments_sanity_check(segments, chr_lengths)
if (length(results) == 0) {
  print("All checks passed")
} else {
  print("Some checks failed:")
  print(results)
}
}

}
\seealso{
Related functions for segment manipulation and validation
}
