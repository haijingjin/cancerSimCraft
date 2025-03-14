% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pattern_extraction.R
\name{subset_interested_profiles}
\alias{subset_interested_profiles}
\title{Subset and Transform Haplotype Copy Number Profiles for Regions of Interest}
\usage{
subset_interested_profiles(
  haplotype_cn,
  haplotype,
  chr,
  interested_region_indices
)
}
\arguments{
\item{haplotype_cn}{A list containing haplotype-specific copy number matrices}

\item{haplotype}{Character string specifying which haplotype to analyze
(must be a name in the haplotype_cn list)}

\item{chr}{Character string specifying the chromosome of interest
(e.g., "chr1", "chrX")}

\item{interested_region_indices}{Numeric vector of indices specifying which
regions within the chromosome to analyze}
}
\value{
A numeric matrix with the same dimensions as the original copy number
matrix, where:
\itemize{
\item Only specified regions of interest contain copy number values
\item All other positions contain NA
\item Column names are preserved from the original matrix
}
}
\description{
Extracts copy number data for specific regions of interest from a haplotype-specific
copy number profile and returns a transformed matrix with the same dimensions as
the original data, where only the regions of interest contain values.
}
\details{
The function performs these steps:
\itemize{
\item Identifies regions in the specified chromosome
\item Extracts copy number data for regions of interest
\item Creates a new matrix with same dimensions as original data
\item Fills in only the specified regions, leaving others as NA
}
}
\examples{
# Create sample haplotype copy number data
sample_cn <- list(
  hap1 = matrix(1:20, nrow = 2, ncol = 10,
                dimnames = list(NULL,
                paste0("chr1_", 1:10))),
  hap2 = matrix(2:21, nrow = 2, ncol = 10,
                dimnames = list(NULL,
                paste0("chr1_", 1:10)))
)

# Extract regions of interest
subset_data <- subset_interested_profiles(
  haplotype_cn = sample_cn,
  haplotype = "hap1",
  chr = "chr1",
  interested_region_indices = 2:4
)

}
\seealso{
Other functions for manipulating copy number profiles in the package
}
