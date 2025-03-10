% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_sim.R
\name{merge_fasta_files}
\alias{merge_fasta_files}
\title{Merge Paternal and Maternal FASTA Files}
\usage{
merge_fasta_files(paternal_fa, maternal_fa, output_fa, tmp_dir = tempdir())
}
\arguments{
\item{paternal_fa}{Character string. Path to the paternal FASTA file.}

\item{maternal_fa}{Character string. Path to the maternal FASTA file.}

\item{output_fa}{Character string. Path where the merged FASTA file will be written.}

\item{tmp_dir}{Character string. Directory for temporary files. Default is tempdir().}
}
\value{
No return value, called for side effect of creating the merged FASTA file.
}
\description{
This function merges two FASTA files (paternal and maternal) into a single output file,
adding suffix identifiers to differentiate the sequences' origins.
}
\details{
The function adds "_paternal" suffix to sequence headers in the paternal file
and "_maternal" suffix to sequence headers in the maternal file before
concatenating them into a single output file.
}
\examples{
\dontrun{
merge_fasta_files(
  paternal_fa = "path/to/paternal.fa",
  maternal_fa = "path/to/maternal.fa",
  output_fa = "path/to/merged.fa"
)
}

}
