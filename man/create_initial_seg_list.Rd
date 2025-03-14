% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segment_relevant_functions.R
\name{create_initial_seg_list}
\alias{create_initial_seg_list}
\title{Create Initial Segment List from Chromosome Arm Table}
\usage{
create_initial_seg_list(chr_arm_table)
}
\arguments{
\item{chr_arm_table}{A data frame containing chromosome arm information with columns:
\itemize{
\item chrom: Chromosome name
\item start: Start position
\item end: End position
\item region_name: Identifier for the chromosome arm
}}
}
\value{
A list with two elements ("maternal" and "paternal"), each containing
a data frame of segments with the following columns:
\itemize{
\item All original columns from chr_arm_table
\item haplotype: "maternal" or "paternal"
\item ref_start: Reference start position
\item ref_end: Reference end position
\item ori_start: Original start position
\item ori_end: Original end position
\item copy_index: Set to 1 for initial segments
\item seg_id: Segment identifier (format: \verb{<region_name>_1})
\item CN_change: Copy number change (initialized to 0)
\item seg_source_edge: Set to "before_root"
\item seg_source_event: Set to "base"
}
}
\description{
Creates initial genomic segments for both maternal and paternal haplotypes from
a chromosome arm table. Handles coordinate system conversion (0-based to 1-based)
and initializes tracking metadata for each segment.
}
\details{
The function performs the following operations:
\enumerate{
\item Checks and converts coordinates from 0-based to 1-based if necessary
\item Creates two copies of the input table (maternal and paternal)
\item Adds metadata columns for tracking segment evolution
\item Uses dplyr for efficient data manipulation
}
}
\note{
Assumes input coordinates are either 0-based or 1-based.
Automatically converts 0-based coordinates to 1-based.
}
\examples{
\dontrun{
chr_arms <- data.frame(
  chrom = c("chr1", "chr1"),
  start = c(0, 1000000),
  end = c(1000000, 2000000),
  region_name = c("chr1p", "chr1q")
)
segments <- create_initial_seg_list(chr_arms)
}

}
