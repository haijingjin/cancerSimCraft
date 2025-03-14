% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/illustration_functions.R
\name{draw_chr_bar}
\alias{draw_chr_bar}
\title{Create Chromosome Bar Annotation for ComplexHeatmap}
\usage{
draw_chr_bar(window_data)
}
\arguments{
\item{window_data}{A data frame containing genomic window information with at least one
column named 'chr' (case-insensitive) containing chromosome identifiers
(e.g., "chr1", "chr2", etc.).}
}
\value{
A HeatmapAnnotation object from the ComplexHeatmap package containing:
\itemize{
\item A text annotation track with chromosome numbers
\item A binary pattern track with alternating black and grey bars
}
}
\description{
Creates a chromosome bar annotation track for ComplexHeatmap visualization, which includes
alternating black and grey bars to distinguish chromosomes and chromosome number labels.
The function handles both standard (22 autosomes + X, Y) and custom chromosome sets.
}
\details{
The function performs the following steps:
\itemize{
\item Standardizes chromosome column names to lowercase
\item Computes chromosome lengths and boundaries
\item Creates alternating binary pattern for visual separation
\item Calculates mean positions for chromosome labels
\item Generates chromosome text annotations
}
}
\examples{
\dontrun{
# Create sample window data
window_data <- data.frame(
  chr = rep(paste0("chr", c(1:22, "X")), each = 100),
  start = 1:2300,
  end = 101:2400
)

# Generate chromosome bar annotation
chr_annotation <- draw_chr_bar(window_data)

# Use in ComplexHeatmap
library(ComplexHeatmap)
Heatmap(matrix_data,
        top_annotation = chr_annotation)
}

}
