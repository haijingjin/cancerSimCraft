% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genome_relevant_functions.R
\name{introduce_snv_sc}
\alias{introduce_snv_sc}
\title{Introduce Single Nucleotide Variants into a Single Cell Genome}
\usage{
introduce_snv_sc(input_genome, sc_mut_table_with_nt)
}
\arguments{
\item{input_genome}{List containing genome sequences to be modified}

\item{sc_mut_table_with_nt}{Data frame containing mutation details with columns:
haplotype, chrom, pos, original_nt, alternative_nt, time}
}
\value{
A list containing the modified genome sequences with all mutations incorporated
}
\description{
Incorporates single nucleotide variants (SNVs) into a genome sequence, handling
both regular mutations and recurrent mutations at the same genomic position.
}
\details{
This function introduces mutations into a genome sequence with special handling
for recurrent mutations (multiple mutations at the same genomic position):
\enumerate{
\item Identifies recurrent mutations by detecting duplicated positions across
haplotype, chromosome, and position
\item Processes non-recurrent mutations first using the insert_mutations() helper function
\item For recurrent mutations:
\itemize{
\item Groups mutations by their genomic coordinates
\item For each position with multiple mutations, creates a consolidated mutation that
applies the cumulative effect (using the original nucleotide from the first mutation
and the alternative nucleotide from the last mutation in the timeline)
\item Applies these consolidated recurrent mutations to the genome
}
}

This approach ensures that the final genome correctly represents the cumulative
effect of sequential mutations at the same position, rather than applying each
mutation independently which could lead to incorrect results.
}
\seealso{
\code{\link{insert_mutations}}, \code{\link{synth_sc_genome}}
}
