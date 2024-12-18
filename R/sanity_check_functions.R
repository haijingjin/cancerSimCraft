# sanity_check_functions_v5.R

segments_sanity_check <- function(segment_list, chr_lengths) {
  # Prepare a list to store the results
  failed_checks <- list()


  for (haplotype in c("maternal", "paternal")) {
    segment_table <- segment_list[[haplotype]]
    chr_names <- names(chr_lengths[[haplotype]])
    for (chr in chr_names) {
      # Extract segments for this chromosome and haplotype, excluding CN_change == -1
      chr_segments <- segment_table[segment_table$chrom == chr & segment_table$CN_change != -1, ]

      # Calculate total segment length
      total_seg_length <- sum(chr_segments$ref_end - chr_segments$ref_start + 1)
      # Compare with the chromosome length
      if (total_seg_length != chr_lengths[[haplotype]][chr]) {
        failed_checks[[paste(haplotype, chr)]] <- paste("Sanity check failed for", haplotype, chr,
                                                        "Expected:", chr_lengths[[haplotype]][chr],
                                                        "Actual:", total_seg_length)
      }
    }
  }

  return(failed_checks)
}

#' Check Nucleotide Content in Lost Segments
#'
#' @description
#' Validates that genomic regions marked as lost (CN_change = -1) contain only 'N'
#' nucleotides in the simulated genome sequence. This ensures proper handling of
#' deletion events.
#'
#' @param clone_genome A nested list containing genome sequences:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: chromosomes with nucleotide sequences
#'   }
#' @param clone_segments A nested list containing segment information:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: data frame with columns:
#'       \itemize{
#'         \item chrom - Chromosome name
#'         \item ori_start - Original start position
#'         \item ori_end - Original end position
#'         \item CN_change - Copy number change (-1 for losses)
#'         \item seg_id - Segment identifier
#'       }
#'   }
#'
#' @return A list where:
#'   \itemize{
#'     \item Names are in format "haplotype_segment_id"
#'     \item Values are nucleotide frequency counts (from Biostrings::alphabetFrequency)
#'     \item Empty list if no loss segments are found
#'   }
#'
#' @details
#' The function:
#' 1. Identifies segments with CN_change = -1
#' 2. For each lost segment:
#'    - Extracts the sequence
#'    - Counts nucleotide frequencies
#'    - Should show only 'N' nucleotides if loss is correctly implemented
#' 3. Prints message if no loss segments found in a haplotype
#'
#' @importFrom Biostrings alphabetFrequency
#' @export
#'
#' @examples
#' \dontrun{
#' # Check lost segments in a clone
#' loss_check <- check_loss_segments(
#'   clone_genome = synthesized_genome,
#'   clone_segments = segment_info
#' )
#' # Results should show only 'N' nucleotides in lost regions
#' }
check_loss_segments <- function(clone_genome, clone_segments){

  results <- list() # Initialize a list to store results
  for(haplotype in c("maternal", "paternal")){
    genome_seg <- clone_segments[[haplotype]]
    loss_segments <- genome_seg[genome_seg$CN_change == -1,]

    # Check if there are no loss segments
    if(nrow(loss_segments) == 0) {
      cat(paste("No loss segments found in", haplotype, "haplotype.\n"))
      next # Skip to the next haplotype if no loss segments found
    }

    for(i in 1:nrow(loss_segments)){
      chrom <- loss_segments$chrom[i]
      ori_start <- loss_segments$ori_start[i]
      ori_end <- loss_segments$ori_end[i]
      loss_seq <- clone_genome[[haplotype]][[chrom]][ori_start:ori_end]
      check_name <- paste(haplotype, loss_segments$seg_id[i], sep = "_")
      results[[check_name]] <- alphabetFrequency(loss_seq, baseOnly = FALSE)

    }

  }

  return(results)
}


#' Check Chromosome Lengths in Simulated Genome
#'
#' @description
#' Validates that the actual chromosome lengths in a simulated genome match their
#' expected lengths. Checks both maternal and paternal copies of each chromosome.
#'
#' @param clone_genome A nested list containing genome sequences:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: chromosomes with nucleotide sequences
#'   }
#' @param expected_chr_lengths A nested list containing expected chromosome lengths:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: named numeric vector of chromosome lengths
#'   }
#'
#' @return A nested list where:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: chromosomes
#'     \item Values are status messages indicating whether lengths match or not,
#'       including both expected and actual lengths when there's a mismatch
#'   }
#'
#' @details
#' The function:
#' 1. Checks each chromosome in both maternal and paternal copies
#' 2. Compares actual sequence length to expected length
#' 3. Generates informative messages about matches and mismatches
#' 4. Useful for validating genome synthesis results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' validation_results <- check_genome_chr_length(
#'   clone_genome = synthesized_genome,
#'   expected_chr_lengths = expected_lengths
#' )
#' # Access results for specific chromosomes
#' validation_results$maternal$chr1  # Check maternal chr1 status
#' }
check_genome_chr_length <- function(clone_genome, expected_chr_lengths) {
  # Initialize a list to hold the results
  results <- list()

  for(haplotype in c("maternal", "paternal")){
    for(chr in names(expected_chr_lengths[[haplotype]])){
      actual_length = length(clone_genome[[haplotype]][[chr]])
      expected_length = expected_chr_lengths[[haplotype]][chr]

      # Compare the actual length to the expected length
      if (actual_length == expected_length) {
        results[[haplotype]][[chr]] <- paste(chr, "in", haplotype, "matches the expected length of", expected_length)
      } else {
        results[[haplotype]][[chr]] <- paste(chr, "in", haplotype, "does NOT match the expected length. Expected:", expected_length, "Actual:", actual_length)
      }

    }

  }

  return(results)
}

#' Collect All Mutations for Each Clone in Tree
#'
#' @description
#' Gathers all mutations that occur along the evolutionary path from root to each clone
#' in the phylogenetic tree. This includes mutations from all ancestral branches.
#'
#' @param mutation_table A data frame containing mutation information with columns:
#'   \itemize{
#'     \item edge_name - Tree edge identifier in format "parent_child"
#'     \item other mutation-specific columns
#'   }
#' @param tree An igraph object representing the phylogenetic tree structure
#'
#' @return A list where:
#'   \itemize{
#'     \item Names are clone names (excluding root)
#'     \item Values are data frames containing all mutations affecting each clone,
#'       including those inherited from ancestors
#'   }
#'
#' @details
#' The function:
#' 1. Identifies the root node and all other clones in the tree
#' 2. For each clone:
#'    - Finds all edges in the path from root to clone
#'    - Collects mutations occurring on these edges
#'    - Returns mutations in their original table format
#'
#' @seealso
#' \code{\link{get_edges_between_clones}}
#'
#' @importFrom igraph V degree
#' @export
#'
#' @examples
#' \dontrun{
#' # For a tree with path Root -> A -> B
#' clone_muts <- acquire_clone_mutations(mutation_table, tree)
#' # clone_muts$B will contain mutations from both Root->A and A->B edges
#' }
acquire_clone_mutations <- function(mutation_table, tree) {

  root_name <- names(V(tree))[degree(tree, mode = "in") == 0]
  clone_names <- names(V(tree))[degree(tree, mode = "in") != 0]
  clone_mutation <- list()
  for(clone in clone_names){
    edges <- get_edges_between_clones(tree, upper_node = root_name, lower_node = clone)
    clone_mutation[[clone]] <- mutation_table[mutation_table$edge_name %in% edges,]
  }
  return(clone_mutation)
}


#' Validate Mutations in Simulated Genome
#'
#' @description
#' Verifies that mutations in the simulated genome match their expected alternative
#' nucleotides, accounting for lost segments that should contain 'N's.
#'
#' @param clone_genome A nested list containing genome sequences:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: chromosomes with nucleotide sequences
#'   }
#' @param clone_mutation_table A data frame containing mutation information:
#'   \itemize{
#'     \item haplotype - Maternal or paternal copy
#'     \item chrom - Chromosome name
#'     \item pos - Mutation position
#'     \item alternative_nt - Expected nucleotide after mutation
#'   }
#' @param clone_seg_list A nested list containing segment information:
#'   \itemize{
#'     \item First level: haplotypes
#'     \item Second level: data frame with segment details including loss information
#'   }
#'
#' @return A nested list where:
#'   \itemize{
#'     \item First level: haplotypes
#'     \item Second level: chromosomes
#'     \item Values are either:
#'       \itemize{
#'         \item Success message if all mutations are correct
#'         \item Data frame of incorrect mutations with details for investigation
#'       }
#'   }
#'
#' @details
#' The function:
#' 1. Processes mutations chromosome by chromosome
#' 2. Checks if mutations fall in lost segments (should be 'N')
#' 3. Compares actual nucleotides to expected ones
#' 4. Provides detailed output only for incorrect mutations
#' 5. Handles special cases:
#'    - Mutations in lost segments
#'    - No mutations in a chromosome
#'
#' @importFrom IRanges IRanges
#' @importFrom Biostrings extractAt
#' @export
#'
#' @examples
#' \dontrun{
#' validation_results <- check_genome_mutations(
#'   clone_genome = synthesized_genome,
#'   clone_mutation_table = mutations,
#'   clone_seg_list = segments
#' )
#' # Check results for specific chromosome
#' validation_results$maternal$chr1
#' }
check_genome_mutations <- function(clone_genome, clone_mutation_table, clone_seg_list) {
  results <- list()

  # Loop through each haplotype
  for (haplotype in c("maternal", "paternal")) {
    chr_names <- names(clone_genome[[haplotype]])
    for(chrom in chr_names){
      # filter mutations in current haplotype and chr
      selected_mutations <- clone_mutation_table[clone_mutation_table$haplotype == haplotype & clone_mutation_table$chrom == chrom,]
      mut_pos <- selected_mutations$pos
      expected_nts <- selected_mutations$alternative_nt

      if(nrow(selected_mutations) > 0){
        for(i in 1:nrow(selected_mutations)){
          mutation <- selected_mutations[i, ]
          genome_seg <- clone_seg_list[[haplotype]]

          all_loss_indices <- which(genome_seg$chrom == chrom &
                                      is.na(genome_seg$start) &
                                      is.na(genome_seg$end) &
                                      genome_seg$CN_change == -1)
          loss_segments <- genome_seg[all_loss_indices,]

          in_loss_segment <- any(
            loss_segments$chrom == mutation$chrom &
              loss_segments$haplotype == mutation$haplotype &
              loss_segments$ori_start <= mutation$pos &
              loss_segments$ori_end >= mutation$pos
          )

          if(in_loss_segment){
            expected_nts[i] <- "N"
          }
        }

        selected_seq <- clone_genome[[haplotype]][[chrom]]
        ir <- IRanges(start = mut_pos, end = mut_pos)
        actual_nts <- extractAt(x = selected_seq, at = ir)
        # Add actual nts to the selected mutations
        selected_mutations$actual_nts <- actual_nts

        # Compare expected vs. actual nucleotides
        correct <- expected_nts == actual_nts

        if (all(correct)) {
          # If all mutations are correct, return a concise message
          results[[haplotype]][[chrom]] <- paste("All mutations in", chrom, "of", haplotype, "are correct.")
        } else {
          # If there are incorrect mutations, return detailed info only for those
          incorrect_indices <- which(!correct)
          incorrect_mutations <- selected_mutations[incorrect_indices, ]
          results[[haplotype]][[chrom]] <- incorrect_mutations

        }
      }
    }
  }

  return(results)

}


validate_mutation_check <- function(check_output) {
  for (haplotype in c("maternal", "paternal")) {
    if (!is.null(check_output[[haplotype]])) {
      for (chr in names(check_output[[haplotype]])) {
        # Expected message
        expected_msg <- paste0("All mutations in ", chr, " of ", haplotype, " are correct.")
        if (check_output[[haplotype]][[chr]] != expected_msg) {
          return(FALSE) # If any entry is incorrect, return FALSE
        }
      }
    }
  }
  return(TRUE) # If all entries are correct, return TRUE
}
