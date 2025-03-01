
#' Calculate Number of SNVs in a Clone
#'
#' @description
#' Calculates the expected number of single nucleotide variants (SNVs) accumulated
#' in a clone over time using a simple linear model where the number of mutations is
#' proportional to time, genome length, and a base mutation rate.
#'
#' @param time Numeric value representing time (units should be consistent with mutation_rate)
#' @param genome_length Numeric value representing the length of the genome in base pairs
#' @param mutation_rate Numeric value representing the mutation rate per base pair per time unit
#'
#' @return Integer number of mutations, rounded to the nearest whole number
#'
#' @details
#' The function implements a simple linear model where the number of mutations is
#' directly proportional to time, genome length, and mutation rate. The result is
#' rounded to the nearest integer.
#'
#' @examples
#' # Calculate mutations for 100 time units, genome length of 1e6, and mutation rate of 1e-8
#' calc_clone_snv_num(100, 1e6, 1e-8)
#'
#' @export
calc_clone_snv_num <- function(time, genome_length, mutation_rate) {
  # A simple model where mutation number is proportional to time, genome length, and a base mutation rate
  num_mutations <- time * genome_length * mutation_rate
  return(round(num_mutations, digits = 0))
}

#' Get Mutations for Sampled Single Cells
#'
#' @description
#' Extracts and organizes mutation information for specified sampled cells, including
#' mutations in their ancestral lineage, and identifies recurrent mutations.
#'
#' @param cell_info Data frame containing cell lineage information with columns:
#'        clone, parent, birth_time, death_time, cell_index
#' @param mutation_info Data frame containing mutation information with columns:
#'        clone, cell_index, haplotype, chrom, pos, time
#' @param sampled_cell_idx Vector of cell indices for which to retrieve mutation information
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item sampled_sc_mutations: List where each element corresponds to a sampled cell,
#'         containing a data frame of mutations for that cell and its ancestors
#'   \item sampled_mutation_table: Data frame containing all unique mutations across all
#'         sampled cells and their ancestors
#'   \item recurrent_mutation_tracker: List tracking recurrent mutations found across
#'         multiple cells or lineages
#' }
#'
#' @details
#' This function reconstructs the complete mutation profile for each sampled cell by:
#'
#' 1. Identifying all ancestral cells for each sampled cell using get_sc_ancestors()
#' 2. Collecting mutations from the sampled cell and all its ancestors
#' 3. Identifying recurrent mutations (mutations that occur multiple times independently)
#' 4. Grouping and processing recurrent mutations for tracking
#' 5. Removing duplicate mutations from the final table
#'
#' The function requires the following helper functions to be defined:
#' - get_sc_ancestors(): To identify ancestral cells
#' - identify_recurrent_mutations(): To find mutations that appear multiple times
#' - group_recurrent_mutations(): To group similar recurrent mutations
#' - process_recurrent_mutation(): To track and analyze recurrent mutations
#'
#' @examples
#' # Create sample cell_info and mutation_info data frames
#' cell_info <- data.frame(
#'   clone = c("A", "A", "B", "A", "B"),
#'   parent = c(NA, 1, 1, 2, 3),
#'   birth_time = c(0, 10, 10, 15, 20),
#'   death_time = c(10, NA, NA, NA, NA),
#'   cell_index = 1:5
#' )
#'
#' mutation_info <- data.frame(
#'   clone = c("A", "A", "B", "A"),
#'   cell_index = c(1, 2, 3, 4),
#'   haplotype = c("hap1", "hap1", "hap2", "hap1"),
#'   chrom = c("chr1", "chr1", "chr2", "chr1"),
#'   pos = c(1000, 2000, 1500, 1000),
#'   time = c(5, 12, 15, 18)
#' )
#'
#' # Get mutations for cells 4 and 5
#' mutations <- get_mutations_sc(cell_info, mutation_info, c(4, 5))
#'
#' # View mutation summary
#' print(mutations$sampled_mutation_table)
#'
#' @seealso
#' \code{\link{simulate_sc_dynamics}}, \code{\link{get_sc_ancestors}}
#'
#' @export
get_mutations_sc <- function(cell_info, mutation_info, sampled_cell_idx = NULL) {

  if (is.null(sampled_cell_idx)) {
    stop("sampled_cell_idx must be provided.")
  }

  # Initialize an empty data frame to store the mutations of all sampled cells
  sampled_sc_mutations <- list()
  sampled_mutation_table <- data.frame()
  recurrent_mutation_tracker <- list()
  # For each sampled cell
  for (i in 1:length(sampled_cell_idx)) {
    cell_index <- sampled_cell_idx[i]
    ancestor_indices <- get_sc_ancestors(cell_info, cell_index)

    # Get all mutations that belong to these cells
    cell_mutations <- mutation_info[mutation_info$cell_index %in% c(cell_index, ancestor_indices), ]
    recurrent_mutations <- identify_recurrent_mutations(cell_mutations)
    grouped_recurrent_mutations <- group_recurrent_mutations(recurrent_mutations)

    # Assuming grouped_recurrent_mutations is a list where each element is a data frame of recurrent mutations
    for (mutation_group_key in names(grouped_recurrent_mutations)) {
      recurrent_mutation_set <- grouped_recurrent_mutations[[mutation_group_key]]
      recurrent_mutation_tracker <- process_recurrent_mutation(mutation = recurrent_mutation_set,
                                                               mutation_key = mutation_group_key,
                                                               tracker = recurrent_mutation_tracker)
    }

    # Add the mutations to the data frame
    sampled_mutation_table <- rbind(sampled_mutation_table, cell_mutations)
    # Add the mutations to the list under the cell's index
    sampled_sc_mutations[[as.character(cell_index)]] <- cell_mutations
  }
  # Remove duplicate mutations
  sampled_mutation_table <- sampled_mutation_table[!duplicated(sampled_mutation_table), ]
  # Return the unique mutation table
  return(list(sampled_sc_mutations = sampled_sc_mutations, sampled_mutation_table = sampled_mutation_table, recurrent_mutation_tracker = recurrent_mutation_tracker))
}

#' Identify Recurrent Mutations in Cell Lineage Data
#'
#' @description
#' Identifies mutations that occur multiple times at the same genomic location
#' across different cells in a lineage.
#'
#' @param cell_mutations Data frame containing mutation information with columns:
#'        clone, cell_index, haplotype, chrom, pos, time
#'
#' @return A data frame containing only the recurrent mutations (mutations that
#'         appear more than once at the same genomic location)
#'
#' @details
#' This function identifies recurrent mutations by:
#'
#' 1. Creating a unique key for each mutation based on its haplotype, chromosome, and position
#' 2. Identifying duplicate keys, which indicate the same mutation occurring multiple times
#' 3. Extracting and returning only the recurrent mutations
#'
#'
#' @seealso
#' \code{\link{get_mutations_sc}}, \code{\link{group_recurrent_mutations}}
#'
#' @export
identify_recurrent_mutations <- function(cell_mutations) {
  # Create a key for each mutation based on haplotype, chromosome, and position
  mutation_keys <- paste(cell_mutations$haplotype, cell_mutations$chrom, cell_mutations$pos, sep = "_")

  # Identify duplicates in these keys, which indicate recurrent mutations
  duplicate_indices <- which(duplicated(mutation_keys) | duplicated(mutation_keys, fromLast = TRUE))

  # Extract the recurrent mutations based on these indices
  recurrent_mutations <- cell_mutations[duplicate_indices, ]

  return(recurrent_mutations)
}

#' Group Recurrent Mutations by Genomic Location
#'
#' @description
#' Organizes recurrent mutations by grouping them according to their genomic location
#' (haplotype, chromosome, and position).
#'
#' @param recurrent_mutations Data frame containing mutation information with columns:
#'        clone, cell_index, haplotype, chrom, pos, time
#'
#' @return A list where each element is a data frame containing all mutations that
#'         occurred at the same genomic location. List names are constructed as
#'         "haplotype_chromosome_position".
#'
#' @details
#' This function takes a data frame of recurrent mutations and organizes them into groups
#' based on their genomic coordinates. It:
#'
#' 1. Creates a unique key for each mutation combining haplotype, chromosome, and position
#' 2. Splits the data frame into a list of smaller data frames, each containing
#'    all mutations that occurred at the same genomic location
#'
#' This grouping is useful for analyzing patterns of mutations at specific sites
#' and for further processing of recurrent mutations in evolutionary analyses.
#'
#'
#' @seealso
#' \code{\link{identify_recurrent_mutations}}, \code{\link{process_recurrent_mutation}}
#'
#' @export
group_recurrent_mutations <- function(recurrent_mutations) {
  # Create a key for grouping mutations
  group_keys <- paste(recurrent_mutations$haplotype, recurrent_mutations$chrom, recurrent_mutations$pos, sep = "_")

  # Split the data frame into a list of data frames, each containing mutations at the same position
  grouped_recurrent_mutations <- split(recurrent_mutations, group_keys)

  return(grouped_recurrent_mutations)
}

#' Process and Track Recurrent Mutations
#'
#' @description
#' Updates a tracker structure that organizes and monitors recurrent mutations
#' across different lineages, categorizing them as identical, inclusive, or branching.
#'
#' @param mutation Data frame containing information about a specific set of recurrent mutations
#'        at the same genomic location
#' @param mutation_key Character string uniquely identifying the genomic location of the mutation
#'        (typically in format "haplotype_chromosome_position")
#' @param tracker List structure that tracks recurrent mutations by organizing them into
#'        sets based on their occurrence patterns
#'
#' @return Updated tracker list with the mutation properly categorized and stored
#'
#' @details
#' This function processes recurrent mutations by comparing them to existing mutation sets in a tracker.
#' It handles three scenarios:
#'
#' 1. If the mutation set already exists and is identical, it returns the tracker unchanged.
#' 2. If the mutation set is a superset of an existing set, it updates the existing set with new mutations.
#' 3. If the mutation set is distinct from existing sets, it creates a new branch in the tracker.
#'
#'
#' The function prints the relationship category for debugging and monitoring purposes.
#'
#' @seealso
#' \code{\link{group_recurrent_mutations}}, \code{\link{get_mutations_sc}}
#'
#' @export
process_recurrent_mutation <- function(mutation, mutation_key, tracker) {
  if (mutation_key %in% names(tracker)) {
    mutation_sets <- tracker[[mutation_key]]

    for (i in 1:length(mutation_sets)) {
      set <- mutation_sets[[i]]
      mutation_times <- mutation$time
      set_times <- set$time

      if (all(set_times == mutation_times) || all(mutation_times %in% set_times)) {
        print("Old & Same")
        return(tracker)
      } else if (!all(set_times == mutation_times) && all(set_times %in% mutation_times)) {
        print("Old & Inclusive")
        new_mut_indices <- which(!mutation_times %in% set_times)
        new_mutations <- mutation[new_mut_indices, ]
        tracker[[mutation_key]][[i]] <- rbind(tracker[[mutation_key]][[i]], new_mutations)
      } else if (!all(mutation_times %in% set_times) && !all(set_times %in% mutation_times) && i == length(mutation_sets)) {
        print("Branching")
        new_idx <- i + 1
        tracker[[mutation_key]][[new_idx]] <- mutation
      }
    }
  } else {
    # Mutation key does not exist in the tracker, create new set
    tracker[[mutation_key]][[1]] <- mutation
  }
  return(tracker)
}

#' Simulate Nucleotide Change for a Single Mutation
#'
#' @description
#' Simulates the specific nucleotide change for a single mutation based on genomic
#' context, handling both regular and recurrent mutations, as well as mutations in
#' lost segments.
#'
#' @param single_mutation Data frame row containing information about a single mutation
#' @param seg_list List structure containing segment information for genome regions
#' @param genome_sequence List of DNA sequences representing the reference genome
#' @param nt_transition_matrix Matrix specifying nucleotide transition probabilities
#' @param recurrent Logical indicating whether this is a recurrent mutation (default: FALSE)
#' @param original_nt Character specifying the original nucleotide for recurrent mutations
#'        (required when recurrent = TRUE, default: NA)
#'
#' @return The updated single_mutation row with additional fields:
#' \itemize{
#'   \item seg_id: Segment identifier for the mutation
#'   \item ref_pos: Reference position within the segment
#'   \item original_nt: Original nucleotide at the mutation site
#'   \item alternative_nt: Mutated nucleotide
#'   \item processed: Set to TRUE indicating the mutation has been processed
#' }
#'
#' @details
#' This function simulates nucleotide-level details for a single mutation by:
#'
#' 1. Identifying the genomic segment containing the mutation
#' 2. Determining if the segment has been lost through deletion or other structural variants
#' 3. For mutations in lost segments:
#'    - Setting both original and alternative nucleotides to "N"
#' 4. For mutations in non-lost segments:
#'    - For regular mutations: Looking up the original nucleotide from the reference genome
#'      and sampling an alternative nucleotide based on the transition matrix
#'    - For recurrent mutations: Using the provided original nucleotide and sampling
#'      an alternative nucleotide based on the transition matrix
#'
#' The function requires the helper function get_segment_info() to map the mutation
#' coordinates to the appropriate genomic segment.
#'
#' @seealso
#' \code{\link{get_segment_info}}, \code{\link{sim_snv_nt_sc}}
#'
#' @export
simulate_single_nt_change <- function(single_mutation, seg_list, genome_sequence, nt_transition_matrix, recurrent = FALSE, original_nt = NA) {
  haplotype <- single_mutation$haplotype
  chrom <- single_mutation$chrom
  clone <- single_mutation$clone
  pos <- single_mutation$pos
  genome_seg <- seg_list[[clone]][[haplotype]]

  # Check if the found mutation segment is lost
  all_loss_indices <- which(is.na(genome_seg$start) &
                              is.na(genome_seg$end) &
                              genome_seg$CN_change == -1)

  loss_segments <- genome_seg[all_loss_indices,]

  in_loss_segment <- any(
    loss_segments$chrom == chrom &
      loss_segments$haplotype == haplotype &
      loss_segments$ori_start <= pos &
      loss_segments$ori_end >= pos
  )

  segment <- get_segment_info(seg_list, single_mutation, in_loss_segment)
  if (is.null(segment)) {
    stop("Error: No segment is found for mutation: ", single_mutation)
  }

  # Check if the segment is lost?
  if(in_loss_segment){
    single_mutation$seg_id <- segment$seg_id
    single_mutation$ref_pos <- pos - segment$ori_start + segment$ref_start
    single_mutation$original_nt <- "N"
    single_mutation$alternative_nt <- "N"
    single_mutation$processed <- TRUE
  }else{
    # Process non-loss segments
    single_mutation$seg_id <- segment$seg_id
    single_mutation$ref_pos <- pos - segment$start + segment$ref_start
    if(recurrent == FALSE){
      # process non-recurrent mutation
      single_mutation$original_nt <- as.character(genome_sequence[[haplotype]][[chrom]][single_mutation$ref_pos])
      single_mutation$alternative_nt <- sample(colnames(nt_transition_matrix), 1, prob = nt_transition_matrix[single_mutation$original_nt, ])
      single_mutation$processed <- TRUE
    }else{
      # process recurrent mutation
      if(is.na(original_nt)){
        stop("original_nt is required for processing recurrent mutations.")
      }
      single_mutation$original_nt <- original_nt
      single_mutation$alternative_nt <- sample(colnames(nt_transition_matrix), 1, prob = nt_transition_matrix[single_mutation$original_nt, ])
      single_mutation$processed <- TRUE
    }
  }

  return(single_mutation)
}

#' Find Index of an Identical Row in a Data Frame
#'
#' @description
#' Searches a data frame for a row that is identical to the provided row
#' and returns its index.
#'
#' @param row A data frame row or vector to search for
#' @param dataframe The data frame to search within
#'
#' @return An integer representing the index of the first identical row found,
#'         or NA if no identical row exists in the data frame
#'
#' @details
#' This function iterates through each row of the provided data frame and compares it
#' with the target row using the all.equal function. The comparison ignores attributes
#' to focus on the content of the data.
#'
#' The function is particularly useful for finding specific mutations or cells in larger
#' data frames when the exact row index is not known but the content is.
#'
#' @seealso
#' \code{\link{sim_snv_nt_sc}}
#'
#' @export
find_identical_row_index <- function(row, dataframe) {
  for (i in 1:nrow(dataframe)) {
    comparison_result <- all.equal(row, dataframe[i, , drop = FALSE], check.attributes = FALSE)
    if (isTRUE(comparison_result)) {
      return(i)
    }
  }
  return(NA)  # Return NA if no identical row is found
}


#' Simulate Single Nucleotide Variants for Single Cell Data
#'
#' @description
#' Simulates nucleotide-level details for single nucleotide variants (SNVs) in single cell
#' mutation data, handling both regular and recurrent mutations with appropriate nucleotide changes.
#'
#' @param genome_sequence List of DNA sequences representing the reference genome
#' @param seg_list Data frame or list containing genomic segment information with mapping between
#'        chromosome coordinates and segment identifiers
#' @param mutation_table Data frame containing mutation information with columns:
#'        clone, cell_index, haplotype, chrom, pos, time
#' @param recurrent_mutation_tracker List structure tracking recurrent mutations, organized by
#'        genomic location and mutation sets (default is NA for no recurrent mutations)
#' @param nt_transition_matrix Matrix specifying nucleotide transition probabilities for
#'        different mutation types
#'
#' @return An extended mutation_table data frame with additional columns:
#' \itemize{
#'   \item seg_id: Segment identifier for the mutation
#'   \item ref_pos: Reference position within the segment
#'   \item original_nt: Original nucleotide at the mutation site
#'   \item alternative_nt: Mutated nucleotide
#'   \item processed: Logical flag indicating whether the mutation has been processed
#' }
#'
#' @details
#' This function simulates the nucleotide-level details of mutations by:
#'
#' 1. Extending the mutation table with columns for segment ID, reference position,
#'    original nucleotide, alternative nucleotide, and processing status
#' 2. Processing recurrent mutations first, maintaining proper nucleotide changes across
#'    mutation sets (later mutations in a set build upon earlier ones)
#' 3. Processing regular (non-recurrent) mutations
#'
#' For recurrent mutations, the function ensures that:
#' - The first mutation in a set is processed normally
#' - Subsequent mutations in the set use the alternative nucleotide from the previous
#'   mutation as their original nucleotide
#'
#' The function relies on helper functions:
#' - find_identical_row_index(): To locate specific mutations in the table
#' - simulate_single_nt_change(): To determine nucleotide changes for individual mutations
#'
#' @seealso
#' \code{\link{simulate_single_nt_change}}, \code{\link{find_identical_row_index}},
#' \code{\link{process_recurrent_mutation}}
#'
#' @export
sim_snv_nt_sc <- function(genome_sequence, seg_list, mutation_table, recurrent_mutation_tracker = NA, nt_transition_matrix){

  original_mut_ncol <- ncol(mutation_table)
  # Initialize vectors
  seg_id <- character(nrow(mutation_table))
  ref_pos <- numeric(nrow(mutation_table))
  original_nt <- character(nrow(mutation_table))
  alternative_nt <- character(nrow(mutation_table))
  processed <- logical(nrow(mutation_table))

  # Create an extended mutation table to fill
  mutation_table <- cbind(mutation_table,
                          seg_id = seg_id,
                          ref_pos = ref_pos,
                          original_nt = original_nt,
                          alternative_nt = alternative_nt,
                          processed = processed)

  # Construct a unique key for each mutation's location (haplotype_chromosome_position)
  mutation_keys <- paste(mutation_table$haplotype, mutation_table$chrom, mutation_table$pos, sep = "_")

  # Get the names of recurrent mutations
  recurrent_keys <- names(recurrent_mutation_tracker)

  # Process reccurent mutations
  for(mutation_key in recurrent_keys){
    mutation_sets <- recurrent_mutation_tracker[[mutation_key]]

    # Process from the first set
    for(i in 1:length(mutation_sets)){
      set <- mutation_sets[[i]]
      for(j in 1:nrow(set)){
        # find the identical row
        recurrent_mut_idx <- find_identical_row_index(set[j, ], mutation_table[,1:original_mut_ncol])
        mut_row <- mutation_table[recurrent_mut_idx,]
        if(j == 1 && mut_row$processed == FALSE){
          # retrieve the remaining info from genome sequence and seg_list
          # simulate original nts
          new_mut_row <- simulate_single_nt_change(single_mutation = mut_row,
                                                      seg_list = seg_list,
                                                      genome_sequence = genome_sequence,
                                                      nt_transition_matrix = nt_transition_matrix)
          mutation_table[recurrent_mut_idx, ] = new_mut_row
        }else if(j > 1 && mut_row$processed == FALSE){
          # for the processed mutation and the second mutation in the set
          # original_nt is the previous alternative_nt
          previous_mut_idx <- find_identical_row_index(set[(j-1), ], mutation_table[,1:original_mut_ncol])
          original_nt <- mutation_table[previous_mut_idx,]$alternative_nt

          new_mut_row <- simulate_single_nt_change(single_mutation = mut_row,
                                                      seg_list = seg_list,
                                                      genome_sequence = genome_sequence,
                                                      nt_transition_matrix = nt_transition_matrix,
                                                      recurrent = TRUE,
                                                      original_nt = original_nt)

          mutation_table[recurrent_mut_idx, ] = new_mut_row
        } # if first or remaining
      } # mutations in a set
    } # mutation set
  } # mutation key # Process remaining non-recurrent mutations

  # Identify mutations that are not processed and are not in the recurrent set
  regular_mut_indices <- which(mutation_table$processed == FALSE & !(mutation_keys %in% recurrent_keys))

  for(regular_mut_idx in regular_mut_indices){
    mut_row <- mutation_table[regular_mut_idx,]
    new_mut_row <- simulate_single_nt_change(single_mutation = mut_row,
                                                seg_list = seg_list,
                                                genome_sequence = genome_sequence,
                                                nt_transition_matrix = nt_transition_matrix)
    mutation_table[regular_mut_idx, ] = new_mut_row
  }

  return(mutation_table)
}


get_mutations <- function(cell_info, mutation_info, sampled_cell_idx = NULL, num_samples = NULL) {
  if (!is.null(num_samples)) {
    # Generate sampled_cells using the sample_cells function if num_samples is provided
    sampled_cell_idx <- sample_cells(cell_info, num_samples)$cell_index
  } else if (is.null(sampled_cell_idx)) {
    stop("Either 'sampled_cells' or 'num_samples' must be provided.")
  }

  # Initialize an empty data frame to store the mutations of all sampled cells
  all_sampled_mutations <- data.frame()

  # For each sampled cell
  for (i in 1:length(sampled_cell_idx)) {
    # Get the cell index
    cell_index <- sampled_cell_idx[i]
    # Get the indices of all ancestor cells
    ancestor_indices <- get_ancestors(cell_info, cell_index)

    # Get all mutations that belong to these cells
    cell_mutations <- mutation_info[mutation_info$cell_index %in% c(cell_index, ancestor_indices), ]

    # Add the mutations to the data frame
    all_sampled_mutations <- rbind(all_sampled_mutations, cell_mutations)
  }

  # Return the data frame
  return(list(all_sampled_mutations = all_sampled_mutations, sampled_cell_idx = sampled_cell_idx))
}


#' Simulate Nucleotide Changes for Clonal Mutations
#'
#' @description
#' Simulates nucleotide changes for mutations in a clonal evolutionary tree, handling both
#' regular and recurrent mutations. Takes into account chromosome segment information and
#' possible loss events.
#'
#' @param genome_sequence A nested list containing reference genome sequences:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: chromosomes
#'     \item Each chromosome contains nucleotide sequence
#'   }
#' @param seg_list A nested list containing segment information for each clone and haplotype,
#'   including coordinates and copy number changes.
#' @param mutation_info A data frame containing mutation locations with columns:
#'   \itemize{
#'     \item clone - Clone name where mutation occurs
#'     \item haplotype - Maternal or paternal copy
#'     \item chrom - Chromosome name
#'     \item pos - Position on chromosome
#'   }
#' @param nt_transition_matrix A 5x5 matrix containing nucleotide transition probabilities
#' @param tree An igraph object representing the phylogenetic tree structure
#'
#' @return An updated mutation_info data frame with additional columns:
#'   \itemize{
#'     \item seg_id - Segment identifier where mutation occurs
#'     \item ref_pos - Position in reference coordinates
#'     \item original_nt - Original nucleotide
#'     \item alternative_nt - Mutated nucleotide
#'   }
#'
#' @details
#' The function processes mutations in two categories:
#' 1. Regular mutations (first occurrence at a position):
#'    - Identifies segment containing the mutation
#'    - Handles lost segments (marked with 'N')
#'    - Samples alternative nucleotide based on transition matrix
#' 2. Recurrent mutations (at previously mutated positions):
#'    - Finds previous mutation in ancestor clones
#'    - Uses previous alternative as new original nucleotide
#'    - Samples new alternative based on transition matrix
#'
#' @seealso
#' \code{\link{get_segment_info}}, \code{\link{get_clone_ancestors}}
#'
#' @importFrom igraph V
#' @export
sim_clonal_mutation_nt <- function(genome_sequence, seg_list, mutation_info, nt_transition_matrix, tree){
  # initialize vectors
  seg_id <- character(nrow(mutation_info))
  ref_pos <- numeric(nrow(mutation_info))
  original_nts <- character(nrow(mutation_info))
  alternative_nts <- character(nrow(mutation_info))

  # Construct a unique key for each mutation's location (haplotype_chromosome_position)
  mutation_keys <- paste(mutation_info$haplotype, mutation_info$chrom, mutation_info$pos, sep = "_")

  # Separate unique and recurrent mutations
  regular_mutations <- !duplicated(mutation_keys) # mutations occured at a position for the first time
  recurrent_mutations <- !regular_mutations

  # Handle non-recurrent mutations
  for (i in which(regular_mutations)) {

    # Get the mutation
    mutation <- mutation_info[i,]
    clone <- mutation$clone
    haplotype <- mutation$haplotype

    genome_seg <- seg_list[[clone]][[haplotype]]

    all_loss_indices <- which(is.na(genome_seg$start) &
                                is.na(genome_seg$end) &
                                genome_seg$CN_change == -1)
    loss_segments <- genome_seg[all_loss_indices,]
    in_loss_segment <- any(
      loss_segments$chrom == mutation$chrom &
        loss_segments$haplotype == mutation$haplotype &
        loss_segments$ori_start <= mutation$pos &
        loss_segments$ori_end >= mutation$pos
    )

    segment <- get_segment_info(seg_list, mutation, in_loss_segment)

    if (is.null(segment)) {
      warning("No segment is found for mutation ", i)
      next
    }
    if(in_loss_segment){
      seg_id[i] <- segment$seg_id
      ref_pos[i] <- mutation$pos - segment$ori_start + segment$ref_start
      original_nts[i] <- "N"
      alternative_nts[i] <- "N"
      next
    }

    seg_id[i] <- segment$seg_id

    ref_pos[i] <- mutation$pos - segment$start + segment$ref_start

    # Get the reference nucleotide from the genome_sequence using the segment_id and position
    original_nts[i] <- as.character(genome_sequence[[mutation$haplotype]][[mutation$chrom]][ref_pos[i]])

    # Sample an alternative nucleotide based on the transition matrix
    alternative_nts[i] <- sample(colnames(nt_transition_matrix), 1, prob = nt_transition_matrix[original_nts[i], ])
  }

  # TODO (Critical): Thoroughly test and validate recurrent mutation handling
  # Issues to check:
  # - Ancestor tracking logic
  # - Nucleotide change chain in multiple ancestors
  # - Edge cases with segment losses
  # - Add specific test cases
  # Handle recurrent mutations
  for (i in which(recurrent_mutations)) {
    # Get the mutation
    mutation <- mutation_info[i,]
    ancestor_clones <- get_clone_ancestors(tree, mutation$clone)
    # Find the previous mutation at this location
    mutation_indices <- which(mutation_info$haplotype == mutation$haplotype &
                                mutation_info$chrom == mutation$chrom &
                                mutation_info$pos == mutation$pos &
                                mutation_info$clone %in% ancestor_clones)

    clone <- mutation$clone
    haplotype <- mutation$haplotype
    genome_seg <- seg_list[[clone]][[haplotype]] # all segments in the haplotype genome the mutation reside

    all_loss_indices <- which(is.na(genome_seg$start) &
                                is.na(genome_seg$end) &
                                genome_seg$CN_change == -1)

    loss_segments <- genome_seg[all_loss_indices,]

    in_loss_segment <- any(
      loss_segments$chrom == mutation$chrom &
        loss_segments$haplotype == mutation$haplotype &
        loss_segments$ori_start <= mutation$pos &
        loss_segments$ori_end >= mutation$pos
    )

    segment <- get_segment_info(seg_list, mutation, in_loss_segment)

    if (is.null(segment)) {
      warning("No segment is found for mutation ", i)
      next
    }

    if(in_loss_segment){
      seg_id[i] <- segment$seg_id
      ref_pos[i] <- mutation$pos - segment$ori_start + segment$ref_start
      original_nts[i] <- "N"
      alternative_nts[i] <- "N"
      next
    }

    mutation_indices <- mutation_indices[mutation_indices < i]
    previous_mutation <- mutation_info[max(mutation_indices),]
    original_nts[i] <- alternative_nts[max(mutation_indices)]
    alternative_nts[i] <- sample(colnames(nt_transition_matrix), 1, prob = nt_transition_matrix[original_nts[i], ])
    seg_id[i] <- segment$seg_id
    ref_pos[i] <- mutation$pos - segment$start + segment$ref_start
  }

  mutation_info$seg_id <- seg_id
  mutation_info$ref_pos <- ref_pos
  mutation_info$original_nt <- original_nts
  mutation_info$alternative_nt <- alternative_nts
  return(mutation_info)
}


#' Simulate Random Mutation Positions in a Clonal Tree
#'
#' @description
#' Simulates random mutation positions along chromosomes for each edge in a phylogenetic tree.
#' Mutations are distributed across chromosomes proportionally to their lengths, and can occur
#' on either maternal or paternal haplotypes.
#'
#' @param tree An igraph object representing the phylogenetic tree structure.
#' @param chr_lengths A nested list containing chromosome lengths for each clone and haplotype:
#'   \itemize{
#'     \item First level: clone names
#'     \item Second level: haplotypes (maternal/paternal)
#'     \item Third level: named numeric vector of chromosome lengths
#'   }
#' @param mutation_number A named numeric vector specifying the number of mutations to simulate
#'   for each edge in the tree. Names should be in format "parent_child".
#'
#' @return A data frame containing simulated mutation information:
#'   \itemize{
#'     \item clone - Name of the clone where mutation occurs
#'     \item edge_name - Tree edge identifier (parent_child)
#'     \item haplotype - Maternal or paternal haplotype
#'     \item chrom - Chromosome where mutation occurs
#'     \item pos - Position of mutation on the chromosome
#'   }
#'
#' @details
#' The function processes the tree in depth-first search order, starting from the first tumor
#' clone (excluding root). For each edge, it simulates the specified number of mutations by:
#' 1. Randomly selecting a haplotype (maternal/paternal)
#' 2. Selecting a chromosome with probability proportional to its length
#' 3. Randomly selecting a position within the chosen chromosome
#'
#' @importFrom igraph V degree dfs get.adjlist
#' @export
sim_clonal_mutation_pos <- function(tree, chr_lengths, mutation_number) {
  mutation_info <- data.frame(clone = character(),
                              edge_name = character(),
                              haplotype = character(),
                              chrom = character(),
                              pos = integer(),
                              stringsAsFactors = FALSE)

  root_name <- names(V(tree))[degree(tree, mode = "in") == 0]
  # Obtain nodes of the tree in Depth-First Search (DFS) order
  dfs_order <- dfs(tree, root = root_name, mode = "out")$order

  # Start from the second node (the first tumor clone)
  for(i in 2:length(dfs_order)) {
    clone <- V(tree)[dfs_order[i]]$name
    parent_node <- V(tree)[get.adjlist(tree, mode = "in")[[clone]]]$name
    edge_name <- paste(parent_node, clone, sep = "_")
    n <- mutation_number[edge_name]

    for (j in 1:n) {
      # Select a position for the mutation
      clone_chr_lengths <- chr_lengths[[clone]]
      mut_hap <- sample(names(clone_chr_lengths), 1)
      mut_chr <- sample(names(clone_chr_lengths[[mut_hap]]), 1, prob = clone_chr_lengths[[mut_hap]])
      mut_pos <- sample(1:clone_chr_lengths[[mut_hap]][mut_chr], 1)

      # Record the mutation
      mutation_info <- rbind(mutation_info, data.frame(clone = clone,
                                                       edge_name = edge_name,
                                                       haplotype = mut_hap,
                                                       chrom = mut_chr,
                                                       pos = mut_pos))
    }
  }

  return(mutation_info)
}

#' Retrieve Genomic Segment Information for a Mutation
#'
#' @description
#' Identifies and returns the genomic segment that contains a specific mutation,
#' handling both regular segments and segments that have been lost through deletions.
#'
#' @param seg_list List structure containing segment information for different
#'        clones and haplotypes
#' @param mutation Data frame row containing information about a single mutation,
#'        with columns: clone, haplotype, chrom, pos
#' @param in_loss_segment Logical indicating whether the mutation is in a segment
#'        that has been lost through deletion
#'
#' @return A data frame row containing the segment information for the mutation,
#'         or NULL if no matching segment is found
#'
#' @details
#' This function locates the appropriate genomic segment for a mutation by:
#'
#' 1. Extracting relevant information (haplotype, chromosome, clone, position) from the mutation
#' 2. Using different search criteria based on whether the segment is lost:
#'    - For lost segments: Matches chromosome and original coordinates (ori_start/ori_end)
#'      where CN_change is -1 and current coordinates (start/end) are NA
#'    - For normal segments: Matches chromosome and current coordinates (start/end)
#' 3. Returns NULL if no matching segment is found
#'
#' This function is primarily used by simulate_single_nt_change() to map mutations
#' to their genomic context before simulating nucleotide changes.
#'
#' @seealso
#' \code{\link{simulate_single_nt_change}}
#'
#' @export
get_segment_info <- function(seg_list, mutation, in_loss_segment) {
  # Extract relevant information from the mutation
  haplotype <- mutation$haplotype
  chrom <- mutation$chrom
  clone <- mutation$clone
  pos <- mutation$pos

  if(in_loss_segment){
    segment_indices <- which(seg_list[[clone]][[haplotype]]$chrom == chrom &
                               seg_list[[clone]][[haplotype]]$ori_start <= pos &
                               seg_list[[clone]][[haplotype]]$ori_end >= pos &
                               seg_list[[clone]][[haplotype]]$CN_change == -1 &
                               is.na(seg_list[[clone]][[haplotype]]$start) &
                               is.na(seg_list[[clone]][[haplotype]]$end))
  }else{
    # Find the segment in seg_list that contains the mutation position
    segment_indices <- which(seg_list[[clone]][[haplotype]]$chrom == chrom &
                               seg_list[[clone]][[haplotype]]$start <= pos &
                               seg_list[[clone]][[haplotype]]$end >= pos)


  }

  if (length(segment_indices) == 0) {
    return(NULL) # Return NULL if no segment is found
  }
  # Return the corresponding segment information
  segment <- seg_list[[clone]][[haplotype]][segment_indices, , drop = FALSE]
  return(segment)
}

