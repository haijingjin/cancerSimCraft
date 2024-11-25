# snv_relevant_functions_v7.R

calc_clone_snv_num <- function(time, genome_length, mutation_rate) {
  # A simple model where mutation number is proportional to time, genome length, and a base mutation rate
  num_mutations <- time * genome_length * mutation_rate
  return(round(num_mutations, digits = 0))
}

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


identify_recurrent_mutations <- function(cell_mutations) {
  # Create a key for each mutation based on haplotype, chromosome, and position
  mutation_keys <- paste(cell_mutations$haplotype, cell_mutations$chrom, cell_mutations$pos, sep = "_")

  # Identify duplicates in these keys, which indicate recurrent mutations
  duplicate_indices <- which(duplicated(mutation_keys) | duplicated(mutation_keys, fromLast = TRUE))

  # Extract the recurrent mutations based on these indices
  recurrent_mutations <- cell_mutations[duplicate_indices, ]

  return(recurrent_mutations)
}

group_recurrent_mutations <- function(recurrent_mutations) {
  # Create a key for grouping mutations
  group_keys <- paste(recurrent_mutations$haplotype, recurrent_mutations$chrom, recurrent_mutations$pos, sep = "_")

  # Split the data frame into a list of data frames, each containing mutations at the same position
  grouped_recurrent_mutations <- split(recurrent_mutations, group_keys)

  return(grouped_recurrent_mutations)
}

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


find_identical_row_index <- function(row, dataframe) {
  for (i in 1:nrow(dataframe)) {
    comparison_result <- all.equal(row, dataframe[i, , drop = FALSE], check.attributes = FALSE)
    if (isTRUE(comparison_result)) {
      return(i)
    }
  }
  return(NA)  # Return NA if no identical row is found
}



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

