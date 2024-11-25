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
        print(actual_nts)
        # Add actual nts to the selected mutations
        selected_mutations$actual_nts <- actual_nts

        # Compare expected vs. actual nucleotides
        correct <- expected_nts == actual_nts

        print(expected_nts)
        print(actual_nts)
        print(correct)
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
