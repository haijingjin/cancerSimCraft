# pattern_extraction_v3.R

subset_interested_profiles <- function(haplotype_cn, haplotype, chr, interested_region_indices) {
  # Extract the chromosome of interest
  all_regions <- colnames(haplotype_cn[[haplotype]])
  chr_regions <- all_regions[grep(chr, all_regions)]
  interested_region <- chr_regions[interested_region_indices]
  subset_data <- haplotype_cn[[haplotype]][, chr_regions]
  all_cn_dat <- haplotype_cn[[haplotype]]

  # Create an empty matrix with appropriate dimensions
  transformed_cn <- matrix(data = NA, nrow = nrow(all_cn_dat), ncol = ncol(all_cn_dat))
  colnames(transformed_cn) <- colnames(all_cn_dat)

  # Fill in the matrix for the regions of interest
  transformed_cn[, interested_region] <- subset_data[, interested_region_indices]

  return(transformed_cn)
}

create_cluster_cytoband_anno <- function(region_clusters, cytoband_mapping) {

  cluster_anno_list <- list()
  # Loop through each chromosome in region_clusters
  for(haplotype in c("maternal", "paternal")){

    for (chr in names(region_clusters[[haplotype]])) {
      cluster_anno_list[[haplotype]][[chr]] <- list()
      # Loop through each cluster in this chromosome
      for (i in 1:length(region_clusters[[haplotype]][[chr]])){
        # Retrieve the regions in this cluster
        regions <- names(region_clusters[[haplotype]][[chr]][[i]])
        # Map each region to its cytoband annotations
        cytoband_annotations <- lapply(regions, function(region) return(cytoband_mapping[[region]]))
        # Store the annotations in cluster_anno_list
        cluster_anno_list[[haplotype]][[chr]][[i]] <- unique(do.call(rbind, cytoband_annotations))
      }
    }

  }
  return(cluster_anno_list)
}


# Define the main function
generate_chr_boundary_data <- function(chr_arm_table, sub_seg_list, chr_names) {
  segment_list <- list()
  windows_list <- list()
  genome_window_vector <- c()  # Initialize an empty vector for the entire genome

  # Supporting function for generating windows
  generate_windows <- function(segments) {
    unique_boundaries <- sort(unique(c(segments$start, segments$end)))
    windows <- data.frame(
      start = head(unique_boundaries, -1),
      end = tail(unique_boundaries, -1)
    )
    return(windows)
  }

  for (chr_name in chr_names) {
    # Extract relevant segments from both tables
    chr_arm_segments <- dplyr::filter(chr_arm_table, chrom == chr_name)[, c('chrom', 'start', 'end', 'region_name')]
    sub_seg_segments <- do.call(rbind, lapply(sub_seg_list, dplyr::filter, chrom == chr_name))[, c('chrom', 'ref_start', 'ref_end', 'region_name')]
    names(sub_seg_segments) <- c('chrom', 'start', 'end', 'region_name')

    # Combine segments and sort
    combined_segments <- dplyr::bind_rows(chr_arm_segments, sub_seg_segments) %>%
      dplyr::arrange(start)

    # Add to segment list
    segment_list[[chr_name]] <- combined_segments

    # Generate new windows
    new_windows <- generate_windows(combined_segments)
    # Add to windows list
    windows_list[[chr_name]] <- new_windows

    # Generate window vectors and append to the genome-wide vector
    window_vectors <- apply(new_windows, 1, function(row) {
      sprintf("%s_%d_%d", chr_name, row["start"], row["end"])
    })
    genome_window_vector <- c(genome_window_vector, window_vectors)
  }

  return(list(segment_info = segment_list, windows_info = windows_list, genome_window_vector = genome_window_vector))
}


generate_blueprint_cn_profiles <- function(clone_events, genome_window_vector, segment_info) {
  # Get a list of all clone names
  clone_names <- names(clone_events)
  # Initialize baseline copy number profiles for paternal and maternal
  initial_values <- rep(1, length(clone_names))
  ground_truth_cn <- list(
    paternal = data.frame(setNames(replicate(length(genome_window_vector), initial_values, simplify = FALSE), genome_window_vector), row.names = clone_names),
    maternal = data.frame(setNames(replicate(length(genome_window_vector), initial_values, simplify = FALSE), genome_window_vector), row.names = clone_names)
  )
  # Function to parse start and end from window vector
  parse_window <- function(window_str) {
    parts <- unlist(strsplit(window_str, "_"))
    start <- as.numeric(parts[2])
    end <- as.numeric(parts[3])
    return(c(start, end))
  }
  # Function to update copy number based on an event
  update_copy_number <- function(event, cn_profile, genome_window_vector, segment_info) {

    region = event$region_name
    cn_change = event$CN_change
    haplotype = event$haplotype

    # Handling Whole Genome Duplication (WGD)
    if (event$event_type == "WGD") {
      cn_profile <- cn_profile * 2
    } else {
      # Check if the event is for a whole chromosome
      if (grepl("^chr[0-9XY]+$", region)) {
        # Update all regions in the chromosome
        affected_regions <- grep(paste0(region, "_"), genome_window_vector)
      } else {
        # Determine affected windows
        affected_seg <- segment_info[segment_info$region_name == region, ]
        affected_regions <- which(genome_window_vector == paste(affected_seg$chrom, affected_seg$start, affected_seg$end, sep = "_"))

      }
      cn_profile[affected_regions] <- cn_profile[affected_regions] + cn_change
    }

    return(cn_profile)
  }

  # Process each event for each clone
  for(i in 1:length(clone_events)) {
    # Check if the clone_events for the clone are empty or NA
    if (all(is.na(clone_events[[i]]))) {
      # If no events, keep the baseline copy number and continue to the next clone
      next
    }else{
      clone <- clone_names[i]
      # If events exist, process them
      for (j in 1:nrow(clone_events[[i]])) {
        event = clone_events[[i]][j, ]
        haplotype = event$haplotype
        ground_truth_cn[[haplotype]][clone,] <- update_copy_number(event, ground_truth_cn[[haplotype]][clone,], genome_window_vector, segment_info)
      }
    }

  }

  return(ground_truth_cn)
}

extend_blueprint_cn_for_equal_bin <- function(seg_dat, window_sizes, bin_unit) {
  replication_factors <- ceiling(window_sizes / bin_unit)
  expanded_df <- data.frame(matrix(nrow = nrow(seg_dat), ncol = 0))

  current_col <- 1
  new_colnames <- character(0)
  for (i in 1:ncol(seg_dat)) {
    cols_to_fill <- current_col:(current_col + replication_factors[i] - 1)
    replicated_cols <- matrix(unlist(rep(seg_dat[, i], times = replication_factors[i])), byrow = FALSE, nrow = nrow(seg_dat))
    expanded_df <- cbind(expanded_df, replicated_cols)
    current_col <- current_col + replication_factors[i]

    # Create new names for the replicated columns
    new_colnames <- c(new_colnames, paste0(colnames(seg_dat)[i], "_", 1:replication_factors[i]))
  }

  colnames(expanded_df) <- new_colnames
  rownames(expanded_df) <- rownames(seg_dat)

  return(expanded_df)
}


calculate_window_sizes <- function(genome_window_vector) {
  # Function to extract the start and end positions from a window vector string
  get_positions <- function(window_str) {
    parts <- unlist(strsplit(window_str, "_"))
    start <- as.numeric(parts[length(parts) - 1])
    end <- as.numeric(parts[length(parts)])
    return(c(start, end))
  }

  # Calculate sizes for each window
  window_sizes <- sapply(genome_window_vector, function(w) {
    positions <- get_positions(w)
    size <- positions[2] - positions[1]
    return(size)
  })

  return(window_sizes)
}
