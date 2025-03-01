# pattern_extraction_v3.R

#' Subset and Transform Haplotype Copy Number Profiles for Regions of Interest
#'
#' @description
#' Extracts copy number data for specific regions of interest from a haplotype-specific
#' copy number profile and returns a transformed matrix with the same dimensions as
#' the original data, where only the regions of interest contain values.
#'
#' @param haplotype_cn A list containing haplotype-specific copy number matrices
#' @param haplotype Character string specifying which haplotype to analyze
#'        (must be a name in the haplotype_cn list)
#' @param chr Character string specifying the chromosome of interest
#'        (e.g., "chr1", "chrX")
#' @param interested_region_indices Numeric vector of indices specifying which
#'        regions within the chromosome to analyze
#'
#' @return A numeric matrix with the same dimensions as the original copy number
#'         matrix, where:
#'         \itemize{
#'           \item Only specified regions of interest contain copy number values
#'           \item All other positions contain NA
#'           \item Column names are preserved from the original matrix
#'         }
#'
#' @details
#' The function performs these steps:
#' \itemize{
#'   \item Identifies regions in the specified chromosome
#'   \item Extracts copy number data for regions of interest
#'   \item Creates a new matrix with same dimensions as original data
#'   \item Fills in only the specified regions, leaving others as NA
#' }
#'
#' @examples
#' # Create sample haplotype copy number data
#' sample_cn <- list(
#'   hap1 = matrix(1:20, nrow = 2, ncol = 10,
#'                 dimnames = list(NULL,
#'                 paste0("chr1_", 1:10))),
#'   hap2 = matrix(2:21, nrow = 2, ncol = 10,
#'                 dimnames = list(NULL,
#'                 paste0("chr1_", 1:10)))
#' )
#'
#' # Extract regions of interest
#' subset_data <- subset_interested_profiles(
#'   haplotype_cn = sample_cn,
#'   haplotype = "hap1",
#'   chr = "chr1",
#'   interested_region_indices = 2:4
#' )
#'
#' @seealso
#' Other functions for manipulating copy number profiles in the package
#'
#' @export
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

#' Create Cytoband Annotations for Region Clusters
#'
#' @description
#' Maps clustered genomic regions to their corresponding cytoband annotations for both
#' maternal and paternal haplotypes. For each cluster of regions, the function collects
#' and combines the cytoband annotations of all regions within that cluster.
#'
#' @param region_clusters A nested list structure containing clustered genomic regions.
#'        The list is organized by haplotype ("maternal"/"paternal"), then by chromosome,
#'        then by cluster number, containing region IDs as the innermost elements.
#' @param cytoband_mapping A list containing cytoband annotations for each genomic region.
#'        Each element should be named with a region ID corresponding to those in region_clusters.
#'
#' @return A nested list with the same structure as region_clusters, where each cluster
#'         contains a data frame of unique cytoband annotations. The structure is:
#'         list[haplotype][chromosome][cluster] -> data frame of cytoband annotations
#'
#' @examples
#' \dontrun{
#' # Assuming region_clusters and cytoband_mapping are already defined:
#' cluster_annotations <- create_cluster_cytoband_anno(
#'   region_clusters = region_clusters,
#'   cytoband_mapping = cytoband_mapping
#' )
#'
#' # Access annotations for a specific cluster:
#' maternal_chr1_cluster1 <- cluster_annotations[["maternal"]][["chr1"]][[1]]
#' }
#'
#' @details
#' The function processes clusters by:
#' 1. Iterating through haplotypes (maternal/paternal)
#' 2. For each chromosome in each haplotype
#' 3. For each cluster in each chromosome
#' 4. Retrieving and combining cytoband annotations for all regions in the cluster
#' 5. Removing duplicate annotations using unique()
#'
#' @seealso
#' Related functions for genomic region analysis and cytoband mapping
#'
#' @export
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


#' Generate Chromosome Boundary Data
#'
#' This function processes chromosome arm and sub-segment data to generate boundary information
#' for specified chromosomes. It combines segments from both input tables, identifies unique
#' boundaries, and creates windows for each chromosome. Additionally, it constructs a genome-wide
#' vector of window identifiers.
#'
#' @param chr_arm_table A data frame containing chromosome arm information:
#'        \itemize{
#'          \item chrom: Chromosome name
#'          \item start: Start position
#'          \item end: End position
#'          \item region_name: Identifier for the chromosome arm
#'        }
#' @param sub_seg_list A list of data frames (by haplotype) containing sub-segment information:
#'        \itemize{
#'          \item chrom: Chromosome name
#'          \item ref_start: Reference start position
#'          \item ref_end: Reference end position
#'          \item region_name: Segment identifier
#'        }
#' @param chr_names A character vector of chromosome names to process
#'
#' @return A list containing three elements:
#'         \itemize{
#'           \item segment_info: A list of data frames, one for each chromosome, containing combined
#'       segments from `chr_arm_table` and `sub_seg_list`.
#'           \item windows_info: A list of data frames, one for each chromosome, containing windows
#'       derived from the unique boundaries of the combined segments.
#'           \item genome_window_vector: Vector of formatted window strings ("chr_start_end")
#'         }
#'
#' @details
#' The function performs these steps for each chromosome:
#' 1. Extracts and combines segments from chromosome arms and sub-segments
#' 2. Generates windows based on unique boundary positions
#' 3. Creates formatted window identifiers
#' 4. Organizes results by chromosome
#'
#' Windows are generated between each unique boundary point found in the
#' combined segment data, ensuring all relevant genomic intervals are captured.
#'
#' @examples
#' \dontrun{
#' chr_arms <- data.frame(
#'   chrom = c("chr1", "chr1"),
#'   start = c(1, 1000000),
#'   end = c(1000000, 2000000),
#'   region_name = c("chr1p", "chr1q")
#' )
#' sub_segs <- list(
#'   maternal = data.frame(...),
#'   paternal = data.frame(...)
#' )
#' chr_names <- c("chr1")
#' boundaries <- generate_chr_boundary_data(chr_arms, sub_segs, chr_names)
#' }
#'
#' @importFrom dplyr filter bind_rows arrange
#' @export
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


#' Generate Blueprint Copy Number Profiles
#'
#' This function generates copy number (CN) profiles for clones based on their events and a genome-wide
#' window vector. It initializes baseline CN profiles for paternal and maternal haplotypes and updates
#' them according to the events (e.g., copy number changes, whole genome duplications) associated with
#' each clone. The function handles both whole chromosome events and specific region events.
#'
#' @param clone_events A named list of data frames, where each data frame contains events for a clone:
#'        \itemize{
#'          \item event_type: Type of event (e.g., "WGD" for whole genome duplication)
#'          \item region_name: Region identifier or chromosome name
#'          \item CN_change: Copy number change value
#'          \item haplotype: "maternal" or "paternal"
#'        }
#' @param genome_window_vector A character vector of window identifiers in format "chr_start_end"
#' @param segment_info A data frame containing segment information:
#'        \itemize{
#'          \item region_name: Region identifier
#'          \item chrom: Chromosome name
#'          \item start: Start position
#'          \item end: End position
#'        }
#'
#' @return A list with two elements ("maternal" and "paternal"), each containing
#'         a data frame where:
#'         \itemize{
#'           \item Rows represent clones
#'           \item Columns represent genomic windows
#'           \item Values represent copy numbers
#'         }
#'
#' @details
#' The function processes events in these steps:
#' 1. Initializes baseline copy numbers (1 for all regions)
#' 2. For each clone and its events:
#'    - Handles WGD by doubling all copy numbers
#'    - Processes chromosome-wide events if region matches "chr[0-9XY]+"
#'    - Updates specific regions based on segment information
#'    - Applies copy number changes progressively
#' 3. Maintains separate profiles for maternal and paternal haplotypes
#'
#' @note
#' - Copy numbers cannot go below 0
#' - WGD events affect all regions simultaneously
#' - Chromosome-wide events are detected by regex pattern
#'
#' @examples
#' \dontrun{
#' clone_events <- list(
#'   clone1 = data.frame(
#'     event_type = c("CN", "WGD"),
#'     region_name = c("chr1_p", "genome"),
#'     CN_change = c(1, 0),
#'     haplotype = c("maternal", "maternal")
#'   ),
#'   clone2 = data.frame(...)
#' )
#' genome_windows <- c("chr1_0_1000000", "chr1_1000000_2000000")
#' segments <- data.frame(
#'   region_name = "chr1_p",
#'   chrom = "chr1",
#'   start = 0,
#'   end = 1000000
#' )
#' profiles <- generate_blueprint_cn_profiles(clone_events, genome_windows, segments)
#' }
#'
#' @export
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


#' Extend Blueprint Copy Number Data for Equal-Sized Bins
#'
#' This function extends a copy number (CN) data frame (`seg_dat`) to ensure that each genomic window
#' is divided into equal-sized bins of a specified size (`bin_unit`). It replicates columns in the
#' input data frame based on the size of each window and the desired bin size, creating a new data
#' frame with expanded columns.
#'
#' @param seg_dat A data frame of copy number segments where:
#'        \itemize{
#'          \item Rows represent clones/samples
#'          \item Columns represent genomic windows
#'          \item Values represent copy numbers
#'        }
#' @param window_sizes A numeric vector containing the size of each window in base pairs,
#'        must match the number of columns in seg_dat
#' @param bin_unit The desired size of each bin in base pairs
#'
#' @return A data frame with:
#'         \itemize{
#'          \item Same number of rows as input
#'          \item Expanded columns based on window sizes and bin unit
#'          \item Column names as "original_window_name_binIndex"
#'          \item Original row names preserved
#'         }
#'
#' @details
#' The function:
#' 1. Calculates how many bins each window needs based on its size
#' 2. Replicates copy number values to fill the required bins
#' 3. Maintains data continuity while creating equal-sized segments
#'
#' @examples
#' \dontrun{
#' seg_data <- data.frame(
#'   chr1_1000_3000 = c(2, 1),
#'   chr1_3000_4000 = c(1, 1),
#'   row.names = c("clone1", "clone2")
#' )
#' window_sizes <- c(2000, 1000)
#' bin_unit <- 1000
#'
#' expanded <- extend_blueprint_cn_for_equal_bin(seg_data, window_sizes, bin_unit)
#' # Returns expanded data frame with bins of 1000bp each
#' }
#'
#' @note
#' Window sizes that aren't exact multiples of bin_unit will be rounded up
#' to ensure complete coverage
#'
#' @export
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


#' Calculate Sizes of Genomic Windows
#'
#' @description
#' Calculates the size (length in base pairs) of each genomic window from a vector
#' of window identifiers in "chr_start_end" format.
#'
#' @param genome_window_vector A character vector of window identifiers in the format
#'        "chr_start_end" (e.g., "chr1_1000_2000")
#'
#' @return A numeric vector containing the size of each window in base pairs,
#'         calculated as (end position - start position)
#'
#' @details
#' The function:
#' 1. Splits each window identifier into its components
#' 2. Extracts start and end positions
#' 3. Calculates window size as the difference between end and start
#'
#' @examples
#' \dontrun{
#' windows <- c("chr1_1000_2000", "chr1_2000_5000", "chr2_1000_3000")
#' sizes <- calculate_window_sizes(windows)
#' # Returns: c(1000, 3000, 2000)
#' }
#'
#' @note
#' Assumes window identifiers are properly formatted with underscore separators
#' and numeric start/end positions
#'
#' @export
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
