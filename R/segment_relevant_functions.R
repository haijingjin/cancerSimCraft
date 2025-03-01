# segment_relevant_functions_v4.R

initiate_sub_seg_list <- function(seg_table) {

  sub_seg_list <- list()

  for (haplotype in c("maternal", "paternal")) {
    one_seg_table <- data.frame(
      haplotype = haplotype,
      chrom = seg_table$chrom,
      ref_start = seg_table$ref_start,
      ref_end = seg_table$ref_end,
      ori_start = seg_table$ori_start,
      ori_end = seg_table$ori_end,
      start = seg_table$start,
      end = seg_table$end,
      region_name = seg_table$region_name,
      copy_index = 0,
      seg_id = 0,
      CN_change = 0,
      seg_source_edge = 'root',
      seg_source_event = "base"
    )

    one_seg_table$seg_id <- paste(one_seg_table$region_name, one_seg_table$copy_index, sep = "_")
    sub_seg_list[[haplotype]] <- one_seg_table
  }

  return(sub_seg_list)
}

#' Identify Clusters of Highly Correlated Genomic Regions
#'
#' @description
#' Identifies clusters of adjacent genomic regions that show high correlation with each other.
#' Uses a sliding window approach to find contiguous regions where the average pairwise
#' correlation exceeds a specified threshold.
#'
#' @param cor_matrix A numeric correlation matrix where rows and columns represent
#'        genomic regions, with column names identifying the regions
#' @param min_cluster_size Integer specifying the minimum number of regions required
#'        to form a cluster (must be >= 2)
#' @param threshold Numeric value between -1 and 1 specifying the minimum average
#'        correlation required for regions to be considered part of the same cluster
#'
#' @return A list where each element represents a cluster. Each cluster contains:
#'         \itemize{
#'           \item Named integer vector with positions of regions in the cluster
#'           \item Names of the vector correspond to the region identifiers from
#'                 the input correlation matrix
#'         }
#'
#' @details
#' The function implements a dynamic window-based clustering algorithm:
#' \itemize{
#'   \item Starts with a window of minimum cluster size
#'   \item Calculates average correlation (excluding diagonal) within the window
#'   \item If average correlation exceeds threshold:
#'     \itemize{
#'       \item Extends window if not at matrix end
#'       \item Finalizes cluster if at matrix end
#'     }
#'   \item If correlation drops below threshold:
#'     \itemize{
#'       \item Finalizes cluster if window size >= minimum
#'       \item Moves to next position if window size < minimum
#'     }
#' }
#'
#' @examples
#' # Create sample correlation matrix
#' set.seed(123)
#' n_regions <- 10
#' cor_mat <- matrix(runif(n_regions^2, -1, 1), nrow = n_regions)
#' cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]
#' diag(cor_mat) <- 1
#' colnames(cor_mat) <- paste0("region_", 1:n_regions)
#'
#' # Find clusters
#' clusters <- identify_region_clusters(
#'   cor_matrix = cor_mat,
#'   min_cluster_size = 3,
#'   threshold = 0.7
#' )
#'
#' @seealso
#' Other clustering functions in the package for genomic analysis
#'
#' @export
identify_region_clusters <- function(cor_matrix, min_cluster_size, threshold) {
  cluster_list <- list()
  window_start <- 1
  window_end <- window_start + min_cluster_size - 1
  region_names <- colnames(cor_matrix)
  new_start = TRUE
  while( new_start == TRUE && window_start <= ncol(cor_matrix) - min_cluster_size + 1){

    extend_window = TRUE
    while( extend_window == TRUE){
      submatrix <- cor_matrix[window_start:window_end, window_start:window_end]
      diag(submatrix) <- diag(submatrix) - 1
      avg_cor <- mean(submatrix[lower.tri(submatrix)])

      if (avg_cor >= threshold) {
        # If we're at the end of the matrix, finalize this as a cluster
        if (window_end == ncol(cor_matrix)) {

          cluster_list[[length(cluster_list) + 1]] <- c(window_start:window_end)
          names(cluster_list[[length(cluster_list)]]) <- region_names[window_start:window_end]
          extend_window <- FALSE
          new_start <- FALSE
        } else {
          # Extend the window to include the next region
          window_end <- window_end + 1
        }
      } else {
        # If the average correlation dropped below the threshold, finalize the cluster
        extend_window <- FALSE

        if ( (window_end-1) >= (window_start + min_cluster_size-1)) {

          cluster_list[[length(cluster_list) + 1]] <- c(window_start:(window_end - 1))
          names(cluster_list[[length(cluster_list)]]) <- region_names[window_start:(window_end - 1)]
          window_start = window_end
        }else{
          window_start = window_start + 1
        }
        window_end = window_start + min_cluster_size - 1
      }
    }
  }

  return(cluster_list)
}

#' Create Genomic Segments from Cytoband Annotations
#'
#' @description
#' Converts clustered cytoband annotations into genomic segments by identifying
#' the spanning regions for each cluster. Creates segments for both maternal and
#' paternal haplotypes with additional metadata for downstream analysis.
#'
#' @param cluster_anno A nested list containing cytoband annotations organized by
#'        haplotype (maternal/paternal), chromosome, and cluster. Each cluster contains
#'        a data frame of cytoband information including chromStart and chromEnd positions.
#'
#' @return A list with two elements (maternal and paternal), each containing a data frame
#'         of genomic segments with the following columns:
#'         \itemize{
#'           \item chrom: Chromosome name
#'           \item start: Starting position of the segment
#'           \item end: Ending position of the segment
#'           \item region_name: Unique identifier (format: "sub_haplotype_chrom_clusterIndex")
#'           \item haplotype: Maternal or paternal
#'           \item ref_start: Reference start position
#'           \item ref_end: Reference end position
#'           \item ori_start: Original start position
#'           \item ori_end: Original end position
#'           \item copy_index: Copy number index (initialized to 0)
#'           \item seg_id: Segment identifier (initialized to 0)
#'           \item CN_change: Copy number change (initialized to 0)
#'           \item seg_source_edge: Source edge info (initialized to "root")
#'           \item seg_source_event: Source event info (initialized to "base")
#'         }
#'
#' @examples
#' \dontrun{
#' # Assuming cluster_anno is already defined:
#' segments <- sub_seg_from_cytoband_anno(cluster_anno)
#'
#' # Access maternal segments:
#' maternal_segs <- segments[["maternal"]]
#' }
#'
#' @details
#' For each cluster in the input annotations, the function:
#' 1. Identifies the minimum start and maximum end positions of all cytobands
#' 2. Creates a unique region name using haplotype, chromosome, and cluster index
#' 3. Initializes metadata fields for downstream analysis
#' 4. Combines all segments into haplotype-specific data frames
#'
#' The resulting segments serve as a baseline for further genomic analyses,
#' with fields prepared for tracking copy number changes and segment evolution.
#'
#' @seealso
#' Related functions for genomic segment analysis and manipulation
#'
#' @export
sub_seg_from_cytoband_anno <- function(cluster_anno) {

  sub_seg_list <- list()

  # Iterate through each chromosome
  for(haplotype in c("maternal", "paternal")){

    for(chr in names(cluster_anno[[haplotype]])){

      for(i in 1:length(cluster_anno[[haplotype]][[chr]])){
        one_cluster_anno <- cluster_anno[[haplotype]][[chr]][[i]]

        # Construct the region_name
        region_name <- paste("sub", haplotype, chr, i, sep = "_")

        # Extract the start and end positions
        start_pos <- min(one_cluster_anno$chromStart)
        end_pos <- max(one_cluster_anno$chromEnd)

        # Create a row for this cluster
        new_row <- data.frame(chrom = chr,
                              start = start_pos,
                              end = end_pos,
                              region_name = region_name,
                              haplotype = haplotype,
                              ref_start = start_pos,
                              ref_end = end_pos,
                              ori_start = start_pos,
                              ori_end = end_pos,
                              copy_index = 0,
                              seg_id = 0,
                              CN_change = 0,
                              seg_source_edge = "root",
                              seg_source_event = "base"
        )

        sub_seg_list[[haplotype]] <- rbind(sub_seg_list[[haplotype]], new_row)
      }
    }
  }
  return(sub_seg_list)
}

cytoband_to_chr_arm  <- function(cytoband_table) {
  # Initialize an empty data frame for the chromosome arm table
  chr_arm_table <- data.frame(chrom = character(),
                              start = integer(),
                              end = integer(),
                              region_name = character(),
                              stringsAsFactors=FALSE)

  # Process each chromosome
  chromosomes <- unique(cytoband_table$chrom)

  for(chrom in chromosomes){
    # subset data for the current chrom
    chrom_dat <- cytoband_table[cytoband_table$chrom == chrom,]
    # find the end of the p arm (start is 0)
    p_arm_end <- chrom_dat[chrom_dat$gieStain == "acen" & grepl(pattern = "^p", chrom_dat$name),]$chromEnd
    q_arm_start <- p_arm_end + 1
    # Find the end of the q arm (maximum chromEnd value of q arm segments)
    q_arm_end <- max(chrom_dat[grepl(pattern = "^q", chrom_dat$name), "chromEnd"])
    # Append the p arm data to the chrArmTable
    chr_arm_table <- rbind(chr_arm_table, data.frame(chromosome = c(chrom, chrom),
                                                     start = c(1, q_arm_start),
                                                     end = c(p_arm_end, q_arm_end),
                                                     region_name = paste0(chrom, c("p", "q"))))
  }
  return(chr_arm_table)
}


#' Create Initial Segment List from Chromosome Arm Table
#'
#' @description
#' Creates initial genomic segments for both maternal and paternal haplotypes from
#' a chromosome arm table. Handles coordinate system conversion (0-based to 1-based)
#' and initializes tracking metadata for each segment.
#'
#' @param chr_arm_table A data frame containing chromosome arm information with columns:
#'        \itemize{
#'          \item chrom: Chromosome name
#'          \item start: Start position
#'          \item end: End position
#'          \item region_name: Identifier for the chromosome arm
#'        }
#'
#' @return A list with two elements ("maternal" and "paternal"), each containing
#'         a data frame of segments with the following columns:
#'         \itemize{
#'           \item All original columns from chr_arm_table
#'           \item haplotype: "maternal" or "paternal"
#'           \item ref_start: Reference start position
#'           \item ref_end: Reference end position
#'           \item ori_start: Original start position
#'           \item ori_end: Original end position
#'           \item copy_index: Set to 1 for initial segments
#'           \item seg_id: Segment identifier (format: `<region_name>_1`)
#'           \item CN_change: Copy number change (initialized to 0)
#'           \item seg_source_edge: Set to "before_root"
#'           \item seg_source_event: Set to "base"
#'         }
#'
#' @details
#' The function performs the following operations:
#' 1. Checks and converts coordinates from 0-based to 1-based if necessary
#' 2. Creates two copies of the input table (maternal and paternal)
#' 3. Adds metadata columns for tracking segment evolution
#' 4. Uses dplyr for efficient data manipulation
#'
#' @note
#' Assumes input coordinates are either 0-based or 1-based.
#' Automatically converts 0-based coordinates to 1-based.
#'
#' @examples
#' \dontrun{
#' chr_arms <- data.frame(
#'   chrom = c("chr1", "chr1"),
#'   start = c(0, 1000000),
#'   end = c(1000000, 2000000),
#'   region_name = c("chr1p", "chr1q")
#' )
#' segments <- create_initial_seg_list(chr_arms)
#' }
#'
#' @importFrom dplyr %>% mutate
#' @export
create_initial_seg_list <- function(chr_arm_table) {

  # If chr_arm_table is 0-based, create 1-based coordinates
  if (any(chr_arm_table$start == 0)) {
    chr_arm_table$start <- chr_arm_table$start + 1
  }


  # Assuming the chr_arm_table has columns: chrom, start, end, region_name
  initial_seg_list <- list()
  for( haplotype in c("maternal", "paternal")){
    initial_seg_list[[haplotype]] = chr_arm_table %>%
      mutate(
        haplotype = haplotype,
        ref_start = start,
        ref_end = end,
        ori_start = start,
        ori_end = end,
        copy_index = 1,
        seg_id = paste0(region_name, "_1"),
        CN_change = 0,
        seg_source_edge = "before_root",
        seg_source_event = "base"
      )
  }

  return(initial_seg_list)
}


#' Initialize Segment Table for Chromosome Regions
#'
#' @description
#' Creates a table of chromosome segments dividing each chromosome into p (short arm),
#' c (centromere), and q (long arm) regions. Each segment is initialized with base
#' state properties including positions, copy numbers, and source information.
#'
#' @param chr_lengths A list containing chromosome lengths for a specific haplotype
#' @param chr_names A vector of chromosome names (e.g., c("chr1", "chr2", "chr3"))
#' @param centromere_length The length of the centromere region
#' @param haplotype Character string specifying the haplotype ("maternal" or "paternal")
#'
#' @return A data frame containing segment information with columns:
#' \itemize{
#'   \item haplotype: Origin of the segment (maternal/paternal)
#'   \item chrom: Chromosome identifier
#'   \item ref_start, ref_end: Reference coordinates
#'   \item ori_start, ori_end: Original coordinates
#'   \item start, end: Current coordinates
#'   \item region_name: Region identifier (e.g., chr1p, chr1c, chr1q)
#'   \item copy_index: Copy number index
#'   \item seg_id: Unique segment identifier
#'   \item CN_change: Copy number change status
#'   \item seg_source_edge: Source edge information
#'   \item seg_source_event: Source event information
#' }
#'
#' @examples
#' chr_lengths <- list(maternal = c(100, 200, 300))
#' chr_names <- c("chr1", "chr2", "chr3")
#' centromere_length <- 10
#' seg_table <- initiate_seg_table(chr_lengths, chr_names, centromere_length, "maternal")
#'
#' @export
initiate_seg_table <- function(chr_lengths, chr_names, centromere_length, haplotype) {
  # Initialize an empty list to store segment data frames
  seg_data_frames <- list()

  # Loop through each chromosome
  for (i in seq_along(chr_names)) {
    # Get the name and length of the chromosome
    chrom_name <- chr_names[i]
    chrom_length <- chr_lengths[[haplotype]][i]

    # Calculate the start and end positions of the p, c, and q regions
    midpoint <- floor(chrom_length / 2)
    c_start <- midpoint - floor(centromere_length / 2)
    c_end <- c_start + centromere_length
    p_end <- c_start - 1
    q_start <- c_end + 1
    q_end <- chrom_length

    # Create a data frame for the p, c, and q regions
    chrom_seg_df <- data.frame(
      haplotype = haplotype,
      chrom = chrom_name,
      ref_start = c(1, c_start, q_start),
      ref_end = c(p_end, c_end, q_end),
      ori_start = c(1, c_start, q_start),
      ori_end = c(p_end, c_end, q_end),
      start = c(1, c_start, q_start),
      end = c(p_end, c_end, q_end),
      region_name = paste0(chrom_name, c("p", "c", "q")),
      copy_index = 1,
      seg_id = paste0(chrom_name, c("p", "c", "q"), "_1"),
      CN_change = 0,
      seg_source_edge = "before_root",
      seg_source_event = "base",
      stringsAsFactors = FALSE
    )

    # Append the data frame to the seg_data_frames list
    seg_data_frames[[i]] <- chrom_seg_df
  }

  # Merge the segment data frames in the list into a single data frame
  seg_table <- do.call(rbind, seg_data_frames)

  return(seg_table)
}

update_wgd_seg <- function(seg_list, chr_lengths, edge_event_table) {
  edge_event_table <- subset(edge_event_table, tolower(region_name) == "wgd")

  for (i in 1:nrow(edge_event_table)) {
    wgd_event_vec <- edge_event_table[i, ]
    haplotype <- wgd_event_vec$haplotype
    edge_name <- paste0(wgd_event_vec$parent, '_', wgd_event_vec$child)

    chr_names <- names(chr_lengths[[haplotype]])

    for (chr_name in chr_names) {
      all_chr_segments <- subset(seg_list[[haplotype]], chrom == chr_name)
      lost_mock_segments <- subset(seg_list[[haplotype]], chrom == chr_name & CN_change == -1)
      remaining_segments <- subset(seg_list[[haplotype]], chrom == chr_name & CN_change != -1)

      if(nrow(remaining_segments) > 0 ){
        # further check if there is any loss segments
        if (nrow(lost_mock_segments) > 0 ) {
          # Detailed processing when there are lost segments and reset the remaining segments
          current_chr_length <- chr_lengths[[haplotype]][chr_name]
          for (l in 1:nrow(lost_mock_segments)) {
            lost_seg <- lost_mock_segments[l,]
            remaining_segments <- remaining_segments[!(remaining_segments$start == lost_seg$ori_start &
                                                         remaining_segments$end == lost_seg$ori_end), ]
          }

          # Create a named vector to store the highest copy index for each region
          highest_copy_index <- tapply(seg_list[[haplotype]]$copy_index, seg_list[[haplotype]]$region_name, max)

          # If there is any remaining segments left, process the data
          if (nrow(remaining_segments) > 0) {
            for (j in 1:nrow(remaining_segments)) {
              segment_length <- remaining_segments$end[j] - remaining_segments$start[j] + 1
              remaining_segments$start[j] <- current_chr_length + 1
              remaining_segments$end[j] <- current_chr_length + segment_length
              current_chr_length <- remaining_segments$end[j]

              # Update copy index for each segment based on its region
              region_name <- remaining_segments$region_name[j]
              if (!is.null(highest_copy_index[region_name])) {
                highest_copy_index[region_name] <- highest_copy_index[region_name] + 1
                remaining_segments$copy_index[j] <- highest_copy_index[region_name]
              }
            }

            chr_lengths[[haplotype]][chr_name] <- current_chr_length
          } else {
            print(paste0("WGD updates: No remaining segments to process for ", haplotype, " ", chr_name, " in ", edge_name, "."))
            next
          }



        }else {
          # Simple duplication when there are no lost segments
          remaining_segments$start <- remaining_segments$start + chr_lengths[[haplotype]][chr_name]
          remaining_segments$end <- remaining_segments$end + chr_lengths[[haplotype]][chr_name]
          chr_lengths[[haplotype]][chr_name] <- chr_lengths[[haplotype]][chr_name]*2
          remaining_segments$copy_index <- ave(remaining_segments$copy_index, remaining_segments$region, FUN = function(x) max(x) + x)

        }

        remaining_segments$seg_source_event <- paste0(remaining_segments$seg_id, ":wgd")
        remaining_segments$seg_source_edge <- edge_name
        remaining_segments$CN_change <- 1
        remaining_segments$seg_id <- paste0(remaining_segments$region_name, "_", remaining_segments$copy_index)

        seg_list[[haplotype]] <- rbind(seg_list[[haplotype]], remaining_segments)

      }else{ # condition of no initial remaining segments
        print("No remaining segments to process.")
        next
      }


    } # for each chr
  } # for each wgd event

  return(list("updated_seg_list" = seg_list, "updated_chr_lengths" = chr_lengths))
}



update_cnv_seg <- function(seg_list, chr_lengths, edge_event_table) {

  # Filter out the non-whole genome duplication events from edge_event_table
  cnv_events <- subset(edge_event_table, region_name != "wgd")

  # For each CNV event
  for (i in 1:nrow(cnv_events)) {

    # Get the current CNV event
    single_cnv_event <- cnv_events[i,]

    # Extract information from the current CNV event
    haplotype <- single_cnv_event$haplotype
    region_name <- single_cnv_event$region_name
    CN_change <- single_cnv_event$CN_change
    copy_index <- single_cnv_event$copy_index

    # Construct edge and seg_id for current CNV event
    edge <- paste0(single_cnv_event$parent, '_', single_cnv_event$child)
    seg_id <- paste0(region_name, "_", copy_index)
    chrom <- gsub("(chr[0-9]+).*", "\\1", region_name)

    # Find the index of the corresponding segment
    segment_indices <- which(seg_list[[haplotype]]$seg_id == seg_id )

    # Get the reference start and end points of the segment
    ref_start <- seg_list[[haplotype]][segment_indices,]$ref_start
    ref_end <- seg_list[[haplotype]][segment_indices,]$ref_end

    # Get the start and end points of the segment
    ori_start <- seg_list[[haplotype]][segment_indices,]$start
    ori_end <- seg_list[[haplotype]][segment_indices,]$end

    # If CN change is -1, mark the segment as inactive
    if(CN_change == -1){
      # Create a mock segment based on the current segment
      mock_segment <- seg_list[[haplotype]][segment_indices,]

      # Update the CN_change to -1 for this mock segment
      mock_segment$CN_change <- -1
      mock_segment$seg_source_edge <- edge
      mock_segment$seg_source_event <- paste0(seg_id, ":", CN_change)

      # update ori_start and ori_end
      new_ori_start <- mock_segment$start
      new_ori_end <- mock_segment$end
      mock_segment$ori_start <- new_ori_start
      mock_segment$ori_end <- new_ori_end

      # turn start and end of loss segment to NA
      mock_segment$start <- mock_segment$end <- NA


      # Append the mock segment to the seg_list
      seg_list[[haplotype]] <- rbind(seg_list[[haplotype]], mock_segment)
    }else{
      # If CN change is positive, add new segments
      seg_length <- ref_end - ref_start + 1
      chrom_length <- chr_lengths[[haplotype]][chrom]
      for (j in 1:CN_change) {
        # Update start and end points for new segment

        new_start <- chrom_length + seg_length*(j-1) + 1
        new_end <- chrom_length + seg_length*j

        # Get new copy index and seg_id for new segment
        new_copy_index <- max(seg_list[[haplotype]]$copy_index[seg_list[[haplotype]]$region_name == region_name]) + 1
        new_seg_id <- paste0(region_name, "_", new_copy_index)
        # Clone new segment from existing one
        new_segment <- seg_list[[haplotype]][segment_indices, ]

        # Update properties of new segment
        new_segment$copy_index <- new_copy_index
        new_segment$start <- new_start
        new_segment$end <- new_end
        new_segment$seg_id <- new_seg_id
        new_segment$CN_change <- CN_change
        new_segment$seg_source_edge <- edge
        new_segment$ori_start <- ori_start
        new_segment$ori_end <- ori_end
        new_segment$seg_source_event <- paste0(seg_id, ":", CN_change)
        # Add the new segment to the list
        seg_list[[haplotype]] <- rbind(seg_list[[haplotype]], new_segment)
      }
      # update chr lengths
      chr_lengths[[haplotype]][chrom] <- chrom_length + CN_change*seg_length
    }

  }
  # Return the updated list of segments and chromosome lengths
  return(list("updated_seg_list" = seg_list, "updated_chr_lengths" = chr_lengths))
}



update_sub_seg <- function(seg_list, sub_seg_list, chr_lengths, one_sub_event) {

  # Extract information from the current sub event
  haplotype <- one_sub_event$haplotype
  region_name <- one_sub_event$region_name  # sub name in this context
  CN_change <- one_sub_event$CN_change
  copy_index <- one_sub_event$copy_index

  # Construct edge and seg_id for current sub event
  edge <- paste0(one_sub_event$parent, '_', one_sub_event$child)
  seg_id <- paste0(region_name, "_", copy_index)
  chrom <- sub_seg_list[[haplotype]]$chrom[sub_seg_list[[haplotype]]$region_name == region_name]  # get chromosome from sub_table
  # Find the index of the corresponding sub segment
  segment_index <- which(sub_seg_list[[haplotype]]$region_name == region_name)
  # Get the reference start and end points of the segment
  ref_start <- sub_seg_list[[haplotype]][segment_index,]$ref_start
  ref_end <- sub_seg_list[[haplotype]][segment_index,]$ref_end

  # Get the original
  ori_start <- sub_seg_list[[haplotype]][segment_index,]$ori_start
  ori_end <- sub_seg_list[[haplotype]][segment_index,]$ori_end
  # If CN change is -1, mark the segment as inactive
  if(CN_change == -1 && copy_index > 0) {
    mock_segment <- seg_list[[haplotype]][segment_index,]
    print("copy_index >0, mock_segment")
    print(mock_segment)
    mock_segment$seg_id = seg_id
    mock_segment$CN_change = -1
    mock_segment$seg_source_edge <- edge
    mock_segment$seg_source_event <- paste0(seg_id, ":", CN_change)

    # update ori_start and ori_end
    new_ori_start <- mock_segment$start
    new_ori_end <- mock_segment$end
    mock_segment$ori_start <- new_ori_start
    mock_segment$ori_end <- new_ori_end

    # turn start and end of loss segment to NA
    mock_segment$start <- mock_segment$end <- NA

    # Append the mock segment to seg_list
    seg_list[[haplotype]] <- rbind(seg_list[[haplotype]], mock_segment)

    # If CN change is -1 and copy_index is 0, clone a segment from sub_seg_list
  }else if(CN_change == -1 && copy_index == 0){

    mock_segment <- sub_seg_list[[haplotype]][segment_index,]
    mock_segment$seg_id = seg_id
    mock_segment$CN_change = -1
    mock_segment$seg_source_edge <- edge
    mock_segment$seg_source_event <- paste0(seg_id, ":", CN_change)

    # update ori_start and ori_end
    new_ori_start <- mock_segment$start
    new_ori_end <- mock_segment$end
    mock_segment$ori_start <- new_ori_start
    mock_segment$ori_end <- new_ori_end

    # turn start and end of loss segment to NA
    mock_segment$start <- mock_segment$end <- NA

    seg_list[[haplotype]] <- rbind(seg_list[[haplotype]], mock_segment)
  }else if(CN_change > 0){ # else case we need to insert sub segments into the seg_list
    # If CN change is positive, add new sub segments
    seg_length <- ref_end - ref_start + 1
    chrom_length <- chr_lengths[[haplotype]][chrom]


    for (j in 1:CN_change) {
      # Update start and end points for new segment
      new_start <- chrom_length + seg_length*(j-1) + 1
      new_end <- chrom_length + seg_length*j

      # Get new copy index and seg_id for new segment
      if(sum(seg_list[[haplotype]]$region_name == region_name) == 0 ){ # add condition to take account into the annotation of loss segments
        new_copy_index = 1
      }else{
        new_copy_index <- max(seg_list[[haplotype]]$copy_index[seg_list[[haplotype]]$region_name == region_name]) + 1
      }
      new_seg_id <- paste0(region_name, "_", new_copy_index)

      # Clone new segment from existing one
      new_segment <- sub_seg_list[[haplotype]][segment_index,]

      # Update properties of new segment
      new_segment$copy_index <- new_copy_index
      new_segment$start <- new_start
      new_segment$end <- new_end
      new_segment$seg_id <- new_seg_id
      new_segment$seg_source_edge <- edge
      new_segment$seg_source_event <- paste0(seg_id, ":", CN_change)
      new_segment$CN_change <- CN_change
      new_segment$ori_start <- ori_start
      new_segment$ori_end <- ori_end

      seg_list[[haplotype]] <- rbind(seg_list[[haplotype]], new_segment)

      # Update chromosome length
      chr_lengths[[haplotype]][chrom] <- chrom_length + CN_change*seg_length
    }
  }
  # Return the updated list of sub segments and chromosome lengths
  return(list("updated_seg_list" = seg_list, "updated_chr_lengths" = chr_lengths))
}
