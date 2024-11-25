# benchmark.R

chisel_read <- function(dat_dir, dat_type = c('barcode', 'clones', 'baf', 'calls', 'rdr', 'combo', 'total')){
  chisel_result <- list()
  for(i in 1:length(dat_type)){

    if(dat_type[i] == 'barcode'){

      tmp <- read.table(paste0(dat_dir, "/barcodedcells.info.tsv"), sep = "\t", header = FALSE)[, 1:2]
      colnames(tmp) <- c("sample", "barcode")
    }else if(dat_type[i] == 'clones'){
      tmp <- read.table(paste0(dat_dir, "/clones/mapping.tsv"), sep = "\t", header = FALSE,
                        col.names = c("barcode", "initial_clone", "final_clone"))
    }else if(dat_type[i] == 'baf'){
      tmp <- read.table(paste0(dat_dir, "/baf/baf.tsv"), sep = "\t", header = FALSE,
                        col.names = c("chr", "pos", "barcode", "a_count", "b_count"))
    }else if(dat_type[i] == 'calls'){
      tmp <- read.table(paste0(dat_dir, "/calls/calls.tsv"), sep = "\t", header = FALSE,
                        col.names = c("chr", "start", "end", "barcode", "normal", "count", "rdr", "a_count", "b_count", "baf", "cluster", "hap_cn", "corrected_hap_cn"))
    }else if(dat_type[i] == 'rdr'){
      tmp <- read.table(paste0(dat_dir, "/rdr/rdr.tsv"), sep = "\t", header = FALSE,
                        col.names = c("chr", "start", "end", "barcode", "normal", "count", "rdr"))

    }else if(dat_type[i] == 'combo'){
      tmp <- read.table(paste0(dat_dir, "/combo/combo.tsv"), sep = "\t", header = FALSE,
                        col.names = c("chr", "start", "end", "barcode", "normal", "count", "rdr", "a_count", "b_count", "baf"))

    }else if(dat_type[i] == 'total'){
      tmp <- read.table(paste0(dat_dir, "/rdr/total.tsv"), sep = "\t", header= FALSE,
                        col.names = c('cell', 'total'))
    }

    chisel_result[[dat_type[i]]] <- tmp
    rm(tmp)
  }

  return(chisel_result)
}

merge_sort_chisel_data <- function(clone_dat, barcode_dat, call_dat) {
  # Merge clone and barcode data
  merged_clone_dat <- merge(clone_dat, barcode_dat, by = "barcode", all.x = TRUE)

  # Set row names for chisel_clone_dat using the 'sample' column
  rownames(merged_clone_dat) <- merged_clone_dat$sample

  # Merge call data with merged_clone_dat
  merged_call_dat <- merge(call_dat, merged_clone_dat, by = "barcode", all.x = TRUE)

  # Derive region_id
  merged_call_dat$region_id <- paste(merged_call_dat$chr, merged_call_dat$start, merged_call_dat$end, sep = "_")

  # Pivot the data to create cn_state_df
  cn_state_df <- merged_call_dat %>%
    dplyr::select(sample, region_id, corrected_hap_cn) %>%
    tidyr::spread(key = region_id, value = corrected_hap_cn)

  # Set the first column as row names
  rownames(cn_state_df) <- cn_state_df[, 1]

  # Remove the first column from the data frame
  cn_state_df <- cn_state_df[, -1]

  # Splitting the region_id into components
  sorted_regions <- region_vector_to_df(colnames(cn_state_df))

  # Reordering the columns of cn_state_df based on sorted regions
  cn_state_df <- cn_state_df[, sorted_regions$region_id]

  # Sort cn_state_df based on clusters identified by chisel
  cn_state_ordered_clone <- merged_clone_dat[rownames(cn_state_df),] %>% dplyr::arrange(final_clone)
  cn_state_df <- cn_state_df[rownames(cn_state_ordered_clone), ]

  # Return the results as a list
  return(list(cn_state_df = cn_state_df, cn_state_ordered_clone = cn_state_ordered_clone, merged_call_dat = merged_call_dat, merged_clone_dat = merged_clone_dat, sorted_regions = sorted_regions))
}


cn_state_split <- function(cn_state_df){
  # separate cn_state_df into two haplotypes
  split_haplotypes <- function(cn_state) {
    haplotypes <- strsplit(as.character(cn_state), split = "\\|")
    sapply(haplotypes, function(x) as.numeric(x))
  }

  # Apply the function to each element of the cn_state_df and transpose the result
  haplotype_cn_arrays <- apply(cn_state_df, c(1:2), split_haplotypes)

  # Creating maternal and paternal matrices
  haplotype_cn <- list()
  haplotype_cn[['hap1']] <- haplotype_cn_arrays[1,,]
  haplotype_cn[['hap2']] <- haplotype_cn_arrays[2,,]
  total_cn <- haplotype_cn[['hap1']] + haplotype_cn[['hap2']]

  # Set the row and column names
  rownames(haplotype_cn[['hap1']]) <- rownames(haplotype_cn[['hap2']]) <- rownames(total_cn) <- rownames(cn_state_df)
  colnames(haplotype_cn[['hap1']]) <- colnames(haplotype_cn[['hap2']]) <- colnames(total_cn) <- colnames(cn_state_df)

  return(list(total_cn = total_cn, haplotype_cn = haplotype_cn))
}


region_vector_to_df <- function(region_vector) {
  # Assuming cn_state_df has a column region_id that contains region identifiers like "chr1_0_5000000"

  # Step 1: Create a dataframe from region_id and separate the CHR, START, and END
  regions <- data.frame(region_id = region_vector, sep_col = region_vector) %>%
    tidyr::separate(col = sep_col, into = c("CHR", "START", "END"), sep = "_") %>%
    dplyr::mutate(
      START = as.numeric(START),
      END = as.numeric(END),
      CHR_idx = as.numeric(gsub("chr", "", CHR))  # Extract numeric part of CHR for sorting
    )

  # Step 2: Sort the regions by chromosome index (CHR_idx) and START position
  sorted_regions <- regions %>%
    dplyr::arrange(CHR_idx, START)

  return(sorted_regions)
}

clone_ground_truth_in_inferred_regions <- function(ground_truth_regions, inferred_regions, ground_truth_cn_mat){
  # Initialize an empty matrix to hold the extended ground truth
  extended_ground_truth <- matrix(0, nrow = nrow(ground_truth_cn_mat), ncol = nrow(inferred_regions))

  # Set the row names and column names of the extended matrix
  rownames(extended_ground_truth) <- rownames(ground_truth_cn_mat)
  colnames(extended_ground_truth) <- inferred_regions$region_id

  # Iterate over each ground truth segment
  for (gt_idx in 1:nrow(ground_truth_regions)) {
    gt_chr <- ground_truth_regions$CHR[gt_idx]
    gt_start <- ground_truth_regions$START[gt_idx]
    gt_end <- ground_truth_regions$END[gt_idx]

    # Find inferred segments within the current ground truth segment
    within_indices <- which(inferred_regions$CHR == gt_chr & inferred_regions$START >= gt_start & inferred_regions$END <= gt_end)

    # Assign the CN value from the ground truth to these fully contained inferred segments
    if (length(within_indices) > 0) {
      extended_ground_truth[, within_indices] <- ground_truth_cn_mat[, gt_idx]
    }

    # Handle inferred segments that cross the boundary of the current ground truth segment
    cross_idx <- which(inferred_regions$CHR == gt_chr & inferred_regions$START > gt_start & inferred_regions$START < gt_end & inferred_regions$END > gt_end)

    if(length(cross_idx) == 1){
      current_start <- inferred_regions$START[cross_idx]
      current_end <- gt_end
      current_gt_idx <- gt_idx

      overlap_length <- current_end - current_start
      inferred_bin_length <- inferred_regions$END[cross_idx] - inferred_regions$START[cross_idx]
      proportion_overlap <- overlap_length / inferred_bin_length
      # Assign proportional CN value for the part within this ground truth segment
      extended_ground_truth[, cross_idx] = extended_ground_truth[, cross_idx] +  proportion_overlap * ground_truth_cn_mat[, current_gt_idx]

      while(inferred_regions$END[cross_idx] > current_end){
        current_gt_idx <- current_gt_idx + 1
        current_start <- ground_truth_regions$START[current_gt_idx]
        current_end <- min(ground_truth_regions$END[current_gt_idx], inferred_regions$END[cross_idx])
        overlap_length <- current_end - current_start
        inferred_bin_length <- inferred_regions$END[cross_idx] - inferred_regions$START[cross_idx]
        proportion_overlap <- overlap_length / inferred_bin_length
        extended_ground_truth[, cross_idx] = extended_ground_truth[, cross_idx] + proportion_overlap * ground_truth_cn_mat[, current_gt_idx]
      }

    }else if(length(cross_idx) > 1){
      warning("More than one overlapping window found for ground truth segment: ", gt_idx)
    }
  }

  extended_ground_truth <- round(extended_ground_truth)
  return(extended_ground_truth)
}

extend_ground_truth_by_clone_identity <- function(clone_ground_truth_cn, sample_clone_identity) {
  # Initialize an empty matrix to hold the extended ground truth
  extended_ground_truth <- matrix(nrow = length(sample_clone_identity), ncol = ncol(clone_ground_truth_cn))

  # Set the column names to match the ground truth CN matrix
  colnames(extended_ground_truth) <- colnames(clone_ground_truth_cn)

  # Traverse through the rownames of the ground truth matrix
  for (clone in rownames(clone_ground_truth_cn)) {
    # Find indices in sample_clone_identity that match the current clone
    clone_indices <- which(sample_clone_identity == clone)

    # Assign the CN values for the current clone to the matching indices in the extended matrix
    if(length(clone_indices) > 0){
      extended_ground_truth[clone_indices, ] <- matrix(rep(clone_ground_truth_cn[clone, ], length(clone_indices)),
                                                       nrow = length(clone_indices), byrow = TRUE)
    }

  }

  rownames(extended_ground_truth) <- names(sample_clone_identity)

  return(extended_ground_truth)
}


# Function to calculate MAE between two matrices
calculate_mae <- function(matrix1, matrix2) {
  # Ensure the matrices have the same dimensions
  if (!all(dim(matrix1) == dim(matrix2))) {
    stop("The two matrices must have the same dimensions")
  }

  # Calculate the absolute differences
  abs_diff <- abs(matrix1 - matrix2)

  # Calculate the mean of the absolute differences
  mae <- mean(abs_diff)

  return(mae)
}


