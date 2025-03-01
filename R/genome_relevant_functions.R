#' Synthesize Founder Genomes for a Cell Lineage Tree
#'
#' @description
#' Synthesizes founder genomes for each clone in a cell lineage tree by traversing
#' the tree from root to leaves, incorporating appropriate mutations and structural
#' variants at each step.
#'
#' @param tree An igraph object representing the cell lineage tree with vertices as clone types
#' @param root_genome List containing the genome sequences for the root clone
#' @param seg_list List structure containing segment information for all clones
#' @param mutation_info Data frame containing mutation details for all cells
#' @param cell_info Data frame containing cell lineage information
#'
#' @return A list containing four elements:
#' \itemize{
#'   \item all_node_genomes: List of genome sequences for each clone in the tree
#'   \item child_node_founders_df: Data frame containing information about the founder cells for each clone
#'   \item child_node_founder_mutations: List of all mutations in founder cells, organized by clone
#'   \item child_node_founder_mutations_filtered: List of mutations that occurred between parent and child nodes
#' }
#'
#' @details
#' The function builds genomes sequentially following the evolutionary relationships
#' defined in the tree structure, ensuring that each clone's genome properly inherits
#' all modifications from its ancestors plus its unique changes.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a cell lineage tree and other required objects
#' founder_genomes <- synth_sc_tree_founder_genomes(
#'   tree = clone_tree,
#'   root_genome = reference_genome,
#'   seg_list = all_segments,
#'   mutation_info = mutations,
#'   cell_info = cell_lineage
#' )}
#'
#'
#' @seealso
#' \code{\link{synth_sc_founder_genome}}, \code{\link{introduce_snv_sc}},
#' \code{\link{get_mutations_sc}}
#'
#' @importFrom igraph V degree dfs get.adjlist
#'
#' @export
synth_sc_tree_founder_genomes <- function(tree, root_genome, seg_list, mutation_info, cell_info) {

  root_name <- names(V(tree))[degree(tree, mode = "in") == 0]

  all_node_genomes <- list()
  all_node_genomes[[root_name]] <- root_genome

  # Create an empty data frame to store child_node_founder information
  child_node_founders_df <- data.frame()

  child_node_founder_mutations <- list()
  child_node_founder_mutations_filtered <- list()

  dfs_order <- dfs(tree, root = root_name, mode = "out")$order

  for(i in 2:length(dfs_order)) {
    child_node_index <- dfs_order[i]
    parent_node_index <- get.adjlist(tree, mode = "in")[[child_node_index]]
    parent_node <- V(tree)[parent_node_index]$name
    child_node <- V(tree)[child_node_index]$name

    # Acquire the first cell of the child_node
    child_node_founder <- cell_info[cell_info$clone == child_node & cell_info$clone_index == 1, ]

    # Append child_node_founder to data frame
    child_node_founders_df <- rbind(child_node_founders_df, child_node_founder)

    # Extract all ancestors of child node founder
    child_node_founder_mut <- get_mutations_sc(cell_info = cell_info,
                                               mutation_info = mutation_info,
                                               sampled_cell_idx = child_node_founder$cell_index)
    founder_index_char <- as.character(child_node_founder$cell_index)
    child_node_founder_mut <- child_node_founder_mut$sampled_sc_mutations[[founder_index_char]]
    child_node_founder_mutations[[child_node]] <- child_node_founder_mut

    # Filter the SNVs that happened between parent node and child node
    parent_node_founder <- cell_info[cell_info$clone == parent_node & cell_info$clone_index == 1, ]
    child_node_founder_mut_filtered <- child_node_founder_mut[child_node_founder_mut$cell_index > parent_node_founder$cell_index &
                                                                child_node_founder_mut$cell_index <= child_node_founder$cell_index, ]
    child_node_founder_mutations_filtered[[child_node]] <- child_node_founder_mut_filtered

    # Start with the parent genome
    child_genome <- all_node_genomes[[parent_node]]

    # Introduce SNVs
    if (!is.null(child_node_founder_mut_filtered)) {
      child_genome <- introduce_snv_sc(input_genome = child_genome, sc_mut_table_with_nt = child_node_founder_mut_filtered)
      print("Pre SNVs introduced to child genome")
    } else {
      print("No new SNVs to introduce")
    }

    all_node_genomes[[child_node]] <- synth_sc_founder_genome(child_node = child_node,
                                                              parent_node = parent_node,
                                                              pre_child_genome = child_genome,
                                                              seg_list = seg_list)

    print(paste0("Synthesis of ", child_node, " genome is finished!"))
  }

  return(list(
    all_node_genomes = all_node_genomes,
    child_node_founders_df = child_node_founders_df,
    child_node_founder_mutations = child_node_founder_mutations,
    child_node_founder_mutations_filtered = child_node_founder_mutations_filtered
  ))
}

#' Synthesize Founder Genome with Structural Variants
#'
#' @description
#' Incorporates structural variants (deletions, duplications) into a genome during
#' the transition from parent to child clone in a cell lineage tree.
#'
#' @param child_node Character string specifying the name of the child clone
#' @param parent_node Character string specifying the name of the parent clone
#' @param pre_child_genome List containing the genome sequences after SNVs have been incorporated
#' @param seg_list List structure containing segment information for all clones
#'
#' @return A list containing the modified genome with structural variants incorporated
#'
#' @details
#' This function applies structural variants to a genome during clone evolution by:
#'
#' 1. Starting with a genome that already has SNVs incorporated (pre_child_genome)
#' 2. Identifying segments that are associated with the specific transition edge
#'    between parent and child clones
#' 3. Processing each relevant segment based on its copy number change (CN_change):
#'    - For deletions (CN_change = -1): Replaces the original sequence with "N" characters
#'    - For duplications (CN_change >= 1): Copies the sequence and adds it to the end of the chromosome
#'
#' The function works on each haplotype (maternal/paternal) separately and modifies
#' only chromosomes that have structural variant events associated with the specific
#' parent-to-child clone transition.
#'
#' This function is typically called after SNVs have been incorporated into the genome
#' but before cell-specific mutations are added, representing the genetic changes that
#' define a new clone's emergence.
#'
#' @seealso
#' \code{\link{synth_sc_tree_founder_genomes}}, \code{\link{introduce_snv_sc}}
#'
#' @importFrom Biostrings subseq
#'
#' @export
synth_sc_founder_genome <- function(child_node, parent_node, pre_child_genome, seg_list){

  # Initialize the target genome as a copy of the nearest genome
  target_genome <- pre_child_genome

  # Get the edge between the target_clone and the parent clone
  edge = paste0(parent_node,"_", child_node)

  for(haplotype in names(seg_list[[child_node]])){
    haplotype_segments <- seg_list[[child_node]][[haplotype]]


    # Filter segments based on the events in the current edge
    relevant_segments <- haplotype_segments[haplotype_segments$seg_source_edge == edge, ]

    if(nrow(relevant_segments) == 0) {
      print(paste("No new segments formed for haplotype", haplotype, "in the edge", edge))
      next
    }

    # Add new segments to the genome
    for(i in 1:nrow(relevant_segments)){
      segment <- relevant_segments[i, ]
      chrom <- segment$chrom
      ori_start <- segment$ori_start
      ori_end <- segment$ori_end
      seg_source_event <- segment$seg_source_event
      CN_change <- segment$CN_change

      # Extract the corresponding sequence from the target genome
      target_seq <- target_genome[[haplotype]][[chrom]]

      if(CN_change == -1){

        target_seq[ori_start:ori_end] <- "N"
        target_genome[[haplotype]][[chrom]] <- target_seq  # Update the genome with the modified sequence

      }else if(CN_change >= 1){

        # Attach the original segment to the new segment
        new_segment_seq <- subseq(target_seq, start = ori_start, end = ori_end)

        # Add the new segment sequence to the genome
        target_genome[[haplotype]][[chrom]] <- c(target_genome[[haplotype]][[chrom]], new_segment_seq)
      }
    }

  }

  return(target_genome)
}

#' Synthesize Single-Cell Genome
#'
#' This function synthesizes a single-cell genome by introducing single nucleotide variants (SNVs)
#' into a backbone genome. It retrieves cell-specific mutations, filters mutations that occurred
#' after the founder cell, and introduces the mutations into the genome.
#'
#' @param cell_index Integer specifying the index of the cell for which to synthesize a genome
#' @param backbone_genome List containing the reference genome sequences of the cell's clone
#' @param sc_mut_list List where each element corresponds to a cell's mutations, with cell indices as names
#' @param mut_table_with_nt Data frame containing mutation details with nucleotide changes
#' @param founder_cell_index Integer specifying the founder cell index of the clone
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item sc_genome: The synthesized single cell genome with mutations incorporated
#'   \item cell_mut_with_nt: Data frame containing the cell-specific mutations with nucleotide details
#' }
#'
#'
#' @seealso
#' \code{\link{find_identical_row_index}}, \code{\link{introduce_snv_sc}},
#' \code{\link{simulate_sc_dynamic_reads_for_batches}}
#'
#' @export
synth_sc_genome <- function(cell_index, backbone_genome, sc_mut_list, mut_table_with_nt, founder_cell_index){
  # Retrieve sc specific mutations for the given cell
  cell_mut <- sc_mut_list[[as.character(cell_index)]]
  cell_mut_colnames <- colnames(cell_mut)

  # Filter mutations that occurred after the founder cell
  cell_mut <- cell_mut[cell_mut$cell_index > founder_cell_index, ]

  # Initialize a data frame to store mutations with nucleotide info
  cell_mut_with_nt <- data.frame()

  # Loop over each mutation in the filtered cell_mut list and retrieve nucleotide info
  if(nrow(cell_mut) > 0){
    for(i in 1:nrow(cell_mut)){
      # Find the mutation nucleotide info in the mut_table_with_nt

      mut_idx <- find_identical_row_index(cell_mut[i, ], mut_table_with_nt[, cell_mut_colnames])
      if (length(mut_idx) > 0) { # Check to avoid missing or unmatched mutations
        cell_mut_with_nt <- rbind(cell_mut_with_nt, mut_table_with_nt[mut_idx, ])
      } else {
        stop(paste("Error: No matching nucleotide information found for mutation at row", cell_mut[i,]))
      }
    }

    # Process and introduce the single nucleotide variants (SNVs)
    sc_genome <- introduce_snv_sc(input_genome = backbone_genome, sc_mut_table_with_nt = cell_mut_with_nt)

  }else{
    print("0 cell-specific mut")
    sc_genome <- backbone_genome
  }
  return(list(sc_genome = sc_genome, cell_mut_with_nt = cell_mut_with_nt))
}


synth_clone_genome_for_sc <- function(target_clone, nearest_genome, nearest_clone, tree, seg_list){
  # Initialize the target genome as a copy of the nearest genome
  target_genome <- nearest_genome

  # Get the edges between the nearest clone and the target clone
  edges <- get_edges_between_clones(tree, upper_node = nearest_clone, lower_node = target_clone)

  # Traverse the edges between the two clones
  for(edge in edges){

    # Traverse the segment list to add new segments based on the events
    for(haplotype in names(seg_list[[target_clone]])){
      haplotype_segments <- seg_list[[target_clone]][[haplotype]]


      # Filter segments based on the events in the current edge
      relevant_segments <- haplotype_segments[haplotype_segments$seg_source_edge == edge, ]

      if(nrow(relevant_segments) == 0) {
        print(paste("No new segments formed for haplotype", haplotype, "in the edge", edge))
        next
      }

      # Add new segments to the genome
      for(i in 1:nrow(relevant_segments)){
        segment <- relevant_segments[i, ]
        chrom <- segment$chrom
        ori_start <- segment$ori_start
        ori_end <- segment$ori_end
        seg_source_event <- segment$seg_source_event
        CN_change <- segment$CN_change
        # Extract the corresponding sequence from the target genome
        target_seq <- target_genome[[haplotype]][[chrom]]

        if(CN_change == -1){
          target_seq[ori_start:ori_end] <- "N"
        }else if(CN_change >= 1){

          # Attach the original segment to the new segment
          new_segment_seq <- subseq(target_seq, start = ori_start, end = ori_end)

          # Add the new segment sequence to the genome
          target_genome[[haplotype]][[chrom]] <- c(target_genome[[haplotype]][[chrom]], new_segment_seq)
        } # if CN_change >= 1
      } # iterate over new segments
    } # iterate over haplotypes
  } # iterate over edges
  return(target_genome)
}

#' Synthesize Clone Genome Based on Events and Mutations
#'
#' @description
#' Constructs a clone's genome by applying sequential genomic alterations (CNVs, WGDs)
#' and mutations (SNVs) starting from a nearest ancestral genome. Processes events along the
#' evolutionary path between the nearest ancestor and target clone.
#'
#' @param target_clone Character. Name of the clone whose genome is to be synthesized.
#' @param nearest_genome Nested list containing the nearest ancestor's genome sequences:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: chromosomes with nucleotide sequences
#'   }
#' @param nearest_clone Character. Name of the nearest ancestor clone.
#' @param tree igraph object. Phylogenetic tree structure.
#' @param seg_list Nested list containing segment information for all clones:
#'   \itemize{
#'     \item First level: clone names
#'     \item Second level: haplotypes
#'     \item Each haplotype contains segment information (data frame)
#'   }
#' @param mut_table Data frame containing mutation information with columns:
#'   \itemize{
#'     \item edge_name - Tree edge identifier
#'     \item other mutation-specific columns required by introduce_snv
#'   }
#'
#' @return A nested list containing the synthesized genome with the same structure
#'   as nearest_genome, but updated with all genomic changes.
#'
#' @details
#' The function:
#' 1. Identifies edges between nearest ancestor and target clone
#' 2. For each edge:
#'    - Processes segment changes (copy number variations)
#'    - Handles segment losses (marks with "N")
#'    - Adds new segments for gains
#'    - Applies SNVs if present
#' 3. Maintains separate tracking for maternal and paternal haplotypes
#'
#' @seealso
#' \code{\link{get_edges_between_clones}}, \code{\link{introduce_snv}}
#'
#' @importFrom Biostrings subseq
#' @export
synth_clone_genome <- function(target_clone, nearest_genome, nearest_clone, tree, seg_list, mut_table, verbose = TRUE){

  # Initialize the target genome as a copy of the nearest genome
  target_genome <- nearest_genome
  rm(nearest_genome)
  # Get the edges between the nearest clone and the target clone
  edges <- get_edges_between_clones(tree, upper_node = nearest_clone, lower_node = target_clone)

  # Traverse the edges between the two clones
  for(edge in edges){
    # Traverse the segment list to add new segments based on the events
    for(haplotype in names(seg_list[[target_clone]])){
      haplotype_segments <- seg_list[[target_clone]][[haplotype]]


      # Filter segments based on the events in the current edge
      relevant_segments <- haplotype_segments[haplotype_segments$seg_source_edge == edge, ]
      if(nrow(relevant_segments) == 0) {
        if(verbose){
          print(paste("No new segments formed for haplotype", haplotype, "in the edge", edge))
        }
        next
      }

      # Add new segments to the genome
      for(i in 1:nrow(relevant_segments)){
        segment <- relevant_segments[i, ]
        chrom <- segment$chrom
        # Add for vignette
        if (!chrom %in% names(target_genome[[haplotype]])) {
          if(verbose){
            print(paste("Chromosome", chrom, "not found in target genome for haplotype", haplotype))
          }
          next
        }

        ori_start <- segment$ori_start
        ori_end <- segment$ori_end
        seg_source_event <- segment$seg_source_event
        CN_change <- segment$CN_change
        # Extract the corresponding sequence from the target genome
        target_seq <- target_genome[[haplotype]][[chrom]]

        if(CN_change == -1){

          target_seq[ori_start:ori_end] <- "N"
          target_genome[[haplotype]][[chrom]] <- target_seq  # Update the genome with the modified sequence

        }else if(CN_change >= 1){
          # Attach the original segment to the new segment
          new_segment_seq <- subseq(target_seq, start = ori_start, end = ori_end)

          # Add the new segment sequence to the genome
          target_genome[[haplotype]][[chrom]] <- c(target_genome[[haplotype]][[chrom]], new_segment_seq)
        }
      }

    }
    edge_mutations <- mut_table[mut_table$edge_name == edge, ]

    if(nrow(edge_mutations) > 0){
      snv_result <- introduce_snv(genome = target_genome, mut_table = edge_mutations, verbose = verbose)
      target_genome <- snv_result$genome
    }
  }

  return(target_genome)
}




#' Introduce Single Nucleotide Variants into a Single Cell Genome
#'
#' @description
#' Incorporates single nucleotide variants (SNVs) into a genome sequence, handling
#' both regular mutations and recurrent mutations at the same genomic position.
#'
#' @param input_genome List containing genome sequences to be modified
#' @param sc_mut_table_with_nt Data frame containing mutation details with columns:
#'        haplotype, chrom, pos, original_nt, alternative_nt, time
#'
#' @return A list containing the modified genome sequences with all mutations incorporated
#'
#' @details
#' This function introduces mutations into a genome sequence with special handling
#' for recurrent mutations (multiple mutations at the same genomic position):
#'
#' 1. Identifies recurrent mutations by detecting duplicated positions across
#'    haplotype, chromosome, and position
#' 2. Processes non-recurrent mutations first using the insert_mutations() helper function
#' 3. For recurrent mutations:
#'    - Groups mutations by their genomic coordinates
#'    - For each position with multiple mutations, creates a consolidated mutation that
#'      applies the cumulative effect (using the original nucleotide from the first mutation
#'      and the alternative nucleotide from the last mutation in the timeline)
#'    - Applies these consolidated recurrent mutations to the genome
#'
#' This approach ensures that the final genome correctly represents the cumulative
#' effect of sequential mutations at the same position, rather than applying each
#' mutation independently which could lead to incorrect results.
#'
#' @seealso
#' \code{\link{insert_mutations}}, \code{\link{synth_sc_genome}}
#'
#' @export
introduce_snv_sc <- function(input_genome, sc_mut_table_with_nt) {
  # Initialize a copy of the genome to store the modified sequences
  modified_genome <- input_genome

  # Identify recurrent mutations
  recurrent_indices <- which(duplicated(sc_mut_table_with_nt[c("haplotype", "chrom", "pos")]) | duplicated(sc_mut_table_with_nt[c("haplotype", "chrom", "pos")], fromLast = TRUE))

  # Process non-recurrent mutations
  if(length(recurrent_indices) > 0){
    non_recurrent_mutations <- sc_mut_table_with_nt[-recurrent_indices, ]
  }else{
    non_recurrent_mutations <- sc_mut_table_with_nt
  }

  # First insert non_recurrent_mutations
  modified_genome <- insert_mutations(modified_genome, non_recurrent_mutations)

  # Process recurrent mutations
  if (length(recurrent_indices) > 0) {
    # Subset the data to only include recurrent mutations
    recurrent_mutations <- sc_mut_table_with_nt[recurrent_indices, ]

    # Sort by coordinates and time
    recurrent_mutations <- recurrent_mutations[order(recurrent_mutations$haplotype, recurrent_mutations$chrom, recurrent_mutations$pos, recurrent_mutations$time),]

    # Initialize an empty dataframe to store the transformed mutations
    transformed_recurrent_mutations <- data.frame()

    # Loop through each unique coordinate
    unique_coords <- unique(paste(recurrent_mutations$haplotype, recurrent_mutations$chrom, recurrent_mutations$pos, sep="_"))
    for(coord in unique_coords) {
      # Subset mutations for this coordinate
      coord_mutations <- subset(recurrent_mutations, paste(haplotype, chrom, pos, sep="_") == coord)

      # Extract the first and last mutations
      first_mutation <- coord_mutations[1, ]
      last_mutation <- tail(coord_mutations, n=1)

      # Create a new mutation record with the original_nts from the first and alternative_nts from the last
      new_mutation <- first_mutation
      new_mutation$alternative_nts <- last_mutation$alternative_nts

      # Append to the transformed mutations dataframe
      transformed_recurrent_mutations <- rbind(transformed_recurrent_mutations, new_mutation)
    }

    # Now, transformed_recurrent_mutations contains the desired data
    modified_genome <- insert_mutations(modified_genome, transformed_recurrent_mutations)
  }

  return(modified_genome)
}


#' Insert Mutations into a Genome
#'
#' This function inserts single nucleotide variants (SNVs) into a genome at specified positions.
#' It groups mutations by haplotype and chromosome for efficient processing and modifies the genome
#' sequence accordingly.
#'
#' @param genome List containing genome sequences organized by haplotype and chromosome
#' @param mutations Data frame containing mutation details with columns:
#'        haplotype, chrom, pos, alternative_nt
#'
#' @return The modified genome list with mutations incorporated
#'
#' @details
#' The function is designed to be efficient when processing large numbers of mutations
#' by minimizing the number of times each sequence needs to be modified.
#'
#' @importFrom Biostrings replaceLetterAt
#'
#' @export
insert_mutations <- function(genome, mutations) {
  # Group mutations by haplotype and chromosome for efficient processing
  grouped_mutations <- split(mutations, list(mutations$haplotype, mutations$chrom))

  # Loop through each group of mutations
  for(group_name in names(grouped_mutations)) {
    mutations <- grouped_mutations[[group_name]]

    # Extract haplotype and chromosome from the group name
    haplotype_chrom <- strsplit(group_name, split = "\\.")
    haplotype <- haplotype_chrom[[1]][1]
    chrom <- haplotype_chrom[[1]][2]

    # Extract the sequence for this haplotype and chromosome
    sequence <- genome[[haplotype]][[chrom]]

    genome[[haplotype]][[chrom]] <- replaceLetterAt(sequence, at = mutations$pos, letter = mutations$alternative_nt)
  }
  return(genome)
}

#' Introduce Single Nucleotide Variants into a Genome
#'
#' @description
#' Introduce single nucleotide variants (SNVs) into a genome sequence by
#' processing mutations grouped by haplotype and chromosome.
#'
#' @param genome A nested list containing genome sequences:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: chromosomes with nucleotide sequences
#'   }
#' @param mut_table A data frame containing mutation information with columns:
#'   \itemize{
#'     \item haplotype - Maternal or paternal copy
#'     \item chrom - Chromosome name
#'     \item pos - Position where mutation occurs
#'     \item alternative_nt - Alternative nucleotide to introduce
#'   }
#'
#' @return A list containing:
#'   \itemize{
#'     \item genome - Modified genome with introduced SNVs
#'   }
#'
#' @details
#' The function:
#' 1. Groups mutations by haplotype and chromosome for efficient processing
#' 2. Processes each group of mutations in a vectorized operation
#' 3. Maintains the original genome structure while updating sequences
#'
#' Note: Uses vectorized operations through replaceLetterAt for efficient mutation introduction
#'
#' @seealso
#' \code{\link{replaceLetterAt}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example mutation table
#' mutations <- data.frame(
#'   haplotype = "maternal",
#'   chrom = "chr1",
#'   pos = 100,
#'   alternative_nt = "A"
#' )
#' modified <- introduce_snv(genome, mutations)
#' }
introduce_snv <- function(genome, mut_table, verbose = TRUE) {
  # Initialize a copy of the genome to store the modified sequences
  modified_genome <- genome

  # Group mutations by haplotype and chromosome for efficient processing
  grouped_mutations <- split(mut_table, list(mut_table$haplotype, mut_table$chrom))

  # Loop through each group of mutations
  for(group_name in names(grouped_mutations)) {
    mutations <- grouped_mutations[[group_name]]
    # Extract haplotype and chromosome from the group name
    haplotype_chrom <- strsplit(group_name, split = "\\.")
    haplotype <- haplotype_chrom[[1]][1]
    chrom <- haplotype_chrom[[1]][2]

    # Check if the chromosome exists in the genome
    if (!chrom %in% names(modified_genome[[haplotype]])) {
      if(verbose){
        print(paste("Chromosome", chrom, "not found in genome for haplotype", haplotype))
      }
      next
    }

    # Extract the sequence for this haplotype and chromosome
    sequence <- modified_genome[[haplotype]][[chrom]]
    # Vectorized operation to introduce SNVs
    # Update the modified genome
    modified_genome[[haplotype]][[chrom]] <- replaceLetterAt(sequence, at = mutations$pos, letter = mutations$alternative_nt)
  }

  return(list("genome" = modified_genome))
}


#' Synthesize Genomes for All Clones in Phylogenetic Tree
#'
#' @description
#' Generates complete genome sequences for all clones in a phylogenetic tree, starting from
#' the root genome and applying sequential genomic alterations following the tree structure.
#'
#' @param tree An igraph object representing the phylogenetic tree structure
#' @param root_genome A nested list containing the initial genome sequences:
#'   \itemize{
#'     \item First level: haplotypes (maternal/paternal)
#'     \item Second level: chromosomes with nucleotide sequences
#'   }
#' @param seg_list A nested list containing segment information for all clones:
#'   \itemize{
#'     \item First level: clone names
#'     \item Second level: haplotypes
#'     \item Each haplotype contains segment information (data frame)
#'   }
#' @param mut_table A data frame containing mutation information for all clones
#'
#' @return A list where:
#'   \itemize{
#'     \item Names are clone names from the tree
#'     \item Values are genome sequences for each clone
#'   }
#'
#' @details
#' The function:
#' 1. Identifies the root node of the tree
#' 2. Processes nodes in depth-first search order
#' 3. For each node:
#'    - Identifies its parent node
#'    - Synthesizes its genome based on parent's genome
#'    - Applies all genomic changes along the branch
#' 4. Tracks progress with print statements
#'
#' @seealso
#' \code{\link{synth_clone_genome}}, \code{\link{dfs}}
#'
#' @importFrom igraph V degree dfs get.adjlist
#' @export
synth_tree_genomes <- function(tree, root_genome, seg_list, mut_table) {

  root_name <- names(V(tree))[degree(tree, mode = "in") == 0]
  all_node_genomes <- list()
  all_node_genomes[[root_name]] <- root_genome
  dfs_order <- dfs(tree, root = root_name, mode = "out")$order

  for(i in 2:length(dfs_order)){
    child_node_index <- dfs_order[i]
    parent_node_index <- get.adjlist(tree, mode = "in")[[child_node_index]]
    parent_node <- V(tree)[parent_node_index]$name
    child_node <- V(tree)[child_node_index]$name
    all_node_genomes[[child_node]] <- synth_clone_genome(target_clone = child_node,
                                                         nearest_genome = all_node_genomes[[parent_node]],
                                                         nearest_clone = parent_node,
                                                         tree = tree,
                                                         seg_list = seg_list,
                                                         mut_table = mut_table)

    print(paste0("Synthesize of ", child_node, " genome is finished!"))
  }

  return(all_node_genomes)
}

#' Get Edge Names Between Two Nodes in Tree
#'
#' @description
#' Finds and returns the names of all edges in the path between two nodes in a phylogenetic tree,
#' formatted as "parent_child" strings.
#'
#' @param tree An igraph object representing the phylogenetic tree structure
#' @param upper_node Character. Name of the starting (ancestor) node
#' @param lower_node Character. Name of the ending (descendant) node
#'
#' @return A character vector containing edge names in the format "parent_child"
#'   for all edges in the path from upper_node to lower_node
#'
#' @details
#' The function:
#' 1. Finds the shortest path between the two nodes
#' 2. Identifies all edges along this path
#' 3. Formats edge names as "parent_child"
#'
#' Note: Assumes that the nodes are connected and upper_node is an ancestor of lower_node
#' in the tree structure.
#'
#' @importFrom igraph shortest_paths E ends
#' @export
#'
#' @examples
#' \dontrun{
#' # For a tree with path A -> B -> C
#' edges <- get_edges_between_clones(tree, "A", "C")
#' # Returns: c("A_B", "B_C")
#' }
get_edges_between_clones <- function(tree, upper_node, lower_node) {
  # Retrieve path from upper_node to lower_node
  path_between_nodes <- shortest_paths(tree, from = upper_node, to = lower_node)$vpath[[1]]

  # Get the edges in the path
  edge_ids <- E(tree)[path_between_nodes[-length(path_between_nodes)] %--% path_between_nodes[-1]]

  # Convert the edges to "node1_node2" format
  edge_endpoints <- ends(tree, edge_ids)
  edge_strings <- paste(edge_endpoints[,1], edge_endpoints[,2], sep = "_")

  return(edge_strings)
}
