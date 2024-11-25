# test_and_check_v2.R

get_downstream_vertices <- function(g, vertex_name) {
  # Check if vertex exists
  if (!(vertex_name %in% V(g)$name)) {
    stop("No vertex found with the given name!")
  }

  # Initialize the list of downstream vertices
  downstream_vertices <- c()

  # Recursive function to traverse downstream
  traverse_downstream <- function(vertex_name) {
    children <- neighbors(g, vertex_name, mode = "out")

    if (length(children) > 0) {
      downstream_vertices <<- c(downstream_vertices, V(g)[children]$name)
      lapply(V(g)[children]$name, traverse_downstream)
    }
  }

  traverse_downstream(vertex_name)
  return(downstream_vertices)
}


event_sanity_check <- function(event_table, tree) {
  # Traverse loss events for sanity check
  conflicting_event_list <- list()
  for(i in 1:nrow(event_table)){
    event <- event_table[i, ]

    downstream_events <- data.frame()
    # check if the event is a loss event
    if(event$CN_change == -1){

      downstream_vertices <- get_downstream_vertices(tree, event$child)
      downstream_events <- rbind(downstream_events, event_table[(i+1):nrow(event_table), ][event_table$child %in% c(event$child, downstream_vertices), ])

      # check for the conflicting events
      conflicting_events <- subset(downstream_events, region_name == event$region_name & haplotype == event$haplotype & copy_index == event$copy_index)

      if(nrow(conflicting_events) > 0){
        conflicting_event_list[[as.character(i)]] <- rbind(event, conflicting_events)
      }
    }
  }

  if(length(conflicting_event_list) == 0) {
    return(print("No conflicting events found."))
  } else {
    warning("Conflicting events detected!")
    return(conflicting_event_list)
  }
}

# Function to process all segments across all clones and haplotypes
all_segments_sanity_check <- function(all_node_segments) {
  sanity_check_results <- list()

  # Iterate over each clone
  for (clone_name in names(all_node_segments)) {
    clone_segments <- all_node_segments[[clone_name]]

    # Iterate over each haplotype within the clone
    for (haplotype_name in names(clone_segments)) {
      segment_table <- clone_segments[[haplotype_name]]
      genome_name <- paste(clone_name, haplotype_name, sep = "_")

      # Perform sanity check and store the results
      result <- segments_sanity_check(segment_table, genome_name)
      sanity_check_results[[genome_name]] <- result
    }
  }

  return(sanity_check_results)
}


segments_sanity_check_old <- function(segment_table, genome_name) {
  # Initialize a list to hold conflicting segments
  conflicting_segment_list <- list()

  # Loop over each segment in the table
  for (i in 1:(nrow(segment_table)-1)) {
    segment <- segment_table[i, ]

    # Check if the segment is a loss segment
    if (segment$CN_change == -1) {
      # Define the range of the lost segment

      # Define the range of the lost segment
      lost_start <- segment$ori_start
      lost_end <- segment$ori_end

      # Find downstream segments that overlap with the lost segment
      downstream_segments <- segment_table[(i+1):nrow(segment_table), ]
      conflicting_segments <- downstream_segments[downstream_segments$chrom == segment$chrom &
                                                    downstream_segments$haplotype == segment$haplotype &
                                                    downstream_segments$ori_start <= lost_end &
                                                    downstream_segments$ori_end >= lost_start, ]
      # Check for conflicting segments
      if (nrow(conflicting_segments) > 0) {
        conflicting_segment_list[[as.character(i)]] <- rbind(segment, conflicting_segments)
      }
    }
  }

  # Check if any conflicting segments were found
  if (length(conflicting_segment_list) == 0) {
    return(paste0("No conflicting segments found on ", genome_name, "."))
  } else {
    warning(paste0("Conflicting segments detected on ", genome_name, "."))
    return(conflicting_segment_list)
  }
}

