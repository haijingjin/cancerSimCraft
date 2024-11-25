# event_relevant_functions_v3.R

process_cnv_events <- function(event_table) {
  new_event_table <- event_table[0, ] # Initialize with the same structure but no rows

  for (i in 1:nrow(event_table)) {
    row <- event_table[i, ]
    # Check if the row is a chromosome-level CNV event
    if (row$event_type == "CNV" && grepl("^chr[0-9]+$", row$region_name)) {
      chr_name <- row$region_name
      # Create new rows for the p and q arms
      new_row_p <- row
      new_row_p$region_name <- paste0(chr_name, "p")
      new_event_table <- rbind(new_event_table, new_row_p)

      new_row_q <- row
      new_row_q$region_name <- paste0(chr_name, "q")
      new_event_table <- rbind(new_event_table, new_row_q)
    } else {
      # Add the row as it is if it's not a chromosome-level CNV event
      new_event_table <- rbind(new_event_table, row)
    }
  }

  return(new_event_table)
}

# Create manual events
create_event <- function(haplotype, parent, child, region_name, CN_change, copy_index, event_type) {
  data.frame(haplotype, parent = parent, child = child, region_name = region_name, CN_change = CN_change, copy_index = copy_index, event_type = event_type)
}

# For hybrid event generation
# Insert manual events among random events
insert_event_at_positions <- function(event_table_x, event_table_y, insert_pos) {
  # Check if the number of positions matches the number of rows in event_table_y
  if(length(insert_pos) != nrow(event_table_y)) {
    stop("The number of positions should match the number of rows in event_table_y.")
  }

  # Create a helper column for position in event_table_x
  event_table_x$insert_pos <- 1:nrow(event_table_x)

  # Adjust positions for event_table_y based on where they need to be inserted
  adjusted_insert_pos <- insert_pos + seq_along(insert_pos)/10  # Add a small fraction to avoid overlap

  # Assign these positions to event_table_y
  event_table_y$insert_pos <- adjusted_insert_pos

  # Combine the tables and arrange by position
  combined_table <- bind_rows(event_table_x, event_table_y) %>%
    arrange(insert_pos) %>%
    select(-insert_pos)  # Drop the helper column

  return(combined_table)
}


insert_event_at_positions <- function(event_table_x, event_table_y, insert_pos) {
  if(length(insert_pos) != nrow(event_table_y)) {
    stop("The number of positions should match the number of rows in event_table_y.")
  }

  # Create a helper column for position in event_table_x
  event_table_x$insert_pos <- 1:nrow(event_table_x)

  # Adjust positions for event_table_y based on where they need to be inserted
  # This will handle duplicate positions by incrementing slightly more for each duplicate
  adjusted_insert_pos <- sapply(insert_pos, function(x) {
    mean(x + which(insert_pos == x)/10)
  })

  # Assign these positions to event_table_y
  event_table_y$insert_pos <- adjusted_insert_pos

  # Combine the tables and arrange by position
  combined_table <- dplyr::bind_rows(event_table_x, event_table_y) %>%
    dplyr::arrange(insert_pos) %>%
    dplyr::select(-insert_pos)  # Drop the helper column

  return(combined_table)
}

create_edge_event_table <- function(event_table, tree, anno_cols = c("haplotype_abbr", "region_name", "CN_change")) {
  tree_table <- data.frame(as_edgelist(tree))
  colnames(tree_table) <- c("parent", "child")
  tree_table$edge_name <- paste(tree_table$parent, tree_table$child, sep = "_")

  anno_cols_symbols <- lapply(anno_cols, sym)

  edge_event_table <- event_table %>%
    mutate(haplotype_abbr = toupper(substr(haplotype, 1, 1))) %>%  # Extract and convert the first letter to uppercase
    mutate(event_name = paste(!!!anno_cols_symbols, sep = ":")) %>%
    mutate(edge_name = paste(parent, child, sep = "_")) %>%  # Create 'edge_name' from 'parent' and 'child'
    group_by(edge_name) %>%
    summarize(edge_label = paste(event_name, collapse = "\n"), n_events = n()) %>%
    arrange(match(edge_name, tree_table$edge_name))  # match returns the position of each element of 'edge_name' in 'edge_order'

  return(edge_event_table)
}


