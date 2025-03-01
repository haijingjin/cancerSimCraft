
#' Process and Expand Chromosome-Level CNV Events
#'
#' @description
#' Processes a table of genomic events, expanding chromosome-level CNV events
#' into separate events for p and q arms while preserving other events.
#'
#' @param event_table A data frame of genomic events containing:
#'        \itemize{
#'          \item event_type: Type of event (e.g., "CNV", "WGD")
#'          \item region_name: Region identifier or chromosome name in format "chrN"
#'                           where N is the chromosome number
#'          \item Other event-specific columns that will be preserved
#'        }
#'
#' @return A data frame with the same structure as input but with chromosome-level
#'         CNV events expanded into p and q arm events. For example:
#'         \itemize{
#'           \item A CNV event for "chr1" becomes two events for "chr1p" and "chr1q"
#'           \item Non-CNV events or non-chromosome-level CNVs remain unchanged
#'         }
#'
#' @details
#' For each row in the input table:
#' 1. If it's a CNV event with region_name matching "chr" followed by numbers:
#'    - Creates two new rows with "p" and "q" suffixes
#'    - Copies all other column values to both new rows
#' 2. Otherwise keeps the original row unchanged
#'
#' @examples
#' \dontrun{
#' events <- data.frame(
#'   event_type = c("CNV", "WGD", "CNV"),
#'   region_name = c("chr1", "genome", "chr1p"),
#'   CN_change = c(1, 0, -1),
#'   stringsAsFactors = FALSE
#' )
#'
#' processed <- process_cnv_events(events)
#' # Returns:
#' #   event_type region_name CN_change
#' # 1      CNV       chr1p        1
#' # 2      CNV       chr1q        1
#' # 3      WGD      genome        0
#' # 4      CNV       chr1p       -1
#' }
#'
#' @export
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

#' Create a Genomic Event Record
#'
#' @description
#' Creates a record of a copy number change event
#' Supports CNV and WGD events
#'
#' @param haplotype Character string specifying the haplotype ("maternal" or "paternal")
#' @param parent Character string identifying the parent cell
#' @param child Character string identifying the child cell
#' @param region_name Character string specifying the chromosome region (e.g., "chr1p", "chr2q")
#' @param CN_change Numeric value indicating the copy number change
#' @param copy_index Integer specifying which copy of the segment is affected
#' @param event_type Character string describing the type of event
#'
#' @return A data frame containing a single row with the event information:
#' \itemize{
#'   \item haplotype: Origin of the segment
#'   \item parent: Parent cell identifier
#'   \item child: Child cell identifier
#'   \item region_name: Name of the affected region
#'   \item CN_change: Copy number change value
#'   \item copy_index: Index of the affected copy
#'   \item event_type: Type of genomic event
#' }
#'
#' @examples
#' event <- create_event(
#'   haplotype = "maternal",
#'   parent = "C1",
#'   child = "C2",
#'   region_name = "chr1p",
#'   CN_change = 1,
#'   copy_index = 1,
#'   event_type = "CNV"
#' )
#'
#' @export
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

#' Create Edge Event Table from Event Table and Tree
#'
#' @description
#' This function generates an edge event table by combining information from an event table
#' and a tree structure. It creates a table that summarizes events associated with each edge
#' in the tree, including event labels and the number of events per edge.
#'
#' @param event_table A data frame containing copy number events with columns:
#'        \itemize{
#'          \item parent: Parent node ID
#'          \item child: Child node ID
#'          \item haplotype: Haplotype information ("maternal" or "paternal")
#'          \item region_name: Region identifier
#'          \item CN_change: Copy number change value
#'          \item Additional annotation columns as specified in anno_cols
#'        }
#' @param tree A phylogenetic tree object that can be converted to an edge list
#'        using as_edgelist()
#' @param anno_cols A character vector specifying which columns to use for creating
#'        event labels. Default is c("haplotype_abbr", "region_name", "CN_change")
#'
#' @return A data frame summarizing events for each edge in the tree. The returned data frame
#'   includes the following columns:
#'   \itemize{
#'     \item `edge_name`: A unique identifier for each edge, constructed as `parent_child`.
#'     \item `edge_label`: A concatenated string of event labels for the edge, separated by newline characters.
#'     \item `n_events`: The number of events associated with the edge.
#'   }
#'
#' @details
#' The function performs these steps:
#' 1. Converts tree to edge list format
#' 2. Creates abbreviated haplotype labels (M/P)
#' 3. Combines annotation columns into event names
#' 4. Groups events by edge
#' 5. Creates multi-line labels for edges with multiple events
#' 6. Orders results to match tree structure
#'
#' @note
#' Event labels are formatted as "HAPLOTYPE:REGION:CN_CHANGE" by default.
#' Multiple events on the same edge are separated by newlines.
#'
#' @examples
#' \dontrun{
#' events <- data.frame(
#'   parent = c("A", "A"),
#'   child = c("B", "C"),
#'   haplotype = c("maternal", "paternal"),
#'   region_name = c("chr1p", "chr2q"),
#'   CN_change = c(1, -1)
#' )
#' tree <- your_tree_object  # Replace with actual tree
#' edge_events <- create_edge_event_table(events, tree)
#' }
#'
#' @importFrom dplyr %>% mutate group_by summarize arrange
#' @importFrom rlang sym
#' @importFrom igraph as_edgelist
#' @export
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


