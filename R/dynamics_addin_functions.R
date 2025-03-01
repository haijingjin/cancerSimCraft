# dynamics_addin_functions.R

#' Prepare Cell Lineage Data for Tree Visualization
#'
#' @description
#' Transforms cell lineage information into a format suitable for phylogenetic tree visualization
#' by creating edge and node tables with appropriate attributes.
#'
#' @param cell_info Data frame containing cell lineage information with columns:
#'        clone, parent, birth_time, death_time, cell_index
#'
#' @return A list containing two data frames:
#' \itemize{
#'   \item edges: Data frame with columns parent, child, edge_label, edge_length, and child_clone
#'   \item nodes: Data frame with columns name, clone_type, and death_time
#' }
#'
#' @details
#' This function processes cell lineage data to create a tree structure representation:
#'
#' 1. Replaces NA parents with "root_node" to establish a common ancestor
#' 2. Assigns theoretical death times to living cells (cells with NA death_time)
#' 3. Calculates edge lengths as the time between birth and death of each cell
#' 4. Creates an edge table that defines connections between cells
#' 5. Creates a node table with information about each cell and the root
#'
#' The resulting data structure is compatible with tree visualization packages
#' such as ggtree and can be used to create cell lineage visualizations.
#'
#' @examples
#' # Example usage:
#' cell_info <- data.frame(
#'   parent = c(NA, 1, 1, 2),
#'   cell_index = c(1, 2, 3, 4),
#'   birth_time = c(0, 1, 1, 2),
#'   death_time = c(5, 3, NA, 4),
#'   clone = c("Clone1", "Clone1", "Clone2", "Clone2")
#' )
#'
#' tree_data <- prepare_tree_data(cell_info)
#' print(tree_data$edges)  # Edge table
#' print(tree_data$nodes)  # Node table
#'
#' @seealso
#' \code{\link{simulate_sc_dynamics}}, \code{\link{plot_cell_tree}}
#'
#' @export
prepare_tree_data <- function(cell_info) {

  # Change NA parents to "root_node"
  cell_info$parent[is.na(cell_info$parent)] <- "root_node"

  # Calculate a theoretical death time for living cells
  max_death_time <- max(cell_info$death_time, na.rm = TRUE)
  cell_info$death_time[is.na(cell_info$death_time)] <- max_death_time + 1

  # Calculate edge lengths as the difference between death_time and birth_time
  edge_table <- data.frame(
    parent = cell_info$parent,
    child = cell_info$cell_index,
    edge_label = paste(cell_info$parent, cell_info$cell_index, sep = "->"),
    edge_length = cell_info$death_time - cell_info$birth_time,
    child_clone = cell_info$clone
  )

  # Create node table
  node_table <- data.frame(
    name = c("root_node", cell_info$cell_index),
    clone_type = c("root", cell_info$clone),  # root_node doesn't belong to any clone
    death_time = c(0, cell_info$death_time)
  )

  return(list(edges = edge_table, nodes = node_table))
}

# Function to create and save tree plot from simulation results
create_sc_tree_plot_archived <- function(dynamic_sim_ob, clone_colors, title = "SC Dynamics Dendrogram") {
  # Extract all cell info from the simulation result
  all_cell_info <- dynamic_sim_ob$cell_info

  # Prepare tree data (assuming `prepare_tree_data()` is a valid function in your environment)
  tree_data <- prepare_tree_data(all_cell_info)

  # Create the graph object using the tree data
  graph_ob <- graph_from_data_frame(d = tree_data$edges, vertices = tree_data$nodes, directed = TRUE)

  # Create the ggraph plot with manual colors
  tree_plot <- ggraph(graph_ob, 'dendrogram', length = edge_length) +
    geom_edge_elbow(aes(color = child_clone)) +
    scale_edge_color_manual(values = clone_colors) +
    theme_minimal() +
    labs(title = title) +  # Use the title parameter here
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("t") +
    xlab(NULL) +  # Removes the x-axis label
    theme(axis.text.x = element_blank(),  # Removes the x-axis tick labels
          axis.ticks.x = element_blank())  # Removes the x-axis tick marks

  return(tree_plot)  # Return the plot object (if you want to display it)
}

#' Create a Cell Lineage Tree Visualization
#'
#' @description
#' Creates a dendrogram visualization of cell lineage data from simulation results,
#' showing the relationships between cells and their clonal evolution over time.
#'
#' @param dynamic_sim_ob List object returned by simulate_sc_dynamics function
#' @param clone_colors Named vector of colors where names match clone types in the simulation
#' @param title Character string for the plot title, defaults to "SC Dynamics Dendrogram"
#'
#' @return A ggraph plot object displaying the cell lineage tree with customized styling
#'
#' @details
#' This function creates a dendrogram visualization of cell lineage data where:
#'
#' 1. The vertical axis represents time
#' 2. Edges represent cells' lifespans
#' 3. Edge colors indicate the clone type of each cell
#' 4. Branch lengths correspond to cell lifetimes
#'
#' The function first prepares the cell information data using the prepare_tree_data function,
#' then creates a graph object and finally generates a dendrogram visualization using ggraph.
#' The resulting plot has publication-ready styling with customized title, legend, and axis labels.
#'
#' @examples
#' # Define colors for each clone type
#' clone_colors <- c(A = "blue", B = "red")
#'
#' # Create and display the lineage tree
#' tree_plot <- create_sc_tree_plot(dynamic_sim_ob, clone_colors, "Cell Lineage Evolution")
#' print(tree_plot)
#'
#' # Save the plot
#' ggsave("cell_lineage_tree.png", tree_plot, width = 10, height = 8)
#'
#' @seealso
#' \code{\link{simulate_sc_dynamics}}, \code{\link{prepare_tree_data}}
#'
#' @import ggraph
#' @import igraph
#' @import ggplot2
#'
#' @export
create_sc_tree_plot <- function(dynamic_sim_ob, clone_colors, title = "SC Dynamics Dendrogram") {
  # Extract all cell info from the simulation result
  all_cell_info <- dynamic_sim_ob$cell_info

  # Prepare tree data (assuming `prepare_tree_data()` is a valid function in your environment)
  tree_data <- prepare_tree_data(all_cell_info)

  # Create the graph object using the tree data
  graph_ob <- graph_from_data_frame(d = tree_data$edges, vertices = tree_data$nodes, directed = TRUE)

  # Create the ggraph plot with manual colors and publication-ready text
  tree_plot <- ggraph(graph_ob, 'dendrogram', length = edge_length) +
    geom_edge_elbow(aes(color = child_clone)) +
    scale_edge_color_manual(values = clone_colors) +
    theme_minimal() +
    labs(title = title, y = "t") +  # Use the title parameter here, set y-axis label
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),    # Centered, larger, bold title
      axis.title.y = element_text(size = 14, face = "bold"),                # Bold y-axis label
      axis.text.y = element_text(size = 12),                                # Larger y-axis tick labels
      axis.text.x = element_blank(),                                        # Remove x-axis text
      axis.ticks.x = element_blank(),                                       # Remove x-axis tick marks
      legend.title = element_text(size = 14, face = "bold"),                # Bold legend title
      legend.text = element_text(size = 12)                                 # Larger legend text
    )


  return(tree_plot)  # Return the plot object (if you want to display it)
}


extract_first_ancestor <- function(cell_info) {
  first_ancestor <- cell_info[cell_info$clone_index == 1, ]
  return(first_ancestor)
}


# Define a function to extract cell indices for each clone, including all live cells
get_all_clone_cell_indices <- function(dynamics_ob) {
  # Filter out dead cells (only keep cells where death_time is NA)
  live_cells <- dynamics_ob$cell_info[is.na(dynamics_ob$cell_info$death_time), ]

  # Extract all unique clone names from the filtered live cells
  clone_names <- unique(live_cells$clone)

  # Create a list to store indices for each clone
  all_clone_cell_indices <- lapply(clone_names, function(clone_name) {
    # Filter the cell indices for the given clone from the live cells
    clone_indices <- live_cells[live_cells$clone == clone_name, ]$cell_index
    # Return the cell indices for this clone
    return(clone_indices)
  })

  # Assign the clone names to the list for easy access
  names(all_clone_cell_indices) <- clone_names

  # Add a panel that contains all live cells
  all_clone_cell_indices[["all_live"]] <- live_cells$cell_index

  # Return the list of cell indices for all clones
  return(all_clone_cell_indices)
}


# Define a function to count the number of times each mutation occurs across all cells
count_mutation_occurrences <- function(mutation_list) {
  # Convert each data frame in the list to a unified format
  combined_mutations <- bind_rows(mutation_list, .id = "cell_id")

  # Count the number of times each mutation occurs
  mutation_counts <- combined_mutations %>%
    group_by(clone, haplotype, chrom, pos) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count))

  # Return the mutation counts data frame
  return(mutation_counts)
}

calculate_clone_mut_summary <- function(all_clone_mut_counts, all_clone_cell_indices) {
  # Initialize an empty data frame to store the mutation and cell counts for each clone
  clone_summary <- data.frame(clone_name = character(), num_mutations = integer(), num_cells = integer(), mutation_ratio = numeric(), shared_mutations = integer(), unique_mutations = integer(), shared_mutation_ratio = numeric(), unique_mutation_ratio = numeric(), stringsAsFactors = FALSE)

  # Extract the "all_live" mutation data from the list
  all_live_mut_data <- all_clone_mut_counts[["all_live"]]

  # Loop through each clone's data in the lists (excluding "all_live")
  for (clone_name in names(all_clone_mut_counts)) {
    if (clone_name == "all_live") {
      next  # Skip the "all_live" panel as we only need specific clone data
    }

    # Get the mutation tibble for the clone
    mut_data <- all_clone_mut_counts[[clone_name]]

    # Get the cell indices for the clone
    cell_indices <- all_clone_cell_indices[[clone_name]]

    # Count the number of cells
    num_cells <- length(cell_indices)

    # Initialize counters for shared and unique mutations
    shared_count <- 0
    unique_count <- 0

    # Loop through each mutation in the clone's mutation data
    for (i in 1:nrow(mut_data)) {
      # Extract mutation information
      clone <- mut_data$clone[i]
      haplotype <- mut_data$haplotype[i]
      chrom <- mut_data$chrom[i]
      pos <- mut_data$pos[i]
      count <- mut_data$count[i]

      # Find the matching mutation in the "all_live" panel using data frame subsetting
      matching_row <- all_live_mut_data[all_live_mut_data$clone == clone &
                                          all_live_mut_data$haplotype == haplotype &
                                          all_live_mut_data$chrom == chrom &
                                          all_live_mut_data$pos == pos, ]

      # Check if the matching mutation was found
      if (nrow(matching_row) == 0) {
        stop(paste("Error: Mutation not found in 'all_live' for clone", clone, ", haplotype", haplotype, ", chromosome", chrom, ", position", pos))
      } else {
        # Extract the count from the matching row
        live_cell_count <- matching_row$count

        # Check if the mutation is unique or shared
        if (count == live_cell_count) {
          unique_count <- unique_count + 1
        } else if (count < live_cell_count) {
          shared_count <- shared_count + 1
        }
      }
    }

    # Calculate the mutation-to-cell ratio
    num_mutations <- nrow(mut_data)
    mutation_ratio <- ifelse(num_cells > 0, num_mutations / num_cells, NA)

    # Calculate shared and unique mutation ratios
    shared_mutation_ratio <- ifelse(num_cells > 0, shared_count / num_cells, NA)
    unique_mutation_ratio <- ifelse(num_cells > 0, unique_count / num_cells, NA)

    # Add the clone's mutation and cell counts to the summary table
    clone_summary <- rbind(clone_summary, data.frame(clone_name = clone_name, num_mutations = num_mutations, num_cells = num_cells, mutation_ratio = mutation_ratio, shared_mutations = shared_count, unique_mutations = unique_count, shared_mutation_ratio = shared_mutation_ratio, unique_mutation_ratio = unique_mutation_ratio))
  }

  return(clone_summary)
}


# Define a function to perform pairwise comparisons of shared mutations across all clones, excluding specific clones
pairwise_shared_mutations <- function(clone_mut_counts, exclude_clones = c()) {
  # Extract all clone names (e.g., MRCA, V, IV, II, etc.)
  clone_names <- names(clone_mut_counts)

  # Remove the clones that need to be excluded
  clone_names <- clone_names[!clone_names %in% exclude_clones]

  # Create an empty list to store pairwise shared mutations
  shared_mutation_list <- list()

  # Loop through each combination of clones to perform pairwise comparisons
  for (i in 1:(length(clone_names) - 1)) {
    for (j in (i + 1):length(clone_names)) {
      # Extract clone names
      clone_name_1 <- clone_names[i]
      clone_name_2 <- clone_names[j]

      # Extract the mutation data frames for each clone
      clone_1_counts <- clone_mut_counts[[clone_name_1]]
      clone_2_counts <- clone_mut_counts[[clone_name_2]]

      # Perform the inner join to find shared mutations
      shared_mutations <- inner_join(
        clone_1_counts,
        clone_2_counts,
        by = c("clone", "haplotype", "chrom", "pos"),
        suffix = c(paste0("_", clone_name_1), paste0("_", clone_name_2))
      )

      # Store the result in the list with a meaningful name
      pair_name <- paste(clone_name_1, clone_name_2, sep = "_vs_")
      shared_mutation_list[[pair_name]] <- shared_mutations
    }
  }

  # Return the list of shared mutations for all pairwise comparisons
  return(shared_mutation_list)
}


# Define a function to count the number of times each mutation occurs across all cells
count_mutation_occurrences <- function(mutation_list) {
  # Convert each data frame in the list to a unified format
  combined_mutations <- bind_rows(mutation_list, .id = "cell_id")

  # Count the number of times each mutation occurs
  # Define unique mutations based on the columns of interest (e.g., `clone`, `haplotype`, `chrom`, `pos`)
  mutation_counts <- combined_mutations %>%
    group_by(clone, haplotype, chrom, pos) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count))

  # Return the mutation counts data frame
  return(mutation_counts)
}
