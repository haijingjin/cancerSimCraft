# tree_relevant_functions_v2.R

generate_random_tree <- function(n, include_MRCA = FALSE) {
  # Start with a single node D
  g <- graph.empty()
  g <- add_vertices(g, nv = 1, name = "D")

  # If MRCA is to be included, add C1 and an edge D -> C1
  if (include_MRCA) {
    if (n < 1) {
      stop("n should be at least 1 to include the MRCA (C1).")
    }
    g <- add_vertices(g, nv = 1, name = "C1")
    g <- add_edges(g, c("D", "C1"))
    start_idx <- 2
  } else {
    start_idx <- 1
  }

  # Add new nodes one by one
  for (i in start_idx:n) {
    # Choose a random node to branch from
    if (include_MRCA) {
      parent_node <- V(g)[name != "D"][sample(length(V(g)[name != "D"]), 1)]
    } else {
      parent_node <- V(g)[sample(vcount(g), 1)]
    }

    # Generate the name for the new node
    new_node_name <- paste0("C", i)

    # Add the new node
    g <- add_vertices(g, nv = 1, name = new_node_name)

    # Add an edge from the parent to the new node
    g <- add_edges(g, c(parent_node$name, new_node_name))
  }

  return(g)
}

#' Get Ancestor Nodes in Phylogenetic Tree
#'
#' This function retrieves the ancestors of a specified node (clone) in a tree structure.
#' It calculates the shortest path from the root node to the target node and returns the names of
#' all ancestor nodes (excluding the target node itself).
#'
#' @param tree An igraph tree object representing the phylogenetic structure
#' @param node The target node for which to find ancestors
#'
#' @return A character vector containing the names of the ancestor nodes, ordered from
#'   the root to the immediate parent of the target node.
#'
#' @details
#' The function:
#' 1. Identifies the root node (assumes first vertex is root)
#' 2. Finds shortest path from root to target node
#' 3. Returns names of all nodes in the path
#'
#' @note
#' Assumes the tree is properly rooted with root node as first vertex
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' # Create a simple tree
#' tree <- make_tree(3, 2)
#' V(tree)$name <- c("root", "A", "B")
#'
#' # Get ancestors of node "B"
#' ancestors <- get_clone_ancestors(tree, "B")
#' # Returns: c("root", "B")
#' }
#'
#' @importFrom igraph V shortest_paths
#' @export
get_clone_ancestors <- function(tree, node){
  # Retrieve path from the root node to the target node
  root_node <- V(tree)[1]  # assuming the root is the first node
  path_to_node <- shortest_paths(tree, from = root_node, to = node)$vpath[[1]]

  # Exclude the last node (itself) to get only ancestors
  ancestors <- V(tree)[path_to_node]$name

  return(ancestors)
}



get_clone_descendants <- function(tree, node) {
  # Ensure that 'node' is a vertex in 'tree'
  if (!node %in% V(tree)$name) {
    stop("Node not found in the tree")
  }

  # Find all vertices that are reachable from the given node
  descendants_idx <- unlist(subcomponent(tree, node, mode = "out"))

  # Get the names of these vertices
  descendants <- V(tree)[descendants_idx]$name

  # Remove the node itself from the list of descendants
  descendants <- descendants[descendants != node]

  return(descendants)
}



#' Identify Ancestral Cells in a Lineage
#'
#' @description
#' Traces and returns all ancestral cell indices for a given cell by recursively
#' following the parent relationships in the cell lineage information.
#'
#' @param cell_info Data frame containing cell lineage information with columns:
#'        cell_index, parent, and other attributes
#' @param cell_index Integer specifying the cell index for which to find ancestors
#'
#' @return A numeric vector containing the cell indices of all ancestors of the specified cell,
#'         ordered from immediate parent to earliest ancestor
#'
#' @details
#' This function traverses the cell lineage tree upward from the specified cell,
#' following the parent relationships until it reaches a cell with no parent (NA in the parent column).
#' It builds and returns a vector of all ancestor cell indices encountered during this traversal.
#'
#' The function is commonly used to reconstruct the complete lineage history of a cell,
#' which is particularly useful for aggregating mutations or other heritable properties
#' that are passed down through cell divisions.
#'
#' @seealso
#' \code{\link{get_mutations_sc}}, \code{\link{simulate_sc_dynamics}}
#'
#' @export
get_sc_ancestors <- function(cell_info, cell_index) {
  # Initialize an empty vector to store the ancestors
  ancestor_indices <- c()

  # While the cell has a parent
  while (!is.na(cell_info[cell_info$cell_index == cell_index, "parent"])) {
    # Add the current cell index to the list of ancestors
    parent_index <- cell_info[cell_info$cell_index == cell_index, "parent"]
    ancestor_indices <- c(ancestor_indices, parent_index)

    # Move up to the parent cell
    cell_index <- parent_index
  }

  return(ancestor_indices)
}
