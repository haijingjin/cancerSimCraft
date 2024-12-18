#' Update Simulation State Based on Event Table
#'
#' @description
#' Updates the simulation state by processing events (WGD, CNV, sub_CNV) along a phylogenetic tree structure.
#' Events are processed in depth-first search order from root to leaves.
#'
#' @param tree An igraph object representing the phylogenetic tree structure
#' @param event_table A data frame containing events with columns:
#'   \itemize{
#'     \item parent - parent node name
#'     \item child - child node name
#'     \item event_type - type of event ("WGD", "CNV", or "sub_CNV")
#'     \item other event-specific columns
#'   }
#' @param initial_chr_arm_seg_list List containing initial chromosome arm segment information
#' @param initial_sub_seg_list List containing initial sub-segment information. Optional,
#'   required only if sub_CNV events are present in event_table
#' @param initial_chr_lengths Vector or list containing initial chromosome lengths
#'
#' @return A list containing:
#'   \itemize{
#'     \item all_node_events - List of events occurring on each tree edge
#'     \item all_node_segments - List of segment states for each node
#'     \item all_node_chr_lengths - List of chromosome lengths for each node
#'   }
#'
#' @details
#' The function traverses the tree in depth-first search order, processing events sequentially.
#' For each edge in the tree, it applies the corresponding events (WGD, CNV, sub_CNV)
#' and updates the chromosome segments and lengths accordingly.
#'
#' @seealso
#' \code{\link{update_wgd_seg}}, \code{\link{update_cnv_seg}}, \code{\link{update_sub_seg}}
#'
#' @importFrom igraph V degree dfs get.adjlist
#' @export
update_sim_from_event_table <- function(tree, event_table,
                                        initial_chr_arm_seg_list,
                                        initial_sub_seg_list = NULL,
                                        initial_chr_lengths) {

  root_name <- names(V(tree))[degree(tree, mode = "in") == 0]
  dfs_order <- dfs(tree, root = root_name, mode = "out")$order

  all_node_events <- list()

  all_node_segments <- list()
  all_node_segments[[root_name]] <- initial_chr_arm_seg_list

  all_node_chr_lengths <- list()
  all_node_chr_lengths[[root_name]] <- initial_chr_lengths

  for (i in 2:length(dfs_order)) {
    child_node_index <- dfs_order[i]

    parent_node_index <- get.adjlist(tree, mode = "in")[[child_node_index]]
    parent_node <- V(tree)[parent_node_index]$name
    child_node <- V(tree)[child_node_index]$name

    tmp_segments <- all_node_segments[[parent_node]]
    tmp_chr_lengths <- all_node_chr_lengths[[parent_node]]
    edge <- paste0(parent_node, sep = "_", child_node)

    relevant_events <- subset(event_table, parent == parent_node & child == child_node)
    for (j in 1:nrow(relevant_events)) {
      event_type <- relevant_events[j, "event_type"]
      if (event_type == "WGD") {
        wgd_event <- relevant_events[j,]
        wgd_updates <- update_wgd_seg(seg_list = tmp_segments,
                                      chr_lengths = tmp_chr_lengths,
                                      edge_event_table = wgd_event)
        tmp_segments <- wgd_updates$updated_seg_list
        tmp_chr_lengths <- wgd_updates$updated_chr_lengths
        all_node_events[[edge]] <- rbind(all_node_events[[edge]], wgd_event)
      } else if (event_type == "CNV") {
        cnv_event <- relevant_events[j,]
        cnv_updates <- update_cnv_seg(seg_list = tmp_segments,
                                      chr_lengths = tmp_chr_lengths,
                                      edge_event_table = cnv_event)
        tmp_segments <- cnv_updates$updated_seg_list
        tmp_chr_lengths <- cnv_updates$updated_chr_lengths
        all_node_events[[edge]] <- rbind(all_node_events[[edge]], cnv_event)
      } else if (event_type == "sub_CNV") {
        sub_event <- relevant_events[j,]
        sub_updates <- update_sub_seg(seg_list = tmp_segments,
                                      sub_seg_list = initial_sub_seg_list,
                                      chr_lengths = tmp_chr_lengths,
                                      one_sub_event = sub_event)
        tmp_segments <- sub_updates$updated_seg_list
        tmp_chr_lengths <- sub_updates$updated_chr_lengths
        all_node_events[[edge]] <- rbind(all_node_events[[edge]], sub_event)
      }
    }
    all_node_segments[[child_node]] <- tmp_segments
    all_node_chr_lengths[[child_node]] <- tmp_chr_lengths
  }

  return(list("all_node_events" = all_node_events,
              "all_node_segments" = all_node_segments,
              "all_node_chr_lengths" = all_node_chr_lengths))
}
