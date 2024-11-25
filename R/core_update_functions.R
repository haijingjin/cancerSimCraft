# core_update_functions_v5.R

update_sim_from_event_table <- function(tree, event_table,
                                        initial_chr_arm_seg_list,
                                        initial_sub_seg_list,
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
