# illustration_functions_v3.R

goodColor <- c("pink", "gold", "darkseagreen", "firebrick4", "lightskyblue", "orange1", "lightyellow", "royalblue", "moccasin","violet", "olivedrab3", "blueviolet", "yellow2", "darkgreen", "royalblue4", "lightsalmon3", "mediumaquamarine", "coral1",  "darkolivegreen1", "cyan4", "tan4", "grey80")
goodColor2 <-  c("black","#DF536B","#61D04F","#2297E6","#28E2E5","#CD0BBC","#F5C710","gray62")
goodColor3 <- c("black", "#CD0BBC", "pink", "gold", "darkseagreen", "firebrick4", "lightskyblue", "orange1", "lightyellow", "royalblue", "moccasin","violet", "olivedrab3", "blueviolet", "yellow2", "darkgreen", "royalblue4", "lightsalmon3", "mediumaquamarine", "coral1",  "darkolivegreen1", "cyan4", "tan4", "grey80")


plot_colors <- function(color_names) {
  n <- length(color_names)
  plot(1:n, rep(1, n), pch = 19, cex = 2, col = color_names, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  text(1:n, rep(1.5, n), labels = color_names, srt = 90, cex = 0.6, adj = c(1, 0.5))
}


assign_colors <- function(df, column, palette) {
  unique_values <- unique(df[, column])
  n <- length(unique_values)
  if (n > length(palette)) {
    warning("Not enough colors in the palette for all unique values; some colors will be reused.")
  }
  # Repeat the palette if there are not enough colors
  color_assignments <- palette[1:n]
  names(color_assignments) <- unique_values
  return(color_assignments)
}


plot_event_tree <- function(edge_event_table, tree) {
  # Generate the tree table
  tree_table <- data.frame(as_edgelist(tree))
  colnames(tree_table) <- c("parent", "child")
  tree_table$edge_name <- paste(tree_table$parent, tree_table$child, sep = "_")

  # Generate the node table
  node_table <- cbind(clone_name = V(tree)$name, clone_type = c("normal", rep("tumor", length(V(tree)) - 1)))

  # Create the graph object
  merge_df_colnames <- c("parent", "child","edge_name", "edge_label", "n_events")
  graph_ob <- graph_from_data_frame(
    merge(tree_table, edge_event_table, by = "edge_name", all.x = TRUE)[,merge_df_colnames],
    node_table,
    directed = TRUE
  )

  tree_plot <- ggraph(graph_ob, layout = 'dendrogram', length = n_events) +
    geom_edge_link(aes(label = edge_label),
                   label_dodge = unit(2.5, 'mm'),
                   color = "grey") +
    geom_node_point(size = 20, aes(colour = clone_type)) +
    geom_node_text(aes(label = name), size = 5) +
    theme_void()

  return(list(graph_ob = graph_ob, tree_plot = tree_plot))
}




archived_draw_chr_bar <- function(window_data){

  chr_ranges <- window_data

  # Standardize column names to lowercase
  colnames(window_data) <- tolower(colnames(window_data))

  chr_lengths <- rle(gsub('chr', '', chr_ranges$chr))$lengths

  chr_n = length(chr_lengths)
  if(chr_n %% 2 == 0){
    times = chr_n/2
    chr_binary <- c(rep(c(2, 1), times))
  }else{
    times = (chr_n - 1)/2
    chr_binary <- c(rep(c(2,1),  times), 2)
  }


  chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
  chr_rl_c <- c(1, cumsum(chr_lengths))

  # creating a data frame to calculate rowMeans
  chr_df <-
    data.frame(
      a = chr_rl_c[1:length(chr_rl_c) - 1],
      b = chr_rl_c[2:length(chr_rl_c)]
    )
  chr_l_means <- round(rowMeans(chr_df))
  chrom.names <- c(1:22, "X", "Y")
  # creating the vector for chr number annotations
  v <- vector(length = sum(chr_lengths), mode = "character") # create vector of the length of windows of ""
  suppressWarnings(v[chr_l_means] <- chrom.names)
  v[is.na(v)] <- ""

  # draw chr bar annotation

  chr_bar <-
    ComplexHeatmap::HeatmapAnnotation(
      chr_text = ComplexHeatmap::anno_text(v[1:nrow(window_data)],
                                           gp = grid::gpar(fontsize = 14)
      ),
      df = as.character(chr[1:nrow(chr), ]),
      show_legend = FALSE,
      show_annotation_name = FALSE,
      which = "column",
      col = list(df = c("1" = "grey88", "2" = "black"))
    )

  return(chr_bar)
}


draw_chr_bar <- function(window_data) {

  # Standardize column names to lowercase
  colnames(window_data) <- tolower(colnames(window_data))

  # Compute chromosome lengths
  chromosome_data <- window_data
  chr_bin_nums <- rle(gsub('chr', '', window_data$chr))$lengths

  # Retrieve exact # of chromosomes, 22 or 23 might be typical for whole genome level data
  n_chr = length(chr_bin_nums)

  # Calc binary patterns for chr bar annotation
  if (n_chr %% 2 == 0) {
    num_pairs = n_chr / 2
    binary_pattern <- c(rep(c(2, 1), num_pairs))
  } else {
    num_pairs = (n_chr - 1) / 2
    binary_pattern <- c(rep(c(2, 1), num_pairs), 2)
  }

  # Create binary pattern data frame, derive exact
  bin_binary_patterns <- rep.int(x = binary_pattern, times = chr_bin_nums)
  bin_indices_at_chr_bounds <- c(1, cumsum(chr_bin_nums))

  # Calculate mean positions for chromosome annotations
  chr_bounds_df <- data.frame(
    start_bin_indices = bin_indices_at_chr_bounds[1:length(bin_indices_at_chr_bounds) - 1],
    end_bin_indices = bin_indices_at_chr_bounds[2:length(bin_indices_at_chr_bounds)]
  )
  chr_segment_mean_pos <- round(rowMeans(chr_bounds_df))

  # Chromosome names for annotation
  chromosome_names <- c(1:22, "X", "Y")[1:n_chr]


  # Prepare vector for chromosome number annotations
  annotation_vector <- vector(length = sum(chr_bin_nums), mode = "character")
  suppressWarnings(annotation_vector[chr_segment_mean_pos] <- chromosome_names)

  # Create chromosome bar annotation for heatmap
  chromosome_bar_annotation <- ComplexHeatmap::HeatmapAnnotation(
    chromosome_text = ComplexHeatmap::anno_text(annotation_vector, gp = grid::gpar(fontsize = 14)),
    binary_pattern = as.character(bin_binary_patterns),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "column",
    col = list(binary_pattern = c("1" = "grey88", "2" = "black"))
  )

  return(chromosome_bar_annotation)
}




plot_truth_heatmap <- function(seg_dat, title, clone_identity_vector, clone_color_anno, chr_bar, integer_col, max_int = 8){
  # set params for complex heatmap
  complex_args <- list(
    use_raster = FALSE,
    raster_quality = 2,
    col = NULL,
    bottom_annotation = NULL,
    right_annotation = NULL,
    left_annotation = rowAnnotation(clones = clone_identity_vector, col = list(clones = clone_color_anno)),
    column_title = "genomic coordinates",
    column_title_gp = grid::gpar(fontsize = 16),
    column_title_side = "bottom",
    row_title = paste0(length(clone_identity_vector), " cells"),
    row_title_gp = grid::gpar(fontsize = 18),
    top_annotation = chr_bar,
    cluster_rows = FALSE,
    border = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_heatmap_legend = TRUE
  )


  if(integer_col == TRUE){

    complex_args <- complex_args[which(names(complex_args) != 'col')]

    seg_dat[seg_dat > max_int] = NA

    color_heat <-
      structure(ocean.balance(length(0:max_int)),
                names = 0:max_int
      )

    heatmap_plot <- do.call(

      ComplexHeatmap::Heatmap,
      c(
        list(
          matrix = seg_dat,
          #   left_annotation = cluster_anno,
          heatmap_legend_param = list(title = "CN"),
          col = color_heat,
          na_col = "yellow"
        ),
        complex_args
      )
    )

    draw(heatmap_plot, column_title = title, column_title_gp = grid::gpar(fontsize = 18))
  }
}

plot_multi_anno_heatmap <- function(seg_dat, title, clone_identity_df, clone_color_anno_list, chr_bar, integer_col, max_int = 8){
  # set params for complex heatmap
  complex_args <- list(
    use_raster = FALSE,
    raster_quality = 2,
    col = NULL,
    bottom_annotation = NULL,
    right_annotation = NULL,
    column_title = "genomic coordinates",
    column_title_gp = grid::gpar(fontsize = 16),
    column_title_side = "bottom",
    row_title = paste0(nrow(clone_identity_df), " cells"),
    row_title_gp = grid::gpar(fontsize = 18),
    top_annotation = chr_bar,
    cluster_rows = FALSE,
    border = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_heatmap_legend = TRUE
  )

  # create annotation bar
  cluster_anno <-
    ComplexHeatmap::rowAnnotation(
      df = clone_identity_df,
      col = clone_color_anno_list,
      show_annotation_name = FALSE
    )

  if(integer_col == TRUE){

    complex_args <- complex_args[which(names(complex_args) != 'col')]

    seg_dat[seg_dat > max_int] = NA

    color_heat <-
      structure(ocean.balance(length(0:max_int)),
                names = 0:max_int
      )
    heatmap_plot <- do.call(

      ComplexHeatmap::Heatmap,
      c(
        list(
          matrix = seg_dat,
          left_annotation = cluster_anno,
          heatmap_legend_param = list(title = "CN"),
          col = color_heat,
          na_col = "yellow"
        ),
        complex_args
      )
    )

    draw(heatmap_plot, column_title = title, column_title_gp = grid::gpar(fontsize = 18))



  }
}

generate_mock_window_data <- function(matrix_colnames) {
  # Extract the chromosome part from the column names
  chr_parts <- sapply(strsplit(matrix_colnames, "_"), `[`, 1)

  # Create a dataframe with a 'chr' column containing these chromosome parts
  window_data <- data.frame(chr = chr_parts)

  return(window_data)
}


