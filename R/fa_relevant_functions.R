# fa_relevant_functions_v1.R

separate_fasta_by_segments <- function(fasta_file_path, output_dir = "segment_fa"){
  # Load the FASTA file
  if (!file.exists(fasta_file_path)) {
    stop("FASTA file does not exist at the specified path: ", fasta_file_path)
  }

  sequences <- readDNAStringSet(fasta_file_path)

  fasta_name_prefix <- tools::file_path_sans_ext(basename(fasta_file_path))

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Iterate through each sequence in the FASTA file
  for (i in 1:length(sequences)) {
    segment_name <- names(sequences[i])
    # Define the output file path
    output_file_path <- file.path(output_dir, paste0(fasta_name_prefix, "_", segment_name, ".fa"))

    # Write the segment to its own FASTA file
    writeXStringSet(sequences[i], filepath = output_file_path)
  }

  cat("Segmentation completed. Files saved in:", output_dir, "\n")

}
