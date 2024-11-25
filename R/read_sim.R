# read_sim_v2.R

simulate_read_bulk <- function(fasta_inputs, output_prefixes, readLen = 100, depth = 30, artPath = "art_illumina", seqSys = "HS25", paired = TRUE, numCores = 2, otherArgs = ""){

  simulateOne <- function(fastaInput, outputPrefix) {
    # Start constructing the basic command
    artCmd <- paste(artPath, "-ss", seqSys, "-i", fastaInput, "-l", readLen, "-f", depth, "-o", outputPrefix)

    # Adjust command based on whether paired-end or single-end reads are desired
    if (paired) {
      # For paired-end reads, add the '-p' flag and specify mean fragment length and its standard deviation
      artCmd <- paste(artCmd, "-p -m 200 -s 10", otherArgs)
    }else{
      artCmd <- paste(artCmd, otherArgs)
    }

    log_file <- paste0(outputPrefix, "_art_output.log")
    artCmd <- paste(artCmd, ">", log_file, "2>&1")
    system(artCmd, intern=FALSE)

    # Return completion message instead of printing directly
    return(paste("Simulation completed for:", outputPrefix))
  }

  # Parallel execution using mclapply
  progressMessages <- mclapply(1:length(fasta_inputs), function(x) simulateOne(fastaInput = fasta_inputs[x], outputPrefix = output_prefixes[x]), mc.cores = numCores)

  # Print all progress messages at once, after all simulations are done
  cat(paste(progressMessages, collapse="\n"), "\n")
}


simulate_read_sc <- function(fasta_input, output_prefixes, readLen=100, depth=30, artPath="art_illumina", seqSys="HS25", paired=TRUE, numCores=2, otherArgs = "") {
  # Define a single simulation function for one cell
  simulateOne <- function(fastaInput, outputPrefix) {
    # Start constructing the basic command
    artCmd <- paste(artPath, "-ss", seqSys, "-i", fastaInput, "-l", readLen, "-f", depth, "-o", outputPrefix)

    # Adjust command based on whether paired-end or single-end reads are desired
    if (paired) {
      # For paired-end reads, add the '-p' flag and specify mean fragment length and its standard deviation
      artCmd <- paste(artCmd, "-p -m 200 -s 10", otherArgs)
    }else{
      artCmd <- paste(artCmd, otherArgs)
    }

    log_file <- paste0(outputPrefix, "_art_output.log")
    artCmd <- paste(artCmd, ">", log_file, "2>&1")
    system(artCmd, intern=FALSE)

    # Return completion message instead of printing directly
    return(paste("Simulation completed for:", outputPrefix))
  }

  # Parallel execution using mclapply
  progressMessages <- mclapply(1:length(output_prefixes), function(x) simulateOne(fastaInput = fasta_input, outputPrefix = output_prefixes[x]), mc.cores = numCores)

  # Print all progress messages at once, after all simulations are done
  cat(paste(progressMessages, collapse="\n"), "\n")
}

merge_fasta_files <- function(paternal_fa, maternal_fa, output_fa, tmp_dir = tempdir()) {
  # Create temporary file names
  tempPaternal <- tempfile(pattern = "temp_paternal", fileext = ".fa", tmpdir = tmp_dir)
  tempMaternal <- tempfile(pattern = "temp_maternal", fileext = ".fa", tmpdir = tmp_dir)

  # Command to add '_paternal' to the headers in the paternal FASTA file
  cmdPaternal <- sprintf("sed '/^>/ s/$/_paternal/' %s > %s", paternal_fa, tempPaternal)

  # Command to add '_maternal' to the headers in the maternal FASTA file
  cmdMaternal <- sprintf("sed '/^>/ s/$/_maternal/' %s > %s", maternal_fa, tempMaternal)

  # Execute the commands
  system(cmdPaternal)
  system(cmdMaternal)

  # Command to merge the temporary files
  cmdMerge <- sprintf("cat %s %s > %s", tempPaternal, tempMaternal, output_fa)

  # Execute the merge command
  system(cmdMerge)

  # Remove the temporary files
  unlink(tempPaternal)
  unlink(tempMaternal)

  cat("Merged file created at:", output_fa, "\n")
}


library(Biostrings)
library(parallel)

simulate_sc_dynamic_reads_for_batches <- function(
    sampled_cell_idx,
    dynamics_ob,
    sc_founder_genomes,
    all_sampled_sc_mutations,
    sc_mut_table_with_nt,
    sim_updates,
    fa_dir,
    output_dir,
    art_path,
    n_cores = 4,
    depth = 30,
    readLen = 150,
    otherArgs = "",
    keep_genome_files = FALSE,
    batch_size = 10
) {
  # List to store sanity check results for each cell
  mutation_checks <- list()

  # Split the sampled cells into batches based on the batch size
  cell_batches <- split(sampled_cell_idx, ceiling(seq_along(sampled_cell_idx) / batch_size))

  # Iterate over each batch
  for (batch in cell_batches) {
    # Lists to store paths to merged genome FASTA files and output prefixes for parallel read simulation
    fasta_inputs <- vector("character", length(batch))
    output_prefixes <- vector("character", length(batch))

    # Step 1: Synthesize genomes and write FASTA files for the current batch
    for (i in seq_along(batch)) {
      cell_index <- batch[i]

      # Extract information about the cell
      sc_info <- dynamics_ob$cell_info[dynamics_ob$cell_info$cell_index == cell_index,]
      clone <- sc_info$clone
      clone_index <- sc_info$clone_index
      founder_cell_index <- sc_founder_genomes$child_node_founders_df[sc_founder_genomes$child_node_founders_df$clone == clone,]$cell_index
      print(founder_cell_index)
      # Synthesize fasta genome for this cell
      tmp_sc_genome <- synth_sc_genome(
        cell_index = cell_index,
        backbone_genome = sc_founder_genomes$all_node_genomes[[clone]],
        sc_mut_list = all_sampled_sc_mutations$sampled_sc_mutations,
        mut_table_with_nt = sc_mut_table_with_nt,
        founder_cell_index = founder_cell_index
      )

      # Perform sanity check on the generated genome
      if(nrow(tmp_sc_genome$cell_mut_with_nt) > 0){
        sc_mut_check <- check_genome_mutations(
          clone_genome = tmp_sc_genome$sc_genome,
          clone_mutation_table = tmp_sc_genome$cell_mut_with_nt,
          clone_seg_list = sim_updates$all_node_segments[[clone]]
        )
      }


      # Store sanity check result
      mutation_checks[[as.character(cell_index)]] <- sc_mut_check

      # Validate the sanity check output
      if (!validate_mutation_check(sc_mut_check)) {
        stop(paste0("Error: Mutation check failed for cell ", cell_index, ". Please investigate."))
      }

      # Initialize an empty DNAStringSet to hold the merged genome for this cell
      merged_genome <- DNAStringSet()

      # Loop through haplotypes (maternal and paternal)
      for (haplotype in c("maternal", "paternal")) {
        # Extract the genome sequence for the current haplotype
        writing_genome <- DNAStringSet(tmp_sc_genome$sc_genome[[haplotype]])

        print(names(writing_genome))
        # Modify the names of the sequences by appending '_maternal' or '_paternal'
        names(writing_genome) <- paste0(names(writing_genome), "_", haplotype)

        # Append the current haplotype sequence to the merged genome for this cell
        merged_genome <- c(merged_genome, writing_genome)

        # Clean up current haplotype's data to free memory
        rm(writing_genome)
      }

      # Define the output file name for the merged haplotypes of this cell
      fa_name <- paste(clone, clone_index, "cell", cell_index, "merged_genome.fa", sep = "_")
      output_fa <- paste0(fa_dir, fa_name)

      # Write the merged genome of this cell to a single FASTA file
      writeXStringSet(merged_genome, output_fa)

      # Output a message indicating that the file has been written
      cat("Merged FASTA file created at:", output_fa, "\n")

      # Add the path to the list of FASTA files and output prefixes for parallel processing later
      fasta_inputs[i] <- output_fa
      output_prefixes[i] <- paste0(output_dir, clone, "_", clone_index, "_cell_", cell_index)

      # Clean up the synthesized genome to save memory
      rm(tmp_sc_genome, merged_genome)
    }

    # Step 2: Simulate Reads in Parallel for the Current Batch
    simulate_read_dynamic_sc(
      fasta_inputs = fasta_inputs,
      output_prefixes = output_prefixes,
      readLen = readLen,
      depth = depth,
      artPath = art_path,
      paired = FALSE,
      numCores = min(n_cores, length(fasta_inputs)),
      otherArgs = otherArgs
    ) # Ensure we don't use more cores than FASTA files

    # Optionally remove genome FASTA files to save disk space
    if (!keep_genome_files) {
      for (fasta in fasta_inputs) {
        if (file.exists(fasta)) {
          file.remove(fasta)
          cat("Genome file removed:", fasta, "\n")
        }
      }
    }
  }

  # Return all mutation check results
  return(mutation_checks)
}

simulate_read_dynamic_sc <- function(fasta_inputs, output_prefixes, readLen=100, depth=30, artPath="art_illumina", seqSys="HS25", paired=TRUE, numCores=2, otherArgs = "") {
  # Define a single simulation function for one cell
  simulateOne <- function(fastaInput, outputPrefix) {
    # Start constructing the basic command
    artCmd <- paste(artPath, "-ss", seqSys, "-i", fastaInput, "-l", readLen, "-f", depth, "-o", outputPrefix)

    # Adjust command based on whether paired-end or single-end reads are desired
    if (paired) {
      # For paired-end reads, add the '-p' flag and specify mean fragment length and its standard deviation
      artCmd <- paste(artCmd, "-p -m 200 -s 10", otherArgs)
    }else{
      artCmd <- paste(artCmd, otherArgs)
    }

    log_file <- paste0(outputPrefix, "_art_output.log")

    # Then, adjust the artCmd to redirect output to this log file
    artCmd <- paste(artCmd, ">", log_file, "2>&1")

    system(artCmd, intern=FALSE)

    # Return completion message instead of printing directly
    return(paste("Simulation completed for:", outputPrefix))
  }

  # Parallel execution using mclapply
  progressMessages <- mclapply(1:length(output_prefixes), function(x) simulateOne(fastaInput = fasta_inputs[x], outputPrefix = output_prefixes[x]), mc.cores = numCores)

  # Print all progress messages at once, after all simulations are done
  cat(paste(progressMessages, collapse="\n"), "\n")
}

