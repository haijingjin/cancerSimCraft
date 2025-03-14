
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

#' Simulate Clonal-level Single-Cell Sequencing Reads
#'
#' This function uses the ART read simulator (Huang et al., Bioinformatics 2012) to generate synthetic sequencing
#' reads from a reference FASTA file for multiple single cells in parallel. It supports
#' both paired-end and single-end reads and allows for parallel execution across multiple cores.
#' The function generates FASTQ files for each cell and logs the ART output.
#'
#' @param fasta_input Character string. Path to the input FASTA file.
#' @param output_prefixes Character vector. Prefixes for output files, one per simulated cell.
#' @param readLen Numeric. Length of the simulated reads in base pairs. Default is 100.
#' @param depth Numeric. Sequencing coverage depth. Default is 1.
#' @param artPath Character string. Path to the ART Illumina executable. Default is "art_illumina".
#' @param seqSys Character string. Sequencing system to simulate. Default is "HS25" (HiSeq 2500).
#' @param paired Logical. Whether to generate paired-end reads (TRUE) or single-end reads (FALSE). Default is TRUE.
#' @param numCores Numeric. Number of cores to use for parallel processing. Default is 2.
#' @param otherArgs Character string. Additional arguments to pass to ART Illumina. Default is "".
#'
#' @details
#' This function creates synthetic sequencing reads by calling ART Illumina for each output prefix
#' in parallel using mclapply. For paired-end reads, it sets a mean fragment length of 200bp with
#' a standard deviation of 10bp. Simulation output and errors are redirected to log files.
#'
#' @return No return value, called for side effects of generating simulated read files and
#' printing completion messages.
#'
#' @examples
#' \dontrun{
#' # Simulate reads for 3 cells with default parameters
#' simulate_read_sc(
#'   fasta_input = "reference.fa",
#'   output_prefixes = c("cell1_", "cell2_", "cell3_")
#' )
#'
#' # Simulate single-end reads with custom parameters
#' simulate_read_sc(
#'   fasta_input = "reference.fa",
#'   output_prefixes = c("cell1_", "cell2_"),
#'   readLen = 150,
#'   depth = 0.01,
#'   paired = FALSE,
#'   ...
#' )
#' }
#'
#' @importFrom parallel mclapply
#' @export
simulate_read_sc <- function(fasta_input, output_prefixes, readLen=100, depth=1, artPath="art_illumina", seqSys="HS25", paired=TRUE, numCores=2, otherArgs = "") {
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

#' Merge Paternal and Maternal FASTA Files
#'
#' This function merges two FASTA files (paternal and maternal) into a single output file,
#' adding suffix identifiers to differentiate the sequences' origins.
#'
#'
#' @param paternal_fa Character string. Path to the paternal FASTA file.
#' @param maternal_fa Character string. Path to the maternal FASTA file.
#' @param output_fa Character string. Path where the merged FASTA file will be written.
#' @param tmp_dir Character string. Directory for temporary files. Default is tempdir().
#'
#' @details
#' The function adds "_paternal" suffix to sequence headers in the paternal file
#' and "_maternal" suffix to sequence headers in the maternal file before
#' concatenating them into a single output file.
#'
#' @return No return value, called for side effect of creating the merged FASTA file.
#'
#' @examples
#' \dontrun{
#' merge_fasta_files(
#'   paternal_fa = "path/to/paternal.fa",
#'   maternal_fa = "path/to/maternal.fa",
#'   output_fa = "path/to/merged.fa"
#' )
#' }
#'
#' @export
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


#' Simulate Sequencing Reads for Single Cells in Batches
#'
#' @description
#' Simulates sequencing reads for multiple single cells in batches by synthesizing
#' cell-specific genomes with mutations and generating sequencing reads using the
#' ART simulator. This function processes cells in batches to manage memory usage.
#'
#' @param sampled_cell_idx Vector of cell indices for which to simulate reads
#' @param dynamics_ob List object returned by simulate_sc_dynamics containing
#'        cell lineage information
#' @param sc_founder_genomes List containing genome information for clone founders
#'        with elements: all_node_genomes, child_node_founders_df
#' @param all_sampled_sc_mutations List containing mutation information for sampled cells
#' @param sc_mut_table_with_nt Data frame containing mutation details with nucleotide changes
#' @param sim_updates List containing updated simulation information including all_node_segments
#' @param fa_dir Character string specifying the directory to store FASTA files
#' @param output_dir Character string specifying the directory for output read files
#' @param art_path Character string specifying the path to the ART sequencing simulator executable
#' @param n_cores Integer specifying the number of CPU cores to use for parallel processing (default: 4)
#' @param depth Numeric value specifying the desired sequencing depth/coverage (default: 30)
#' @param readLen Integer specifying the length of simulated reads in base pairs (default: 150)
#' @param otherArgs Character string with additional arguments to pass to ART (default: "")
#' @param keep_genome_files Logical indicating whether to retain the temporary FASTA
#'        genome files after simulation (default: FALSE)
#' @param batch_size Integer specifying the number of cells to process in each batch (default: 10)
#'
#' @return A list where each element contains mutation check results for a cell, with the
#'         cell index as the name of each element
#'
#' @details
#' This function performs the following steps for each batch of cells:
#'
#' 1. For each cell in the batch:
#'    - Extracts cell lineage information from dynamics_ob
#'    - Synthesizes a cell-specific genome with mutations using synth_sc_genome()
#'    - Performs checks on the generated genome using check_genome_mutations()
#'    - Validates the sanity check results with validate_mutation_check()
#'    - Combines maternal and paternal haplotypes into a single FASTA file
#'    - Writes the merged genome to disk
#'
#' 2. Simulates sequencing reads for all cells in the batch in parallel using simulate_read_dynamic_sc()
#'
#' 3. Optionally removes temporary genome files to save disk space
#'
#' Processing cells in batches helps manage memory usage when dealing with large numbers of cells.
#' The function relies on external functions like synth_sc_genome(), check_genome_mutations(),
#' validate_mutation_check(), and simulate_read_dynamic_sc().
#'
#' @examples
#' \dontrun{
#' # Assuming you have already run a simulation and have necessary objects
#' mutation_checks <- simulate_sc_dynamic_reads_for_batches(
#'   sampled_cell_idx = c(50, 51, 52, 53, 54),
#'   dynamics_ob = sc_dynamics_result,
#'   sc_founder_genomes = founder_genomes,
#'   all_sampled_sc_mutations = sampled_mutations,
#'   sc_mut_table_with_nt = mutations_with_nt,
#'   sim_updates = simulation_updates,
#'   fa_dir = "data/genomes/",
#'   output_dir = "results/reads/",
#'   art_path = "/usr/local/bin/art_illumina",
#'   n_cores = 4,
#'   depth = 30,
#'   readLen = 150,
#'   keep_genome_files = FALSE,
#'   batch_size = 2
#' )
#'
#' # Check results for a specific cell
#' print(mutation_checks[["50"]])
#' }
#'
#' @seealso
#' \code{\link{synth_sc_genome}}, \code{\link{check_genome_mutations}},
#' \code{\link{validate_mutation_check}}, \code{\link{simulate_read_dynamic_sc}}
#'
#' @importFrom Biostrings DNAStringSet writeXStringSet
#'
#' @export
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

#' Simulate Single-cell Level (dynamics incorporated) Sequencing Reads Using ART in Parallel
#'
#' @description
#' Simulates next-generation sequencing reads from input FASTA files in parallel
#' using the ART sequencing simulator (Huang et al., Bioinformatics 2012). This function
#' supports both single-end and paired-end read generation.
#'
#' @param fasta_inputs Character vector of paths to input FASTA files
#' @param output_prefixes Character vector of output prefixes for generated read files
#' @param readLen Integer specifying the length of simulated reads in base pairs (default: 100)
#' @param depth Numeric value specifying the desired sequencing depth/coverage (default: 30)
#' @param artPath Character string specifying the path to the ART executable (default: "art_illumina")
#' @param seqSys Character string specifying the sequencing system to simulate (default: "HS25" for HiSeq 2500)
#' @param paired Logical indicating whether to generate paired-end reads (default: TRUE)
#' @param numCores Integer specifying the number of CPU cores to use for parallel processing (default: 2)
#' @param otherArgs Character string with additional arguments to pass to ART (default: "")
#'
#' @return None. The function generates sequencing read files at the locations specified
#'         by output_prefixes and prints completion messages.
#'
#' @details
#' This function provides a parallel interface to the ART sequencing simulator by:
#'
#' 1. Defining an internal function `simulateOne()` that constructs and executes
#'    the ART command for a single FASTA input
#' 2. Using `mclapply()` to run multiple simulations in parallel
#' 3. Redirecting ART's output to log files
#'
#' For paired-end reads (when paired=TRUE), the function sets a mean fragment length
#' of 200bp with a standard deviation of 10bp. For single-end reads, these parameters
#' are omitted.
#'
#' The function requires that the ART executable is installed and accessible.
#'
#' @examples
#' \dontrun{
#' # Simulate single-end reads for 3 genomes
#' simulate_read_dynamic_sc(
#'   fasta_inputs = c("data/genome1.fa", "data/genome2.fa", "data/genome3.fa"),
#'   output_prefixes = c("results/sim1", "results/sim2", "results/sim3"),
#'   readLen = 150,
#'   depth = 30,
#'   artPath = "/usr/local/bin/art_illumina",
#'   seqSys = "HS25",
#'   paired = FALSE,
#'   numCores = 3
#' )
#'
#' # Simulate paired-end reads with custom ART arguments
#' simulate_read_dynamic_sc(
#'   fasta_inputs = c("data/genome1.fa", "data/genome2.fa"),
#'   output_prefixes = c("results/paired1", "results/paired2"),
#'   paired = TRUE,
#'   otherArgs = "--noALN --rndSeed 42"
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_sc_dynamic_reads_for_batches}}
#'
#' @importFrom parallel mclapply
#'
#' @export
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

