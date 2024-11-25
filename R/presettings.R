# presettings_v3.R
# Set genome-relevant variables
chr_indices <- c(1:22, 'X')
chr_names <- paste0("chr", chr_indices)

goodColor <- c("pink", "gold", "darkseagreen", "firebrick4", "lightskyblue", "orange1", "lightyellow", "royalblue", "moccasin","violet", "olivedrab3", "blueviolet", "yellow2", "darkgreen", "royalblue4", "lightsalmon3", "mediumaquamarine", "coral1",  "darkolivegreen1", "cyan4", "tan4", "grey80")

nt_transition_matrix <- matrix(c(
  0, 0.1, 0.7, 0.2, 0, # Probabilities for A -> C, A -> G, A -> T
  0.1, 0, 0.2, 0.7, 0, # Probabilities for C -> A, C -> G, C -> T
  0.7, 0.2, 0, 0.1, 0, # Probabilities for G -> A, G -> C, G -> T
  0.2, 0.7, 0.1, 0, 0, # Probabilities for T -> A, T -> C, T -> G
  0, 0, 0, 0, 1

), nrow = 5, byrow = TRUE)
rownames(nt_transition_matrix) <- colnames(nt_transition_matrix) <- c("A", "C", "G", "T", "N")

