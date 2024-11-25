# params_setting_v4.R

generate_JC_matrix <- function() {
  # Initialize matrix with 1/3 for all off-diagonal elements
  JC_matrix <- matrix(1/3, nrow = 5, ncol = 5)
  diag(JC_matrix) <- 0  # Diagonal elements are 0 for A, C, G, T
  JC_matrix[5, ] <- 0   # Transitions to 'N' are not allowed
  JC_matrix[, 5] <- 0   # Transitions from 'N' are not allowed
  JC_matrix[5, 5] <- 1  # 'N' transitions to itself with probability 1

  rownames(JC_matrix) <- colnames(JC_matrix) <- c("A", "C", "G", "T", "N")

  return(JC_matrix)
}


generate_K2P_matrix <- function(alpha, beta) {
  # Initialize a 5x5 matrix with specified transition and transversion rates
  K2P_matrix <- matrix(beta/(2*beta + alpha), nrow = 5, ncol = 5)

  # Set transition rates for A, C, G, T
  K2P_matrix[1, 3] <- alpha/(2*beta + alpha)  # A->G
  K2P_matrix[3, 1] <- alpha/(2*beta + alpha)  # G->A
  K2P_matrix[2, 4] <- alpha/(2*beta + alpha)  # C->T
  K2P_matrix[4, 2] <- alpha/(2*beta + alpha)  # T->C

  # Set diagonal to 0 for A, C, G, T
  diag(K2P_matrix) <- 0

  # Add the N state transitions
  K2P_matrix[5, ] <- c(0, 0, 0, 0, 1)  # N can only transition to N
  K2P_matrix[, 5] <- c(0, 0, 0, 0, 1)  # No other state can transition to N

  rownames(K2P_matrix) <- colnames(K2P_matrix) <- c("A", "C", "G", "T", "N")

  return(K2P_matrix)
}

