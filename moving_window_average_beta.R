## Moving Window averaging over distance matrices
## Matt Rees
## mathew.rees@ed.ac.uk
## 17/11/2024

# In this script, I want to calculate a moving window average based on 2
# distance matrices. 
# The first one will be a floristic distance matrix
# The second one will be a geographical distance matrix.

# First create an example with 2 matrices of 10 x 10

# Initialize a 10x10 matrix with zeros
distance_matrix <- matrix(0, nrow = 10, ncol = 10)

# Fill the matrix with increasing values by 0.1
for (i in 1:10) {
  for (j in 1:10) {
    if (i != j) {
      distance_matrix[i, j] <- abs(i - j) * 0.1
    }
  }
}

# give the row and col names

rownames(distance_matrix) <- colnames(distance_matrix) <- paste0("n", 1:10)

# Print the distance matrix
distance_matrix

# Create a second matrix of geographical distances

geographical_matrix <- matrix(0, nrow = 10, ncol = 10)

for (i in 1:10) {
  for (j in 1:10) {
    if (i != j) {
      geographical_matrix[i, j] <- abs(i - j) * 10
    }
  }
}

rownames(geographical_matrix) <- colnames(geographical_matrix) <- paste0("n", 1:10)

geographical_matrix

# create an empty object to store the data

moving_avg <- data.frame(matrix(nrow=10, ncol=0))
rownames(moving_avg) <- paste0("n", 1:10)
moving_avg

# Now let's write a function to calculate a moving window average of the 
# floristic distance matrix based on the geographical distance matrix 
# The function can either select a number of adjacent cells, or it can select
# all cells that are closer than a certain threshold in km for example

# n is the number of adjacent cells
# window_size is the minimum distance

moving_window_avg <- function(distance_matrix, geographical_matrix, n, window_size) {
  
  if (!all(dim(distance_matrix) == dim(geographical_matrix))) { 
    stop("Error: The dimensions of distance_matrix and geographical_matrix must be the same.") 
  }
  
  if (!all(colnames(distance_matrix) == colnames(geographical_matrix) )) {
    stop("Error: the names of columns and rows must be the same in both matrices.")
  }
  
  avg_value <- vector()  
  for (i in 1:nrow(distance_matrix)) {
    if(!is.null(n)) {
      x <- geographical_matrix[i,-i] %>% sort()
      y <- x[1:(n)] %>% names()
      z <- mean(distance_matrix[i, y])
      avg_value[i] <- z
    }
    else {
      valid_columns <- which(geographical_matrix[i,-i] <= window_size)
      avg_value[i] <- mean(distance_matrix[i, names(valid_columns)])
    }
  }
  return(avg_value)
  
}

moving_avg$neighbours_beta <- moving_window_avg(distance_matrix = distance_matrix, 
                  geographical_matrix = geographical_matrix, 
                  n = 4, window_size = NULL)

moving_avg$window_beta <- moving_window_avg(distance_matrix = distance_matrix, 
                  geographical_matrix = geographical_matrix, 
                  n = NULL, window_size = 40)

moving_avg
