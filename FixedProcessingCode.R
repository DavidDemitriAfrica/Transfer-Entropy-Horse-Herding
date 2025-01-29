## Clean Horse Paper 
# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(RTransferEntropy)
library(igraph)
library(cluster)

# Set working directory
setwd("C:/Users/David Africa/Masters Code Stuff/Horses with JP")

# Function to rename columns for consistent structure
rename_columns <- function(df) {
  col_names <- names(df)
  new_col_names <- col_names
  
  # Always rename the first column to `Time (s)`
  new_col_names[1] <- "Time (s)"
  
  # Loop through the column names and rename based on assumed x/y coordinate pattern
  for (i in seq_along(col_names)) {
    if (i %% 2 == 0 && i > 1) {
      new_col_names[i] <- paste0("Horse", i/2, "_x")
    } else if (i %% 2 == 1 && i > 1) {
      new_col_names[i] <- paste0("Horse", (i+1)/2, "_y")
    }
  }
  names(df) <- new_col_names
  return(df)
}

# Load data and apply renaming
horse_data <- read_csv("HorseMovementData.csv")
horse_data <- rename_columns(horse_data)

# Calculate Transfer Entropy for pairwise horse movements
calculate_transfer_entropy <- function(df) {
  entropy_results <- data.frame()
  
  # Loop through each pair of horses to compute Transfer Entropy
  num_horses <- (ncol(df) - 1) / 2
  for (i in 1:num_horses) {
    for (j in (i + 1):num_horses) {
      x_data <- df[[paste0("Horse", i, "_x")]]
      y_data <- df[[paste0("Horse", j, "_x")]]
      
      te_result <- transfer_entropy(x_data, y_data, lx = 1, ly = 1)
      entropy_results <- rbind(entropy_results, 
                               data.frame(Source = paste0("Horse", i),
                                          Target = paste0("Horse", j),
                                          TransferEntropy = te_result$coef))
    }
  }
  return(entropy_results)
}

# Calculate entropy results and save
entropy_results <- calculate_transfer_entropy(horse_data)
write.csv(entropy_results, "TransferEntropyResults.csv", row.names = FALSE)

# Cluster analysis based on Transfer Entropy
perform_clustering <- function(entropy_df, num_clusters = 3) {
  dist_matrix <- as.dist(1 - as.matrix(entropy_df$TransferEntropy))
  clustering <- hclust(dist_matrix, method = "ward.D2")
  cluster_assignments <- cutree(clustering, k = num_clusters)
  
  return(cluster_assignments)
}

# Perform clustering on transfer entropy results
cluster_results <- perform_clustering(entropy_results)

# Visualization of clustering results
plot_clusters <- function(cluster_assignments, df) {
  df$Cluster <- factor(cluster_assignments)
  ggplot(df, aes(x = Source, y = Target, fill = Cluster)) +
    geom_tile() +
    theme_minimal() +
    labs(title = "Herding Clusters Based on Transfer Entropy",
         x = "Source Horse", y = "Target Horse") +
    scale_fill_brewer(palette = "Set3")
}

# Plot the clustering results
plot_clusters(cluster_results, entropy_results)

