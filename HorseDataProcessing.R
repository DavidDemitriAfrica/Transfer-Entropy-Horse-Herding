library('tidyverse')
library('ggplot2')
library("readxl")
library("gganimate")
library("lmtest")
library("qpcR")
library(RTransferEntropy)
library(lsa)
library(igraph)
library(dplyr)
library(ggpubr)
library("rinform")
library(reshape2)
library("patchwork")

setwd("C:/Users/David Africa/Masters Code Stuff/Horses with JP")

rename_columns <- function(df) {
  col_names <- names(df)
  new_col_names <- col_names
  
  # Always rename the first column to `Time (s)`
  new_col_names[1] <- "Time (s)"
  
  # Loop through the column names and apply renaming logic
  for (i in seq_along(col_names)) {
    if (i %% 2 == 0 && i > 1) { # For even columns (assumed to be x-coordinates)
      base_name <- gsub(" ", "_", tolower(col_names[i]))
      new_col_names[i] <- paste0(base_name, "_x")
    } else if (i %% 2 == 1 && i > 2) { # For odd columns (assumed to be y-coordinates)
      base_name <- gsub(" ", "_", tolower(col_names[i-1]))
      new_col_names[i] <- paste0(base_name, "_y")
    }
  }
  
  # Set the new column names
  names(df) <- new_col_names
  
  return(df)
}

folder_path <- "data/herding_only"
grazing_folder_path <- "data/filtered_fullsampled/filtered_fullsampled"
# List all Excel files in the directory
excel_files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
grazing_excel_files <- list.files(path = grazing_folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Read each Excel file into a list of data frames
list_of_dfs <- lapply(excel_files, read_excel)
grazing_list_of_dfs <- lapply(grazing_excel_files, read_excel)

# Optionally, if you want to name the list elements based on the file names (without path and extension)
names(list_of_dfs) <- sapply(excel_files, function(x) tools::file_path_sans_ext(basename(x)))
names(grazing_list_of_dfs) <- sapply(grazing_excel_files, function(x) tools::file_path_sans_ext(basename(x)))
list_of_dfs_renamed <- lapply(list_of_dfs, rename_columns)
grazing_list_of_dfs_renamed <- lapply(grazing_list_of_dfs, rename_columns)

bin_and_summarize <- function(df) {
  df %>%
    mutate(discretized_time = cut(`Time (s)`, breaks = seq(floor(min(`Time (s)`)), ceiling(max(`Time (s)`)), by = 0.25), include.lowest = TRUE, labels = FALSE)) %>%
    group_by(discretized_time) %>%
    summarize(across(everything(), mean, na.rm = TRUE), .groups = "drop")
}

list_of_dfs_binned <- lapply(list_of_dfs_renamed, bin_and_summarize)
grazing_list_of_dfs_binned <- lapply(grazing_list_of_dfs_renamed, bin_and_summarize)
grazing_list_of_dfs_binned
plot_movement <- function(df, entities) {
  p <- ggplot() +
    labs(title = "Horses' Movements",
         x = "X-coordinate",
         y = "Y-coordinate") +
    theme_minimal()
  
  for (entity in entities) {
    x_col <- paste0(entity, "_x")
    y_col <- paste0(entity, "_y")
    p <- p + geom_point(data = df, aes_string(x = x_col, y = y_col, color = sprintf('"%s"', entity)), size = 2)
  }
  
  return(p)
}

# Specify the entities you want to plot
entities <- c("kobe", "tarumi", "himeji", "maiko", "kakogawa", "akashi", "uji")

plot_movement(list_of_dfs_binned[[4]], entities)

calculate_velocity <- function(df, dt = 0.25) {
  velocity_df <- data.frame(matrix(NA, ncol = 0, nrow = nrow(df)))
  
  for (i in seq(1, ncol(df), by = 2)) {
    x_col <- i
    y_col <- i + 1
    
    dx <- c(0, diff(df[[x_col]]))
    dy <- c(0, diff(df[[y_col]]))
    
    velocity_df[paste0("velocity_", names(df)[x_col])] <- dx / dt
    velocity_df[paste0("velocity_", names(df)[y_col])] <- dy / dt
  }
  
  return(velocity_df)
}

list_of_dfs_velocity <- lapply(list_of_dfs_binned, calculate_velocity)
grazing_list_of_dfs_velocity <-lapply(grazing_list_of_dfs_binned, calculate_velocity)
compute_speed <- function(velocity_df) {
  # Initialize an empty list to store speed calculations
  speed_list <- list()
  
  # Iterate over each pair of velocity columns (x and y) to calculate speed
  for (i in seq(1, ncol(velocity_df), by = 2)) {
    # Extract the column names for x and y velocities
    x_vel_col_name <- names(velocity_df)[i]
    y_vel_col_name <- names(velocity_df)[i + 1]
    
    # Calculate speed from x and y velocities
    speed <- sqrt(velocity_df[[x_vel_col_name]]^2 + velocity_df[[y_vel_col_name]]^2)
    
    # Generate a new column name for speed based on the entity's name
    # Assumes velocity column names are in the format "velocity_ENTITY_x" and "velocity_ENTITY_y"
    entity_name <- gsub("velocity_(.*)_x", "\\1", x_vel_col_name)
    speed_col_name <- paste0(entity_name, "_speed")
    
    # Add the speed calculation to the list with the appropriate column name
    speed_list[[speed_col_name]] <- speed
  }
  
  # Convert the list of speed calculations into a dataframe
  speed_df <- as.data.frame(speed_list)
  
  return(speed_df)
}

list_of_dfs_speed <- lapply(list_of_dfs_velocity, compute_speed)
grazing_list_of_dfs_speed <- lapply(grazing_list_of_dfs_velocity, compute_speed)
compute_transfer_entropy_matrix <- function(speed_df) {
  n <- ncol(speed_df)
  matrixinp <- matrix(data = NA, nrow = n, ncol = n) 
  colnames(matrixinp) <- colnames(speed_df) 
  rownames(matrixinp) <- colnames(speed_df)
  
  for (i in 1:n) { 
    for (j in 1:n) { 
      matrixinp[i, j] <- RTransferEntropy::transfer_entropy(speed_df[[i]], speed_df[[j]], lx = 2, ly = 2, nboot = 0)$coef[1]
    } 
  } 
  
  return(round(matrixinp, 2))
}


list_of_temats <- lapply(head(list_of_dfs_speed,-1), compute_transfer_entropy_matrix)

list_of_temats
calculate_transfer_entropy <- function(df) {
  # Columns to exclude
  exclude_columns <- c("velocity_discretized_time_speed", "bato_1_speed", "bato_2_speed")
  
  # Filter columns
  filtered_df <- df[, !(names(df) %in% exclude_columns)]
  
  # Initialize list to store transfer entropy results
  te_results <- list()
  
  # Calculate transfer entropy for each pair of columns
  for (col_x in names(filtered_df)) {
    for (col_y in names(filtered_df)) {
      if (col_x != col_y) {
        # Use the transfer_entropy function from rinform
        te_value <- rinform::transfer_entropy(as.integer(filtered_df[[col_x]]), as.integer(filtered_df[[col_y]]),k = 2, local = TRUE)
        # Store results in a named list
        te_results[[paste(col_x, "to", col_y)]] <- c(te_value, tail(filtered_df[[col_x]],1), tail(filtered_df[[col_y]],1))
      }
    }
  }
  return(te_results)
}
list_te_results <- lapply(list_of_dfs_speed, calculate_transfer_entropy)

length(list_te_results)
count_results = data.frame()
for (i in 1:length(list_te_results)){
  for (j in 1:length(list_te_results[[i]])){
    threshold = 0.05
    result = sum(head(list_te_results[[i]][[j]],-2) > threshold)
    count_results = rbind(count_results, c(result, tail(list_te_results[[i]][[j]],1), head(tail(list_te_results[[i]][[j]],2),1)))
  }
}


to_long <- list_te_results$`filtered_fullsampled_inverness_0011-1`

filtered_lists <- to_long[grep("^tarumi_speed to .*_speed$", names(to_long))]

# Create a dataframe
df <- do.call(rbind, lapply(names(filtered_lists), function(name) {
  data.frame(
    timestep = seq_along(filtered_lists[[name]]) - 1,
    value = filtered_lists[[name]],
    variable = sub("tarumi_speed to (.*)_speed", "\\1", name)  # Extract X from "tarumi_speed to X_speed"
  )
}))

# Convert to long format
base_size <- 12
df_long <- df %>%
  mutate(variable = factor(variable, levels = unique(variable))) %>%
  gather(key = "key", value = "value", -timestep, -variable)
df_long
plot <- ggplot(df_long, aes(x = timestep / 2, y = value)) +
  geom_line(color = "lightblue", size = 1.5) +
  labs(
    x = "Time (s)",
    y = "TE") +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = base_size) +  # Use a minimal theme
  theme(
    text = element_text(size = base_size),  # Increase text size
    plot.title = element_text(size = base_size + 2, face = "bold"),  # Title text size and bold
    axis.title = element_text(size = base_size + 1),  # Axis title size
    strip.text = element_text(size = base_size + 1),  # Facet label size
    panel.border = element_rect(color = "black", fill = NA),
    panel.spacing = unit(1.5, "lines")
  )
plot

#### Heatmap / N/A Counts ####

# Function to process a single DataFrame
process_df <- function(df) {
  # Get all horse columns in the DataFrame, skipping the first three non-horse columns
  cols <- colnames(df)[-(1:3)]
  n_cols <- length(cols)
  
  # Compute average speed over the first eight timesteps for each horse
  avg_speeds <- sapply(cols, function(col) {
    mean(df[[col]][1:8], na.rm = TRUE)
  })
  
  # Rank the horses based on average speed (descending order)
  # Assign ordinal positions: 1 = fastest, 2 = second fastest, etc.
  ranking <- order(avg_speeds, decreasing = TRUE)
  ordinal_positions <- setNames(seq_along(ranking), cols[ranking])
  
  # Initialize an empty data frame to store results
  results <- data.frame(
    leader = character(),
    follower = character(),
    leader_ordinal = integer(),
    follower_ordinal = integer(),
    relationship_found = logical(),
    na_count = integer(),
    total_computations = integer(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each pair of columns
  for (i in 1:n_cols) {
    for (j in 1:n_cols) {
      if (i != j) {
        col_x <- df[[cols[i]]]
        col_y <- df[[cols[j]]]
        
        # Get ordinal positions
        leader_ordinal <- ordinal_positions[cols[i]]
        follower_ordinal <- ordinal_positions[cols[j]]
        
        # Initialize counters
        na_count <- 0
        total_computations <- 0
        relationship_found <- FALSE
        
        # Iterate over possible lags to find the lag with maximum TE
        max_te <- 0
        optimal_lag <- NA
        for (lag in 1:10) {  # Adjust lag range as needed
          total_computations <- total_computations + 1
          te_result <- tryCatch({
            RTransferEntropy::transfer_entropy(col_x, col_y, lx = lag, ly = lag, entropy = "Shannon")
          }, error = function(e) NA)
          
          te_xy <- tryCatch({te_result$coef[1]}, error = function(e) NA)
          
          if (is.na(te_xy)) {
            na_count <- na_count + 1
          } else {
            # Check if te_xy exceeds the threshold
            if (te_xy > 0.1 && te_xy > max_te) {
              max_te <- te_xy
              optimal_lag <- lag
              relationship_found <- TRUE
            }
          }
        }
        
        # Append the result for this pair
        results <- rbind(results, data.frame(
          leader = cols[i],
          follower = cols[j],
          leader_ordinal = leader_ordinal,
          follower_ordinal = follower_ordinal,
          relationship_found = relationship_found,
          na_count = na_count,
          total_computations = total_computations,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  return(results)
}

# Apply the function across all DataFrames in the list
# For both regular and grazing periods (assuming you have both lists)
results_list_speed <- lapply(list_of_dfs_speed, process_df)
results_list_grazing <- lapply(grazing_list_of_dfs_speed, process_df)  # Assuming you have this list



# Save the list object to a file
saveRDS(results_list_speed, file = "results_list_speed.rds")
saveRDS(results_list_grazing, file = "results_list_grazing.rds")

# To load it back
# results_list_grazing <- readRDS("results_list_grazing.rds")

# Combine results for speed data
results_combined_speed <- do.call(rbind, results_list_speed)

# Combine results for grazing data
results_combined_grazing <- do.call(rbind, results_list_grazing)

# Function to summarize results and create matrices for heatmaps
create_heatmap_data <- function(results_combined) {
  # Summarize results by horse names
  results_summary_names <- results_combined %>%
    group_by(leader, follower) %>%
    summarise(
      total_relationships = sum(relationship_found),
      total_na = sum(na_count),
      total_computations = sum(total_computations),
      .groups = "drop"  # Add this line
    )
  
  # Compute proportions
  results_summary_names <- results_summary_names %>%
    mutate(
      proportion_relationships = total_relationships / (total_computations - total_na),
      proportion_na = total_na / total_computations
    )
  
  # Create proportion matrix for horse names
  na_proportion_matrix_names <- results_summary_names %>%
    dplyr::select(leader, follower, proportion_na) %>%
    pivot_wider(names_from = follower, values_from = proportion_na, values_fill = NA)
  
  # Summarize results by ordinal positions
  results_summary_ordinals <- results_combined %>%
    group_by(leader_ordinal, follower_ordinal) %>%
    summarise(
      total_relationships = sum(relationship_found),
      total_na = sum(na_count),
      total_computations = sum(total_computations),
      .groups = "drop"  # Add this line
    )
  
  # Compute proportions
  results_summary_ordinals <- results_summary_ordinals %>%
    mutate(
      proportion_relationships = total_relationships / (total_computations - total_na),
      proportion_na = total_na / total_computations
    )
  
  # Create proportion matrix for ordinal positions
  na_proportion_matrix_ordinals <- results_summary_ordinals %>%
    dplyr::select(leader_ordinal, follower_ordinal, proportion_na) %>%
    pivot_wider(names_from = follower_ordinal, values_from = proportion_na, values_fill = NA)
  
  return(list(
    names_matrix = na_proportion_matrix_names,
    ordinals_matrix = na_proportion_matrix_ordinals
  ))
}

# Create heatmap data for speed data
heatmap_data_speed <- create_heatmap_data(results_combined_speed)

# Create heatmap data for grazing data
heatmap_data_grazing <- create_heatmap_data(results_combined_grazing)

plot_heatmap <- function(matrix_data, title) {
  # Convert the matrix to long format for ggplot
  matrix_long <- matrix_data %>%
    pivot_longer(-leader, names_to = "follower", values_to = "proportion_na")
  
  # Plot heatmap with consistent color scale
  ggplot(matrix_long, aes(x = follower, y = leader, fill = proportion_na)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "blue", high = "red",
      na.value = "grey50", limits = c(0, 1)
    ) +
    labs(
      title = title, x = "Follower", y = "Leader",
      fill = "N/A Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1, size = 20
      ),
      axis.text.y = element_text(size = 20),
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 25
      ),
      plot.title.position = "plot",
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-10, 0, 0, 0)
    )
}

# Generate the heatmaps
heatmap_speed <- plot_heatmap(
  heatmap_data_speed$names_matrix,
  "Herding Leadership Heatmap"
)
heatmap_grazing <- plot_heatmap(
  heatmap_data_grazing$names_matrix,
  "Grazing Leadership Heatmap"
)

# Combine heatmaps into one horizontal layout with a shared legend at the bottom
combined_heatmap <- (heatmap_speed + heatmap_grazing) +
  plot_layout(
    ncol = 2, widths = c(1, 1), guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 0, 0, 0)
  )

# Display the combined heatmap
print(combined_heatmap)

# Plot heatmaps for ordinal positions
plot_heatmap_ordinals <- function(matrix_data, title) {
  # Convert leader_ordinal and follower_ordinal to factors
  matrix_data$leader_ordinal <- factor(matrix_data$leader_ordinal)
  matrix_long <- matrix_data %>%
    pivot_longer(-leader_ordinal, names_to = "follower_ordinal", values_to = "proportion_na")
  matrix_long$follower_ordinal <- factor(matrix_long$follower_ordinal)
  
  # Plot heatmap
  ggplot(matrix_long, aes(x = follower_ordinal, y = leader_ordinal, fill = proportion_na)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "blue", high = "red", na.value = "grey50") +
    labs(title = title, x = "Follower Ordinal", y = "Leader Ordinal", fill = "N/A Proportion") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1, size = 20
      ),
      axis.text.y = element_text(size = 20),
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 25
      ),
      plot.title.position = "plot",
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-10, 0, 0, 0)
    )
}

ordinalheatmap_speed <- plot_heatmap_ordinals(heatmap_data_speed$ordinals_matrix, "Herding - Ordinal Positions")
ordinalheatmap_grazing <- plot_heatmap_ordinals(heatmap_data_grazing$ordinals_matrix, "Grazing - Ordinal Positions")

combined_heatmap <- (ordinalheatmap_speed + ordinalheatmap_grazing) +
  plot_layout(
    ncol = 2, widths = c(1, 1), guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 0, 0, 0)
  )

print(combined_heatmap)

# Output overall N/A counts and relationship counts
overall_summary <- function(results_combined) {
  total_pairs <- nrow(results_combined)
  total_relationships <- sum(results_combined$relationship_found)
  total_na <- sum(results_combined$na_count)
  total_computations <- sum(results_combined$total_computations)
  
  cat("Total pairs analyzed:", total_pairs, "\n")
  cat("Total relationships found:", total_relationships, "\n")
  cat("Total N/A computations:", total_na, "\n")
  cat("Total computations attempted:", total_computations, "\n")
  cat("Overall N/A proportion:", total_na / total_computations, "\n")
}

cat("Speed Data Summary:\n")
overall_summary(results_combined_speed)

cat("\nGrazing Data Summary:\n")
overall_summary(results_combined_grazing)

##### Average Time Lag #####

# Basic parameters for TE estimation
n_lags <- 10  # number of time lags to test
max_te <- 0
optimal_lag <- 0

# Loop over possible time lags
for (lag in 1:n_lags) {
  # Calculate TE from X to Y with a given time lag
  te_result <- transfer_entropy(data_X, data_Y, lx = lag, ly = 1, entropy = "Shannon")
  
  # Extract the TE value from X to Y
  te_xy <- te_result$tx2y
  
  # Check if this is the maximum TE observed so far
  if (te_xy > max_te) {
    max_te <- te_xy
    optimal_lag <- lag
  }
}

# Output the optimal time lag
cat("Optimal time lag between X and Y:", optimal_lag, "\n")
cat("Maximum Transfer Entropy:", max_te, "\n")

average_lag_by_event <- function(df) {
  # Get all pairs of columns in the DataFrame
  cols <- colnames(df)[-(1:3)]
  n_cols <- length(cols)
  te_lags <- c()  # Vector to store lags with TE > 0.1
  
  # Iterate over each pair of columns
  for (i in 1:(n_cols - 1)) {
    for (j in (i + 1):n_cols) {
      col_x <- df[[cols[i]]]
      col_y <- df[[cols[j]]]
      
      # Iterate over possible lags to find the lag with maximum TE
      max_te <- 0
      optimal_lag <- NA
      for (lag in 1:10) {  # Adjust lag range as needed
        te_result <- tryCatch({
          RTransferEntropy::transfer_entropy(col_x, col_y, lx = lag, ly = lag, entropy = "Shannon")
        }, error = function(e) NA)
        
        te_xy <- tryCatch({te_result$coef[1]}, error = function(e) NA)
        
        # Check if te_xy is not NA and exceeds the threshold
        if (!is.na(te_xy) && te_xy > 0.1 && te_xy > max_te) {
          max_te <- te_xy
          optimal_lag <- lag
        }
      }
      
      # Store the lag if TE > 0.1 was found
      if (!is.na(optimal_lag)) {
        te_lags <- c(te_lags, optimal_lag)
      }
    }
  }
  
  # Calculate the average lag for the horse if lags are found
  if (length(te_lags) > 0) {
    return(mean(te_lags))
  } else {
    return(NA)  # Return NA if no pairs met the TE > 0.1 threshold
  }
}
# Apply the function across all DataFrames and calculate the total average lag
average_lags_by_event <- sapply(list_of_dfs_speed, average_lag_by_event)
average_lag_total <- mean(average_lags_by_event, na.rm = TRUE)
average_lags_by_event
average_lag_total

# Output results
cat("Average lag per horse:\n")
print(average_lags_by_horse)
cat("\nTotal average lag across all horses:", average_lag_total, "\n")

#### Aggregated by Horse Across Datasets ####

average_lag_by_horse <- function(df) {
  # Get all horse columns in the DataFrame, skipping the first three non-horse columns
  cols <- colnames(df)[-(1:3)]
  n_cols <- length(cols)
  
  # Initialize lists to store lags for each horse as leader and follower
  leader_lags <- setNames(vector("list", n_cols), cols)
  follower_lags <- setNames(vector("list", n_cols), cols)
  
  # Iterate over each pair of columns
  for (i in 1:(n_cols - 1)) {
    for (j in (i + 1):n_cols) {
      col_x <- df[[cols[i]]]
      col_y <- df[[cols[j]]]
      
      # Find the lag with the maximum transfer entropy from col_x to col_y
      max_te_xy <- 0
      optimal_lag_xy <- NA
      for (lag in 1:10) {
        te_result_xy <- tryCatch({
          RTransferEntropy::transfer_entropy(col_x, col_y, lx = lag, ly = lag, entropy = "Shannon")
        }, error = function(e) NA)
        
        te_xy <- tryCatch({te_result_xy$coef[1]}, error = function(e) NA)
        
        if (!is.na(te_xy) && te_xy > 0.1 && te_xy > max_te_xy) {
          max_te_xy <- te_xy
          optimal_lag_xy <- lag
        }
      }
      
      # Store the lag where col_x is the leader and col_y is the follower
      if (!is.na(optimal_lag_xy)) {
        leader_lags[[cols[i]]] <- c(leader_lags[[cols[i]]], optimal_lag_xy)
        follower_lags[[cols[j]]] <- c(follower_lags[[cols[j]]], optimal_lag_xy)
      }
      
      # Now, find the lag with the maximum transfer entropy from col_y to col_x
      max_te_yx <- 0
      optimal_lag_yx <- NA
      for (lag in 1:10) {
        te_result_yx <- tryCatch({
          RTransferEntropy::transfer_entropy(col_y, col_x, lx = lag, ly = lag, entropy = "Shannon")
        }, error = function(e) NA)
        
        te_yx <- tryCatch({te_result_yx$coef[1]}, error = function(e) NA)
        
        if (!is.na(te_yx) && te_yx > 0.1 && te_yx > max_te_yx) {
          max_te_yx <- te_yx
          optimal_lag_yx <- lag
        }
      }
      
      # Store the lag where col_y is the leader and col_x is the follower
      if (!is.na(optimal_lag_yx)) {
        leader_lags[[cols[j]]] <- c(leader_lags[[cols[j]]], optimal_lag_yx)
        follower_lags[[cols[i]]] <- c(follower_lags[[cols[i]]], optimal_lag_yx)
      }
    }
  }
  
  # Calculate the average lag per horse as leader and follower
  average_leader_lag <- sapply(leader_lags, function(lags) if (length(lags) > 0) mean(lags) else NA)
  average_follower_lag <- sapply(follower_lags, function(lags) if (length(lags) > 0) mean(lags) else NA)
  
  return(list(leader = average_leader_lag, follower = average_follower_lag))
}

# Apply the function across all DataFrames in the list
average_lags_per_horse <- lapply(list_of_dfs_speed, average_lag_by_horse)

# Combine results and average across all DataFrames for each horse
combine_role_lags <- function(role) {
  role_lags <- lapply(average_lags_per_horse, `[[`, role)
  combined_lags <- do.call(rbind, role_lags)
  colMeans(combined_lags, na.rm = TRUE)
}

average_leader_lag_by_horse <- combine_role_lags("leader")
average_follower_lag_by_horse <- combine_role_lags("follower")

# Output the results
cat("Average lag per horse as Leader across all DataFrames:\n")
print(average_leader_lag_by_horse)
cat("\nAverage lag per horse as Follower across all DataFrames:\n")
print(average_follower_lag_by_horse)


#### Aggregated by First Mover Across Datasets ####

# Assuming you have already computed 'average_lags_per_horse' for each DataFrame
# 'average_lags_per_horse' is a list where each element corresponds to a DataFrame
# Each element is a list with 'leader' and 'follower' lags per horse

# Initialize lists to collect lags per ordinal position across all DataFrames
leader_lags_by_ordinal <- list()
follower_lags_by_ordinal <- list()

# Iterate over each DataFrame and its corresponding lags
for (idx in seq_along(list_of_dfs_speed)) {
  df <- list_of_dfs_speed[[idx]]
  lags <- average_lags_per_horse[[idx]]
  
  # Get all horse columns, skipping the first three non-horse columns
  cols <- colnames(df)[-(1:3)]
  
  # Compute average speed over the first eight timesteps for each horse
  avg_speeds <- sapply(cols, function(col) {
    mean(df[[col]][1:8], na.rm = TRUE)
  })
  
  # Rank the horses based on average speed (descending order)
  # Assign ordinal positions: 1 = fastest, 2 = second fastest, etc.
  ranking <- order(avg_speeds, decreasing = TRUE)
  ordinal_positions <- setNames(seq_along(ranking), cols[ranking])
  
  # Map lags from horse names to ordinal positions
  leader_lags <- lags$leader
  follower_lags <- lags$follower
  
  # Map leader lags to ordinal positions
  for (horse_name in names(leader_lags)) {
    ordinal_pos <- as.character(ordinal_positions[horse_name])
    lag_value <- leader_lags[horse_name]
    if (!is.na(lag_value)) {
      leader_lags_by_ordinal[[ordinal_pos]] <- c(leader_lags_by_ordinal[[ordinal_pos]], lag_value)
    }
  }
  
  # Map follower lags to ordinal positions
  for (horse_name in names(follower_lags)) {
    ordinal_pos <- as.character(ordinal_positions[horse_name])
    lag_value <- follower_lags[horse_name]
    if (!is.na(lag_value)) {
      follower_lags_by_ordinal[[ordinal_pos]] <- c(follower_lags_by_ordinal[[ordinal_pos]], lag_value)
    }
  }
}

# Identify all ordinal positions across DataFrames
all_ordinals <- unique(c(names(leader_lags_by_ordinal), names(follower_lags_by_ordinal)))

# Calculate the average lag per ordinal position as Leader
average_leader_lag_by_ordinal <- sapply(all_ordinals, function(pos) {
  lags <- leader_lags_by_ordinal[[pos]]
  if (!is.null(lags) && length(lags) > 0) {
    mean(unlist(lags), na.rm = TRUE)
  } else {
    NA
  }
})

# Calculate the average lag per ordinal position as Follower
average_follower_lag_by_ordinal <- sapply(all_ordinals, function(pos) {
  lags <- follower_lags_by_ordinal[[pos]]
  if (!is.null(lags) && length(lags) > 0) {
    mean(unlist(lags), na.rm = TRUE)
  } else {
    NA
  }
})

# Sort the results by ordinal position
average_leader_lag_by_ordinal <- average_leader_lag_by_ordinal[order(as.numeric(names(average_leader_lag_by_ordinal)))]
average_follower_lag_by_ordinal <- average_follower_lag_by_ordinal[order(as.numeric(names(average_follower_lag_by_ordinal)))]

# Output the results
cat("Average lag per ordinal position as Leader across all DataFrames:\n")
print(average_leader_lag_by_ordinal)

cat("\nAverage lag per ordinal position as Follower across all DataFrames:\n")
print(average_follower_lag_by_ordinal)


##### Duration #####

# Load necessary library
library(rinform)

# Initialize storage for durations
total_durations <- c()
durations_per_event <- list()
durations_per_horse <- list()
durations_per_ordinal <- list()

# Iterate over each event (DataFrame in list_of_dfs_speed)
for (event_idx in 1:length(list_of_dfs_speed)) {
  df <- list_of_dfs_speed[[event_idx]]
  
  # Get all horse columns, skipping the first three non-horse columns
  cols <- colnames(df)[-(1:3)]
  
  # Compute average speed over the first eight timesteps for each horse
  avg_speeds <- sapply(cols, function(col) {
    mean(df[[col]][1:8], na.rm = TRUE)
  })
  
  # Rank the horses based on average speed (descending order)
  # Assign ordinal positions: 1 = fastest, 2 = second fastest, etc.
  ranking <- order(avg_speeds, decreasing = TRUE)
  ordinal_positions <- setNames(seq_along(ranking), cols[ranking])
  
  # Initialize durations for this event
  event_durations <- c()
  
  # Iterate over each pair of horses (x, y)
  for (x_idx in 1:length(cols)) {
    for (y_idx in 1:length(cols)) {
      if (x_idx != y_idx) {
        x_col <- cols[x_idx]
        y_col <- cols[y_idx]
        
        # Get time series for horses x and y
        ts_x <- df[[x_col]]
        ts_y <- df[[y_col]]
        
        # Remove any NA values and ensure the time series are of equal length
        min_length <- min(length(ts_x), length(ts_y))
        ts_x <- ts_x[1:min_length]
        ts_y <- ts_y[1:min_length]
        
        # Remove any NA values from both time series
        valid_indices <- which(!is.na(ts_x) & !is.na(ts_y))
        ts_x <- ts_x[valid_indices]
        ts_y <- ts_y[valid_indices]
        
        # Ensure that the time series are non-empty after removing NAs
        if (length(ts_x) == 0 || length(ts_y) == 0) {
          next  # Skip this pair if data is insufficient
        }
        
        # Ensure that the time series are integer sequences for rinform
        if (!is.integer(ts_x)) {
          ts_x <- as.integer(ts_x)
        }
        if (!is.integer(ts_y)) {
          ts_y <- as.integer(ts_y)
        }
        
        # Compute local transfer entropy from x to y
        te_result <- tryCatch({
          rinform::transfer_entropy(ts_x, ts_y, k = 1, local = TRUE)
        }, error = function(e) NA)
        
        if (sum(is.na(te_result)) > 0) {
          next  # Skip if there's an error in computation
        }
        
        # Extract local TE values and identify where TE > 0.01
        local_te <- as.vector(te_result)
        is_high_te <- local_te > 0.01
        
        # Identify contiguous sequences where local TE > 0.01
        rle_te <- rle(is_high_te)
        
        # Extract durations (number of contiguous time steps with high TE)
        durations <- rle_te$lengths[rle_te$values == TRUE]
        
        # Convert durations to seconds (since timestep = 0.25s)
        durations_sec <- durations * 0.25  # Multiply by 0.25 to convert to seconds
        
        # Store total durations
        total_durations <- c(total_durations, durations_sec)
        
        # Store durations for the current event
        event_durations <- c(event_durations, durations_sec)
        
        # Store durations per horse
        # x_col is the leader, y_col is the follower
        # Initialize lists if they don't exist
        if (is.null(durations_per_horse[[x_col]]$leader)) {
          durations_per_horse[[x_col]]$leader <- c()
        }
        if (is.null(durations_per_horse[[y_col]]$follower)) {
          durations_per_horse[[y_col]]$follower <- c()
        }
        durations_per_horse[[x_col]]$leader <- c(durations_per_horse[[x_col]]$leader, durations_sec)
        durations_per_horse[[y_col]]$follower <- c(durations_per_horse[[y_col]]$follower, durations_sec)
        
        # Store durations per ordinal position
        ordinal_x <- as.character(ordinal_positions[x_col])
        ordinal_y <- as.character(ordinal_positions[y_col])
        
        # Initialize lists if they don't exist
        if (is.null(durations_per_ordinal[[ordinal_x]]$leader)) {
          durations_per_ordinal[[ordinal_x]]$leader <- c()
        }
        if (is.null(durations_per_ordinal[[ordinal_y]]$follower)) {
          durations_per_ordinal[[ordinal_y]]$follower <- c()
        }
        durations_per_ordinal[[ordinal_x]]$leader <- c(durations_per_ordinal[[ordinal_x]]$leader, durations_sec)
        durations_per_ordinal[[ordinal_y]]$follower <- c(durations_per_ordinal[[ordinal_y]]$follower, durations_sec)
      }
    }
  }
  
  # Store the durations for the current event
  durations_per_event[[event_idx]] <- event_durations
}
total_durations
# Create and save histogram of leader-follower durations
png(filename="leader_follower_durations_histogram.png", width=900, height=600)
par(mar=c(5, 6, 4, 2)) # Bottom, left, top, right margins

# Setting large font sizes for the plot
hist(total_durations, breaks=30, 
     main="Histogram of Leader-Follower Durations", 
     xlab="Duration (seconds)", ylab="Frequency", 
     col="skyblue", border="white", 
     cex.main=2.5, # Font size for title
     cex.lab=2, # Font size for axis labels
     cex.axis=2) # Font size for axis tick labels

dev.off()


# Create a data frame for the durations
duration_df <- data.frame(Duration = total_durations)

# Create and save the enhanced histogram
ggplot(duration_df, aes(x = Duration)) +
  geom_histogram(
    binwidth = 0.5,
    fill = "#69b3a2",
    color = "white",
    alpha = 0.9,
    boundary = 0
  ) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Distribution of Leader-Follower Relationship Durations",
    x = "Duration (seconds)",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(0, max(duration_df$Duration, na.rm = TRUE), by = 1)) +
  ggsave("enhanced_leader_follower_durations_histogram.png", width = 12, height = 8, dpi = 300)

  # Calculate the average total duration
average_total_duration <- mean(total_durations, na.rm = TRUE)

# Calculate the average duration per event
average_duration_per_event <- sapply(durations_per_event, function(d) {
  if (length(d) > 0) mean(d, na.rm = TRUE) else NA
})

# Calculate the average duration per horse (as leader and follower)
horses <- names(durations_per_horse)

average_duration_per_horse_leader <- sapply(horses, function(horse) {
  durations <- durations_per_horse[[horse]]$leader
  if (length(durations) > 0) mean(durations, na.rm = TRUE) else NA
})

average_duration_per_horse_follower <- sapply(horses, function(horse) {
  durations <- durations_per_horse[[horse]]$follower
  if (length(durations) > 0) mean(durations, na.rm = TRUE) else NA
})

# Calculate the average duration per ordinal position (as leader and follower)
ordinals <- sort(unique(names(durations_per_ordinal)))

average_duration_per_ordinal_leader <- sapply(ordinals, function(ord) {
  durations <- durations_per_ordinal[[ord]]$leader
  if (length(durations) > 0) mean(durations, na.rm = TRUE) else NA
})

average_duration_per_ordinal_follower <- sapply(ordinals, function(ord) {
  durations <- durations_per_ordinal[[ord]]$follower
  if (length(durations) > 0) mean(durations, na.rm = TRUE) else NA
})

# Sort the average durations per ordinal position by ordinal number
average_duration_per_ordinal_leader <- average_duration_per_ordinal_leader[order(as.numeric(names(average_duration_per_ordinal_leader)))]
average_duration_per_ordinal_follower <- average_duration_per_ordinal_follower[order(as.numeric(names(average_duration_per_ordinal_follower)))]

# Output the results
cat("Average Total Duration (seconds):\n")
print(average_total_duration)

cat("\nAverage Duration per Event (seconds):\n")
print(average_duration_per_event)

cat("\nAverage Duration per Horse as Leader (seconds):\n")
print(average_duration_per_horse_leader)

cat("\nAverage Duration per Horse as Follower (seconds):\n")
print(average_duration_per_horse_follower)

cat("\nAverage Duration per Ordinal Position as Leader (seconds):\n")
print(average_duration_per_ordinal_leader)

cat("\nAverage Duration per Ordinal Position as Follower (seconds):\n")
print(average_duration_per_ordinal_follower)


#### Duration Burst ####

library(rinform)
library(reshape2)
library(ggplot2)
library(dplyr)

# Select a single herding event
event_idx <- 3  # Change this to select a different event
df <- list_of_dfs_speed[[event_idx]]

# Get all horse columns, skipping the first three non-horse columns
cols <- colnames(df)[-(1:3)]
n_horses <- length(cols)
n_timepoints <- nrow(df)

# Initialize matrices to store leader and follower statuses
leader_matrix <- matrix(FALSE, nrow = n_horses, ncol = n_timepoints)
follower_matrix <- matrix(FALSE, nrow = n_horses, ncol = n_timepoints)
rownames(leader_matrix) <- cols
rownames(follower_matrix) <- cols

# Iterate over each pair of horses (x, y)
for (x_idx in 1:length(cols)) {
  for (y_idx in 1:length(cols)) {
    if (x_idx != y_idx) {
      x_col <- cols[x_idx]
      y_col <- cols[y_idx]
      
      # Get time series for horses x and y
      ts_x <- df[[x_col]]
      ts_y <- df[[y_col]]
      
      # Ensure the time series are of equal length
      min_length <- min(length(ts_x), length(ts_y))
      ts_x <- ts_x[1:min_length]
      ts_y <- ts_y[1:min_length]
      
      # Replace NA values with zeros (or interpolate if appropriate)
      ts_x[is.na(ts_x)] <- 0
      ts_y[is.na(ts_y)] <- 0
      
      # Ensure that the time series are integer sequences for rinform
      ts_x <- as.integer(ts_x)
      ts_y <- as.integer(ts_y)
      
      # Compute local transfer entropy from x to y
      te_result <- tryCatch({
        rinform::transfer_entropy(ts_x, ts_y, k = 1, local = TRUE)
      }, error = function(e) NA)
      
      if (is.na(any(te_result))) {
        next  # Skip if there's an error in computation
      }
      
      # Extract local TE values and identify where TE > 0.01
      local_te <- as.vector(te_result)
      is_high_te <- local_te > 0.1
      
      # Identify indices where TE > 0.01
      high_te_indices <- which(is_high_te)
      
      # Update leader_matrix and follower_matrix
      # For horse x (leader)
      leader_matrix[x_col, high_te_indices] <- TRUE
      # For horse y (follower)
      follower_matrix[y_col, high_te_indices] <- TRUE
    }
  }
}

# Combine leader_matrix and follower_matrix to create status_matrix
# Status codes: 0 = None, 1 = Leader, 2 = Follower, 3 = Both
status_matrix <- matrix(0, nrow = n_horses, ncol = n_timepoints)
rownames(status_matrix) <- cols

for (i in 1:n_horses) {
  for (t in 1:n_timepoints) {
    leader_status <- leader_matrix[i, t]
    follower_status <- follower_matrix[i, t]
    
    status_code <- 0
    if (leader_status && !follower_status) {
      status_code <- 1  # Leader only
    } else if (!leader_status && follower_status) {
      status_code <- 2  # Follower only
    } else if (leader_status && follower_status) {
      status_code <- 3  # Both leader and follower
    }
    # If neither leader nor follower, status_code remains 0
    
    status_matrix[i, t] <- status_code
  }
}

# Melt the status_matrix into a long-format data frame suitable for ggplot2
status_df <- melt(status_matrix)
colnames(status_df) <- c("Horse", "Time", "Status")

# Add Time in seconds (assuming time steps of 0.25s)
status_df$TimeSec <- (status_df$Time - 1) * 0.25

# Convert Status to a factor with meaningful labels
status_df$StatusFactor <- factor(
  status_df$Status,
  levels = c(0, 1, 2, 3),
  labels = c("None", "Leader", "Follower", "Both")
)

# Adjust the ordering of horses for consistent display
status_df$Horse <- factor(status_df$Horse, levels = unique(status_df$Horse))

# Create a custom color palette
status_colors <- c(
  "None" = "#f0f0f0",
  "Leader" = "#377eb8",
  "Follower" = "#e41a1c",
  "Both" = "#984ea3"
)

ggplot(status_df, aes(x = TimeSec, y = Horse, color = StatusFactor)) +
  geom_point(size = 6) +
  scale_color_manual(values = status_colors) +
  scale_x_continuous(expand = expansion(mult = 0.05)) + # Adds padding to x-axis
  scale_y_discrete(expand = expansion(mult = 0.05)) + # Adds padding to y-axis
  labs(
    title = paste("Leader-Follower Status Over Time"),
    x = "Time (seconds)",
    y = "Horse",
    color = "Status"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 40),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)) +
  ggsave(
    filename = paste0("corrected_leader_follower_status_event", event_idx, ".png"),
    width = 14,
    height = 8,
    dpi = 300
  )


#### #####


count_results
nonzero_count_results = count_results[count_results$X1 != 0,]
nonzero_count_results$X3

nonzero_count_results
mean(nonzero_count_results[nonzero_count_results[2] == 1])
mean(nonzero_count_results[nonzero_count_results[2] != 1])
mean(nonzero_count_results[nonzero_count_results[3] == 1])
mean(nonzero_count_results[nonzero_count_results[3] != 1])

mean(nonzero_count_results[nonzero_count_results[2] == 7])
mean(nonzero_count_results[nonzero_count_results[2] != 7])
mean(nonzero_count_results[nonzero_count_results[3] == 7])
mean(nonzero_count_results[nonzero_count_results[3] != 7])

# Recursively extract all matrices from a potentially nested list
extract_matrices <- function(lst) {
  if (is.matrix(lst)) return(list(lst))
  else if (is.list(lst)) return(unlist(lapply(lst, extract_matrices), recursive = FALSE))
  else return(NULL)
}

# Extract all matrices
all_matrices <- extract_matrices(list_te_results)

# Count values greater than a threshold in a matrix
count_values_above_threshold <- function(mat, threshold = 0.05) {
  sum(mat > threshold)
}

# Apply the function to each matrix and collect results
counts <- sapply(all_matrices, count_values_above_threshold)

# Calculate the average, excluding matrices with no values above the threshold
average_count <- median(counts[counts > 0])
print(average_count/2)

counts
# Create a histogram
hist_data <- hist(result_df$Duration/2, breaks = num_bins, plot = FALSE)
hist_table <- data.frame(table(counts))
hist_table <- hist_table[hist_table$counts != 0,]
expanded_data <- rep(hist_table$counts, hist_table$Freq)
expanded_data <- as.numeric(as.character(expanded_data))
hist(expanded_data, breaks=seq(0, max(expanded_data) + 2, by=1), main="Histogram with 2-Unit Bins", xlab="Values", ylab="Frequency")
expanded_data
axis(1, at=seq(0, max(expanded_data) + 2, by=))
hist_data <- hist(expanded_data, breaks=seq(0, max(expanded_data) + 1, by=1), plot=FALSE)

# Step 6: Create the barplot with specified formatting and labels
hist_plot <- barplot(hist_data$counts, col = "light blue", border = "black", 
                     main = "Duration Histogram", 
                     xlab = "Duration of Leader-Follower Relationship (s)", 
                     ylab = "Frequency", 
                     names.arg = hist_data$mids, 
                     las = 1)  # las=1 for horizontal axis labels

plot_local_te <- function(df) {
  # Get column names
  col_names <- names(df)
  
  # Initialize a list to store plots
  plot_list <- list()
  plot_index <- 1
  
  # Nested loops to calculate TE from each column to every other column
  for (source in col_names) {
    for (target in col_names) {
      if (source != target) {  # Avoid self-transfer entropy
        # Convert data to integers (assuming the data is already discretized appropriately)
        source_data <- as.integer(df[[source]])
        target_data <- as.integer(df[[target]])
        
        # Calculate local transfer entropy
        te_values <- rinform::transfer_entropy(source_data, target_data, k = 2, local = TRUE)
        
        # Create a dataframe for plotting
        te_df <- data.frame(Time = 1:length(te_values), TE = te_values, Source = source, Target = target)
        
        # Generate plot
        p <- ggplot(te_df, aes(x = Time, y = TE)) +
          geom_line() +
          labs(title = paste("Transfer Entropy from", source, "to", target),
               x = "Time", y = "Transfer Entropy") +
          theme_minimal()
        
        # Store the plot in the list
        plot_list[[paste(source, "to", target)]] <- p
        plot_index <- plot_index + 1
      }
    }
  }
  
  return(plot_list)
}

# Apply the function to each dataframe in the list and store the result
list_of_plots <- lapply(list_of_dfs_speed, plot_local_te)

generate_transfer_entropy_heatmap <- function(temat, title = "Transfer Entropy Heatmap") {
  # Exclude specific rows and columns
  exclude_vars <- c("bato_2_speed", "bato_1_speed", "velocity_discretized_time_speed")
  filtered_temat <- temat[!(rownames(temat) %in% exclude_vars), !(colnames(temat) %in% exclude_vars)]
  
  melted_temat <- melt(filtered_temat)
  
  ggheatmap <- ggplot(data = melted_temat, aes(x = Var2, y = Var1, fill = value)) +
    ggtitle(title) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                         midpoint = 0, limit = c(-0.35, 0.35), space = "Lab", 
                         name = "Transfer Entropy") +
    theme_minimal() + # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1),
          axis.text.y = element_text(size = 9)) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = sprintf("%.2f", value)), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "vertical"
    ) +
    guides(fill = guide_colorbar(barwidth = 1, barheight = 5,
                                 title.position = "right", title.hjust = 0.5))
  
  return(ggheatmap)
}

#list_of_heatmap_plots <- lapply(list_of_temats, generate_transfer_entropy_heatmap, names(list_of_dfs))


list_of_heatmap_plots <- generate_transfer_entropy_heatmap(list_of_temats[[1]], names(list_of_temats)[[1]])

for (idx in 2:length(list_of_temats)) {
  list_of_heatmap_plots[[idx]] = generate_transfer_entropy_heatmap(list_of_temats[[idx]], names(list_of_temats)[[idx]])
}

for (plot in list_of_heatmap_plots){
  print(plot)
}


library(igraph)

graph <- graph_from_data_frame(melted_temat, directed = TRUE)
filtered_edges <- get.data.frame(graph, what = "edges")[melted_temat$value >= 0.1, ]
filtered_graph <- graph_from_data_frame(filtered_edges, directed = TRUE)
layout <- layout_with_fr(filtered_graph)
vertex_color <- "lightblue"
edge_color <- "gray"
vertex_size <- 20
edge_width <- 2


# Plot the graph
plot(
  filtered_graph,
  layout = layout,
  vertex.color = vertex_color,
  vertex.size = vertex_size,
  edge.color = edge_color,
  edge.width = edge_width,
  edge.label = filtered_edges$value,
  main = "Filtered Network Graph",
  margin = 0.1
)

png("data/figures/santiago-herding-network.png")
dev.off()

###########

matrixinp = matrix(data=NA, nrow=5, ncol=5) 
colnames(matrixinp) <- colnames(speed_herding_df) 
rownames(matrixinp) <- colnames(speed_herding_df)
for(i in 1:5){ 
  for(j in 1:5){ 
    test <- rinform::transfer_entropy(speed_herding_df[[i]], speed_herding_df[[j]], k=1, local=TRUE)
    
    rle_result <- rle(abs(c(test)) > 0.05)
    
    max_true_length <- max(max(rle_result$lengths[rle_result$values]),0)
    
    matrixinp[i,j] = max_true_length/2
  } 
} 

matrixinp
lagmat <- round(matrixinp,2)
melted_lagmat <- melt(lagmat)

ggheatmap <- ggplot(data = melted_lagmat, aes(x=Var2, y=Var1, fill=value)) +
  ggtitle("Duration of Leadership Influence During Herding Periods") + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "blue", 
                       midpoint = 0, limit = c(0,11), space = "Lab", 
                       name="Duration (s)") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1),
        axis.text.y = element_text(size = 9))+
  coord_fixed()
to_save <- ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "vertical")+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5,
                               title.position = "right", title.hjust = 0.5))
to_save
ggsave("data/figures/yangon-duration.png")



#rinform::transfer_entropy(ys, xs, k = 1, local = TRUE)

acf(list_of_dfs_speed$`filtered_fullsampled_chiangmai_160617-9.1-1`$kobe_speed,list_of_dfs_speed$`filtered_fullsampled_chiangmai_160617-9.1-1`$tarumi_speed)


names(list_of_dfs_speed$`filtered_fullsampled_chiangmai_160617-9.1-1`)

list_of_dfs_speed[[1]][,1]
acf(list_of_dfs_speed[[1]][,4],list_of_dfs_speed[[1]][,5])

#### Duration ####

results = data.frame()
#acf_result$acf[,1,data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAIAAAAUCAYAAACnOeyiAAAAXklEQVR42mNgAIL///8zMYSGhjIDGYIMIiIMvECGMwMDN4M4kFEDUqIIZKwDMdSBjAsghj6Q8QPEMAAy/lOBoQekv4AYKkDGfgZeXl4RICOLQUtLiw3IUAJJMQIZ7AC2tU2tXJxOYgAAAABJRU5ErkJggg==1]
threshold <- (exp(2*1.96/sqrt(N-3)-1))/(exp(2*1.96/sqrt(N-3)+1))
threshold

transfer_entropy(list_of_dfs_speed[[1]][[4]], list_of_dfs_speed[[1]][[5]], k= 1, local = TRUE)

for (event in 1:length(list_of_dfs_speed)) {
  for (x in 4:ncol(list_of_dfs_speed[[event]])) {
    for (y in 4:ncol(list_of_dfs_speed[[event]])){
    if (x != y){
    acf_result <- acf(list_of_dfs_speed[[event]][,c(x,y)], plot = FALSE, na.action = na.pass)
    acf_result_series <- acf_result$acf[,1,1]
    N <- acf_result$n.used
    threshold <- (exp(2*1.96/sqrt(N-3)-1))/(exp(2*1.96/sqrt(N-3)+1)) + 0.25
    # Initialize the significant_lag variable
    significant_lag <- 0
    
    # Find the indices of significant lags
    significant_indices <- which(abs(acf_result_series) > threshold)
    
    # Check if there are any significant lags
    if (length(significant_indices) > 0) {
      # Find the first non-significant lag
      first_non_significant <- min(which(abs(acf_result_series) <= threshold))
      
      # Check if a non-significant lag exists
      if (!is.na(first_non_significant)) {
        # Select the continuous significant lags up to the first non-significant lag
        continuous_significant <- significant_indices[significant_indices <= first_non_significant]
        
        # Find the highest index of continuous significant lag
        significant_lag <- max(continuous_significant)
      } else {
        # If there are no non-significant lags, set to the highest significant index
        significant_lag <- max(significant_indices)
      }
    }
    
    # Store the highest significant lag for the current vehicle
    results = rbind(results, c(significant_lag))    
    }
    }
  }
}
results$X6L

table(results$X6L/4)
lag_table <- data.frame(table(results$X7L/2))
lag_table
expanded_data <- rep(lag_table$Var1, lag_table$Freq)
expanded_data <- as.numeric(as.character(expanded_data))
hist_data <- hist(expanded_data, breaks=seq(0.5, max(expanded_data) + 1, by=1), plot=FALSE)

# Step 6: Create the barplot with specified formatting and labels
hist_plot <- barplot(hist_data$counts, col = "light blue", border = "black", 
                     main = "Lag Histogram", 
                     xlab = "Lag of Leader-Follower Relationship (s)", 
                     ylab = "Frequency", 
                     names.arg = hist_data$mids, 
                     las = 1)  # las=1 for horizontal axis labels

mean(expanded_data)
median(expanded_data)

## Order analysis

#Function to tag dataframes based on the second row values in descending order
tag_dataframes <- function(df_list) {
  for (i in seq_along(df_list)) {
    df <- df_list[[i]]
    # Extract the second row values, excluding the first three columns
    second_row_values <- as.numeric(df[2, -c(1:3)])
    # Sort column indices by the second row values in descending order
    col_order <- order(-second_row_values)
    # Create a vector of tags based on the sorted order
    tags <- seq_along(col_order)
    # Map the tags to the corresponding column names
    col_names <- colnames(df)[-c(1:3)]
    names(tags) <- col_names[col_order]
    # Add the tags to the dataframe as a new row
    new_row <- c(rep(NA, 3), tags[match(colnames(df)[-c(1:3)], names(tags))])
    df <- rbind(df, new_row)
    # Update the dataframe in the list
    df_list[[i]] <- df[,-c(1:3)]
  }
  return(df_list)
}
# Tag the dataframes in the list
tagged_list_of_dfs <- tag_dataframes(list_of_dfs_speed)

print(tagged_list_of_dfs[[3]])

results = data.frame()

for (event in 1:length(tagged_list_of_dfs)) {
  for (x in 1:ncol(tagged_list_of_dfs[[event]])) {
    for (y in 1:ncol(tagged_list_of_dfs[[event]])){
      if (x != y){
        acf_result <- acf(head(tagged_list_of_dfs[[event]][,c(x,y)],-1), plot = FALSE, na.action = na.pass)
        acf_result_series <- acf_result$acf[,1,1]
        N <- acf_result$n.used
        threshold <- (exp(2*1.96/sqrt(N-3)-1))/(exp(2*1.96/sqrt(N-3)+1)) + 0.25
        # Initialize the significant_lag variable
        significant_lag <- 0
        
        # Find the indices of significant lags
        significant_indices <- which(abs(acf_result_series) > threshold)
        
        # Check if there are any significant lags
        if (length(significant_indices) > 0) {
          # Find the first non-significant lag
          first_non_significant <- min(which(abs(acf_result_series) <= threshold))
          
          # Check if a non-significant lag exists
          if (!is.na(first_non_significant)) {
            # Select the continuous significant lags up to the first non-significant lag
            continuous_significant <- significant_indices[significant_indices <= first_non_significant]
            
            # Find the highest index of continuous significant lag
            significant_lag <- max(continuous_significant)
          } else {
            # If there are no non-significant lags, set to the highest significant index
            significant_lag <- max(significant_indices)
          }
        }
        
        # Store the highest significant lag for the current vehicle
        results = rbind(results, c(significant_lag, tail(tagged_list_of_dfs[[event]][,x],1), tail(tagged_list_of_dfs[[event]][,y],1)))    
      }
    }
  }
}

df = data.frame(results)
mean(df[df[2] == 1])
mean(df[df[2] != 1])
mean(df[df[3] == 1])
mean(df[df[3] != 1])

mean(df[df[2] == 7])
mean(df[df[2] != 7])
mean(df[df[3] == 7])
mean(df[df[3] != 7])

plot(x = 1,
     type = "n",
     xlim = c(0, 100), 
     ylim = c(0, 5),
     pch = 16,
     xlab = "Index", 
     ylab = "Speed",
     main = "Speed Trajectories of the First Horse")
for (i in 1:length(tagged_list_of_dfs)){
  lines(head(tagged_list_of_dfs[[i]][,tail(tagged_list_of_dfs[[i]],1) == 1],-1))
}

plot(x = 1,
     type = "n",
     xlim = c(0, 100), 
     ylim = c(0, 5),
     pch = 16,
     xlab = "Index", 
     ylab = "Speed",
     main = "Speed Trajectories of the Last Horse")
for (i in 1:length(tagged_list_of_dfs)){
  lines(head(tagged_list_of_dfs[[i]][,tail(tagged_list_of_dfs[[i]],1) == 5],-1), col = 'red')
}
