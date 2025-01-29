# R script for simulating horse herding with specified coefficients and grazing periods.

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Constants for distances in Body Lengths (BL)
d_repulsion_stallion <- 9    # Repulsion distance from the stallion
d_short_range <- 4           # Short-range repulsion

# Initialize position and direction of mares and stallion
set.seed(42)
initial_positions <- data.frame(x = runif(6, 4, 6), y = runif(6, 4, 6))  # Initial positions for 6 mares
stallion_position <- data.frame(x = 5, y = 5)  # Example fixed position for the stallion

# Provided coefficients for each horse
coefficients <- data.frame(
  inertia = c(0.600, 0.510, 0.355, 0.608, 0.912, 0.675),
  repulsion = c(0.029, 0.193, 0.403, 0.365, 0.398, 0.281),
  short_range = c(0.0, 0.030, 0.209, 0.205, 0.037, 0.0),
  com_attraction = c(0.058, 0.037, 0.620, 0.426, 0.287, 0.330)
)

# Function to calculate center of mass (COM) for the group
calculate_COM <- function(mares) {
  com_x <- mean(mares$x)
  com_y <- mean(mares$y)
  return(data.frame(x = com_x, y = com_y))
}

# Function to calculate forces for each mare in the simplified model
calculate_forces <- function(mares, stallion, coeffs) {
  n <- nrow(mares)
  total_forces <- data.frame(x = rep(0, n), y = rep(0, n))
  
  # Calculate the center of mass
  COM <- calculate_COM(mares)
  
  # Calculate each force for each mare
  for (i in 1:n) {
    # Inertia (keep moving in the same direction)
    inertia_x <- coeffs$inertia[i] * rnorm(1, mean = 0, sd = 0.1)
    inertia_y <- coeffs$inertia[i] * rnorm(1, mean = 0, sd = 0.1)
    
    # Repulsion from stallion
    dist_to_stallion <- sqrt((mares$x[i] - stallion$x)^2 + (mares$y[i] - stallion$y)^2)
    repulsion_x <- 0
    repulsion_y <- 0
    if (dist_to_stallion < d_repulsion_stallion) {
      repulsion_x <- coeffs$repulsion[i] * (mares$x[i] - stallion$x) / dist_to_stallion
      repulsion_y <- coeffs$repulsion[i] * (mares$y[i] - stallion$y) / dist_to_stallion
    }
    
    # Short-range repulsion with other mares
    short_repulsion_x <- 0
    short_repulsion_y <- 0
    for (j in 1:n) {
      if (i != j) {
        dist_to_mare <- sqrt((mares$x[i] - mares$x[j])^2 + (mares$y[i] - mares$y[j])^2)
        if (dist_to_mare < d_short_range) {
          short_repulsion_x <- short_repulsion_x + coeffs$short_range[i] * (mares$x[i] - mares$x[j]) / dist_to_mare
          short_repulsion_y <- short_repulsion_y + coeffs$short_range[i] * (mares$y[i] - mares$y[j]) / dist_to_mare
        }
      }
    }
    
    # Attraction to the center of mass (COM)
    dist_to_com <- sqrt((mares$x[i] - COM$x)^2 + (mares$y[i] - COM$y)^2)
    com_attraction_x <- coeffs$com_attraction[i] * (COM$x - mares$x[i]) / dist_to_com
    com_attraction_y <- coeffs$com_attraction[i] * (COM$y - mares$y[i]) / dist_to_com
    
    # Total force acting on each mare
    total_forces$x[i] <- inertia_x + repulsion_x + short_repulsion_x + com_attraction_x
    total_forces$y[i] <- inertia_y + repulsion_y + short_repulsion_y + com_attraction_y
  }
  
  return(total_forces)
}

# Function to generate a random target point at least 5 units away from the stallion's current position
generate_random_target <- function(stallion, min_distance = 30) {
  repeat {
    target_x <- runif(1, stallion$x - 10, stallion$x + 10)
    target_y <- runif(1, stallion$y - 10, stallion$y + 10)
    distance <- sqrt((target_x - stallion$x)^2 + (target_y - stallion$y)^2)
    
    if (distance >= min_distance) {
      return(data.frame(x = target_x, y = target_y))
    }
  }
}

simulate_herding <- function(mares, stallion, grazing_steps, herding_steps, coeffs) {
  history <- list()
  target_position <- stallion  # Initialize target_position for stallion
  
  # Grazing period at the beginning
  for (t in 1:grazing_steps) {
    # Mares graze undirectedly
    mares$x <- mares$x + rnorm(nrow(mares), mean = 0, sd = 0.5)
    mares$y <- mares$y + rnorm(nrow(mares), mean = 0, sd = 0.5)
    
    # Stallion also grazes undirectedly
    stallion$x <- stallion$x + rnorm(1, mean = 0, sd = 0.5)
    stallion$y <- stallion$y + rnorm(1, mean = 0, sd = 0.5)
    
    # Capture current positions with separate columns for each mare
    mare_positions <- as.vector(t(mares))
    names(mare_positions) <- paste0("mare", rep(1:nrow(mares), each = 2), "_", c("x", "y"))
    
    history[[t]] <- data.frame(time = t, state = "Grazing", stallion_x = stallion$x, stallion_y = stallion$y, mare_positions)
  }
  
  # Set a random target point for the stallion at the start of herding
  target_position <- generate_random_target(stallion, min_distance = 5)
  
  # Herding period
  for (t in (grazing_steps + 1):(grazing_steps + herding_steps)) {
    # Calculate forces for mares
    forces <- calculate_forces(mares, stallion, coeffs)
    
    # Update positions of mares based on forces
    mares$x <- mares$x + forces$x
    mares$y <- mares$y + forces$y
    
    # Move the stallion towards the target point
    dist_to_target <- sqrt((target_position$x - stallion$x)^2 + (target_position$y - stallion$y)^2)
    if (dist_to_target > 0) {
      stallion$x <- stallion$x + (target_position$x - stallion$x) / dist_to_target * 0.1  # Step size towards target
      stallion$y <- stallion$y + (target_position$y - stallion$y) / dist_to_target * 0.1
    }
    
    # Capture current positions with separate columns for each mare
    mare_positions <- as.vector(t(mares))
    names(mare_positions) <- paste0("mare", rep(1:nrow(mares), each = 2), "_", c("x", "y"))
    
    history[[t]] <- data.frame(time = t, state = "Herding", stallion_x = stallion$x, stallion_y = stallion$y, mare_positions)
  }
  
  # Grazing period at the end
  for (t in (grazing_steps + herding_steps + 1):(grazing_steps * 2 + herding_steps)) {
    # Mares graze undirectedly
    mares$x <- mares$x + rnorm(nrow(mares), mean = 0, sd = 0.5)
    mares$y <- mares$y + rnorm(nrow(mares), mean = 0, sd = 0.5)
    
    # Stallion also grazes undirectedly
    stallion$x <- stallion$x + rnorm(1, mean = 0, sd = 0.5)
    stallion$y <- stallion$y + rnorm(1, mean = 0, sd = 0.5)
    
    # Capture current positions with separate columns for each mare
    mare_positions <- as.vector(t(mares))
    names(mare_positions) <- paste0("mare", rep(1:nrow(mares), each = 2), "_", c("x", "y"))
    
    history[[t]] <- data.frame(time = t, state = "Grazing", stallion_x = stallion$x, stallion_y = stallion$y, mare_positions)
  }
  
  return(do.call(rbind, history))
}

# Run the simulation
grazing_steps <- 1000
herding_steps <- 100
simulation_data <- simulate_herding(initial_positions, stallion_position, grazing_steps, herding_steps, coefficients)

# Reshape the data to long format for easier plotting with ggplot2
simulation_long <- simulation_data %>%
  pivot_longer(cols = starts_with("mare"), 
               names_to = c("mare", "position"),
               names_sep = "_") %>%
  pivot_wider(names_from = position, values_from = value)

# Plot the data
ggplot(simulation_long, aes(x = x, y = y, color = mare)) +
  geom_path(aes(group = mare), alpha = 0.7) +
  geom_path(data = simulation_data, aes(x = stallion_x, y = stallion_y), 
            color = "black", linetype = "dashed", size = 1.2, inherit.aes = FALSE) +
  labs(title = "Simulated Herding and Grazing Behavior of Mares and Stallion Over Time",
       x = "X Position", y = "Y Position") +
  theme_minimal() +
  theme(legend.position = "right")

