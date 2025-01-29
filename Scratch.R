library('tidyverse')
library('ggplot2')
library("readxl")
library("gganimate")
library("gifski")
library("lmtest")
library("hexbin")
library("qpcR")
library(RTransferEntropy)
library(lsa)
library(igraph)
library(dplyr)
library(ggpubr)
setwd("C:/Users/David Africa/Masters Code Stuff/Horses with JP")
# Set plot theme
windowsFonts(name=windowsFont("Times New Roman"))
par(family = "Times New Roman")
theme_set(theme_minimal(base_size = 14) +
            theme(text = element_text(family = "Times New Roman"),
                  plot.title = element_text(size = 18, face = "bold"),
                  axis.title = element_text(size = 20, face = "bold"),
                  axis.text = element_text(size = 17),
                  axis.line = element_line(size = 1.5),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 16)))

# Read file to get data
dataset <- read_excel("data/filtered_fullsampled/filtered_fullsampled/filtered_fullsampled_chiangmai_160617-9.1.xlsx")

dataset <- dataset %>%
  rename(`Time (s)` = ...1,
         `bato 1_x` = `bato 1`,
         `bato 1_y` = ...3,
         `bato 2_x` = `bato 2`,
         `bato 2_y` = ...5,
         `kobe_x` = kobe,
         `kobe_y` = ...7,
         `tarumi_x` = tarumi,
         `tarumi_y` = ...9,
         `maiko_x` = maiko,
         `maiko_y` = ...11,
         `kakogawa_x` = kakogawa,
         `kakogawa_y` = ...13, 
         `himeji_x` = himeji, 
         `himeji_y` = ...15)#, 
         #`akashi_x` = akashi,
         #`akashi_y` = ...17, 
         #`uji_x` = uji,
         #`uji_y` = ...19)

dataset$`Time (s)`

dataset_binned <- dataset %>%
  mutate(discretized_time = cut(`Time (s)`, breaks = seq(-0.5, max(dataset$`Time (s)`)+0.5, by = 0.5), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(discretized_time) %>%
  summarize_all(mean)  # You can change 'mean' to other aggregation functions if needed
dataset_binned <- dataset_binned %>%
  mutate(in_time_range = (`Time (s)` >= 5 & `Time (s)` <= 27))
# Display the updated tibble
print(dataset_binned)

dataset_binned

############################
# MESSING AROUND #
ggplot(dataset_binned, aes(x = `kobe_x`, y = `kobe_y`)) +
  geom_point(aes(color = "Kobe"), size = 2) +
  geom_point(aes(x = `tarumi_x`, y = `tarumi_y`, color = "Tarumi"), size = 2) +
  #geom_point(aes(x = `bato 1_x`, y = `bato 1_y`, color = "Bato 1"), size = 2) +
  #geom_point(aes(x = `bato 2_x`, y = `bato 2_y`, color = "Bato 2"), size = 2) +
  geom_point(aes(x = `himeji_x`, y = `himeji_y`, color = "Himeji"), size = 2) +
  geom_point(aes(x = `maiko_x`, y = `maiko_y`, color = "Maiko"), size = 2) +
  geom_point(aes(x = `kakogawa_x`, y = `kakogawa_y`, color = "Kakogawa"), size = 2) +
  #geom_point(aes(x = `akashi_x`, y = `akashi_y`, color = "Akashi"), size = 2) +
  #geom_point(aes(x = `uji_x`, y = `uji_y`, color = "Uji"), size = 2) +
  # Add black outlines for points in the time range
  geom_point(aes(x = kobe_x, y = kobe_y, group = in_time_range),
             data = subset(dataset_binned, in_time_range),
             color = "black", size = 2, shape = 1) +
  geom_point(aes(x = tarumi_x, y = tarumi_y, group = in_time_range),
             data = subset(dataset_binned, in_time_range),
             color = "black", size = 2, shape = 1) +
  geom_point(aes(x = himeji_x, y = himeji_y, group = in_time_range),
             data = subset(dataset_binned, in_time_range),
             color = "black", size = 2, shape = 1) +
  geom_point(aes(x = maiko_x, y = maiko_y, group = in_time_range),
             data = subset(dataset_binned, in_time_range),
             color = "black", size = 2, shape = 1) +
  geom_point(aes(x = kakogawa_x, y = kakogawa_y, group = in_time_range),
             data = subset(dataset_binned, in_time_range),
             color = "black", size = 2, shape = 1) +
  labs(title = "Horses' Movements",
       x = "X-coordinate",
       y = "Y-coordinate") +
  theme_minimal()

dataset
############################

times <- function(df, vn) {
  index <- na.omit(df[c(which(df["Vehicle Number"]==vn)),]['Time (sec)'])
  return(index)
}
#1b. Get all other cars that share x timeslices with that car

firsts <- dataset %>% group_by(`Vehicle Number`)
firsts <- firsts %>% summarize(`Vehicle Number`, first = min(`Time (sec)`)) %>% unique()

vehicles <- function(df, vn) {
  last <- as.numeric(tail(times(df,vn),n=1))
  result <- which(firsts$`first`-last < -10 & firsts$`Vehicle Number` > vn)
  return(result)
}
#2. Run Granger tests on all reasonable orders for them, aligned to time
#2a. Make a vector for each and put dummy entries at the end

maketables <- function(df, vn){
  if(length(vehicles(df,vn))==0){
    print("Vehicle has no other eligible vehicles in sliding window",)
  }
  else{
  data <- data.frame(df$`Lat Distance (m)`[df$`Vehicle Number` == vn])
  names <- c(paste0("Vehicle", vn))
  for(i in vehicles(df, vn)) {                                   # Head of for-loop
    new <- df$`Lat Distance (m)`[df$`Vehicle Number` == i]                      # Create new column
    data <- qpcR:::cbind.na(data, new)               # Append new column
    names <- append(names, paste0("Vehicle", i))
  }
  colnames(data) <- names
  return(data)
  }
}
dataset

######################
# Transfer Entropy #

test <- maketables(testset, 1002)
test
nrow(test)

entropy_results <- data.frame(timestep = 1:50)
for (n in 2:ncol(test)){
  entropy <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("timestep", "TE"))))
  if(length(test[,n]>2)){
    for(i in 2:min(length(test[,n]),length(test[,1]))){
      entropy <- rbind(entropy, c(i, coef(RTransferEntropy::transfer_entropy(test[1:i,1], test[1:i,n],nboot=100),quiet = T)[1,2]))
    }
    col_name <- paste("TE -> ", colnames(test)[n],  sep = "")
    colnames(entropy) <- c("timestep", col_name)
    entropy_results <- merge(entropy_results, entropy, by = "timestep", all = T)
  }
}

df_filtered <- entropy_results[!is.na(entropy_results$timestep), ]

# Create a line plot
line_colors <- c("red", "blue", "green", "purple", "orange", "pink", "brown", "gray")
df_filtered

# Create a line plot with a for loop
p <- ggplot(df_filtered, aes(x = timestep / 2)) +
  labs(title = "Line Plot of TE vs. timestep",
       x = "Time (s)",
       y = "TE")

for (i in 2:ncol(df_filtered)) {
  col_name <- colnames(df_filtered)[i]
  p <- p + geom_line(aes_string(y = as.name(col_name)), color = line_colors[i - 1], size = 1.5)
}

p

df_filtered

df_long <- df_filtered %>%
  pivot_longer(cols = -timestep, names_to = "variable", values_to = "value")

# Create a faceted line plot
ggplot(df_long, aes(x = timestep / 2, y = value)) +
  geom_line(color = "lightblue", size = 1.5) +
  labs(title = "Line Plot of TE vs. timestep",
       x = "Time (s)",
       y = "TE") +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.spacing = unit(1.5, "lines"))

base_size <- 12  # Adjust this as needed

# Increase the size of the text elements
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

# Print the plot
ggsave("TEfacet.png", plot, width = 8, height = 6, dpi = 300)

############### VELOCITY IFY

df <- dataset_binned[, -c(1,2,3,4,5,6)]
df
velocity_df <- data.frame(
  matrix(NA, ncol = 0, nrow = 57)
)

velocity_df

for (i in seq(1, ncol(df), by = 2)) {
  x_col <- i
  y_col <- i + 1
  
  # Calculate differences in positions
  dx <- c(0, diff(df[[x_col]]))
  dy <- c(0, diff(df[[y_col]]))
  
  # Calculate velocity as the rate of change of position
  dt <- 0.25
  velocity_df[paste0("velocity_", names(df)[x_col])] <- dx / dt
  velocity_df[paste0("velocity_", names(df)[y_col])] <- dy / dt
}

# Print the resulting data frame with velocity data
print(df)
velocity_df
############### SPEEDIFY (S = sqrt(V_X^2 + V_Y^2))
library(dplyr)

velocity_df

bato_1_speed <- sqrt(velocity_df$velocity_bato_1_x^2 + velocity_df$velocity_bato_1_y^2)
bato_2_speed <- sqrt(velocity_df$velocity_bato_2_x^2 + velocity_df$velocity_bato_2_y^2)
kobe_speed <- sqrt(velocity_df$velocity_kobe_x^2 + velocity_df$velocity_kobe_y^2)
tarumi_speed <- sqrt(velocity_df$velocity_tarumi_x^2 + velocity_df$velocity_tarumi_y^2)
maiko_speed <- sqrt(velocity_df$velocity_maiko_x^2 + velocity_df$velocity_maiko_y^2)
kakogawa_speed <- sqrt(velocity_df$velocity_kakogawa_x^2 + velocity_df$velocity_kakogawa_y^2)
himeji_speed <- sqrt(velocity_df$velocity_himeji_x^2 + velocity_df$velocity_himeji_y^2)
akashi_speed <- sqrt(velocity_df$velocity_akashi_x^2 + velocity_df$velocity_akashi_y^2)
uji_speed <- sqrt(velocity_df$velocity_uji_x^2 + velocity_df$velocity_uji_y^2)

kobe_speed

# Step 2: Create a new dataframe with the calculated speeds
speed_df <- data.frame(
  kobe_speed,
  tarumi_speed,
  maiko_speed,
  kakogawa_speed,
  himeji_speed,
  akashi_speed,
  uji_speed
)

print(speed_df)

speed

############### LOCAL TRANSFER ENTROPY USING DEPRECATED RINFORM LIBRARY
library(devtools)
install_github("ELIFE-ASU/rinform")
library("rinform")
rinform::transfer_entropy(
  na.omit(test$Vehicle1006[1:28]),
  na.omit(test$Vehicle1002[1:28]),
  k = 1,
  local = TRUE
)

dataset
dataset[,6]

colnames(dataset)

xx <- rinform::transfer_entropy(
  dataset$kobe_x,
  dataset$tarumi_x,
  k = 1,
  local = TRUE
)

xy <- rinform::transfer_entropy(
  dataset$kobe_x,
  dataset$tarumi_y,
  k = 1,
  local = TRUE
)

yx <- rinform::transfer_entropy(
  dataset$kobe_y,
  dataset$tarumi_x,
  k = 1,
  local = TRUE
)

yy <- rinform::transfer_entropy(
  dataset$kobe_y,
  dataset$tarumi_y,
  k = 1,
  local = TRUE
)

df_plot <- data.frame(
  Time = head(dataset$`Time (s)`, -1),
  KobeXtoTarumiX = xx,
  KobeXtoTarumiY = xy,
  KobeYtoTarumiX = yx,
  KobeYtoTarumiY = yy
)

# Reshape the data for ggplot
df_plot_long <- tidyr::gather(df_plot, key = "Series", value = "TransferEntropy", -Time)

# Round up the minimum and maximum values to the nearest whole number
min_time <- ceiling(min(df_plot$Time))
max_time <- ceiling(max(df_plot$Time))

# Plot the data using ggplot with increased line size and removed grid lines
ggplot(df_plot_long, aes(x = Time, y = TransferEntropy, color = Series)) +
  geom_line(size = 1.35, alpha = 0.5) +  # Increase line size
  labs(title = "Transfer Entropy from Kobe to Tarumi",
       x = "Time (s)",
       y = "Transfer Entropy") +
  scale_x_continuous(breaks = seq(min_time, max_time, by = 2)) +  # Set x-axis ticks at every 0.5 seconds
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


xx <- rinform::transfer_entropy(
  dataset$kakogawa_x,
  dataset$tarumi_x,
  k = 1,
  local = TRUE
)

xy <- rinform::transfer_entropy(
  dataset$kakogawa_x,
  dataset$tarumi_y,
  k = 1,
  local = TRUE
)

yx <- rinform::transfer_entropy(
  dataset$kakogawa_y,
  dataset$tarumi_x,
  k = 1,
  local = TRUE
)

yy <- rinform::transfer_entropy(
  dataset$kakogawa_y,
  dataset$tarumi_y,
  k = 1,
  local = TRUE
)

df_plot <- data.frame(
  Time = head(dataset$`Time (s)`, -1),
  KakogawaXtoTarumiX = xx,
  KakogawaXtoTarumiY = xy,
  KakogawaYtoTarumiX = yx,
  KakogawaYtoTarumiY = yy
)

# Reshape the data for ggplot
df_plot_long <- tidyr::gather(df_plot, key = "Series", value = "TransferEntropy", -Time)

ggplot(df_plot_long, aes(x = Time, y = TransferEntropy, color = Series)) +
  geom_line(size = 1.35, alpha = 0.5) +  # Increase line size
  labs(title = "Transfer Entropy from Kakogawa to Tarumi",
       x = "Time (s)",
       y = "Transfer Entropy") +
  scale_x_continuous(breaks = seq(min_time, max_time, by = 2)) +  # Set x-axis ticks at every 0.5 seconds
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

rinform::transfer_entropy(speed_df[[2]], speed_df[[4]], k = 1, local = TRUE)


-#####

rinform::transfer_entropy()

horse_pairs <- c("kobe", "tarumi", "maiko", "kakogawa", "himeji", "akashi", "uji")
plot_list <- c()
plot_names <- c()

for (i in seq_along(horse_pairs)){
  for (j in seq_along(horse_pairs)) {
    ix_col <- paste(horse_pairs[i], "_x", sep = "")
    iy_col <- paste(horse_pairs[i], "_y", sep = "")
    jx_col <- paste(horse_pairs[j], "_x", sep = "")
    jy_col <- paste(horse_pairs[j], "_y", sep = "")
    
    plot_name <- paste("plot_", horse_pairs[i], horse_pairs[j], sep = "")
    
    # Create a plot for each pair
    xx <- rinform::transfer_entropy(
      dataset[[ix_col]],
      dataset[[jx_col]],
      k = 1,
      local = TRUE
    )
    
    xy <- rinform::transfer_entropy(
      dataset[[ix_col]],
      dataset[[jy_col]],
      k = 1,
      local = TRUE
    )
    
    yx <- rinform::transfer_entropy(
      dataset[[iy_col]],
      dataset[[jx_col]],
      k = 1,
      local = TRUE
    )
    
    yy <- rinform::transfer_entropy(
      dataset[[iy_col]],
      dataset[[jy_col]],
      k = 1,
      local = TRUE
    )
    
    df_plot <- data.frame(
      Time = head(dataset$`Time (s)`, -1),
      XtoX = xx,
      XtoY = xy,
      YtoX = yx,
      YtoY = yy
    )
    
    df_plot_long <- tidyr::gather(df_plot, key = "Series", value = "TransferEntropy", -Time)
    
    min_time <- ceiling(min(df_plot$Time))
    max_time <- ceiling(max(df_plot$Time))
    
    plot <- ggplot(df_plot_long, aes(x = Time, y = TransferEntropy, color = Series)) +
      geom_line(size = 1.35, alpha = 0.5) +
      labs(title = plot_name,
           x = "Time (s)",
           y = "Transfer Entropy") +
      scale_x_continuous(breaks = seq(min_time, max_time, by = 2)) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    assign(plot_name, plot)
    plot_list <- append(plot_list, plot)
    plot_names <- append(plot_names, plot_name)
  }
}

plot_kobetarumi

n <- length(horse_pairs)

horse_pairs
allplots <- ggarrange(plot_kobekobe, plot_kobetarumi, plot_kobemaiko, plot_kobekakogawa, plot_kobehimeji,plot_kobeakashi,plot_kobeuji,
                      plot_tarumikobe,plot_tarumitarumi,plot_tarumimaiko,plot_tarumikakogawa, plot_tarumihimeji,plot_tarumiakashi,plot_tarumiuji,
                      ncol = 7, nrow = 7)

allplots

#####

# HEAT MAP #

dataset
for_cor <- dataset[ -c(1:5) ]

cormat <- round(cor(for_cor),2)
cormat
head(cormat)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
ggheatmap <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggheatmap + 
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



# Getting a TE matrix WITH SPEED #

test <- rinform::transfer_entropy(diff_speed_df$tarumi_speed, diff_speed_df$kobe_speed, k =1, local = TRUE)
plot(seq(0,length(test),length=length(test)), test, type = 'l')
length(test)
speed_df

matrixinp = matrix(data=NA, nrow=7, ncol=7) 
colnames(matrixinp) <- colnames(speed_df) 
rownames(matrixinp) <- colnames(speed_df)
for(i in 1:7){ 
  for(j in 1:7){ 
    matrixinp[i,j] = RTransferEntropy::transfer_entropy(speed_df[[i]],speed_df[[j]], lx = 1, ly = 1,nboot = 300)$coef[1]
  } 
} 
print(matrixinp)

temat <- round(matrixinp,2)
temat


library(reshape2)
melted_temat <- melt(temat)
melted_temat
ggheatmap <- ggplot(data = melted_temat, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0, limit = c(-0.2,0.2), space = "Lab", 
                       name="Transfer Entropy") +
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


diff_speed_df <- as.data.frame(lapply(speed_df, diff))

# Print the new data frame
print(diff_speed_df)
matrixinp2 = matrix(data=NA, nrow=7, ncol=7) 
colnames(matrixinp2) <- colnames(diff_speed_df) 
rownames(matrixinp2) <- colnames(diff_speed_df)
for(i in 1:7){ 
  for(j in 1:7){ 
    matrixinp2[i,j] = RTransferEntropy::transfer_entropy(diff_speed_df[[i]],diff_speed_df[[j]], lx = 1, ly = 1,nboot = 100)$coef[1]
  } 
} 
print(matrixinp2)

temat2 <- round(matrixinp2,2)
temat2
library(reshape2)
melted_temat2 <- melt(temat2)
melted_temat2
ggheatmap <- ggplot(data = melted_temat2, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0, limit = c(-0.2,0.2), space = "Lab", 
                       name="Transfer Entropy") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1),
        axis.text.y = element_text(size = 9))+
  coord_fixed()
to_save2 <- ggheatmap + 
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
to_save2

library(igraph)

graph <- graph_from_data_frame(melted_temat, directed = TRUE)
filtered_edges <- get.data.frame(graph, what = "edges")[melted_temat$value >= 0.05, ]
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
######


p <- ggplot(entropy_df, aes(x = Index, y = Entropy)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_smooth(color = "light blue", size = 1.25, se = FALSE) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  labs(x = "Time (s)", y = "Local Transfer Entropy from Vehicle 1002 to Candidate Followers") +
  theme(
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Adjust size here
    panel.spacing = unit(1, "lines")
  ) + 
  scale_x_continuous(breaks = c(0,2,4,6,8,10))

runLocalEtest <- function(df, vn){
  entropy_df <- data.frame(Vehicle = character(), Variable = character(), Entropy = numeric())
  test <- maketables(df, vn)
  main_name <- paste('Vehicle',vn,sep ="")
  for (n in 2:ncol(test)) {
  col_name <- names(test)[n]  # Get the column name
  len <- min(length(na.omit(test[,1])), length(na.omit(test[, n])))
    # Calculate entropy
    entropy <- rinform::transfer_entropy(
      na.omit(test[,1][1:len]),
      na.omit(test[,n][1:len]),
      k = 1,
      local = TRUE
    )
    
    # Create a data frame for the current variable
    var_df <- data.frame(
      Vehicle = main_name,
      Variable = col_name,
      Entropy = entropy
    )
    
    # Append to the main data frame
    entropy_df <- bind_rows(entropy_df, var_df)
  }
  return(entropy_df)
  }


runLocalEtestAll <- function(df) {
  entropy_df <- data.frame(Vehicle = character(), Variable = character(), Entropy = numeric())
  
  for (vn in unique(df$`Vehicle Number`)) {
    test <- maketables(df, vn)
    main_name <- paste('Vehicle', vn, sep = "")
    if(is.data.frame(test)){
      for (n in 2:ncol(test)) {
        col_name <- names(test)[n]  # Get the column name
        len <- min(length(na.omit(test[, 1])), length(na.omit(test[, n])))
        
        # Calculate entropy
        entropy <- rinform::transfer_entropy(
          na.omit(test[, 1][1:len]),
          na.omit(test[, n][1:len]),
          k = 1,
          local = TRUE
        )
        
        # Create a data frame for the current variable
        var_df <- data.frame(
          Vehicle = main_name,
          Variable = col_name,
          Entropy = entropy
        )
        
        # Append to the main data frame
        entropy_df <- rbind(entropy_df, var_df)
      }
    }
  }
  
  return(entropy_df)
}

result <- runLocalEtestAll(testset)

filtered_result <- result[result$Entropy > 0.1, ]

# Get all unique pairs of Vehicle values
pairs <- unique(cbind(filtered_result$Vehicle, filtered_result$Variable))

# Calculate the duration of each pair by counting the rows
pair_durations <- sapply(1:nrow(pairs), function(i) {
  pair <- pairs[i, ]
  duration <- nrow(filtered_result[filtered_result$Vehicle == pair[1] & filtered_result$Variable == pair[2], ])
  return(duration)
})

# Create a data frame with pairs and their durations
result_df <- data.frame(Pair = sapply(1:nrow(pairs), function(i) {
  pair <- pairs[i, ]
  return(paste(pair[1], "-", pair[2]))
}), Duration = pair_durations)

mean(result_df$Duration)

num_bins <-  max(result_df$Duration)/2
x_range <- c(0, max(result_df$Duration)/2)

breaks <- seq(min(result_df$Duration), max(result_df$Duration), length.out = num_bins + 1)

# Create a histogram
hist_data <- hist(result_df$Duration, breaks = breaks, main = "Duration Histogram", xlab = "Duration", ylab = "Frequency", col = "blue", border = "black")

# Plot the histogram
durationhistplot<- plot(hist_data, col = "light blue", border = "black", main = "Histogram of Duration of Leader-Follower Relationships", xlab = "Duration of Leader-Follower Relationship (s)", ylab = "Frequency")
  
  x_ticks <- seq(0, max(result_df$Duration), 1)
  axis(1, at = x_ticks, labels = x_ticks)

  ggsave("duration_histogram.png", plot = durationhistplot, width = 8, height = 6, dpi = 800)
  

# Create a histogram
hist_data <- hist(result_df$Duration/2, breaks = num_bins, plot = FALSE)

# Plot the histogram with more frequent x-axis ticks and labels
hist_plot <- barplot(hist_data$counts, xlim = x_range, col = "light blue", border = "black", 
                     main = "Duration Histogram", xlab = "Duration of Leader-Follower Relationship (s)", ylab = "Frequency")
# Define the positions for x-axis ticks and labels
x_ticks <- seq(0, max(result_df$Duration), 1)  # You can adjust the "by" value for the desired frequency


# Add x-axis ticks and labels
axis(1, at = x_ticks, labels = x_ticks)

# Add a grid
grid()

test <- maketables(dataset, 1002)

testset

# Create an empty data frame to store the entropy values
entropy_df <- data.frame(Vehicle = character(), Variable = character(), Entropy = numeric())

# Iterate over columns of the test dataframe (starting from column 2)
for (n in 2:ncol(test)) {
  col_name <- names(test)[n]  # Get the column name
  len <- min(length(na.omit(test$Vehicle1002)), length(na.omit(test[, n])))
  
  # Calculate entropy
  entropy <- rinform::transfer_entropy(
    na.omit(test$Vehicle1002)[1:len],
    na.omit(test[, n])[1:len],
    k = 1,
    local = TRUE
  )
  
  # Create a data frame for the current variable
  var_df <- data.frame(
    Vehicle = "Vehicle1002",
    Variable = col_name,
    Entropy = entropy
  )
  
  # Append to the main data frame
  entropy_df <- bind_rows(entropy_df, var_df)
}

entropy_df <- entropy_df %>%
  group_by(Vehicle) %>%
  mutate(Index = row_number()%%11)

# Create a ggplot object
p <- ggplot(entropy_df, aes(x = Index, y = Entropy)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_smooth(color = "light blue", size = 1.25, se = FALSE) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  labs(x = "Time (s)", y = "Local Transfer Entropy from Vehicle 1002 to Candidate Followers") +
  theme(
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Adjust size here
    panel.spacing = unit(1, "lines")
  ) + 
  scale_x_continuous(breaks = c(0,2,4,6,8,10))
p
entropy_df

p

shaded_regions_df <- entropy_df %>%
  group_by(Variable) %>%
  filter(any(Entropy > 0.1))

# Create the plot with shaded region

p <- ggplot(entropy_df, aes(x = Index, y = Entropy)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_smooth(color = "light blue", size = 1.25, se = FALSE) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  labs(x = "Time (s)", y = "Local Transfer Entropy from Vehicle 1002 to Candidate Followers") +
  theme(
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(1, "lines")
  ) + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  geom_rect(
    aes(xmin = 2, xmax = 6.5, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1003")
  ) +
  geom_rect(
    aes(xmin = 2, xmax = 6.5, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1004")
  ) +
  geom_rect(
    aes(xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1005")
  ) +
  geom_rect(
    aes(xmin = 0.5, xmax = 6, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1006")
  ) +
  geom_rect(
    aes(xmin = 2, xmax = 9, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1007")
  ) +
  geom_rect(
    aes(xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1008")
  ) +
  geom_rect(
    aes(xmin = 9.5, xmax = 10, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1008")
  ) +
  geom_rect(
    aes(xmin = 0, xmax = 4.5, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1009")
  ) +
  geom_rect(
    aes(xmin = 7, xmax = 10, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1009")
  ) +
  coord_cartesian(ylim = c(-0.05, 0.4))  # Set the y-axis limits
  
# Print the plot
print(p)
ggsave("TEfacet.png", p, width = 8, height = 8, dpi = 800)


##################

library('tidyverse')
library('ggplot2')
library("readxl")
library("gganimate")
library("gifski")
library("lmtest")
library("hexbin")
library("qpcR")
library(RTransferEntropy)
library(lsa)
library(igraph)
library(dplyr)
library(deSolve)
library(viridis)
setwd("~/Thesis/chennai data from Kanagaraj")

