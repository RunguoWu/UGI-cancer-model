library(ggplot2)

# Create assumptional data structure
# You can replace these values with your actual data
# Only include region-cancer combinations that actually exist
data <- data.frame(
  cancer_type = c(
    # Gastric cancer (total: 33)
    rep("Gastric", 5),
    # Oesophageal cancer (total: 23)
    rep("Oesophageal", 5),
    # Pancreatic cancer (total: 5)
    rep("Pancreatic", 4),
    # Multi-cancer (total: 4)
    rep("Multi-cancer", 3)
  ),
  region = c(
    # Gastric regions
    "North America", "East Asia", "Europe", "Latin America", "Africa",
    # Oesophageal regions
    "North America", "East Asia", "Europe", "Latin America", "Oceania",
    # Pancreatic regions
    "North America", "East Asia", "Europe", "Latin America",
    # Multi-cancer regions
    "North America", "East Asia", "Europe"
  ),
  count = c(
    # Gastric cancer by region
    10, 15, 5, 2, 1,
    # Oesophageal cancer by region
    8, 10, 3, 1, 1,
    # Pancreatic cancer by region
    2, 2, 0.5, 0.5,
    # Multi-cancer by region
    1, 2, 1
  )
)

# Calculate totals for cancer types
cancer_totals <- aggregate(count ~ cancer_type, data, sum)
cancer_totals$percentage <- round(cancer_totals$count / sum(cancer_totals$count) * 100, 1)

# Create the nested pie chart
# Inner ring: Cancer types
# Outer ring: Regions within each cancer type

# Prepare data for plotting
data$cancer_type <- factor(data$cancer_type, 
                           levels = c("Gastric", "Oesophageal", "Pancreatic", "Multi-cancer"))

# Calculate cumulative positions for inner pie
cancer_totals$fraction <- cancer_totals$count / sum(cancer_totals$count)
cancer_totals$ymax <- cumsum(cancer_totals$fraction)
cancer_totals$ymin <- c(0, head(cancer_totals$ymax, -1))

# Calculate positions for outer pie
# Must align with inner ring by cancer type
data <- data[order(data$cancer_type, data$region), ]

# Calculate positions within each cancer type
data_list <- list()
for (cancer in unique(data$cancer_type)) {
  # Get data for this cancer type
  cancer_data <- data[data$cancer_type == cancer, ]
  
  # Get the boundaries for this cancer type from cancer_totals
  cancer_info <- cancer_totals[cancer_totals$cancer_type == cancer, ]
  ymin_start <- cancer_info$ymin
  ymax_end <- cancer_info$ymax
  cancer_range <- ymax_end - ymin_start
  
  # Calculate fraction of each region within this cancer type
  cancer_data$fraction_within <- cancer_data$count / sum(cancer_data$count)
  
  # Calculate cumulative positions within this cancer's section
  cancer_data$ymin <- ymin_start + c(0, cumsum(cancer_data$fraction_within[-nrow(cancer_data)])) * cancer_range
  cancer_data$ymax <- ymin_start + cumsum(cancer_data$fraction_within) * cancer_range
  
  data_list[[cancer]] <- cancer_data
}

# Combine back into one dataframe
data <- do.call(rbind, data_list)
rownames(data) <- NULL

# Define colors - only for regions
region_colors <- c("North America" = "#D55E00", "East Asia" = "#CC79A7",
                   "Europe" = "#0072B2", "Latin America" = "#999999", 
                   "Africa" = "#E6AB02", "Oceania" = "#66A61E")

# Create the plot
p <- ggplot() +
  # Outer ring (regions)
  geom_rect(data = data, 
            aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = region),
            color = "white", size = 0.5) +
  # Inner ring (cancer types) - no fill, just white background
  geom_rect(data = cancer_totals,
            aes(xmin = 2, xmax = 3, ymin = ymin, ymax = ymax),
            fill = "white", color = "gray30", size = 0.5) +
  # Add labels for cancer types with line breaks for long names
  geom_text(data = cancer_totals,
            aes(x = 2.5, y = (ymin + ymax)/2, 
                label = paste0(gsub("-", "-\n", cancer_type), "\n(n=", count, ")")),
            size = 3, fontface = "bold") +
  # Convert to polar coordinates
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  scale_fill_manual(values = region_colors, name = "Region") +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  labs(title = "Literature Review: Studies by Cancer Type and Region (N=65)",
       fill = "Category")

print(p)

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("\nBy Cancer Type:\n")
print(cancer_totals[, c("cancer_type", "count", "percentage")])

cat("\nBy Region (across all cancers):\n")
region_totals <- aggregate(count ~ region, data, sum)
region_totals$percentage <- round(region_totals$count / sum(region_totals$count) * 100, 1)
print(region_totals)

cat("\nDetailed breakdown:\n")
print(data[, c("cancer_type", "region", "count")])