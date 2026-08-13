library(openxlsx)
library(reshape2)
library(ggplot2)

# Load the Excel file
file_path <- "~/Desktop/PLOS_PATHOGENS_NEW2.xlsx"
data2 <- read.xlsx(file_path, sheet = "UPDATED FIG8B-MVLG COMMON")

# Extract the Functions column
function_data <- data2$Functions

# Remove Functions column and re-attach as factor
plot_data <- data2[, !(names(data2) %in% "Functions")]
plot_data$Functions <- factor(function_data, levels = unique(function_data))

# Melt data to long format
melted_data <- melt(plot_data, id.vars = "Functions", 
                    variable.name = "Stage", 
                    value.name = "Frequency")

# Check unique Stage values
print(unique(melted_data$Stage))  # Very important!

ggplot(melted_data, aes(x = Functions, y = Frequency, fill = Stage)) + 
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.5) +  # <-- Adjusted bar width
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "", x = "", y = "Frequency of editing sites") +
  theme(
    text = element_text(size = 14, family = "Arial", face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 24, face = "bold", color = "black"),
    axis.text.y = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 23, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 0, face = "bold")
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black")))
