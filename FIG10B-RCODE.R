# Load necessary libraries
library(ggplot2)
library(readxl)

# Read the data from the Excel file
data <- read_excel("~/Desktop/NEW2.xlsx")

# Plotting the data with increased text size
ggplot(data, aes(x = log2FoldChange, y = reorder(Functions, log2FoldChange))) + 
  geom_bar(stat = "identity", fill = "blue") +
  labs(y = "Functions", x = "log2 Fold Change") +
  ggtitle("") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),  # Increase size of axis titles
    axis.text = element_text(size = 18),   # Increase size of axis text
    plot.title = element_text(size = 16)   # Increase size of plot title
  )



