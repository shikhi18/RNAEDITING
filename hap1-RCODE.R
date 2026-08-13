# Load necessary libraries
library(readxl)
library(ggplot2)

# Read the Excel file
fig1ef <- read_excel("Desktop/Upadated figures APRIL2025/fig1efmay16th.xlsx")

# Ensure the 'Species' factor is ordered as desired
fig1ef$Species <- factor(fig1ef$Species, levels = c("p1A1", "p1A2", "6P", "6D", "MIA1", "MIA2"))

# Create the plot with Y-axis range limited
p <- ggplot(data = fig1ef, aes(x = edit, y = frequency, fill = Species)) +
  geom_bar(position = "dodge", stat = "identity", colour = "black") +
  scale_fill_manual(values = c(
    "p1A1" = "red", "p1A2" = "red",
    "6P" = "green", "6D" = "green",
    "MIA1" = "blue", "MIA2" = "blue"
  )) +
  coord_cartesian(ylim = c(0, 30)) +  # Limit y-axis to make smaller bars visible
  facet_wrap(~haploid, ncol = 1, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 1, hjust = 1, color = "black",
      size = 22, face = "bold", family = "Arial"
    ),
    axis.text.y = element_text(
      color = "black", size = 24, face = "bold", family = "Arial"
    ),
    axis.ticks.x = element_line(color = "black", linewidth = 1.5),
    axis.title.x = element_text(
      color = "black", size = 34, face = "bold", family = "Arial"
    ),
    axis.title.y = element_text(
      color = "black", size = 18, face = "bold", family = "Arial"
    ),
    strip.text = element_text(
      color = "black", size = 24, face = "bold", family = "Arial"
    ),
    legend.text = element_text(
      color = "black", size = 14, face = "bold", family = "Arial"
    ),
    legend.title = element_text(
      color = "black", size = 20, face = "bold", family = "Arial"
    ),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  xlab("") +
  ylab("Frequency of amino acid substitution")

# Print the plot
print(p)

# Save the plot
ggsave("fig1ef1_capped.png", plot = p, width = 8.58, height = 10.4, dpi = 300)
