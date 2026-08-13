library(wordcloud)
library(RColorBrewer)

# Data — MvSup, 56 genes common to all four stages (6P, 6D, Mated, Infection)
# prediction (19 "hypothetical protein / no confident annotation" genes excluded,
# "domain-containing protein" genes excluded, and the miscellaneous).
#
#   Kinase (4): g1196, g8780, g8794, g13895
#   RNA-binding protein (3): g479, g11705, g12441
#   DNA-binding protein (3): g1225, g8179, g14676
#   Dehydrogenase (2): g1165, g13095
#   Isomerase (2): g1106, g1222
#   Polymerase (2): g3766, g13846
#   Transporter (2): g547, g13842
#   Ligase (1): g14321
#   Acyltransferase (1): g12565
#   Dehydratase (1): g1189

categories <- c(
  "Kinase",
  "RNA-binding protein",
  "DNA-binding protein",
  "Dehydrogenase",
  "Isomerase",
  "Polymerase",
  "Transporter",
  "Ligase",
  "Acyltransferase",
  "Dehydratase"
)
counts <- c(4, 3, 3, 2, 2, 2, 2, 1, 1, 1)

# Add counts to labels
labels <- paste0(categories, " (", counts, ")")

# Remove plot margins
par(mar = c(0, 0, 0, 0))

# Reproducible layout
set.seed(123)

# Create word cloud
wordcloud(
  words = labels,
  freq = counts,
  scale = c(3.5, 1),
  min.freq = 1,
  random.order = FALSE,
  rot.per = 0,
  colors = c(
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#984ea3",
    "#a65628",
    "#66c2a5",
    "#fc8d62"
  ),
  ordered.colors = TRUE,
  font = 2
)
