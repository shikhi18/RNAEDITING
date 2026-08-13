library(wordcloud)
library(RColorBrewer)

# Data — MVLG genes uniquely edited only during the Mitate stage
# Categories built from the 9 genes (of 31 unique gene IDs) with a confident
# or GO-term-level PANNZER prediction (22 "UNKNOWN" genes excluded).
#
#   Transporter (2): MVLG_05383, MVLG_01436 (Major Facilitator Superfamily)
#   Kinase (1): MVLG_03327 (protein phosphorylation)
#   Protease (1): MVLG_07239 (proteolysis)
#   Helicase (1): MVLG_07364 (RNA helicase)
#   Hydrolase (1): MVLG_07084 (hydrolase activity)
#   Oxidoreductase (1): MVLG_03008 (Pyridine nucleotide-disulphide oxidoreductase)
#   DNA-binding protein (1): MVLG_07223 (DNA binding)
#   Membrane trafficking protein (1): MVLG_06006 (protein involved in membrane traffic)
categories <- c(
  "Major Facilitator Superfamily",
  "Kinase",
  "Protease",
  "Helicase",
  "Hydrolase",
  "Oxidoreductase",
  "DNA binding",
  "Membrane trafficking"
)
counts <- c(2, 1, 1, 1, 1, 1, 1, 1)

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
    "#a65628"
  ),
  ordered.colors = TRUE,
  font = 2
)

