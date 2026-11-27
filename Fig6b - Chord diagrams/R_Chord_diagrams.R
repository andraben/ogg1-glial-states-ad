# ------------------------------------------------------------------------------
# Fig. 6b - Chord diagram of GO terms per data set and genotype
#
# Input:  GO_terms.xlsx
# Output: Chord_diagrams_all.svg
# ------------------------------------------------------------------------------

# Load libraries
library(circlize)
library(readxl)
library(scales)

# Read input data
data <- read_excel("GO_terms.xlsx")

# Trim whitespace in column names
colnames(data) <- gsub("^\\s+|\\s+$", "", colnames(data))

# Convert log10.P.value to numeric and set missing values to 0
data$log10.P.value <- as.numeric(data$log10.P.value)
data$log10.P.value[is.na(data$log10.P.value)] <- 0

# Category colors
category_colors <- c(
  "RNA" = "#25516a",
  "ATAC" = "#02323b",
  "Microglia" = "#5da8a3",
  "Astrocytes" = "#3f8ab2",
  "Oligodendrocytes" = "#1c503a"
)

# Pathway colors
pathway_colors <- setNames(
  c("#ff6347",  "#ff6347",  "#ff6347",  "#ff6347",  "#ff6347",  "#ff6347",  "#ff6347",  # immune (7)
    "#99c1df",  # neuron projection (1)
    "#4682b4", "#4682b4", "#4682b4", "#4682b4", "#4682b4","#4682b4",  # synaptic signaling (6)
    "#20B2AA", "#20B2AA", "#20B2AA"),  # myelination (3)
  unique(data$Pathways)
)

# Genotype colors
genotype_colors <- c(
  "Ogg1" = rgb(84/255, 60/255, 108/255, alpha = 0.75),
  "APP" = rgb(5/255, 100/255, 118/255, alpha = 0.75),
  "AO" = rgb(244/255, 122/255, 11/255, alpha = 0.75)
)

# Sector colors
grid_colors <- c(category_colors, pathway_colors)

# Set category order
desired_order <- c("RNA", "ATAC", "Microglia", "Astrocytes", "Oligodendrocytes")
data$Category <- factor(data$Category, levels = desired_order)

# Links between categories and pathways
links <- data.frame(
  from = data$Category,
  to = data$Pathways,
  value = data$log10.P.value,
  genotype = data$Genotype
)

# Link colors by genotype
link_colors <- sapply(links$genotype, function(g) genotype_colors[g])

# Open SVG device
svg(filename = "Chord_diagrams_all.svg", width = 14, height = 14)

# Create chord diagram
circos.clear()
circos.par(clock.wise = TRUE)
chordDiagram(
  x = links,
  grid.col = grid_colors,
  transparency = 0.4,
  annotationTrack = "grid",
  annotationTrackHeight = c(0.05, 0.1),
  preAllocateTracks = 1,
  col = link_colors,
  link.lty = 1,
  link.border = NA,
  link.sort = TRUE,
  link.lwd = links$value
)

# Add sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(
    CELL_META$xcenter, CELL_META$ylim[1],
    sector.name,
    facing = "clockwise",
    niceFacing = TRUE,
    adj = c(0, 0.5))
}, bg.border = NA)

# Close device
dev.off()
