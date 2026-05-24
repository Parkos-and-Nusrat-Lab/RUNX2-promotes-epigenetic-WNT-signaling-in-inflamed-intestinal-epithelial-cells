# Doing motifs heatmap
library(tidyverse)
library(ggrepel)
library(grid)

# File obtained from the function plotEnrichHeatmap() with returnMatrix = TRUE

mat <- read_rds("C:/Users/ricsilva/Downloads/heatmap_motifs.rds")

# Transforming the dataframe

df <- mat %>%  reshape2::melt()


# Adding colnames
colnames(df) <- c("Cell", "TF","value")
df$TF <- str_remove(df$TF, "_.*")
# Renaming cell anno


new_anno <- list(
  Differentiated_2D = "Differentiated (2D)",
  Transitional_2D = "Transitional (2D)",
  Undifferentiated_2D = "Undifferentiated (2D)",
  Differentiated_3D = "Late G1 (3D)",
  Progenitor_2D = "BMI+Stem (2D)",
  UC_specific = "UC-specific (3D)",
  M_like = "M-like (2D)",
  Intermediate_2D = "Intermediate (2D)",
  HES1_progenitor = "Hes+Stem (2D)",
  Progenitor_3D = "LGR5+Stem (3D)",
  Intermediate_3D = "Early G1 (3D)"
)

cell_anno <- c()
for (x in df$Cell){
  cell_anno <- c(cell_anno,new_anno[[x]])
  }


df$Cell <- cell_anno


# Reordering the values based on celltypes

celltypes <- rev(c("Hes+Stem (2D)",
                   "Early G1 (3D)",
                   "LGR5+Stem (3D)",
                   "Undifferentiated (2D)",
                   "Intermediate (2D)",
                   "M-like (2D)",
                   "UC-specific (3D)",
                   "Late G1 (3D)",
                   "BMI+Stem (2D)",
                   "Transitional (2D)",
                   "Differentiated (2D)"))

# Loop to obtain the new order in which to place the Transcription Factor motifs
tf_order <- c()
for(cell in celltypes){
  tf_names <- df[df$Cell == cell & df$value == 100, "TF"] %>% as.vector()

  tf_order <- c(tf_order, tf_names)
}

# Adding new order by making the variable a factor
df$TF <- factor(df$TF, levels = tf_order)

# Ensuring adaquate order of cell populations
df$Cell <- factor(df$Cell, levels = c("Hes+Stem (2D)",
                                        "Early G1 (3D)",
                                        "LGR5+Stem (3D)",
                                        "Undifferentiated (2D)",
                                        "Intermediate (2D)",
                                        "M-like (2D)",
                                        "UC-specific (3D)",
                                        "Late G1 (3D)",
                                        "BMI+Stem (2D)",
                                        "Transitional (2D)",
                                        "Differentiated (2D)"))

# Custom color vector
custom_colors <- c(
  "grey93", "#E0F0FF", "#C0E0FF", "#A0D0FF", "#80C0FF",
  "#60A0FF", "#4080FF", "#4060E0", "#4040C0", "#5020A0",
  "#601080", "#700060", "#800040", "#700030", "#600020",
  "#500010", "#400008", "#300004", "#200002", "#000000"
)


# Removing all the names but the runx names for labeling purpose
# Due to overlap labels are removed and added outside of R
uc_motifs <- c("RUNX1", "RUNX2", "RUNX3", "TCF3", "TCF4", "TCF12")


# Plotting
# Due to overlap Runx and TCF x ticks were re-placed with text and line Grobs
# To visualize the position based on the df un-comment the sections
# in the plot code below

font_size = 7

gg <- ggplot(df, aes(x = TF, y = Cell, fill = value)) +
  geom_tile(color = "darkgrey") +
  scale_fill_gradientn(
    colors = custom_colors,
    limits = c(0, 100),
    breaks = c(0, 100),
    name = "Normalized Enrichment",
    guide = guide_colourbar(title.position = "top", title.hjust = 0.5)
  ) +
  scale_y_discrete(position = "right")+
  #scale_x_discrete(labels = function(x) ifelse(x %in% uc_motifs, x, "")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "bold", size = 9),
    axis.text.x = element_blank(), # Comment this one if uncommenting the one below
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.35, face = "bold", size = 5),
    legend.title = element_text(size = 10, vjust = 1.8, hjust = 5, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(t = 5, r = 5, b = 30, l = 5)
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  coord_cartesian(clip = "off") +
  # Comment everything under here to remove manual re-labeling
  annotation_custom(grob = grid::textGrob("RUNX1", # Adding text
                                          x = 0.595,
                                          y = -0.08,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = font_size))) +
  annotation_custom(grob = grid::linesGrob(x = c(0.586, 0.5959), # Adding line
                                           y  =c(0.01, -0.009),
                                           gp = gpar(lwd= unit(1, "pt"),  col = "black"))) +
  annotation_custom(grob = grid::textGrob("RUNX2", # Adding text
                                          x = 0.65,
                                          y = -0.08,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = font_size))) +
  annotation_custom(grob = grid::linesGrob(x = 0.6495, # Adding line
                                           y  =c(-0.009, 0.01),
                                           gp = gpar(lwd= unit(1, "pt"),  col = "black"))) +
  annotation_custom(grob = grid::textGrob("RUNX3", # Adding text
                                          x = 0.57,
                                          y = -0.08,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = font_size))) +
  annotation_custom(grob = grid::linesGrob(x = c(0.57,0.5818), # Adding line
                                           y  =c(-0.009, 0.01),
                                           gp = gpar(lwd= unit(1, "pt"),  col = "black"))) +
  annotation_custom(grob = grid::textGrob("TCF3", # Adding text
                                          x = 0.6933,
                                          y = -0.065,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = font_size))) +
  annotation_custom(grob = grid::linesGrob(x = 0.6933, # Adding line
                                           y  =c(-0.009, 0.01),
                                           gp = gpar(lwd= unit(1, "pt"),  col = "black"))) +
  annotation_custom(grob = grid::textGrob("TCF12", # Adding text
                                          x = 0.6155,
                                          y = -0.075,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = font_size))) +
  annotation_custom(grob = grid::linesGrob(x = c(0.6155,0.6255), # Adding line
                                           y  =c(-0.009, 0.01),
                                           gp = gpar(lwd= unit(1, "pt"),  col = "black"))) +
  annotation_custom(grob = grid::textGrob("TCF4", # Adding text
                                          x = 0.6325,
                                          y = -0.065,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = font_size))) +
  annotation_custom(grob = grid::linesGrob(x = 0.6325, # Adding line
                                           y  =c(-0.009, 0.01),
                                           gp = gpar(lwd= unit(1, "pt"),  col = "black")))



# Saving plot
ggsave("Figure3E.png",
       plot = gg,
       device = "png",
       dpi = 600,
       height = 4,
       width = 6,
       units = "in",
       bg = "white")





