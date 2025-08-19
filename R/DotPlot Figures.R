
# Loading Librarires
library(Seurat)
library(tidyverse)


# Loading File

seu <- readRDS("Seurat_file.rds")


# Updating levels to ensure celltype placement in figures

seu@active.ident <- factor(seu@active.ident, levels = c("LGR5+Stem_3D",
                                               "Early_G1_3D",
                                               "Late_G1_3D",
                                               "UC-specific_3D",
                                               "Hes+Stem_2D",
                                               "BMI+Stem_2D",
                                               "Undifferentiated_2D",
                                               "Intermediate_2D",
                                               "Transitional_2D",
                                               "Differentiated_2D",
                                               "M-like_2D")
                        )


#Figure 2 E and F


p <- DotPlot(seu, features =c("HNF4A","REG4","GUCA2A","KRT18","KRT20","GPA33","CDH17"),
        dot.scale = 12, assay = "RNA")  +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") + ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

cowplot::ggsave2(filename = "Fig_2E.png", dpi = 300, plot=p, device = "png", bg ="white")





p <- DotPlot(seu, features = c("LGR5","LRIG1","SMOC2","EPHB2","BMI1","ASCL2","OLFM4", "MSI1", "HOPX"), dot.scale = 12)  +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") + ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


cowplot::ggsave2(filename = "Fig_2F.png", dpi = 300, plot=p, device = "png", bg = "white")




# Supplementary Figure 2

genes <- c("EXT2","EHBP1","CKAP5","MIS18BP1","TTK","LSAMP","ANLN","CDC25C","KIF14","SQSTM1","CASC19","UST","NFAT5","FUT8",
           "ROCK2","PTEN","ETV6","MGRN1","LINC01088","SLC11A2","AC002460.2","DDB2","EGLN3","MGAT5","EEF1A1","MT-ND4L","MT-ND2",
           "LGALS3BP","RCAN3","ANKRD37","GPBP1","ERO1A","HSPH1","FAM13A","PDK1","ELOVL5","HIST1H2AC","PGK1","ZC3H3","NR1D1",
           "KLF10","GADD45B","ZBTB10","HES1","JUN","TXNIP","IER2","GDAP1","EGR1","GDF15","TENT5A","IGFBP3","LYZ","CXCL8",
           "GAPDH","NRIP1","BRAP","HSPA1A","DENND4A","GSDMB","TMEM170A","RBM23","RPL36A","TMEM120B",
           "RPS26","DST","SH3RF1","SLC20A2","PPARG","SLC49A4","NRXN3","MECOM","NTN4","SGMS2","FLNB","PHGR1","FTH1","SLC26A3",
           "CKB","TSPAN8","CA1","COBLL1","PI3","SELENOP","SLC5A3","FABP1","TFF1","B2M","TFF2","REG4","AGR2","S100A6","LGALS4",
           "GUCA2A","S100A4","AC083837.1","CCL20","SORBS2","MMP7","MYOF","ZBTB20","TRIO","ITGBL1","ADAMTS9","HSP90B1")

p <-DotPlot(seu, features = genes, dot.scale = 5)  +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") + ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


cowplot::ggsave2(filename = "Suppl_Fig_2.png", dpi = 300, plot=p, device = "png", width = 24, height = 20, units = "in", bg = "white")


# Supplementary Figure 3

genes <- c("CCL20","IFNA1","IFNB1","IFNG","IL1B","IL17A","IL17F","IL25","IL33","CXCL8","CCL3","TSLP","IL6","TNF")

p <-DotPlot(seu, features = genes, dot.scale = 11)  +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") + ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

cowplot::ggsave2(filename = "Suppl_Fig_3_A.png", dpi = 300, plot=p, device = "png", bg = "white")




genes <- c("CCL20","IFNA1","IFNB1","IFNG","IL1B","IL17A","IL17F","IL25","IL33","CXCL8","CCL3","TSLP","IL6","TNF")

p <-DotPlot(seu, features = genes, dot.scale = 11)  +
  ggtitle("Cytokines") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12),
        plot.title = element_text(vjust = 0.5, size = 15, hjust = 0.5)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") +
  ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  coord_cartesian(clip = "off") + # Removing the clipping so I can add things outside of the plot
  annotation_custom(grob = grid::linesGrob(x = -0.325, # Adding the first vertical line
                                           y  =c(0.38, 0.99),
                                           gp = gpar(lwd= unit(2, "pt")))) +
  annotation_custom(grob = grid::linesGrob(x = -0.325, # Adding the second vertical line
                                           y  =c(0.01, 0.35),
                                           gp = gpar(lwd= unit(2, "pt")))) +
  annotation_custom(grob = grid::textGrob("2D", # Adding the text
                                          x = -0.355,
                                          y = 0.67,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = 15))) +
  annotation_custom(grob = grid::textGrob("3D", # Adding the second text
                                          x =  -0.355,
                                          y = 0.2,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = 15)))

p


cowplot::ggsave2(filename = "Suppl_Fig_3_A.png", dpi = 300, plot=p, device = "png", bg = "white")





genes <- c("CCL20","TNFAIP2","ANXA5","MARCKSL1","GP2","CCL15")




p<-DotPlot(seu, features = genes, dot.scale = 11)  +
  ggtitle("M-cell markers") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12),
        plot.title = element_text(vjust = 0.5, size = 15, hjust = 0.5),
        axis.text.y = element_blank()) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") +
  ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
coord_cartesian(clip = "off") + # Removing the clipping so I can add things outside of the plot
  annotation_custom(grob = grid::linesGrob(x = -0.022, # Adding the first vertical line
                                           y  =c(0.38, 0.99),
                                           gp = gpar(lwd= unit(2, "pt")))) +
  annotation_custom(grob = grid::linesGrob(x = -0.022, # Adding the second vertical line
                                           y  =c(0.01, 0.35),
                                           gp = gpar(lwd= unit(2, "pt")))) +
  annotation_custom(grob = grid::textGrob("2D", # Adding the text
                                          x = -0.040,
                                          y = 0.67,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = 15))) +
  annotation_custom(grob = grid::textGrob("3D", # Adding the second text
                                          x =  -0.040,
                                          y = 0.2,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = 15)))



cowplot::ggsave2(filename = "Suppl_Fig_3_B.png", dpi = 300, plot=p, device = "png", bg = "white")




genes <- c("TEAD4","CLU","YAP1","ANXA5","ANXA1","ANXA3","TNFRSF12A","LAMC2","CD44","LY6D","ANXA8","SPRR1A")

p <-DotPlot(seu, features = genes, dot.scale = 11)  +
  ggtitle("Fetal regenerative markers") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12),
        plot.title = element_text(vjust = 0.5, size = 15, hjust = 0.5)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") +
  ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  coord_cartesian(clip = "off") + # Removing the clipping so I can add things outside of the plot
  annotation_custom(grob = grid::linesGrob(x = -0.325, # Adding the first vertical line
                                           y  =c(0.38, 0.99),
                                           gp = gpar(lwd= unit(2, "pt")))) +
  annotation_custom(grob = grid::linesGrob(x = -0.325, # Adding the second vertical line
                                           y  =c(0.01, 0.35),
                                           gp = gpar(lwd= unit(2, "pt")))) +
  annotation_custom(grob = grid::textGrob("2D", # Adding the text
                                          x = -0.355,
                                          y = 0.67,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = 15))) +
  annotation_custom(grob = grid::textGrob("3D", # Adding the second text
                                          x =  -0.355,
                                          y = 0.2,
                                          rot = 90,
                                          gp = gpar(fontface = "bold",
                                                    fontsize = 15)))

p


cowplot::ggsave2(filename = "Suppl_Fig_3_C.png", dpi = 300, plot=p, device = "png", bg = "white")


