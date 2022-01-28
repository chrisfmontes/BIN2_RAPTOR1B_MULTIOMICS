library(readxl)
library(tidyverse)
library(ComplexHeatmap)

#The following libraries and code lines are useful to create colorblind-friendly palettes
library(RColorBrewer)
library(viridis)
library(circlize)
col_fun <- colorRamp2(c(-0.5, 0, 5),plasma(3))

# For our heatmaps, we first read out value tables

### First, hypocotyl response to BL
phenoBL <- read_excel("../../DE_tables/BL_GLM_results.xlsx", sheet = 12)
# We make a table with hypocotyl lenght by genotype
phenoBL %>%
  filter(data.Treatment == "BL0") %>%
  select(paper.name,`mut/wt logFC`) %>%
  column_to_rownames(var = "paper.name") %>% as.matrix ->BLgenoHypoDat
# We make a table with hypocotyl response to BL treatment (log2FC of BL/mock hypocotyl length)
phenoBL %>%
  filter(data.Treatment == "BL100") %>%
  select(paper.name, `mut/wt logFC`) %>%
  column_to_rownames(var = "paper.name") %>%
  as.matrix -> BLHypoDat
#We make a table with genotype response to BL treatment (log2FC of mutant hypocotyl response over log2FC of WT hypocotyl response)
phenoBL %>%
  filter(data.Treatment == "BL0") %>%
  mutate(`BL_resp_logFC` = log2(`mut[BL/mock] / WT[BL/mock] FC`)) %>%
  select(paper.name, BL_resp_logFC) %>%
  column_to_rownames(var = "paper.name") %>%
  as.matrix -> BLrespDat

# We now build the heatmaps
BLgenoHypo <- Heatmap(matrix = t(BLgenoHypoDat),
              name = "Basal Hypocotyl length\n mutant vs WT log2FC",
              border = T,
#              col = viridis(100),
#              col = c("#FAFE8D", "black", "#E17CCA"),
              col = inferno(100),
              row_names_side = "left",
              row_labels = "mut/WT mock",
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_gp = gpar(fontface = "italic", fontsize = 10),
              heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
              cell_fun = function(j, i, x, y, width, h, fill) {
                grid.text(ifelse(na.omit(phenoBL[j+43,13] <= 0.054), yes = "*", no = ""), x = x, y-h*0.15,
                          gp = gpar(fontsize = 19, fontface = "bold", col = "black"))
        })

#draw(BLgenoHypo)

BLHypo <- Heatmap(matrix = t(BLHypoDat),
                      name = "Hypocotyl length on 100nM BL\n mutant vs WT log2FC",
                      border = T,
                      col = viridis(100),
                      #              col = c("#FAFE8D", "black", "#E17CCA"),
                      #              col = inferno(100),
                      row_names_side = "left",
                      row_labels = "mut/wt BL",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      rect_gp = gpar(col = "white", lwd = 1),
                      column_names_gp = gpar(fontface = "italic", fontsize = 10),
                      heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
                      cell_fun = function(j, i, x, y, width, h, fill) {
                        grid.text(ifelse(na.omit(phenoBL[j,13] <= 0.054), yes = "*", no = ""), x = x, y-h*0.15,
                                  gp = gpar(fontsize = 19, fontface = "bold", col = "black"))
                      })
#draw(BLHypo)

BLresp <- Heatmap(matrix = t(BLrespDat),
                  name = "WT-normalized hypocotyl response to BL\nlog2[mut. hypocotyl (BL/mock) FC / WT hypocotyl (BL/mock) FC]",
                  border = T,
                  col = plasma(100),
#                  col = viridis(100),
                  #              col = c("#FAFE8D", "black", "#E17CCA"),
                  #              col = inferno(100),
                  row_names_side = "left",
                  row_labels = "WT-norm. hypocotyl\nresponse to BL",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  rect_gp = gpar(col = "white", lwd = 1),
                  column_names_gp = gpar(fontface = "italic", fontsize = 10),
                  heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
                  cell_fun = function(j, i, x, y, width, h, fill) {
                    grid.text(ifelse(na.omit(phenoBL[j+43,14] <= 0.1), yes = "*", no = ""), x = x, y-h*0.15,
                              gp = gpar(fontsize = 19, fontface = "bold", col = "black"))
                  })

#draw(BLresp)

# We finally create and export the heatmap as a PDF filw

BL <- BLgenoHypo %v% BLHypo %v% BLresp

pdf(file = "BL_response_colors.pdf", width = 16, height = 4, fonts = "sans")
draw(BL,ht_gap = unit(0.7, "cm"), heatmap_legend_side = "bottom")
dev.off()

## We make tables for ATG8e fluorescence (autophagy levels) on protoplast
phenoGFP <- read_excel("../../DE_tables/FP-ATG8e.summary.stats.xlsx", sheet = 1) %>%
  filter(in.figure == "yes") %>%
  mutate(across(c(8:19), as.numeric))

# Table with basal autophagy levels per genotype
phenoGFP %>%
  filter(data.Treatment == "control") %>%
  select(mutantID,`basal log2FC`) %>%
  column_to_rownames(var = "mutantID") %>%
  as.matrix ->GFPdat
#Table with autophagy levels after sucrose starvation
phenoGFP %>%
  filter(data.Treatment == "suc.starv") %>%
  select(mutantID,`basal log2FC`) %>%
  column_to_rownames(var = "mutantID") %>%
  as.matrix -> GFPstressdat
#Table wirh WT-normalized autophagy increase upon sucrose starvation
phenoGFP %>%
  filter(data.Treatment == "control") %>%
  select(mutantID,`effect logFC`) %>%
  column_to_rownames(var = "mutantID") %>%
  as.matrix -> GFPrespDat

#Now make the heatmaps

GFP <- Heatmap(matrix = t(GFPdat),
               name = "% Protoplasts with active autophagy\nlog2([mut/WT] FC)",
#               col = plasma(100),
               col = inferno(100),
               border = T,
               row_names_side = "left",
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               rect_gp = gpar(col = "white", lwd = 1),
               column_names_gp = gpar(fontface = "italic", fontsize = 10),
               heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
               cell_fun = function(j, i, x, y, width, h, fill) {
                 grid.text(ifelse(na.omit(phenoGFP[j,13] <= 0.054), yes = "*", no = ""), x, y-h*0.15,
                           gp = gpar(fontsize = 19, fontface = "bold", col = "black"))
               })

#draw(GFP)
#col_suc <- colorRamp2(c(-2, 0, 2),inferno(3))

GFPsuc <- Heatmap(matrix = t(GFPstressdat),
               name = "% Protoplasts with active autophagy\nafter sucrose starvation\nlog2([mut/WT] FC)",
#               col = plasma(100),
               col = viridis(100),
               border = T,
               row_names_side = "left",
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               rect_gp = gpar(col = "white", lwd = 1),
               column_names_gp = gpar(fontface = "italic", fontsize = 10),
               heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
               cell_fun = function(j, i, x, y, width, h, fill) {
                 grid.text(ifelse(na.omit(phenoGFP[j+30,13] <= 0.054), yes = "*", no = ""), x, y-h*0.15,
                           gp = gpar(fontsize = 19, fontface = "bold", col = "black"))
               })

GFPresp <- Heatmap(matrix = t(GFPrespDat),
                  name = "WT-normalized autophagy increase\nafter sucrose starvation",
                  col = plasma(100),
                  border = T,
                  row_names_side = "left",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  rect_gp = gpar(col = "white", lwd = 1),
                  column_names_gp = gpar(fontface = "italic", fontsize = 10),
                  heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
                  cell_fun = function(j, i, x, y, width, h, fill) {
                    grid.text(ifelse(na.omit(phenoGFP[j,14] <= 0.054), yes = "*", no = ""), x, y-h*0.15,
                              gp = gpar(fontsize = 19, fontface = "bold", col = "black"))
                  })
#draw(GFPsuc)

GFP_list <- GFP %v% GFPsuc %v% GFPresp
draw(GFP_list, ht_gap = unit(0.7, "cm"), heatmap_legend_side = "bottom")

pdf(file = "GFP_autophagy_colors.pdf", width = 16, height = 5, fonts = "sans")
draw(GFP_list, ht_gap = unit(0.7, "cm"), heatmap_legend_side = "bottom")
dev.off()

# Same for MDC staining measurements

phenoMDC <- read_excel("../figure_ideas/table for MDC_spp figure.xlsx", sheet = 1)
phenoMDC %>%
  select(1,10) %>% column_to_rownames(var = "mutant line") %>% as.matrix -> MDCdat

MDC <- Heatmap(matrix = t(MDCdat),
               name = "MDC staining\nAutophagosomes per frame\n[mut/WT] FC",
               col = plasma(100),
               border = T,
               row_names_side = "left",
               row_labels = "MDC",
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               rect_gp = gpar(col = "white", lwd = 1),
               column_names_gp = gpar(fontface = "italic", fontsize = 10),
               heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
               cell_fun = function(j, i, x, y, width, h, fill) {
                 grid.text(ifelse(na.omit(pheno1[j,7] <= 0.05), yes = "*", no = ""), x, y-h*0.15,
                           gp = gpar(fontsize = 18, fontface = "bold", col = "#7F878F"))
               })

pdf(file = "MDC_autophagy.pdf", width = 16, height = 3, fonts = "sans")
draw(MDC, heatmap_legend_side = "bottom")
dev.off()

# Finally, we make the plots for MDC under sucrose starvation

phenoMDCsuc <- read_excel("../figure_ideas/table for MDC_spp figure.xlsx", sheet = 2, skip = 1)
phenoMDCsuc %>%
  filter(treatment == "control") %>%
  select(1,12) %>%
  column_to_rownames(var = "mutant line") %>%
  as.matrix -> MDCbasalDat

phenoMDCsuc %>%
  filter(treatment == "control") %>%
  select(1,13) %>%
  column_to_rownames(var = "mutant line") %>%
  as.matrix -> MDCsucDat

phenoMDCsuc %>%
  filter(treatment == "control") %>%
  select(1,14) %>%
  column_to_rownames(var = "mutant line") %>%
  as.matrix -> MDCeffectDat

MDCbasal <- Heatmap(matrix = t(MDCbasalDat),
                  name = "MDC staining control media [mutant/WT] log2FC",
                  col = rev(plasma(100)),
                  border = T,
                  row_names_side = "left",
                  row_labels = "MDC basal",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  rect_gp = gpar(col = "white", lwd = 1),
                  column_names_gp = gpar(fontface = "italic", fontsize = 10),
                  heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
                  cell_fun = function(j, i, x, y, width, h, fill) {
                    grid.text(ifelse(na.omit(pheno[j,9] <= 0.05), yes = "*", no = ""), x, y-h*0.15,
                              gp = gpar(fontsize = 18, fontface = "bold", col = "#7F878F"))
                  })

MDCsuc <- Heatmap(matrix = t(MDCsucDat),
                  name = "MDC staining under sucrose starvation [mutant/WT] log2FC",
                  col = rev(plasma(100)),
                  border = T,
                  row_names_side = "left",
                  row_labels = "MDC stress",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  rect_gp = gpar(col = "white", lwd = 1),
                  column_names_gp = gpar(fontface = "italic", fontsize = 10),
                  heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
                  cell_fun = function(j, i, x, y, width, h, fill) {
                    grid.text(ifelse(na.omit(pheno[j,9] <= 0.05), yes = "*", no = ""), x, y-h*0.15,
                              gp = gpar(fontsize = 18, fontface = "bold", col = "#7F878F"))
                  })

MDCeffect <- Heatmap(matrix = t(MDCeffectDat),
                  name = "MDC staining (sucrose starvation)\n WT-nomalized\nAutophagosomes per frame\n[suc-/suc+] FC",
                  col = rev(plasma(100)),
                  border = T,
                  row_names_side = "left",
                  row_labels = "MDC effect",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  rect_gp = gpar(col = "white", lwd = 1),
                  column_names_gp = gpar(fontface = "italic", fontsize = 10),
                  heatmap_legend_param = list(legend_width = unit(6, "cm"),direction = "horizontal", title_position = c("topcenter")),
                  cell_fun = function(j, i, x, y, width, h, fill) {
                    grid.text(ifelse(na.omit(pheno[j,11] <= 0.05), yes = "*", no = ""), x, y-h*0.15,
                              gp = gpar(fontsize = 18, fontface = "bold", col = "#7F878F"))
                  })
MDC_list <- MDCbasal %v% MDCsuc %v% MDCeffect

pdf(file = "MDC_starvation_autophagy.pdf", width = 10, height = 4, fonts = "sans")
draw(MDC_list, ht_gap = unit(0.7, "cm"))
dev.off()