library(tidyverse)
library(UpSetR)
library(readxl)

##Let's make separate UpSet plots for each data type. First start with protein data
#First read all sheets in your excel file
sets <- excel_sheets(path = "Prot_DE_up_down_table.xlsx")

#create a list with each separate dataset
fullUpDnDE <- list()
for (i in 1:length(sets)) {
  fullUpDnDE[[i]] <- read_excel(path = "Prot_DE_up_down_table.xlsx", sheet = i)
}
#join all the datasets into one binary matrix
fullUpDnDE <- reduce(fullUpDnDE, full_join, by = "Accession")
fullUpDnDE[is.na(fullUpDnDE)] <- 0
fullUpDnDE <- as.data.frame(fullUpDnDE)
fullUpDnDE[2:7] <- lapply(X = fullUpDnDE[2:7], function(x) as.integer(x))
colnames(fullUpDnDE) <- c("Accession", sets)
fullUpDnDE$Accession <- as.factor(fullUpDnDE$Accession)

#create UpSet plots
setEPS()
postscript("proteinUpSet.eps", width = 25, height = 15)
upset(fullUpDnDE, sets = c(colnames(select(fullUpDnDE,matches("prot")))), point.size = 8, line.size = 3.5, keep.order = T,
      mainbar.y.label = "Shared proteins", sets.x.label = "Total DE proteins",
      mb.ratio = c(0.55,0.45), text.scale = c(5,3.5,3,3.5,4.5,3), number.angles = 30)
dev.off()

svg("proteinUpSet.svg", width = 25, height = 15)
upset(fullUpDnDE, sets = c(colnames(select(fullUpDnDE,matches("prot")))), point.size = 8, line.size = 3.5, keep.order = T,
      mainbar.y.label = "Shared proteins", sets.x.label = "Total DE proteins",
      mb.ratio = c(0.55,0.45), text.scale = c(5,3.5,3,3.5,4.5,3), number.angles = 30)
dev.off()

##Now make UpSet of transcript data
#First read all sheets in your excel file
sets <- excel_sheets(path = "Transc_DE_up_down_table.xlsx")

#create a list with each separate dataset
fullUpDnDE <- list()
for (i in 1:length(sets)) {
  fullUpDnDE[[i]] <- read_excel(path = "Transc_DE_up_down_table.xlsx", sheet = i)
}
#join all the datasets into one binary matrix
fullUpDnDE <- reduce(fullUpDnDE, full_join, by = "Accession")
fullUpDnDE[is.na(fullUpDnDE)] <- 0
fullUpDnDE <- as.data.frame(fullUpDnDE)
fullUpDnDE[2:7] <- lapply(X = fullUpDnDE[2:7], function(x) as.integer(x))
colnames(fullUpDnDE) <- c("Accession", sets)
fullUpDnDE$Accession <- as.factor(fullUpDnDE$Accession)

setEPS()
postscript("upset_transcripts.eps", width = 25, height = 15)
upset(fullUpDnDE, sets = c(colnames(select(fullUpDnDE,matches("transc")))), point.size = 8, line.size = 3.5, keep.order = T,
      mainbar.y.label = "Shared transcripts", sets.x.label = "Total DE transcripts",
      mb.ratio = c(0.55,0.45), text.scale = c(5,3.5,3,3.5,4.5,3), number.angles = 30)
dev.off()

pdf(file = "upset_transcripts.pdf", width = 25, height = 15, fonts = "sans")
upset(fullUpDnDE, sets = c(colnames(select(fullUpDnDE,matches("transc")))), point.size = 8, line.size = 3.5, keep.order = T,
      mainbar.y.label = "Shared transcripts", sets.x.label = "Total DE transcripts",
      mb.ratio = c(0.55,0.45), text.scale = c(5,3.5,3,3.5,4.5,3), number.angles = 30)
dev.off()

svg("transcriptUpSet.svg", width = 25, height = 15)
upset(fullUpDnDE, sets = c(colnames(select(fullUpDnDE,matches("transc")))), point.size = 8, line.size = 3.5, keep.order = T,
      mainbar.y.label = "Shared transcripts", sets.x.label = "Total DE transcripts",
      mb.ratio = c(0.55,0.45), text.scale = c(5,3.5,3,3.5,4.5,3), number.angles = 30)
dev.off()

#Finally, we make UpSet for phosphosite data
#First read all sheets in your excel file
sets <- excel_sheets(path = "phosite_DE_up_down_table.xlsx")

#create a list with each separate dataset
fullUpDnDE <- list()
for (i in 1:length(sets)) {
  fullUpDnDE[[i]] <- read_excel(path = "phosite_DE_up_down_table.xlsx", sheet = i)
}
#join all the datasets into one binary matrix
fullUpDnDE <- reduce(fullUpDnDE, full_join, by = "UID")
fullUpDnDE[is.na(fullUpDnDE)] <- 0
fullUpDnDE <- as.data.frame(fullUpDnDE)
fullUpDnDE[2:7] <- lapply(X = fullUpDnDE[2:7], function(x) as.integer(x))
colnames(fullUpDnDE) <- c("UID", sets)
fullUpDnDE$UID <- as.factor(fullUpDnDE$UID)

setEPS()
postscript("phosphoUpSet.eps", width = 25, height = 15)
upset(fullUpDnDE, sets = c(colnames(select(fullUpDnDE,matches("phos")))), point.size = 8, line.size = 3.5, keep.order = T,
      mainbar.y.label = "Shared phosphosites", sets.x.label = "Total DE phosphosites",
      mb.ratio = c(0.55,0.45), text.scale = c(5,3.5,3,3.5,4.5,3), number.angles = 30)
dev.off()

svg("phosphoUpSet.svg", width = 25, height = 15)
upset(fullUpDnDE, sets = c(colnames(select(fullUpDnDE,matches("phos")))), point.size = 8, line.size = 3.5, keep.order = T,
      mainbar.y.label = "Shared phosphosites", sets.x.label = "Total DE phosphosites",
      mb.ratio = c(0.55,0.45), text.scale = c(5,3.5,3,3.5,4.5,3), number.angles = 30)
dev.off()