library(influential)
library(readxl)
library(tidyverse)
for (i in list.files(pattern = "*.sif$",full.names = T)) {
  #we read the network data
  myData <- read_tsv(file = i)
  #then we make the network reconstruction
  my_graph <- graph_from_data_frame(d=myData[c(1,3)],directed = T)
  GraphVertices <- V(my_graph)
  #let's calculate centrality measurements
  list(DC = as.data.frame(degree(my_graph, v = GraphVertices, normalized = FALSE)),
       OD = as.data.frame(degree(my_graph, v = GraphVertices, normalized = FALSE, mode = "out")),
       ID = as.data.frame(degree(my_graph, v = GraphVertices, normalized = FALSE, mode = "in")),
       BC = as.data.frame(betweenness(my_graph, v = GraphVertices, normalized = FALSE, directed = T)),
       NC = as.data.frame(neighborhood.connectivity(my_graph, vertices = GraphVertices)),
       LH_index = as.data.frame(lh_index(my_graph, vertices = GraphVertices)),
       CI = as.data.frame(collective.influence(my_graph, vertices = GraphVertices)),
       CR = as.data.frame(clusterRank(my_graph, vids = GraphVertices, directed = T,loops = T))) %>%
    lapply(function(x) rownames_to_column(x, var = "GeneID")) %>%
    reduce(left_join, by = "GeneID") %>%
    set_names(c("GeneID","DC","OD","ID","BC","NC","LH_index","CI","CR")) -> centrality
  
  as.data.frame(ivi(graph = my_graph, vertices = GraphVertices, 
      directed = TRUE, mode = "all")) %>%
  
#  as.data.frame(ivi.from.indices(DC = centrality$DC,
#                   CR = centrality$CR,
#                   NC = centrality$NC,
#                   LH_index = centrality$LH_index,
#                   BC = centrality$BC,
#                   CI = centrality$CI)) %>%
    rownames_to_column( var = "GeneID") %>%
    rename(IVI = 2) %>%
    inner_join(centrality, by = "GeneID") %>%
    write_csv(path = paste0(i,"_IVI.csv")) -> myIVIscore_table
  
}

myIVIfiles <- list.files(pattern = "*_IVI.csv$", full.names = T)
myGeneAnnot <- read_csv(file = "../Ath_symbol_file.csv")

for (i in myIVIfiles) {
  read_csv(file = i) %>%
    arrange(desc(IVI)) %>%
    inner_join(myGeneAnnot, by = c("GeneID" = "AGI")) %>%
    write_csv(path = i)
}
