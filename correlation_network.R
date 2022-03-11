NatCorrNetwork <- function(target_list = "target_list.csv",
                           reg_list = "reg_list.csv",
                           data_matrix = "phospho_norm.csv",
                           out_pref = "",
                           pears = 0.5,
                           spear = 0.6) {
  #correlation network
  #read in regulator and target lists
  target_genes = read.csv(target_list,stringsAsFactors=FALSE)
  reg_genes = read.csv(reg_list,stringsAsFactors=FALSE)
  
  #read in the data matrix
  reg_data = read.csv(data_matrix,row.names=1)
  
  #get data for regulators and targets
  myregdata = reg_data[row.names(reg_data)%in%reg_genes[,1],]
  mytargetdata = reg_data[row.names(reg_data)%in%target_genes[,1],]
  
  #for each regulator, calculate spearman and pearson correlation with each target
  #if spearman>0.6 OR pearson>0.5 we keep the edge (from Walley PNAS paper)
  networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Pearson=double(), Spearman=double(),
                            stringsAsFactors=FALSE)
  row = 1
  for (i in 1:(dim(myregdata)[1])){
    print(i)
    for (j in 1:dim(mytargetdata)[1]){
      #skip correlation for the same gene
      if (row.names(myregdata)[i]==row.names(mytargetdata)[j]){
        next
      }
      reg = myregdata[i,]
      tar = mytargetdata[j,]
      pearson = cor(t(reg),t(tar),method="pearson")
      spearman = cor(t(reg),t(tar),method="spearman")
      if (pearson>=pears || spearman>=spear){
        networktable[row,]=c(row.names(myregdata)[i],"regulates",row.names(mytargetdata)[j],pearson,spearman)
        row=row+1
      }
    }
  }
  
  #if there are . in the networktable, replace them with - (for phospho)
  if (length(grepl('\\.',networktable$Regulator))>0){
    networktable$Regulator <- sapply(networktable$Regulator, function(x) gsub("\\.", "-", x))
  }
  
  #write table to file
  dir.create("output")
  write.csv(networktable, paste("output/", out_pref, "_correlation_network.csv", sep = ""),row.names=FALSE, quote = FALSE)
}
myTargets <- list.files(path = "input_tables/", pattern = "Target_list",full.names = T)
myRegulators <- list.files(path = "input_tables/", pattern = "Regulator", full.names = T)
myMatrix <- list.files(path = "input_tables/",pattern = "matrix", full.names = T)
for (f in seq(myTargets)) {
  NatCorrNetwork(target_list = myTargets[f], reg_list = myRegulators[f],data_matrix = myMatrix[f],out_pref = f)
}
target_list = "formatted_tables_kinase/TargetList_all_kinase.csv"
reg_list = "formatted_tables_kinase/RegulatorList_all_DEkinases.csv"
data_matrix = "formatted_tables_kinase/Matrix_kinase_all.csv"
out_pref = "allDE"
pears = 0.5
spear = 0.6
NatCorrNetwork(target_list = target_list,reg_list = reg_list, data_matrix = data_matrix, out_pref = out_pref)
