#MLM script
#first we need to read in the file and get the data for each line
library(openxlsx)
library(MASS)
library(pals)

data = read.xlsx('measurements.xlsx')
colnames(data) <- c("LineNumber", "Scan","Plate","Treatment","Rep","Seed","Hypocotyl.Length")
data <- data[1:7]

lines = unique(data$LineNumber)
lines = lines[!lines%in%"WT"] #get rid of WT
treats = unique(data$Treatment)

WT = data[data$LineNumber=="WT",] #get WT data for all comps

#run model for each line and each treatment
#compare to WT and treatment 0
for (i in 1:length(lines)){
  for(j in seq(2,length(treats),by=2)){
    #get data for this line and treatment
    mydata = data[data$LineNumber==lines[i],]
    mydata = mydata[mydata$Treatment==treats[j-1]| mydata$Treatment==treats[j],]
    testdata = rbind(WT[WT$Treatment==treats[j-1] | WT$Treatment==treats[j],],mydata)
    
    #fit the model
    #primary effects genotype, treatment, and interaction
    #random effects plate, replicate, and experiment
    mymodel = glmmPQL(Hypocotyl.Length ~ LineNumber + Treatment + LineNumber*Treatment, ~1 | Rep/Plate/Seed, 
                      family = gaussian,data=testdata, verbose=FALSE)
    myresults = summary(mymodel)
    #save pvalues for later
    pvals = myresults$tTable
    if (i ==1 & j==2){
      genop = data.frame(pvals[2,5])
      treatp = data.frame(pvals[3,5])
      genoxtreatp = data.frame(pvals[4,5])
      myrownames = paste(lines[i],treats[j],sep="_")
    } else{
      genop = data.frame(genop,pvals[2,5])
      treatp = data.frame(treatp,pvals[3,5])
      genoxtreatp = data.frame(genoxtreatp,pvals[4,5])
      myrownames = data.frame(myrownames,paste(lines[i],treats[j],sep="_"))
    }
  }
}

#get corrected pvals
allpvals = data.frame(t(genop),t(treatp),t(genoxtreatp))
rownames(allpvals)=t(myrownames)
colnames(allpvals) = c("Genotype","Treatment","Genotype.x.Treatment")
for (i in 1:3){
  if (i==1){
    myqvals = data.frame(p.adjust(allpvals[,i],method="fdr"))
  } else{
    myqvals = data.frame(myqvals,p.adjust(allpvals[,i],method="fdr"))
  }
}
colnames(myqvals) = c("Genotype_FDR","Treatment_FDR","Genotype.x.Treatment_FDR")
finalvals = data.frame(allpvals,myqvals)
write.csv(finalvals,'GLM_results_BL_CMS01_5_FDR.csv')
