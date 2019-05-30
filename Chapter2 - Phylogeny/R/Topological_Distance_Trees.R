###############################################################################
######################## Tree Topology Comparison #############################
################################## TUTORIAL ###################################
###############################################################################

###############################################################################
######################### JERONYMO DALAPICOLLA, 2019 ##########################
###################### LABORATORIO DE MAMIFEROS - LAMA ########################
####################### ESALQ - UNIVERSIDADE DE SAO PAULO #####################
###############################################################################

##1. FOR OTHER STEPS, PLEASE GO TO: https://github.com/jdalapicolla/

##2.INPUTS FOR THIS ANALYSIS:
#A. ".tre" TREE FILES FROM CONCATENATED DATA IN RAxML.
#B. ".nex" TREE FILES FROM COALESCENT ANALYSES IN SVDquartets.

##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("tcltk","pacman", "ape")) #tcltk is for Linux!
#load the packages
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("tcltk", "ape")
#OR
install.packages("ape") # Install ape package if you don't have it already
library(ape) # Load the ape package


##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())


##5. RAxML/Concatenated trees:
#read trees do you want to compare as unroot trees:
pc30 = unroot(read.tree("RAxML_bipartitions.proe_73_30pc"))
pc40 = unroot(read.tree("RAxML_bipartitions.proe_73_40pc"))
pc50 = unroot(read.tree("RAxML_bipartitions.proe_73_50pc"))
pc60 = unroot(read.tree("RAxML_bipartitions.proe_73_60pc"))
pc70 = unroot(read.tree("RAxML_bipartitions.proe_73_70pc"))
pc80 = unroot(read.tree("RAxML_bipartitions.proe_73_80pc"))
pc90 = unroot(read.tree("RAxML_bipartitions.proe_73_90pc"))

#create a object with all trees for comparison
phylogenies = c(pc30, pc40, pc50, pc60, pc70, pc80, pc90) 

#calculate the topological distance:
distance = round(as.matrix(dist.topo(phylogenies, method = "PH85")),3)
rownames(distance) = c("30%", "40%", "50%", "60%", "70%", "80%", "90%")
colnames(distance) = c("30%", "40%", "50%", "60%", "70%", "80%", "90%")

##6. Calculate bootstrap indices with RAxML/Concatenated trees:
loop = 0
tab = matrix(NA, 7, 4) # you will save the results in this "tab" matrix. Number of columns refers to arguments in cbind functions (nodes, mean, sd_boots, good); rows the trees number, in my case 7 trees. 

for (i in phylogenies) {
  loop = loop +1
  nodes = length(i$edge.length)
  z = as.numeric(i$node.label)
  mean = round(mean(z, na.rm = T),2)
  z1 = z < 50 #bootstrap value to be considered as "good"
  good = length(z1[z1 == TRUE])
  sd_boots = round(sd(as.numeric(i$node.label), na.rm=T),2)
  tab[loop,] = cbind(nodes, mean, sd_boots, good)
}

#rename rows and colunms
rownames(tab) = c("30%", "40%", "50%", "60%", "70%", "80%", "90%")
colnames(tab) = c("Nodes", "Bootstrap (Mean)", "Bootstrap (SD)", "Clades with Low Support (<50%)")

#convert matrix to data frame to save the results
distance = as.data.frame(distance)
tab = as.data.frame(tab)

#save the results
write.csv(distance, "distance_topology_concatenaded.csv", row.names = T)
write.csv(tab, "bootstrap_concatenaded.csv", row.names = T)


##7. SVDquartets/Coalescent trees:
#read trees as unroot trees:
pc30_svd = unroot(read.nexus("proe_73.u.snps_30pc.tre"))
pc40_svd = unroot(read.nexus("proe_73.u.snps_40pc.tre"))
pc50_svd = unroot(read.nexus("proe_73.u.snps_50pc.tre"))
pc60_svd = unroot(read.nexus("proe_73.u.snps_60pc.tre"))
pc70_svd = unroot(read.nexus("proe_73.u.snps_70pc.tre"))
pc80_svd = unroot(read.nexus("proe_73.u.snps_80pc.tre"))
pc90_svd = unroot(read.nexus("proe_73.u.snps_90pc.tre"))

phylogenies_svd = c(pc30_svd, pc40_svd, pc50_svd, pc60_svd, pc70_svd, pc80_svd, pc90_svd) 

distance_svd = round(as.matrix(dist.topo(phylogenies_svd, method = "PH85")),3)
rownames(distance_svd) = c("30%", "40%", "50%", "60%", "70%", "80%", "90%")
colnames(distance_svd) = c("30%", "40%", "50%", "60%", "70%", "80%", "90%")


##6. Calculate bootstrap indices with RAxML/Concatenated trees:
loop = 0
tab_svd = matrix(NA, 7, 4)

for (i in phylogenies_svd) {
  loop = loop +1
  nodes = length(i$edge.length)
  z = as.numeric(i$edge.length)
  mean = round(mean(z, na.rm = T),2)
  z1 = z < 50 #bootstrap value to be considered as "good"
  good = length(z1[z1 == TRUE])
  sd_boots = round(sd(as.numeric(i$edge.length), na.rm=T),2)
  tab_svd[loop,] = cbind(nodes, mean, sd_boots, good)
}

rownames(tab_svd) = c("30%", "40%", "50%", "60%", "70%", "80%", "90%")
colnames(tab_svd) = c("Nodes", "Bootstrap (Mean)", "Bootstrap (SD)", "Clades with Low Support (<50%)")

distance_svd = as.data.frame(distance_svd)
tab_svd = as.data.frame(tab_svd)

write.csv(distance_svd, "distance_topology_coalescent.csv", row.names = T)
write.csv(tab_svd, "bootstrap_coalescent.csv", row.names = T)