###############################################################################
############################# HYPERVOLUME ANALYSIS ############################
################################## TUTORIAL ###################################
###############################################################################

###############################################################################
######################## Based on Mariah Kenney scripts #######################
###################### KNOWLES LAB, UNIVERSITY OF MICHIGAN ####################
###############################################################################

###############################################################################
######################### JERONYMO DALAPICOLLA, 2019 ##########################
###################### LABORATORIO DE MAMIFEROS - LAMA ########################
####################### ESALQ - UNIVERSIDADE DE SAO PAULO #####################
###############################################################################

##1.FOR OTHERS SCRIPTS, PLEASE GO TO: https://github.com/jdalapicolla/

##2.INPUTS FOR THIS ANALYSIS:
#A. MORPHOMETRIC DATA OF THE TARGET SPECIES FOR THE MORPHOLOGICAL HYPERVOLUME

##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("pacman","RStoolbox", "raster", "rgdal", "dismo", "hypervolume","alphahull", "tcltk", "vegan", "adegenet"))
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("RStoolbox", "raster", "rgdal", "dismo", "hypervolume", "alphahull", "tcltk", "vegan", "adegenet")


################################################################################
####################### MORPHOLOGICAL HYPERVOLUME ##############################
################################################################################


##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())
#Any computer
setwd('C:/Users/Johnny/Downloads/bio_10m_bil') # change the path


##5. LOAD THE MORPHOMETRIC DATA:
data=read.csv("WA.csv")
#data without missing data
#do not need to eliminate co-variables
#do not need to eliminate non-normal distribution variables

#select the adults
#data = data[data$age_patton == "8" | data$age_patton == "9" | data$age_patton == "10", ]
summary(data)

#create a file for each species
bre = data[data$species == "BRE", ]
sim = data[data$species == "SIM", ]
ste = data[data$species == "STE", ]

#save tables
write.csv(bre, "BRE_all.csv")
write.csv(sim, "SIM_all.csv")
write.csv(ste, "STE_all.csv")


#minimum number of samples in all species: 32
bre = bre[sample(nrow(bre), 36), ]
sim = sim[sample(nrow(sim), 36), ]
ste = ste[sample(nrow(ste), 36), ]


##6. CHOOSE THE MORPHOLOGICAL VARIABLES MORE IMPORTANTS FOR THE VARIATION WITHOUT REMOVE THE SIZE IN A PCA:
#join all tables together
dados_all = rbind(bre, sim, ste)
head(dados_all)
tail(dados_all)

#PCA analyses
dados = dados_all[,11:32]
sapply(dados, class)
head(dados)
summary(dados)
#number of variables
nvariable = length(dados)
nind = nrow(dados)
#normalize and stardard the morphological data
#log
dados = log(dados)
#standard
dados = decostand (dados, method="standardize")
summary(dados)
#center = TRUE means centring by the mean, in this case is 0
#scannf	= FALSE, indicates the screeplot should not be displayed
pca_input = dudi.pca(dados, center = TRUE , scannf = FALSE, nf = nvariable)

#% contribution of each pc
pcs = pca_input$eig/sum(pca_input$eig)
summary(pca_input) #must be the same values

#save results. PCA do not need to be performed with normal distribution variables or non-correlated variables
write.csv(pca_input$co, "contribution_variables_all.csv")
write.csv(pca_input$li, "contribution_individuals_all.csv")
write.csv(pcs, "contribution_pcs_all.csv")

##sum of PCs to discovery the number of PC that explain 90% of the variation
n.pc = length(pcs)
sum_pca = matrix(NA, n.pc, 2)
for (i in 1:n.pc) {sum_pca[i,] = c(i, sum(pcs[1:i]))}
colnames(sum_pca) = c("PCs","% Variance")
sum_pca
write.csv(sum_pca, file = "PC_variance_sum_all.csv", row.names = F)

#how many dimensions: log(sample size)> dimensions for KDE. Blonder 2016 (Holes in Hypervolumes)
log(36)
pc = 3 #number of PC with 91,6% of the variance to KDE
pc = 5 # PC = 96,4% to SVM


##7.CHOOSE FOR THE PC WITH AROUND 95% OF VARIANCE THE VARIABLES THAT HAVE VALUES BIGGER >0.32 - 10% OF THE VARIANCE (Tabachnick and Fidell 1989). THEY WILL FORM THE CLUSTERS
contr= as.data.frame(pca_input$co)
tab = abs(contr) #don't matter the direction
lista = list()

for (i in 1:pc)
{
  linhas=tab[tab[, i] > 0.32, ]
  lista[[i]]=row.names(linhas)
}

linhas.resultado=unlist(lista)
variables_morpho=unique(linhas.resultado)
write.table(variables_morpho, file = "Variables_10%_PCA.csv", sep = "\n", row.names = F, col.names = F)
variables_morpho #see the number of variables, good number < or equal 6; more than it you need to perform a new PCA to reduce the dimensionality or use the PCs instead of variables to create the hypervolume. 

##8. CREATE A NEW DATASET FOR THE IMPORTANT VARIABLES YOU WILL USE
dados_all_pc = dados_all[c(1:14,15:31,34)] #numbers are the positions of the columns in the initial table, 1:14 are information about the samples.

#in this case there are a lot of importants variables, so to reduce the dimensionality of data we will use the PCA to analysis:
dados_all_pc = as.data.frame(pca_input$li[,1:3]) #7 PCs to 95,5% of variance

#create PCs per species
values_bre = dados_all_pc[1:36,]
values_sim = dados_all_pc[37:72,]
values_ste = dados_all_pc[73:108,]


#hypervolume calculation:
hv_bre = hypervolume_gaussian(values_bre, name = "BRE", 
                              weight = NULL,
                              samples.per.point = ceiling((10^(3 + sqrt(ncol(values_bre))))/nrow(values_bre)),
                              kde.bandwidth = estimate_bandwidth(values_bre), 
                              sd.count = 3, 
                              quantile.requested = 0.95,
                              quantile.requested.type = "probability", 
                              chunk.size = 1000,
                              verbose = TRUE)

hv_sim = hypervolume_gaussian(values_sim, name = "SIM", 
                              weight = NULL,
                              samples.per.point = ceiling((10^(3 + sqrt(ncol(values_sim))))/nrow(values_sim)),
                              kde.bandwidth = estimate_bandwidth(values_sim), 
                              sd.count = 3, 
                              quantile.requested = 0.95,
                              quantile.requested.type = "probability", 
                              chunk.size = 1000,
                              verbose = TRUE)

hv_ste = hypervolume_gaussian(values_ste, name = "STE", 
                              weight = NULL,
                              samples.per.point = ceiling((10^(3 + sqrt(ncol(values_ste))))/nrow(values_ste)),
                              kde.bandwidth = estimate_bandwidth(values_ste), 
                              sd.count = 3, 
                              quantile.requested = 0.95,
                              quantile.requested.type = "probability", 
                              chunk.size = 1000,
                              verbose = TRUE)

##10. SAVE THE VALUES OF VOLUME AND CENTROID
hv_all = hypervolume_join(hv_bre, hv_sim, hv_ste)
centroid = get_centroid(hv_all)
vol = get_volume(hv_all)

write.csv(centroid, "centroid_WA_KDE.csv")
write.csv(vol, "volume_WA_KDE.csv")


##11. OVERLAP AMONG HYPERVOLUMES
#for two hypervolumes:

set_bre_sim = hypervolume_set(hv_bre, hv_sim, check.memory=FALSE)
bre_sim = hypervolume_overlap_statistics(set_bre_sim)

set_bre_ste = hypervolume_set(hv_bre, hv_ste, check.memory=FALSE)
bre_ste = hypervolume_overlap_statistics(set_bre_ste)

set_sim_ste = hypervolume_set(hv_sim, hv_ste, check.memory=FALSE)
sim_ste = hypervolume_overlap_statistics(set_sim_ste)


#creating a data frame with a results:
stat = rbind(bre_sim, bre_ste,
             sim_ste
             )

write.csv(stat, "comparison_morpho_hypervolumes_WA_KDE.csv")


##12. SAVE A FIGURE AS SHAPE
#for save the results
pdf("./Morphological_hypervolume_WA_KDE.pdf", onefile = F ) #change the name of pdf
plot.new()

plot(hv_all,
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=F, show.density=F,show.data=F,
     names=NULL, show.legend=F, limits=NULL,
     show.contour=T, contour.lwd=2,
     contour.type="kde",
     #contour.alphahull.alpha=0.25,
     #contour.ball.radius.factor=1,
     contour.kde.level=0.01,
     #contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2, pt.bg.centroid="black",
     colors= c('slategray4', 'blue', 'red'), #choose the colors
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1.5,cex.axis=0.9,cex.names=2,cex.legend=2,
     num.points.max.data = 1000, num.points.max.random = 2000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE)

#create a legend box for your graph
legend =c("P. brevicauda", "P. simonsi", "P. steerei")
col = c('slategray4', 'blue', 'red')
legend("bottomleft", legend = legend, cex = 1, bty="n",  col = "white", pt.bg = col, pch= 21, pt.cex=1.5)

dev.off()


#######################################################################
#######################################################################
#############################  SVM  ###################################
#######################################################################
#######################################################################

#in this case there are a lot of importants variables, so to reduce the dimensionality of data we will use the PCA to analysis:
dados_all_pc = as.data.frame(pca_input$li[,1:10]) #5 PCs to 95% of variance

#create PCs per species
values_bre = dados_all_pc[1:36,]
values_sim = dados_all_pc[37:72,]
values_ste = dados_all_pc[73:108,]

#hypervolume calculation:
hv_bre = hypervolume_svm(values_bre, name = "BRE", 
                            #  weight = NULL,
                              samples.per.point = ceiling((10^(3 + sqrt(ncol(values_bre))))/nrow(values_bre)),
                           #   kde.bandwidth = estimate_bandwidth(values_bre), 
                            #  sd.count = 3, 
                            #  quantile.requested = 0.95,
                             # quantile.requested.type = "probability", 
                         svm.nu = 0.01,
                         svm.gamma = 0.5, 
                         scale.factor = 1,
                              chunk.size = 1000,
                              verbose = TRUE)

hv_sim = hypervolume_gaussian(values_sim, name = "SIM", 
                             # weight = NULL,
                              samples.per.point = ceiling((10^(3 + sqrt(ncol(values_sim))))/nrow(values_sim)),
                            #  kde.bandwidth = estimate_bandwidth(values_sim), 
                             # sd.count = 3, 
                              #quantile.requested = 0.95,
                            #  quantile.requested.type = "probability", 
                            svm.nu = 0.01,
                            svm.gamma = 0.5, 
                            scale.factor = 1,
                              chunk.size = 1000,
                              verbose = TRUE)

hv_ste = hypervolume_gaussian(values_ste, name = "STE", 
                             # weight = NULL,
                              samples.per.point = ceiling((10^(3 + sqrt(ncol(values_ste))))/nrow(values_ste)),
                            #  kde.bandwidth = estimate_bandwidth(values_ste), 
                             # sd.count = 3, 
                            #  quantile.requested = 0.95,
                           #   quantile.requested.type = "probability", 
                              svm.nu = 0.01,
                              svm.gamma = 0.5, 
                              scale.factor = 1,
                              chunk.size = 1000,
                              verbose = TRUE)

##10. SAVE THE VALUES OF VOLUME AND CENTROID
hv_all = hypervolume_join(hv_bre, hv_sim, hv_ste)
centroid = get_centroid(hv_all)
vol = get_volume(hv_all)

write.csv(centroid, "centroid_WA_SVM.csv")
write.csv(vol, "volume_WA_SVM.csv")


##11. OVERLAP AMONG HYPERVOLUMES
#for two hypervolumes:
set_bre_sim = hypervolume_set(hv_bre, hv_sim, check.memory=FALSE)
bre_sim = hypervolume_overlap_statistics(set_bre_sim)

set_bre_ste = hypervolume_set(hv_bre, hv_ste, check.memory=FALSE)
bre_ste = hypervolume_overlap_statistics(set_bre_ste)

set_sim_ste = hypervolume_set(hv_sim, hv_ste, check.memory=FALSE)
sim_ste = hypervolume_overlap_statistics(set_sim_ste)


#creating a data frame with a results:
stat = rbind(bre_sim, bre_ste,
             sim_ste
)

write.csv(stat, "comparison_morpho_hypervolumes_WA_SVM.csv")


##12. SAVE A FIGURE AS SHAPE
#for save the results
pdf("./Morphological_hypervolume_WA_SVM.pdf", onefile = F ) #change the name of pdf
plot.new()

plot(hv_all,
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=F, show.density=F,show.data=F,
     names=NULL, show.legend=F, limits=NULL,
     show.contour=T, contour.lwd=2,
     contour.type="kde",
     #contour.alphahull.alpha=0.25,
     #contour.ball.radius.factor=1,
     contour.kde.level=0.01,
     #contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2, pt.bg.centroid="black",
     colors= c('slategray4', 'blue', 'red'), #choose the colors
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1.5,cex.axis=0.9,cex.names=2,cex.legend=2,
     num.points.max.data = 1000, num.points.max.random = 2000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE)

#create a legend box for your graph
legend =c("P. brevicauda", "P. simonsi", "P. steerei")
col = c('slategray4', 'blue', 'red')
legend("bottomleft", legend = legend, cex = 1, bty="n",  col = "white", pt.bg = col, pch= 21, pt.cex=1.5)

dev.off()