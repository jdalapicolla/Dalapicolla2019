###############################################################################
############################# HYPERVOLUME ANALYSES ############################
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
#A. BINARY MAP REPRESENTING A SPECIES DISTRIBUTION MODEL FOR PRESENT-DAY 
#B. ENVIRONMENTAL LAYERS FOR CLIMATE HYPERVOLUME

##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("pacman","RStoolbox", "raster", "rgdal", "dismo", "hypervolume","alphahull", "tcltk", "vegan", "adegenet"))
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("RStoolbox", "raster", "rgdal", "dismo", "hypervolume", "alphahull", "tcltk", "vegan", "adegenet")


################################################################################
####################### ENVIRONMENTAL HYPERVOLUME ##############################
################################################################################


##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())
#Any computer
setwd('C:/Users/Johnny/Downloads/bio_10m_bil') # change the path

##5. LOAD THE ENVIRONMENTAL LAYERS
layers = stack(sapply(list.files(pattern='asc$', recursive = F), raster)) #load all 19 variables or more. Variables must be clipped and resampled for the study area to use in this analysis!!!
crs(layers) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers

##6. PCA FOR THE GEOGRAPHICAL AREA
#spca = standard; nSamples = pixels number used, maskCheck = avoid NA; nComp = number of PC to use 
PCA_climate_WA = rasterPCA(layers, nSamples = NULL, nComp = 19, spca = TRUE, maskCheck = TRUE)
#plot the map with the 3 first PC
plotRGB(PCA_climate_WA$map,r = 1, g = 2, b = 3, alpha =120, stretch="lin", axes = F)
#create a new object for the PCA map
map_PCA_WA = PCA_climate_WA$map
#% of variation for each PC
summary(PCA_climate_WA$model) #choose the number of PCs that explain 90% of variance or more

#save the results
#save the PCA map in raster
writeRaster(PCA_climate_WA$map, "PCA_WA_19PC_RASTER.tif", driver='GTiff', bylayer=F)

#save the % of variance of each PC
std_dev = PCA_climate_WA$model$sdev
pr_var = std_dev^2
prop_varex = pr_var/sum(pr_var)
prop_varex = prop_varex*100 # must be the same values than summary(PCA_climate_WA$model)
write.csv(prop_varex, file = "contribution_pc_eig_WA.csv")

#save the variables contribution to the PC's
write.csv(PCA_climate_WA$model$loadings, "contribuition_variables_WA.csv")

##7. LOAD THE BINARY MAPS FOR THE TARGET SPECIES:
binary_bre = raster("Binary_BRE_PRES.asc")
binary_sim = raster("Binary_SIM_PRES.asc")
binary_ste = raster("Binary_STE_PRES.asc")
#plot to verify
plot(binary_bre)
plot(binary_sim)
plot(binary_ste)

##8. CREATE 1,000 RANDOM POINTS INSIDE THE AREA OF BINARY MAPS, PLOT AND SAVE THEM
randomPoints_bre = as.data.frame(randomPoints(binary_bre, 1000))
plot(randomPoints_bre, cex = 0.5)
colnames(randomPoints_bre) = c("long", "lat")
write.csv(randomPoints_bre, "BRE_randomPoints.csv", row.names = F)

randomPoints_sim = as.data.frame(randomPoints(binary_sim, 1000))
plot(randomPoints_sim, cex = 0.5)
colnames(randomPoints_sim) = c("long", "lat")
write.csv(randomPoints_sim, "SIM_randomPoints.csv", row.names = F)

randomPoints_ste = as.data.frame(randomPoints(binary_ste, 1000))
plot(randomPoints_ste, cex = 0.5)
colnames(randomPoints_ste) = c("long", "lat")
write.csv(randomPoints_ste, "STE_randomPoints.csv", row.names = F)

#load random points saved if you need to:
randomPoints_bre=read.csv("BRE_randomPoints.csv")
randomPoints_sim=read.csv("SIM_randomPoints.csv")
randomPoints_ste=read.csv("STE_randomPoints.csv")


##9. CREATE THE HYPERVOLUME WITH THE NUMBER OF PC EQUAL OR MORE THAN 90% OF VARIATION
#extracting values of variables for the Random Points per species 
values_bre = extract(map_PCA_WA, randomPoints_bre)
values_sim = extract(map_PCA_WA, randomPoints_sim)
values_ste = extract(map_PCA_WA, randomPoints_ste)

#Use only the first 5 PC for analysis: ~95,15% of variation
values_bre=as.data.frame(values_bre[,1:5]) #change the number 5 for the number of PC you need to use!
values_sim=as.data.frame(values_sim[,1:5])
values_ste=as.data.frame(values_ste[,1:5])

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
hv_wa = hypervolume_join(hv_bre, hv_sim, hv_ste)
centroid = get_centroid(hv_wa)
vol = get_volume(hv_wa)

write.csv(centroid, "centroid_WA.csv")
write.csv(vol, "volume_WA.csv")

##11. COMPARING HYPERVOLUMES AMONG THE TARGET SPECIES
#3 species, pair-to-pair, are 3 combinations:

set_bre_sim = hypervolume_set(hv_bre, hv_sim, check.memory=FALSE)
bre_sim = hypervolume_overlap_statistics(set_bre_sim)

set_bre_ste = hypervolume_set(hv_bre, hv_ste, check.memory=FALSE)
bre_ste = hypervolume_overlap_statistics(set_bre_ste)

set_sim_ste = hypervolume_set(hv_sim, hv_ste, check.memory=FALSE)
sim_ste = hypervolume_overlap_statistics(set_sim_ste)

#creating a data frame with a results:
stat = cbind(bre_sim, bre_ste, sim_ste)
write.csv(stat, "comparison_enriv_hypervolumes_WA.csv")


##12. SAVE A FIGURE AS SHAPE
#for save the results
pdf("./Environmental_hypervolumes_WA.pdf", onefile = F ) #change the name of pdf
plot.new()

plot(hv_wa,
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