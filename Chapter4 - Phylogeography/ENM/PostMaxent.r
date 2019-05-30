###############################################################################
######################## SPECIES DISTRIBUTION MODELING ########################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
############################ STEP 09: POST MAXENT #############################

###############################################################################
########### Based on Mariah Kenney and Christopher Tracey scripts #############
###################### KNOWLES LAB, UNIVERSITY OF MICHIGAN ####################
###############################################################################

###############################################################################
######################### JERONYMO DALAPICOLLA, 2019 ##########################
###################### LABORATORIO DE MAMIFEROS - LAMA ########################
####################### ESALQ - UNIVERSIDADE DE SAO PAULO #####################
###############################################################################


##1.FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/


##2.INPUTS FOR THIS STEP:
#A. MODELS FROM STEP 08.


##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("raster", "rgdal", "rgeos", "pacman","tcltk", "ROCR", "vcd", "boot"))
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("raster", "rgdal", "rgeos", "pacman", "tcltk", "ROCR", "vcd", "boot")


##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())
#Any computer
setwd('C:/Users/Johnny/Downloads/bio_10m_bil') # change the path


## 5. LOAD THE OCCURRENCE POINTS USED TO BUILD THE MODELS
#unbiased points
occ = read.csv("ROB_A_UNBIAS.csv")
# if you have more that long/lat columns - we need to remove all other columns
occ = occ[,2:3]
head(occ)


##6. LOAD THE MODELS FOR THE SAME SPECIES
species = "ROB_A" #species name - the same used in MaxEnt
folder = c("PRES", "HOL", "LGM", "LIG") #folder names with the projections
ending = "_avg.asc" #final name of the average models

models = c()
layers = list()
for (i in 1:4) #number of models, same as folder object
{models[i] = paste0(species,"_",folder[i],ending)
layers[[i]] = stack(sapply(list.files(pattern=models[i], recursive = TRUE), raster))}
model = stack(layers[[1]], layers[[2]], layers[[3]], layers[[4]]) #number of models


##7. CHOOSE THE THRESHOLD
#extract values at all occurrence points as use top 90% as presence to breate binary. 10% threshold.
suit = extract(model[[1]], occ, latlong=FALSE) #ATTENTION! Values from the present-day model, in the case the first one in the model object
thresh = quantile(suit, 0.1)
thresh #value of the threshold


##8. CREATE A BINARY MAP
#reclassify the model
pxReclass = list()
for (i in 1:4) {pxReclass[[i]] = reclassify(model[[i]], rcl = c(0, thresh, NA, thresh, 1, 1))}
plot(pxReclass[[1]])
pxReclass = stack(pxReclass[[1]], pxReclass[[2]], pxReclass[[3]], pxReclass[[4]])

#save binary maps
dir.create("Binary_Maps")
setwd ("./Binary_Maps")
writeRaster(pxReclass, paste0(models,".asc"), driver='ascii', overwrite=TRUE, bylayer=T)
setwd("..")


##9. CUT THE MODEL BY THE THRESHOLD FOR PUBLICATION
cortados = list()
for (i in 1:4) {cortados[[i]] = mask(model[[i]], pxReclass[[i]])}
plot(cortados[[1]])
cortados = stack(cortados[[1]], cortados[[2]], cortados[[3]], cortados[[4]])

#save the raster
dir.create("Final_Models")
setwd ("./Final_Models")
writeRaster(cortados, paste0(models,".asc"), driver='ascii', overwrite=TRUE, bylayer=T)
setwd("..")


##################################################################################
######################### STABLE AREAS THROUGH TIME ##############################
##################################################################################

##10. LOAD THE BIANRY MAPS
setwd ("./Binary_Maps")
bin_maps = stack(sapply(list.files(pattern='asc$', recursive = F), raster))
crs(bin_maps) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
setwd ("..")


##11. CONVERT THE RASTER INTO POLYGON
pol = list()
for (i in 1:length(names(bin_maps))) {pol[[i]] = rasterToPolygons(bin_maps[[i]], fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=T)
crs(pol[[i]]) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"}
plot(pol[[1]])

##12. CALCULATE THE INTERSECTION AMONG ALL MAPS
stable = pol[[1]]
for (i in 2:length(names(bin_maps))) {stable = intersect(stable, pol[[i]])}
plot(stable)
class(stable)

#save the polygon
dir.create("Stable_Areas")
setwd ("./Stable_Areas")
writeOGR(stable, dsn = ".", layer=paste0("Stable",species), driver='ESRI Shapefile', overwrite=TRUE)
setwd ("..")

##################################################################################
################### CALCULATE THE MODEL AREA IN KM2 ##############################
##################################################################################

##13. CALCULATE THE MODEL AREAS
pol[[5]] = stable
names_pol = c(names(bin_maps), "stable")
area_km2 = matrix(NA,1,5)
row.names(area_km2) = species
colnames(area_km2) = names_pol

for (i in 1:length(pol)) {area_km2[[i]] = area(pol[[i]]) / 1000000}
area_km2

write.csv(area_km2, paste0("areas_model_", species, ".csv"))

##################################################################################
#################### CALCULATE THE TSS FROM MAXENT OUTPUT ########################
##################################################################################
#---------------------------------------------------------------------------------------------
# Name: TSS.R
# Purpose: calculate True Skill Statistic (TSS) for maxent models
# Author: Christopher Tracey
# Created: 2015-03-05
# Updated: 2015-03-05
#
# Updates:
# * 2015-03-06 - added possible values for the threshold rule
#
# To Do List/Future Ideas:
# * put the threshold rule into a variable
#---------------------------------------------------------------------------------------------

# based on the script by  referenced here:
# https://groups.google.com/forum/#!msg/maxent/CUeI5xT9wTI/i3aibLdDOEYJ
# as well as the discussion located here:
# https://groups.google.com/forum/#!topic/maxent/eCgJ_0vdOb0

#ATTENTION! "Write background predictions" has to be enabled in Maxent

# set the working directory to where the maxent files are at
setwd("C:/Users/ctracey/Dropbox/_PNHP_ConservationPlanning/WRCP SWAP Predictive Modeling/ModelEvaluationR/Swainsons")

list.files(pattern="_samplePredictions.csv")->listaoutput
sub("_samplePredictions.csv","",listaoutput)->listaoutput

tss_general <- as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(tss_general)<-c("Replica","TSS","Bin.Prob")

lapply(listaoutput, function(nome){
  
  nome->sp
  presence <- read.csv(paste(sp,"_samplePredictions.csv",sep=""))
  background <- read.csv(paste(sp,"_backgroundPredictions.csv",sep=""))
  pp <- presence$Logistic.prediction                # get the column of predictions
  testpp <- pp[presence$Test.or.train=="test"]      # select only test points
  trainpp <- pp[presence$Test.or.train=="train"]    # select only test points
  bb <- background$Logistic
  combined <- c(testpp, bb)                         # combine into a single vector
  
  read.csv("maxentResults.csv")->maxres
  
  # Set the threshold rule here
  threshold<-maxres[maxres[,1]==nome,"X10.percentile.training.presence.Logistic.threshold"]
  bin.prob<-maxres[maxres[,1]==nome,"X10.percentile.training.presence.binomial.probability"]
  # other possible choices are (cut and paste as needed:
  # "Fixed.cumulative.value.1.logistic.threshold"                                       
  # "Fixed.cumulative.value.1.binomial.probability"                                     
  # "Fixed.cumulative.value.5.logistic.threshold"                                       
  # "Fixed.cumulative.value.5.binomial.probability"                                                                    
  # "Fixed.cumulative.value.10.logistic.threshold"                                                                                                                                                                        
  # "Fixed.cumulative.value.10.binomial.probability"                                                                       
  # "Minimum.training.presence.logistic.threshold"                                                                                                                                                                          
  # "Minimum.training.presence.binomial.probability"                                                                
  # "X10.percentile.training.presence.logistic.threshold"                                                                                                                                             
  # "X10.percentile.training.presence.binomial.probability"                                               
  # "Equal.training.sensitivity.and.specificity.logistic.threshold"                                                                                                   
  # "Equal.training.sensitivity.and.specificity.binomial.probability"                                 
  # "Maximum.training.sensitivity.plus.specificity.logistic.threshold"                                                                                       
  # "Maximum.training.sensitivity.plus.specificity.binomial.probability"                                     
  # "Equal.test.sensitivity.and.specificity.logistic.threshold"                                                                                                                     
  # "Equal.test.sensitivity.and.specificity.binomial.probability"                                       
  # "Maximum.test.sensitivity.plus.specificity.logistic.threshold"                                                
  # "Maximum.test.sensitivity.plus.specificity.binomial.probability"                    
  # "Balance.training.omission..predicted.area.and.threshold.value.logistic.threshold"       
  # "Balance.training.omission..predicted.area.and.threshold.value.binomial.probability"   
  # "Equate.entropy.of.thresholded.and.original.distributions.logistic.threshold"               
  # "Equate.entropy.of.thresholded.and.original.distributions.binomial.probability" 
  
  # Number of values greater and less than the threshold test
  sum(testpp > threshold) -> majortest
  sum(testpp < threshold) -> minortest
  
  # Number of high and low values that the threshold in the background
  sum(bb > threshold) -> majorbb
  sum(bb < threshold) -> minorbb
  
  ### Calculate sensitivity and specificity
  sensitivity <- (majortest) / (majortest+minortest)
  specificity <- (minorbb) / (majorbb+minorbb)
  
  tss <- sensitivity + specificity - 1
  
  tsssp <- as.data.frame(sp)
  tsssp[2] <- tss
  tsssp[3] <- bin.prob
  colnames(tsssp) <- c("Replica","TSS","Bin.Prob")
  
  rbind(tss_general,tsssp)->>tss_general 
})

# write the final results
dir.create("../TSS")
setwd("../TSS")
write.csv(tss_general, paste0("TSS_general_", species, ".csv"))

