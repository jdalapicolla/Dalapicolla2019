###############################################################################
######################## SPECIES DISTRIBUTION MODELING ########################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
################## STEP 04: EDITING ENVIRONMENTAL LAYERS ######################

###############################################################################
######################## Based on Mariah Kenney scripts #######################
###################### KNOWLES LAB, UNIVERSITY OF MICHIGAN ####################
###############################################################################

###############################################################################
######################### JERONYMO DALAPICOLLA, 2019 ##########################
###################### LABORATORIO DE MAMIFEROS - LAMA ########################
####################### ESALQ - UNIVERSIDADE DE SAO PAULO #####################
###############################################################################


##1.FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/


##2.INPUTS FOR THIS STEP:
#A. ENVIRONMENTAL LAYERS
#B. MASK REPRESENTING THE STUDY AREA.
#C. CHOOSE THE RESOLUTION OF THE ENVIROMENTAL LAYERS FOR THE WORK.


##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("raster", "rgdal", "pacman","tcltk"))
#load the packages
library(raster)
library(rgdal)
library(tcltk)
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("raster", "rgdal", "tcltk")


##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())
#Any computer
setwd('C:/Users/Johnny/Downloads/bio_10m_bil') # change the path to the folder 


##5. SELECT THE LAYERS WITH THE EXTENSION YOU WANT TO. MY CASE IS '.bil' or '.tif':
#load an example! A layer with the resulution you want to. In my case one with 2,5 arc-seconds:
example = raster("bio1.bil")
crs(example) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#do a list of all the .bil files and stack them together. Recursive = TRUE include all subfolders. I put layers from LGM, Holocene and present in subfolders in the same folder. Layers must be the same resolution to became a stack.
layers = stack(sapply(list.files(pattern='bil$', recursive = F), raster))
crs(layers) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#Verify if all your raster are in 'layers' object.
layers

#do a list of all the .tif files and stack them together
layers = stack(sapply(list.files(pattern='tif$', recursive = TRUE), raster))
crs(layers) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers

#do a list of all the .asc files and stack them together
layers = stack(sapply(list.files(pattern='asc$', recursive = TRUE), raster))
crs(layers) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers


##6. LOAD THE MASK FOR THE STUDY AREA 
#load the mask, if you need to
study_area = shapefile("WA+CER_projection_all.shp")
#define a projection, the same one for all!!!!
crs(study_area) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(study_area)


##7. RESAMPLING AND CLIPPING THE LAYERS WITH THE EXAMPLE
#reduce the area for the resample. Make the analysis faster.
example2 = crop(example, extent(study_area), snap="out") #retangular area
plot(example2)

#Resample all layers. This step will take a long time.
layers_res = resample(layers, example2, method="bilinear", bylayer=TRUE, progress='text', snap="out")
layers_res
plot(layers_res[[1]])

#Clip the resampled layers with the study area:
layers_res = mask(layers_res, study_area, bylayer=TRUE) #exactly the same area than mask
plot(layers_res[[1]])


##8. EXPORT IN .ASC
#create a new folder and move to there:
dir.create("Layers_Resampled_2.5")
setwd(file.path(getwd(),"Layers_Resampled_2.5"))

#Save as .ascii in the new folder
writeRaster(layers_res, paste0(names(layers_res),".asc"), driver='ascii', bylayer=TRUE)
#come back to former folder
setwd('..')