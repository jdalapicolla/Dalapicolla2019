###############################################################################
######################## SPECIES DISTRIBUTION MODELING ########################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
#################### STEP 05: FILTERING OCCURRENCE POINTS #####################

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
#A. OCCURRENCE POINTS BY SPECIES. SEE STEP 1.
#B. ONE ENVIRONMENTAL LAYER EDITED FOLLOWING THE STEP 4.


##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("raster", "rgdal", "pacman","tcltk","dismo", "rgeos", "MASS"))
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("raster", "rgdal", "tcltk", "dismo", "rgeos", "MASS")


##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())
#Any computer
setwd('C:/Users/Johnny/Downloads/bio_10m_bil') # change the path


##5. LOAD THE OCCURRENCE POINTS AND A LAYER YOU WILL USE IN THE MODELS:
#load and verify the table of points
points_raw = read.table ("ROB_C.csv", header=TRUE, sep=",")
head(points_raw)
tail(points_raw)
str(points_raw)

#load a layer
variavel = raster("bio1.asc")
crs(variavel) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
variavel
plot(variavel)


##6. REMOVING DUPLICATE RECORDS AND POINTS OUT OF THE STUDY AREA
#duplicate data
points_raw = points_raw [!duplicated(points_raw[c("longitude","latitude")]), ]

#remove points out of the mask
ocorrencia = points_raw [,2:3] #only lon/lat columns
head(ocorrencia)
tail(ocorrencia)
str(ocorrencia)
#add projection in table
coordinates(ocorrencia) = ~longitude+latitude #name of the columns.
crs(ocorrencia) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#identify points out of mask
ovr = extract(variavel, ocorrencia)
head(ovr)
i = which(ovr != "-9999") #replace the value if the NoData value is different of -9999
i #lines in the point_raw

#update the points raw and create a SpatialPoints for ocorrence points.
points_raw = points_raw[i,]
ocorrencia = points_raw[,2:3]
coordinates(ocorrencia) = ~longitude+latitude #name of the columns.
crs(ocorrencia) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


##7. ELIMINATING THE SAMPLING BIAS:
#create a buffer of 10Km around the points
buffer = circles(ocorrencia, d = 10000, lonlat=TRUE) #d is radius of each circle in meters
plot(buffer)
class(buffer)

#convert circles in polygon
buffer = polygons(buffer)
#rasterize the buffer as the layer
buffer= rasterize(buffer, variavel)

#select 1 point in each circle
sel = gridSample(ocorrencia, buffer, n=1)


##8. SAVE THE RESULTS IN A CSV FILE
sel = as.data.frame(sel)
write.csv(sel, "ROB_C_UNBIAS.csv", row.names = FALSE)


##9.CREATE A BIAS GRID/FILE:
#create the Two-Dimensional Kernel Density Estimation with same proprieties than the example raster
dens = kde2d(sel[,1], sel[,2], n = c(ncol(variavel), nrow(variavel)), lims = c(xmin(variavel), xmax(variavel), ymin(variavel), ymax(variavel)))
#convert the density in raster             
dens.ras = raster(dens)

#resample:
dens.ras = resample(dens.ras, variavel, method="bilinear", snap="out")

#compare the nrow, ncol, ncell, the resolution and the extent of the dens.ras to the variavel, Everything must to be the same! 
dens.ras
variavel

#plot the density map:
plot(dens.ras)

#save the bias files as .ascii
writeRaster(dens.ras, filename= "biasfile_ROB_C.asc", driver='ascii')

#################################################################################################
#################################################################################################
## RESAMPLE FUNCTION DO NOT CORRECT SOME DECIMALS NUMBERS OF THE BIAS FILE, AFTER THE 5TH NUMBER AFTER COMMA. SO YOU HAVE TO OPEN THE ASC FILES (BIAS FILE AND THE VARIABLES) IN A TEXT EDITOR AND CORRECT MANUALLY THE NUMBER OF THE BIAS FILE. IN MY CASE WERE IN THE 3 LAST NUMBERS. AFTER THAT THE BIAS FILE WORKED IN THE MAXENT SOFTWARE.