###############################################################################
####################### PHYLOGEOGRAPHICAL ANALYSES ############################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
############       STEP 03: COMPARATIVE PHYLOGEOGRAPHY       ##################
##################  IBD IBE AMONG DIFERENT SPECIES  ###########################

###############################################################################
################ Based on Joyce R. Prado and J. P. Huang scripts ##############
###################### KNOWLES LAB, UNIVERSITY OF MICHIGAN ####################
###############################################################################

###############################################################################
######################### JERONYMO DALAPICOLLA, 2019 ##########################
###################### LABORATORIO DE MAMIFEROS - LAMA ########################
####################### ESALQ - UNIVERSIDADE DE SAO PAULO #####################
###############################################################################


##1.FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/


##2.INPUTS FOR THIS STEP:
#A. THE FILE ".fst_summary.tsv" CLEANED AFTER POPULATIONS STEP IN STACKS, ORGANIZED BY POPULATIONS/GROUPS.
#B. GEOGRAPHICAL COORDENATES FOR INDIVIDUALS IN DECIMAL DEGREES.


##3. INSTALL AND LOAD THE PACKAGES#
#install the packages
#install.packages(c("ade4","adegenet", "plyr", "fossil", "vegan", "qiimer" , "geosphere", "gdistance", 'RStoolbox', "tcltk"))
#load the packages
library("ade4")
library("adegenet")
library("plyr")
library("fossil")
library("vegan")
library("qiimer")
library ("geosphere")
library("gdistance")
library('RStoolbox')
library("tcltk") #for Linux

#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("ade4","adegenet", "plyr", "fossil", "vegan", "qiimer" ,"geosphere" , "gdistance", 'RStoolbox', "tcltk")


##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())


####################################################
####### 5. CREATE COORDENATES BY POPULATIONS #######
####################################################

###################
####### SP1 #######
###################

# IMPORT THE LON-LAT MATRIX FOR THE SAMPLES
coord_mat = read.table('coord_bre.txt', header = F, sep = "\t")
head(coord_mat)
tail(coord_mat)
str(coord_mat)
coord_mat2= coord_mat
coord_mat = coord_mat[ ,c(14:15)]

#SPLIT COORDENATES BY POPULATIONS:
coord_GALVEZ = coord_mat[1:3,]
coord_MADRE = coord_mat[4:5,]
coord_BOLIVIA = coord_mat[6:8,]
coord_AMAZONAS = coord_mat[9:10,]
coord_ACRE = coord_mat[11:13,]
coord_IQUITOS = coord_mat[14:16,]

#DELIMIT A CENTROID OR MEAN FOR VARIOUS COORDENATES BY POPULATION:
coord_GALVEZ_centro = coord_GALVEZ[1,]
rownames(coord_GALVEZ_centro) = "GALVEZ"
colnames(coord_GALVEZ_centro) = c("lon", "lat")

coord_MADRE_centro = data.frame(mean(coord_MADRE[,1]), mean(coord_MADRE[,2])) # two
rownames(coord_MADRE_centro) = "MADRE"
colnames(coord_MADRE_centro) = c("lon", "lat")

coord_BOLIVIA_centro = data.frame(mean(coord_BOLIVIA[,1]), mean(coord_BOLIVIA[,2])) # two samples used mean
rownames(coord_BOLIVIA_centro) = "BOLIVIA"
colnames(coord_BOLIVIA_centro) = c("lon", "lat")

coord_AMAZONAS_centro = data.frame(mean(coord_AMAZONAS[,1]), mean(coord_AMAZONAS[,2]))
rownames(coord_AMAZONAS_centro) = "AMAZONAS"
colnames(coord_AMAZONAS_centro) = c("lon", "lat")
coord_AMAZONAS_centro[1,1] = -77.76667 #mean fall out of the area. You need to verify points

coord_ACRE_centro = data.frame(mean(coord_ACRE[,1]), mean(coord_ACRE[,2])) # two
rownames(coord_ACRE_centro) = "ACRE"
colnames(coord_ACRE_centro) = c("lon", "lat")

coord_IQUITOS_centro = coord_IQUITOS[1,]
rownames(coord_IQUITOS_centro) = "IQUITOS"
colnames(coord_IQUITOS_centro) = c("lon", "lat")


#CREATE A MATRIX WITH COORDENATES PER POPULATION:
pop_coord_SP1 = rbind(coord_GALVEZ_centro, coord_MADRE_centro, coord_BOLIVIA_centro,
                  coord_AMAZONAS_centro, coord_ACRE_centro, coord_IQUITOS_centro)

write.csv(pop_coord_SP1, file = "coord_bypop_bre.csv", row.names = T)


###################
####### SP2 #######
###################

# IMPORT THE LON-LAT MATRIX FOR THE SAMPLES
coord_mat = read.table('coord_sim.txt', header = F, sep = "\t")
head(coord_mat)
tail(coord_mat)
str(coord_mat)
coord_mat2= coord_mat
coord_mat = coord_mat[ ,c(14:15)]

#SPLIT COORDENATES BY POPULATIONS:
coord_SOLIMOES = coord_mat[1:2,]
coord_CENTRAL = coord_mat[3:5,]
coord_MADEIRA = coord_mat[6:7,]
coord_GALVEZ2 = coord_mat[8:9,]
coord_YUNGAS = coord_mat[10:11,]
coord_MADRE2 = coord_mat[12:13,]
coord_JURUA = coord_mat[14:15,]

#DELIMIT A CENTROID OR MEAN FOR VARIOUS COORDENATES BY POPULATION:
coord_SOLIMOES_centro = data.frame(mean(coord_SOLIMOES[,1]), mean(coord_SOLIMOES[,2]))
rownames(coord_SOLIMOES_centro) = "SOLIMOES"
colnames(coord_SOLIMOES_centro) = c("lon", "lat")

coord_CENTRAL_centro = centroid(coord_CENTRAL)
rownames(coord_CENTRAL_centro) = "CENTRAL"
colnames(coord_CENTRAL_centro) = c("lon", "lat")
coord_CENTRAL_centro[1,] = c(-68.76670, -13.58330) #mean fall out of the area. You need to verify points

coord_MADEIRA_centro = data.frame(mean(coord_MADEIRA[,1]), mean(coord_MADEIRA[,2]))
rownames(coord_MADEIRA_centro) = "MADEIRA"
colnames(coord_MADEIRA_centro) = c("lon", "lat")

coord_GALVEZ2_centro =  coord_GALVEZ2[1,]
rownames(coord_GALVEZ2_centro) = "GALVEZ"
colnames(coord_GALVEZ2_centro) = c("lon", "lat")

coord_YUNGAS_centro = data.frame(mean(coord_YUNGAS[,1]), mean(coord_YUNGAS[,2]))
rownames(coord_YUNGAS_centro) = "YUNGAS"
colnames(coord_YUNGAS_centro) = c("lon", "lat")

coord_MADRE2_centro = data.frame(mean(coord_MADRE2[,1]), mean(coord_MADRE2[,2]))
rownames(coord_MADRE2_centro) = "MADRE2"
colnames(coord_MADRE2_centro) = c("lon", "lat")

coord_JURUA_centro = data.frame(mean(coord_JURUA[,1]), mean(coord_JURUA[,2]))
rownames(coord_JURUA_centro) = "JURUA"
colnames(coord_JURUA_centro) = c("lon", "lat")


#CREATE A DISTANCE MATRIX FOR GEOGRAPHIC LOCALITIES:
pop_coord_SP2 = rbind(coord_SOLIMOES_centro, coord_CENTRAL_centro, coord_MADEIRA_centro, coord_GALVEZ2_centro, coord_YUNGAS_centro, coord_MADRE2_centro, coord_JURUA_centro)

write.csv(pop_coord_SP2, file = "coord_bypop_sim.csv", row.names = T)


###################
####### SP3 #######
###################

# IMPORT THE LON-LAT MATRIX FOR THE SAMPLES
coord_mat = read.table('coord_ste.txt', header = F, sep = "\t")
head(coord_mat)
tail(coord_mat)
str(coord_mat)
coord_mat2= coord_mat
coord_mat = coord_mat[ ,c(13:14)]

#SPLIT COORDENATES BY POPULATIONS:
coord_BENI = coord_mat[1:2,]
coord_PANDO = coord_mat[3:4,]
coord_SOLIMOES2 = coord_mat[5:6,]
coord_ACRE2 = coord_mat[7:8,]
coord_JURUA2 = coord_mat[9:12,]
coord_JAINU = coord_mat[13:14,]

#DELIMIT A CENTROID OR MEAN FOR VARIOUS COORDENATES BY POPULATION:
coord_BENI_centro = data.frame(mean(coord_BENI[,1]), mean(coord_BENI[,2])) # with two different samples you must use the mean
rownames(coord_BENI_centro) = "MADRE-BENI"
colnames(coord_BENI_centro) = c("lon", "lat")

coord_PANDO_centro = data.frame(mean(coord_PANDO[,1]), mean(coord_PANDO[,2])) #more than two, use centroid.
rownames(coord_PANDO_centro) = "PANDO"
colnames(coord_PANDO_centro) = c("lon", "lat")

coord_SOLIMOES2_centro = data.frame(mean(coord_SOLIMOES2[,1]), mean(coord_SOLIMOES2[,2]))
rownames(coord_SOLIMOES2_centro) = "SOLIMOES"
colnames(coord_SOLIMOES2_centro) = c("lon", "lat")

coord_ACRE2_centro =  coord_ACRE2[1,]
rownames(coord_ACRE2_centro) = "ACRE2"
colnames(coord_ACRE2_centro) = c("lon", "lat")

coord_JURUA2_centro = centroid(coord_JURUA2)
rownames(coord_JURUA2_centro) = "JURUA"
colnames(coord_JURUA2_centro) = c("lon", "lat")

coord_JAINU_centro =  coord_JAINU[1,]
rownames(coord_JAINU_centro) = "JAINU"
colnames(coord_JAINU_centro) = c("lon", "lat")

#CREATE A DISTANCE MATRIX FOR GEOGRAPHIC LOCALITIES:
pop_coord_SP3 = rbind(coord_BENI_centro, coord_PANDO_centro, coord_SOLIMOES2_centro, coord_ACRE2_centro, coord_JURUA2_centro, coord_JAINU_centro)

write.csv(pop_coord_SP3, file = "coord_bypop_sim.ste", row.names = T)

##Add the number of species you need to.


#####################################
## 6. GEOGRAPHICAL DISTANCE MATRIX ##
#####################################

###################
####### SP1 #######
###################

geoDIST_SP1 = earth.dist(pop_coord_SP1, dist = TRUE)
geoDIST_SP1 = as.dist(geoDIST_SP1)
#geoDIST_SP1 = pcnm(geoDIST_SP1)

###################
####### SP2 #######
###################

geoDIST_SP2 = earth.dist(pop_coord_SP2, dist = TRUE)
geoDIST_SP2 = as.dist(geoDIST_SP2)
#geoDIST_SP2 = pcnm(geoDIST_SP2)

###################
####### SP3 #######
###################

geoDIST_SP3 = earth.dist(pop_coord_SP3, dist = TRUE)
geoDIST_SP3 = as.dist(geoDIST_SP3)
#geoDIST_SP3 = pcnm(geoDIST_SP3)


################################
## 7. GENETIC DISTANCE MATRIX ##
################################

###################
####### SP1 #######
###################

#based on pairwise FST by population:
fst_pop = read.table("batch_1.fst_summary_bre.tsv", header = T, row.names = 1, sep="\t",  na.strings = "")
fst_pop = t(fst_pop)
fst_pop_sp1 = fst_pop/(1-fst_pop)
fst_pop = as.data.frame(fst_pop_sp1)
fst_pop$EM = c("NA", "NA", "NA", "NA", "NA", "NA")
fst_pop = as.matrix(fst_pop)
fst_DIST_SP1 = as.dist(fst_pop) #warnings will be create about the empty cells.


###################
####### SP2 #######
###################

#based on pairwise FST by population:
fst_pop = read.table("batch_1.fst_summary_sim.tsv", header = T, row.names = 1, sep="\t",  na.strings = "")
fst_pop = t(fst_pop)
fst_pop_sp2 = fst_pop/(1-fst_pop)
fst_pop = as.data.frame(fst_pop_sp2)
fst_pop$EM = c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
fst_pop = as.matrix(fst_pop)
fst_DIST_SP2 = as.dist(fst_pop)


###################
####### SP3 #######
###################

#based on pairwise FST by population:
fst_pop = read.table("batch_1.fst_summary_ste.tsv", header = T, row.names = 1, sep="\t",  na.strings = "")
fst_pop = t(fst_pop)
fst_pop_sp3 = fst_pop/(1-fst_pop)
fst_pop = as.data.frame(fst_pop_sp3)
fst_pop$EM = c("NA", "NA", "NA", "NA", "NA", "NA")
fst_pop = as.matrix(fst_pop)
fst_DIST_SP3 = as.dist(fst_pop)


######################################
## 8. ENVIRONMENTAL DISTANCE MATRIX ##
######################################

# LOAD THE ENVIRONMENTAL LAYERS
layers = stack(sapply(list.files(pattern='asc$', recursive = F), raster)) #load all 19 variables or more. Variables must be clipped and resampled for the study area to use in this analysis!!!
crs(layers) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers

# PCA FOR THE GEOGRAPHICAL AREA
#spca = standard; nSamples = pixels number used, maskCheck = avoid NA; nComp = number of PC to use 
PCA_climate_WA = rasterPCA(layers, nSamples = NULL, nComp = 5, spca = TRUE, maskCheck = TRUE)
#plot the map with the 3 first PC
plotRGB(PCA_climate_WA$map,r = 1, g = 2, b = 3, alpha =120, stretch="lin", axes = F)
#create a new object for the PCA map
map_PCA_WA = PCA_climate_WA$map
#% of variation for each PC
summary(PCA_climate_WA$model) #analyze % of variance explained by PC1.

# EXTRACT THE PC LOADINGS PER POPULATION:
#extracting values of variables for the coordenates per populations in each species. Save the PC1 results:
PCs_SP1 = extract(map_PCA_WA, pop_coord_SP1)
PC1_SP1 = PCs_SP1[,1]

PCs_SP2 = extract(map_PCA_WA, pop_coord_SP2)
PC1_SP2 = PCs_SP2[,1]

PCs_SP3 = extract(map_PCA_WA, pop_coord_SP3)
PC1_SP3 = PCs_SP3[,1]


#############################
## 9. COST DISTANCE MATRIX ##
#############################

###################
####### SP1 #######
###################

# ENM FOR PRESENT-DAY
enm_cur_SP1 = raster("BRE_avg.asc")
crs(enm_cur_SP1) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_cur_SP1)
#invert values
enm_cur_SP1 = 1-enm_cur_SP1
plot(enm_cur_SP1)

#transition matrix
enm_cur_SP1_tr = transition(enm_cur_SP1, transitionFunction=mean, 16)

#correct the transition matrix
#The argument scl is set to TRUE to scale the transition values to a reasonable range. If the transition values are too large, commute distance and randomized shortest path functions will not work well. No scaling should be done if the user wants to obtain absolute distance values as output.
enm_cur_SP1_trC = geoCorrection(enm_cur_SP1_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR HOLOCENE
enm_hol_SP1 = raster("BRE_HOL_avg.asc")
crs(enm_hol_SP1) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_hol_SP1)
enm_hol_SP1 = 1-enm_hol_SP1
plot(enm_hol_SP1)
enm_hol_SP1_tr = transition(enm_hol_SP1, transitionFunction=mean, 16)
enm_hol_SP1_trC = geoCorrection(enm_hol_SP1_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR LGM
enm_lgm_SP1 = raster("BRE_LGM_avg.asc")
crs(enm_lgm_SP1) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_lgm_SP1)
enm_lgm_SP1 = 1-enm_lgm_SP1
plot(enm_lgm_SP1)
enm_lgm_SP1_tr = transition(enm_lgm_SP1, transitionFunction=mean, 16)
enm_lgm_SP1_trC = geoCorrection(enm_lgm_SP1_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR LIG
enm_lig_SP1 = raster("BRE_LIG_avg.asc")
crs(enm_lig_SP1) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_lig_SP1)
enm_lig_SP1 = 1-enm_lig_SP1
plot(enm_lig_SP1)
enm_lig_SP1_tr = transition(enm_lig_SP1, transitionFunction=mean, 16)
enm_lig_SP1_trC = geoCorrection(enm_lig_SP1_tr, type="c", multpl=FALSE, scl=TRUE)

#Need to remove the row and col names
points_SP1 = as.matrix(pop_coord_SP1) #points
rownames(points_SP1) = NULL
colnames(points_SP1) = NULL

#verify if all points are inside the map
plot(enm_cur_SP1)
points(points_SP1[,1], points_SP1[,2])

#calculate the cost distance for SP1
cosDIST_cur_SP1 = costDistance(enm_cur_SP1_trC, points_SP1)
cosDIST_hol_SP1 = costDistance(enm_hol_SP1_trC, points_SP1)
cosDIST_lgm_SP1 = costDistance(enm_lgm_SP1_trC, points_SP1)
cosDIST_lig_SP1 = costDistance(enm_lig_SP1_trC, points_SP1)

#verify if distance are different:
cosDIST_cur_SP1
cosDIST_hol_SP1
cosDIST_lgm_SP1
cosDIST_lig_SP1


###################
####### SP2 #######
###################

# ENM FOR PRESENT-DAY
enm_cur_SP2 = raster("SIM_avg.asc")
crs(enm_cur_SP2) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_cur_SP2)
enm_cur_SP2 = 1-enm_cur_SP2
plot(enm_cur_SP2)
enm_cur_SP2_tr = transition(enm_cur_SP2, transitionFunction=mean, 16)
enm_cur_SP2_trC = geoCorrection(enm_cur_SP2_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR HOLOCENE
enm_hol_SP2 = raster("SIM_HOL_avg.asc")
crs(enm_hol_SP2) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_hol_SP2)
enm_hol_SP2 = 1-enm_hol_SP2
plot(enm_hol_SP2)
enm_hol_SP2_tr = transition(enm_hol_SP2, transitionFunction=mean, 16)
enm_hol_SP2_trC = geoCorrection(enm_hol_SP2_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR LGM
enm_lgm_SP2 = raster("SIM_LGM_avg.asc")
crs(enm_lgm_SP2) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_lgm_SP2)
enm_lgm_SP2 = 1-enm_lgm_SP2
plot(enm_lgm_SP2)
enm_lgm_SP2_tr = transition(enm_lgm_SP2, transitionFunction=mean, 16)
enm_lgm_SP2_trC = geoCorrection(enm_lgm_SP2_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR LIG
enm_lig_SP2 = raster("SIM_LIG_avg.asc")
crs(enm_lig_SP2) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_lig_SP2)
enm_lig_SP2 = 1-enm_lig_SP2
plot(enm_lig_SP2)
enm_lig_SP2_tr = transition(enm_lig_SP2, transitionFunction=mean, 16)
enm_lig_SP2_trC = geoCorrection(enm_lig_SP2_tr, type="c", multpl=FALSE, scl=TRUE)

#Need to remove the row and col names
points_SP2 = as.matrix(pop_coord_SP2) #points
rownames(points_SP2) = NULL
colnames(points_SP2) = NULL

#verify if all points are inside the map
plot(enm_cur_SP2)
points(points_SP2[,1], points_SP2[,2])

#calculate the cost distance for SP1
cosDIST_cur_SP2 = costDistance(enm_cur_SP2_trC, points_SP2)
cosDIST_hol_SP2 = costDistance(enm_hol_SP2_trC, points_SP2)
cosDIST_lgm_SP2 = costDistance(enm_lgm_SP2_trC, points_SP2)
cosDIST_lig_SP2 = costDistance(enm_lig_SP2_trC, points_SP2)

#verify if distance are different:
cosDIST_cur_SP2
cosDIST_hol_SP2
cosDIST_lgm_SP2
cosDIST_lig_SP2


###################
####### SP3 #######
###################

# ENM FOR PRESENT-DAY
enm_cur_SP3 = raster("STE_avg.asc")
crs(enm_cur_SP3) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_cur_SP3)
enm_cur_SP3 = 1-enm_cur_SP3
plot(enm_cur_SP3)
enm_cur_SP3_tr = transition(enm_cur_SP3, transitionFunction=mean, 16)
enm_cur_SP3_trC = geoCorrection(enm_cur_SP3_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR HOLOCENE
enm_hol_SP3 = raster("STE_HOL_avg.asc")
crs(enm_hol_SP3) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_hol_SP3)
enm_hol_SP3 = 1-enm_hol_SP3
plot(enm_hol_SP3)
enm_hol_SP3_tr = transition(enm_hol_SP3, transitionFunction=mean, 16)
enm_hol_SP3_trC = geoCorrection(enm_hol_SP3_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR LGM
enm_lgm_SP3 = raster("STE_LGM_avg.asc")
crs(enm_lgm_SP3) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_lgm_SP3)
enm_lgm_SP3 = 1-enm_lgm_SP3
plot(enm_lgm_SP3)
enm_lgm_SP3_tr = transition(enm_lgm_SP3, transitionFunction=mean, 16)
enm_lgm_SP3_trC = geoCorrection(enm_lgm_SP3_tr, type="c", multpl=FALSE, scl=TRUE)

# ENM FOR LIG
enm_lig_SP3 = raster("STE_LIG_avg.asc")
crs(enm_lig_SP3) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(enm_lig_SP3)
enm_lig_SP3 = 1-enm_lig_SP3
plot(enm_lig_SP3)
enm_lig_SP3_tr = transition(enm_lig_SP3, transitionFunction=mean, 16)
enm_lig_SP3_trC = geoCorrection(enm_lig_SP3_tr, type="c", multpl=FALSE, scl=TRUE)

#Need to remove the row and col names
points_SP3 = as.matrix(pop_coord_SP3) #points
rownames(points_SP3) = NULL
colnames(points_SP3) = NULL

#verify if all points are inside the map
plot(enm_cur_SP3)
points(points_SP3[,1], points_SP3[,2])

#calculate the cost distance for SP1
cosDIST_cur_SP3 = costDistance(enm_cur_SP3_trC, points_SP3)
cosDIST_hol_SP3 = costDistance(enm_hol_SP3_trC, points_SP3)
cosDIST_lgm_SP3 = costDistance(enm_lgm_SP3_trC, points_SP3)
cosDIST_lig_SP3 = costDistance(enm_lig_SP3_trC, points_SP3)

#verify if distance are different:
cosDIST_cur_SP3
cosDIST_hol_SP3
cosDIST_lgm_SP3
cosDIST_lig_SP3


#####################################################
## 10. MANTEL TEST WITH GEOGRAPHIC DISTANCE MATRIX ##
#####################################################

#CALCULATE MANTEL TEST FOR IBD:
mantel_sp1 = mantel.rtest(fst_DIST_SP1, log(geoDIST_SP1), 10000)
mantel_sp2 = mantel.rtest(fst_DIST_SP2, log(geoDIST_SP2), 10000)
mantel_sp3 = mantel.rtest(fst_DIST_SP3, log(geoDIST_SP3), 10000)


#CREATE TABLES TO SAVE THE RESULTS
sp1_mantel = cbind (mantel_sp1$obs, mantel_sp1$pvalue)
sp2_mantel = cbind (mantel_sp2$obs, mantel_sp2$pvalue)
sp3_mantel = cbind (mantel_sp3$obs, mantel_sp3$pvalue)


res_mantel = as.data.frame(rbind(sp1_mantel, sp2_mantel, sp3_mantel))
colnames(res_mantel) = c("r", "p-value")
rownames(res_mantel) = c("P. brevicauda", "P. simonsi", "P. steerei")

write.csv(res_mantel, "MATEL_geoDIST.csv", row.names = T)

#Save results as PDF
pdf("./FST_IBD_geoDIST.pdf", onefile=FALSE)
plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(0,15), ylim=c(0, 1));

summary(fit_sp1)
fit_sp1 = lm(as.vector(fst_DIST_SP1) ~ as.vector(log(geoDIST_SP1)))
abline(fit_sp1, lwd=3, col='slategray4')
fit_sp2 = lm(as.vector(fst_DIST_SP2) ~ as.vector(log(geoDIST_SP2)))
abline(fit_sp2, lwd=3, col='red')
fit_sp3 = lm(as.vector(fst_DIST_SP3) ~ as.vector(log(geoDIST_SP3)))
abline(fit_sp3, lwd=3, col='blue')

points(log(geoDIST_SP1), fst_DIST_SP1,  col = 'black', bg = 'slategray4', cex = 2.5, pch=21)
points(geoDIST_SP2, fst_DIST_SP2, col = 'black', bg = 'red', cex = 2.5, pch=21)
points(geoDIST_SP3, fst_DIST_SP3, col = 'black', bg = 'blue', cex = 2.5, pch=21)

axis(1, at=seq(0, 15, by=5), cex.axis=1.5);
axis(2, at=seq(0, 1, by=0.2), cex.axis=1.5, las=1);
mtext(side=1, text='Geographic Distances (ln)', line=2.8, cex=1.5)
mtext(side=2, text='FST/(1-FST)', line=3.3 , cex=1.5)
legend= c("P. brevicauda", "P. simonsi", "P. steerei")
col=c("slategray4",'red','blue')
legend("topleft", legend = legend, cex=1.5, bty="n",  col = "black", pt.bg = col , pch= 21)
dev.off()


###################################################################
## 11. MANTEL TEST WITH GEOGRAPHIC DISTANCE MATRIX - PERMUTATION ##
###################################################################

###################
####### SP1 #######
###################

pop_counter = 6

sp1_mantel_perm = matrix(NA,pop_counter,2)

for(i in 1:pop_counter) {
  temp_coord = pop_coord_SP1[-i,]
  mDIST_temp = earth.dist(temp_coord, dist = TRUE)
  mDIST_temp = as.dist(mDIST_temp)
  fst_cor_temp = fst_pop_sp1[-i,-i]
  fst_cor_temp = as.dist(fst_cor_temp)
  mt = mantel.rtest(fst_cor_temp, log(mDIST_temp), 10000)
  sp1_mantel_perm[i, ] = cbind (mt$obs, mt$pvalue)
}

sp1_mantel_perm
colnames(sp1_mantel_perm) = c("r", "p-value")
rownames(sp1_mantel_perm) = rownames(fst_pop_sp1)

write.csv(sp1_mantel_perm, "MATEL_geoDIST_PERM_bre.csv", row.names = T)


###################
####### SP2 #######
###################

pop_counter = 7

sp2_mantel_perm = matrix(NA,pop_counter,2)

for(i in 1:pop_counter) {
  temp_coord = pop_coord_SP2[-i,]
  mDIST_temp = earth.dist(temp_coord, dist = TRUE)
  mDIST_temp = as.dist(mDIST_temp)
  fst_cor_temp = fst_pop_sp2[-i,-i]
  fst_cor_temp = as.dist(fst_cor_temp)
  mt = mantel.rtest(fst_cor_temp, log(mDIST_temp), 10000)
  sp2_mantel_perm[i, ] = cbind (mt$obs, mt$pvalue)
}

sp2_mantel_perm
colnames(sp2_mantel_perm) = c("r", "p-value")
rownames(sp2_mantel_perm) = rownames(fst_pop_sp2)

write.csv(sp2_mantel_perm, "MATEL_geoDIST_PERM_sim.csv", row.names = T)


###################
####### SP3 #######
###################

pop_counter = 6

sp3_mantel_perm = matrix(NA,pop_counter,2)

for(i in 1:pop_counter) {
  temp_coord = pop_coord_SP3[-i,]
  mDIST_temp = earth.dist(temp_coord, dist = TRUE)
  mDIST_temp = as.dist(mDIST_temp)
  fst_cor_temp = fst_pop_sp3[-i,-i]
  fst_cor_temp = as.dist(fst_cor_temp)
  mt = mantel.rtest(fst_cor_temp, log(mDIST_temp), 10000)
  sp3_mantel_perm[i, ] = cbind (mt$obs, mt$pvalue)
}

sp3_mantel_perm
colnames(sp3_mantel_perm) = c("r", "p-value")
rownames(sp3_mantel_perm) = rownames(fst_pop_sp3)

write.csv(sp3_mantel_perm, "MATEL_geoDIST_PERM_ste.csv", row.names = T)

mantel.rtest

######################################################################################
## 12. PARTIAL MANTEL TEST WITH GEOGRAPHIC DISTANCE MATRIX AND COST DISTANCE MATRIX ##
######################################################################################

###################
####### SP1 #######
###################

mantel_partial_cur_SP1 = mantel.partial(cosDIST_cur_SP1, fst_DIST_SP1, log(geoDIST_SP1), method="pearson", permutations=999)
cur_sp1 = c(mantel_partial_cur_SP1$statistic, mantel_partial_cur_SP1$signif, mantel_partial_cur_SP1$permutations)

###################
####### SP2 #######
###################

mantel_partial_cur_SP2 = mantel.partial(cosDIST_cur_SP2, fst_DIST_SP2, log(geoDIST_SP2), method="pearson", permutations=10000)
cur_sp2 = c(mantel_partial_cur_SP2$statistic, mantel_partial_cur_SP2$signif, mantel_partial_cur_SP2$permutations)

###################
####### SP3 #######
###################

mantel_partial_cur_SP3 = mantel.partial(cosDIST_cur_SP3, fst_DIST_SP3, log(geoDIST_SP3), method="pearson", permutations=10000)
cur_sp3 = c(mantel_partial_cur_SP3$statistic, mantel_partial_cur_SP3$signif, mantel_partial_cur_SP3$permutations)

#CREATE TABLES TO SAVE THE RESULTS
sp_partial_mantel = rbind (cur_sp1, cur_sp2, cur_sp3)
colnames(sp_partial_mantel) = c("r", "p-value", "Permutations")
rownames(sp_partial_mantel) = c("P. brevicauda", "P. simonsi", "P. steerei")

write.csv(sp_partial_mantel, "PARTIAL_MATEL_costDIST.csv", row.names = T)




###################################################################
## 13. MANTEL TEST WITH GEOGRAPHIC DISTANCE MATRIX - PERMUTATION ##
###################################################################

###################
####### SP1 #######
###################

pop_counter = 6

sp1_partial_mantel_perm = matrix(NA,pop_counter,2)

for(i in 1:pop_counter) {
  temp_coord = pop_coord_SP1[-i,]
  mDIST_temp = earth.dist(temp_coord, dist = TRUE)
  mDIST_temp = as.dist(mDIST_temp)
  fst_cor_temp = fst_pop_sp1[-i,-i]
  fst_cor_temp = as.dist(fst_cor_temp)
  cosDIST_cur_SP1_temp = costDistance(enm_cur_SP1_trC, points_SP1[-i, ])
  mt = mantel.partial(cosDIST_cur_SP1_temp, fst_cor_temp, log(mDIST_temp), method="pearson", permutations=999)
  sp1_partial_mantel_perm[i, ] = cbind (mt$statistic, mt$signif)
}

sp1_partial_mantel_perm
colnames(sp1_partial_mantel_perm) = c("r", "p-value")
rownames(sp1_partial_mantel_perm) = rownames(fst_pop_sp1)

write.csv(sp1_partial_mantel_perm, "Partial_MATEL_geoDIST_PERM_bre.csv", row.names = T)


###################
####### SP2 #######
###################

pop_counter = 7

sp2_partial_mantel_perm = matrix(NA,pop_counter,2)

for(i in 1:pop_counter) {
  temp_coord = pop_coord_SP2[-i,]
  mDIST_temp = earth.dist(temp_coord, dist = TRUE)
  mDIST_temp = as.dist(mDIST_temp)
  fst_cor_temp = fst_pop_sp2[-i,-i]
  fst_cor_temp = as.dist(fst_cor_temp)
  cosDIST_cur_SP2_temp = costDistance(enm_cur_SP2_trC, points_SP2[-i, ])
  mt = mantel.partial(cosDIST_cur_SP2_temp, fst_cor_temp, log(mDIST_temp), method="pearson", permutations=999)
  sp2_partial_mantel_perm[i, ] = cbind (mt$statistic, mt$signif)
}

sp2_partial_mantel_perm
colnames(sp2_partial_mantel_perm) = c("r", "p-value")
rownames(sp2_partial_mantel_perm) = rownames(fst_pop_sp2)

write.csv(sp2_partial_mantel_perm, "Partial_MATEL_geoDIST_PERM_sim.csv", row.names = T)


###################
####### SP3 #######
###################

pop_counter = 6

sp3_partial_mantel_perm = matrix(NA,pop_counter,2)

for(i in 1:pop_counter) {
  temp_coord = pop_coord_SP3[-i,]
  mDIST_temp = earth.dist(temp_coord, dist = TRUE)
  mDIST_temp = as.dist(mDIST_temp)
  fst_cor_temp = fst_pop_sp3[-i,-i]
  fst_cor_temp = as.dist(fst_cor_temp)
  cosDIST_cur_SP3_temp = costDistance(enm_cur_SP3_trC, points_SP3[-i, ])
  mt = mantel.partial(cosDIST_cur_SP3_temp, fst_cor_temp, log(mDIST_temp), method="pearson", permutations=999)
  sp3_partial_mantel_perm[i, ] = cbind (mt$statistic, mt$signif)
}

sp3_partial_mantel_perm
colnames(sp3_partial_mantel_perm) = c("r", "p-value")
rownames(sp3_partial_mantel_perm) = rownames(fst_pop_sp3)

write.csv(sp3_partial_mantel_perm, "Partial_MATEL_geoDIST_PERM_ste.csv", row.names = T)



######################################################################################
## 13. dbRDA ANALYSIS TEST WITH GEOGRAPHIC DISTANCE MATRIX AND COST DISTANCE MATRIX ##
######################################################################################

#transform geographic distance matrix in CCA scores
geoDIST_SP1_pcnm = pcnm(log(geoDIST_SP1))
geoDIST_SP2_pcnm = pcnm(geoDIST_SP2)
geoDIST_SP3_pcnm = pcnm(geoDIST_SP3)

###################
####### SP1 #######
###################

dbrda_Env_Geo_SP1 = capscale(fst_DIST_SP1 ~ PC1_SP1 + Condition(scores(geoDIST_SP1_pcnm)))
sig_Env_Geo_SP1 =  anova(dbrda_Env_Geo_SP1, by="terms")
RES_Env_Geo_SP1 = cbind(sig_Env_Geo_SP1$Df, sig_Env_Geo_SP1$SumOfSqs, sig_Env_Geo_SP1$F, sig_Env_Geo_SP1$`Pr(>F)`)
colnames(RES_Env_Geo_SP1) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Env_Geo_SP1) = c("Model FST ~ Env|Geo", "Residual")

dbrda_Env_SP1 = capscale(fst_DIST_SP1 ~ PC1_SP1)
sig_Env_SP1 =  anova(dbrda_Env_SP1)
RES_Env_SP1 = cbind(sig_Env_SP1$Df, sig_Env_SP1$SumOfSqs, sig_Env_SP1$F, sig_Env_SP1$`Pr(>F)`)
colnames(RES_Env_SP1) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Env_SP1) = c("Model FST ~ Env", "Residual")

dbrda_Geo_SP1 = capscale(fst_DIST_SP1 ~ scores(geoDIST_SP1_pcnm))
sig_Geo_SP1 =  anova(dbrda_Geo_SP1)
RES_Geo_SP1 = cbind(sig_Geo_SP1$Df, sig_Geo_SP1$SumOfSqs, sig_Geo_SP1$F, sig_Geo_SP1$`Pr(>F)`)
colnames(RES_Geo_SP1) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Geo_SP1) = c("Model FST ~ Geo", "Residual")

RES_dbRDA_SP1 = rbind (RES_Geo_SP1, RES_Env_SP1, RES_Env_Geo_SP1)
write.csv(RES_dbRDA_SP1, "RES_dbRDA_SP1_bre.csv", row.names = T)


###################
####### SP2 #######
###################

dbrda_Env_Geo_SP2 = capscale(fst_DIST_SP2 ~ PC1_SP2 + Condition(scores(geoDIST_SP2_pcnm)))
sig_Env_Geo_SP2 =  anova(dbrda_Env_Geo_SP2)
RES_Env_Geo_SP2 = cbind(sig_Env_Geo_SP2$Df, sig_Env_Geo_SP2$SumOfSqs, sig_Env_Geo_SP2$F, sig_Env_Geo_SP2$`Pr(>F)`)
colnames(RES_Env_Geo_SP2) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Env_Geo_SP2) = c("Model FST ~ Env|Geo", "Residual")

dbrda_Env_SP2 = capscale(fst_DIST_SP2 ~ PC1_SP2)
sig_Env_SP2 =  anova(dbrda_Env_SP2)
RES_Env_SP2 = cbind(sig_Env_SP2$Df, sig_Env_SP2$SumOfSqs, sig_Env_SP2$F, sig_Env_SP2$`Pr(>F)`)
colnames(RES_Env_SP2) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Env_SP2) = c("Model FST ~ Env", "Residual")

dbrda_Geo_SP2 = capscale(fst_DIST_SP2 ~ scores(geoDIST_SP2_pcnm))
sig_Geo_SP2 =  anova(dbrda_Geo_SP2)
RES_Geo_SP2 = cbind(sig_Geo_SP2$Df, sig_Geo_SP2$SumOfSqs, sig_Geo_SP2$F, sig_Geo_SP2$`Pr(>F)`)
colnames(RES_Geo_SP2) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Geo_SP2) = c("Model FST ~ Geo", "Residual")

RES_dbRDA_SP2 = rbind (RES_Geo_SP2, RES_Env_SP2, RES_Env_Geo_SP2)
write.csv(RES_dbRDA_SP2, "RES_dbRDA_SP2_sim.csv", row.names = T)


###################
####### SP3 #######
###################

dbrda_Env_Geo_SP3 = capscale(fst_DIST_SP3 ~ PC1_SP3 + Condition(scores(geoDIST_SP3_pcnm)))
sig_Env_Geo_SP3 =  anova(dbrda_Env_Geo_SP3)
RES_Env_Geo_SP3 = cbind(sig_Env_Geo_SP3$Df, sig_Env_Geo_SP3$SumOfSqs, sig_Env_Geo_SP3$F, sig_Env_Geo_SP3$`Pr(>F)`)
colnames(RES_Env_Geo_SP3) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Env_Geo_SP3) = c("Model FST ~ Env|Geo", "Residual")

dbrda_Env_SP3 = capscale(fst_DIST_SP3 ~ PC1_SP3)
sig_Env_SP3 =  anova(dbrda_Env_SP3)
RES_Env_SP3 = cbind(sig_Env_SP3$Df, sig_Env_SP3$SumOfSqs, sig_Env_SP3$F, sig_Env_SP3$`Pr(>F)`)
colnames(RES_Env_SP3) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Env_SP3) = c("Model FST ~ Env", "Residual")

dbrda_Geo_SP3 = capscale(fst_DIST_SP3 ~ scores(geoDIST_SP3_pcnm))
sig_Geo_SP3 =  anova(dbrda_Geo_SP3)
RES_Geo_SP3 = cbind(sig_Geo_SP3$Df, sig_Geo_SP3$SumOfSqs, sig_Geo_SP3$F, sig_Geo_SP3$`Pr(>F)`)
colnames(RES_Geo_SP3) = c("D.f.", "Sum Of Squares", "F", "p-value")
rownames(RES_Geo_SP3) = c("Model FST ~ Geo", "Residual")

RES_dbRDA_SP3 = rbind (RES_Geo_SP3, RES_Env_SP3, RES_Env_Geo_SP3)
write.csv(RES_dbRDA_SP3, "RES_dbRDA_SP3_ste.csv", row.names = T)