###############################################################################
######################## SPECIES DISTRIBUTION MODELING ########################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
######################### STEP 07: MODELS COMPLEXITY ##########################

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
#A. ENVIRONMENTAL LAYERS SELECTED IN THE STEP 06.
#B. UNBIASED OCCURRENCE POINTS IN THE STEP 05.


##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("rJava", "dismo", "ENMeval", "plyr" , "SDMTools", "pacman","tcltk"))
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("rJava", "dismo", "ENMeval", "plyr", "SDMTools","tcltk")


##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())
#Any computer
setwd('C:/Users/Johnny/Downloads/bio_10m_bil') # change the path


## 5. YOU NEED TO COPY THE FILE maxent.jar IN THE FOLDER R/3.4/dismo/java/. VERIFY IF YOU DID IT CORRECTLY.
jar = paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar) & require(rJava)){
  print("maxent.jar is in the correct place")
}else{
  print("maxent.jar file needs to be moved to dismo/java folder in R")
}


## 6. LOAD THE OCCURRENCE POINTS AND THE LAYERS SELECTED
#unbiased points
occ = read.csv("ROB_A_UNBIAS.csv")
# if you have more that long/lat columns - we need to remove all other columns
occ = occ[,2:3]
head(occ)

#unbiased layers
predictors = stack(sapply(list.files(pattern='asc$', recursive = F), raster))
predictors

##7. USE ENMeval TO CHOOSE THE FEATURE CLASS AND BETA VALUES
#increase the memory
options(java.parameters = "-Xmx1g" )

#create background points
bg_pnts = randomPoints(predictors, n = 10000)

#testing the parameters
EVALme = ENMevaluate (occ, predictors, bg.coords = bg_pnts, RMvalues = seq(0.5, 3, 0.5), fc = c("L", "LQ", "LQP", "H", "T", "LQH", "LQHP", "LQHPT"), method = "block", clamp = TRUE, overlap= FALSE, rasterPreds = TRUE, parallel = TRUE, numCores = 3)

dir.create("ENMeval") 
setwd ("./ENMeval")

#save the results
saveRDS(EVALme, file = paste("Results_eval_", "ROB_A", ".rds", sep= ""))

#read the RDS file
#EVALme = readRDS(paste("Results_eval_", species, ".rds", sep= ""))

#Best parameters
res = arrange(EVALme@results, EVALme@results$AICc)
best_res = res[1,c(2:3, 14)]
best_res

#save the table
write.csv(res, file="Best_Parameters_ROB_A.csv", row.names = F)
setwd ("..")
