###############################################################################
######################## SPECIES DELIMITATION ANALYSES ########################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
################################     iBPP    ##################################
######################## STEP 01: MORPHOLOGICAL DATA ##########################

###############################################################################
######################### JERONYMO DALAPICOLLA, 2019 ##########################
###################### LABORATÓRIO DE MAMÍFEROS - LAMA ########################
####################### ESALQ - UNIVERSIDADE DE SAO PAULO #####################
###############################################################################

##1. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/

##2. INPUTS FOR THIS STEP:
#A. MORPHOMETRIC DATA
#B. FUNCTIONS .R IN THE SAME FOLDER. GO TO THE GITHUB TO DOWNLOAD THEM

##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("tcltk","pacman", "dplyr", "adegenet"))
#load the packages
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("tcltk", "dplyr", "adegenet")

##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILE MUST BE THERE! 
#Windows
choose.dir()
#Linux. You will need the tcltk package for this command:
setwd(tk_choose.dir())

##5. LOAD THE FUNCTIONS FOR THIS SCRIPTS. THE 'R' FILES MUST BE IN THE FOLDER CHOSE ON THE PREVIOUS STEP.
source("id.outliers.R")
source("contar.NA.R")
source("normalidade.R")

##6. LOAD THE FILE WITH MORPHOMETRICS DATA. IN MY FILE THERE ARE 39 COLUNMS, 1:17 CATEGORICAL VARIABLES AND 18:39 THE MORPHOMETRICS MEASUREMENTS.  
tab = read.table (file="ibbp_raw_morpho.csv", header = T, sep=",", as.is = TRUE)
#verify if the load of the file is ok.
head(tab)
tail(tab)
summary(tab)

##7. CHOOSE ONLY ADULTS FOR THE ANALYSIS. AGE_PATTON EQUALS TO 8,9, AND 10 FOR PROECHIMYS.
#copy the original data, for security proposals
tab_save = tab

#selecting the adults
tab = tab[tab$age_patton=="8" | tab$age_patton=="9" | tab$age_patton=="10", ]


##8. IDENTIFY THE OUTLIERS INSIDE THE MORPHOGROUPS/SPECIES/SPECIES GROUP AND IN ALL SAMPLES TOGETHER:
#all samples together
out = id.outliers(x=tab, quant=12:33, group=0, id="z", NUMBER=100, visual="biplot", res="MED", csv=T)

#########################################################
#########################################################
#verify if the samples indicated by the test are really outliers, if yes, replace them by "NA". Replace in the original file too. Remember to load the file again, and redo the steps 6 and 7!
#########################################################
#########################################################


##9. TRANSFORM THE DATA IN LOG. DO NOT STANDARD THE DATA. iBPP WILL DO IT LATTER.
#save a copy of raw data
tab_raw = tab

#save a table with log transformated data
tab_log = tab
tab_log[,12:33] = log(tab[,12:33])

dados = tab_log[,12:33]


##10. PCA WITH SIZE
#center = TRUE means centring by the mean, in this case is 0
#scannf	= FALSE, indicates the screeplot should not be displayed
#nf = 3, if scannf = FALSE, an integer must be indicated. It will be the number of kept axes in the PCA
pca_input = dudi.pca(dados, center = TRUE , scannf = FALSE, nf = 20)


##11. REGRESSION TO REMOVE THE SIZE EFFECT:
residuo1 = list(NA)

for(i in 1:length(dados)){
  residuo1[i] = lm(pca_input$li[,1] ~ dados[,i])[2]
}

output1 = matrix(unlist(residuo1), nrow = 479, byrow = FALSE)

dados2 = as.data.frame(output1)

tab_log[,12:33] = dados2



#####################################################
#################### NORMALITY ######################
#####################################################

##12. CALCULATE THE NORMALITY FOR EACH VARIABLE FOR ALL SAMPLES TOGETHER. For Shapiro Test sample size must be between 3 and 5000.

norm = normalidade (x=tab_log, quant=12:33, group=0)

# verify the graphs if they are really not normal.
# I deleteD the variable that are non-normal by Shapiro-test and by qqplot

#verify the column position
colnames(tab_log)

#remove the 
tab_red = tab_log
tab_red = tab_red [,c(-13, -19, -23, -24, -25, -26, -28, -29, -32, -33)]
head (tab_red)


#####################################################
######## CORRELATIONS BETWEEN VARIABLES #############
#####################################################

#13.CALCULATE WHICH VARIABLE HAVE HIGHER CORRELATIONS.
#use only the morphometrics variables
tab2 = tab_red[,12:23]

#calculate the correlation between
resu = round (cor(tab2, use="pairwise.complete.obs", method ="pearson"),2)
write.table(abs(resu), "cor_pearson.csv", sep = ",", quote = F)

#Assumptions For the Pearson r correlation, both variables should be normally distributed (normally distributed variables have a bell-shaped curve).  Other assumptions include linearity and homoscedasticity. Linearity assumes a straight line relationship between each of the two variables and homoscedasticity assumes that data is equally distributed about the regression line.

#open the file and verify the variables with r value greater than 0.5:

#verify the column position
colnames(tab_red)

#remove them
tab_red = tab_red [,c(-16, -17)]
head (tab_red)


#####################################################
################## MISSING DATA #####################
#####################################################

##14.CALCULATE THE NUMBER OF NA'S FOR EACH VARIABLE AND EACH INDIVIDUALS. 
#all samples together
na.group = contar.NA(x=tab_red, var=18:37, group=0)

#verify which variables and individual have more NA's. Delete the individual if it has more replicates, last case, delete a variable.

#verify the column position
colnames(tab_red)

#remove
tab_red = tab_red [-19,] #rows
tab_red = tab_red [,c(-22, -23, -24, -26, -27, -28, -29,-30,-37)] #columns
head (tab_red)

##15.EXPORT THE RESULTS:
#save data:
write.table(tab_red, "morpho_clean_log.csv", sep = ",", quote = F, row.names = F)

#create the morphological input
Individual = tab_red$sample #ID for specimens
Species = dados2$clade_code #Code for putative species, should be the same for genetic data
input = cbind(Individual, Species, tab_red[,18:28]) # columns with morphometric data

#save
write.table(input, file = "IBPP_MORPHO_AN1.txt", sep="\t", quote = FALSE, row.names = FALSE)
