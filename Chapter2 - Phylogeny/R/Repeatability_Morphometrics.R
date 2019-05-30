############################################################################################
######################### Repeatability for Morphometric Data  #############################
############################## Jeronymo Dalapicolla, 2019 ##################################
############################################################################################

############################################################################################
################################## JERONYMO DALAPICOLLA, 2019 ##############################
############################### LABORATORIO DE MAMIFEROS - LAMA ############################
############################## ESALQ - UNIVERSIDADE DE SAO PAULO ###########################
############################################################################################

##1. FOR OTHER STEPS, PLEASE GO TO: https://github.com/jdalapicolla/

##2.INPUTS FOR THIS ANALYSIS:
#A. None in the first step to figure out the number of repeats.
#B. ".csv" file with measurements.

##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("tcltk","pacman", "ICC")) #tcltk is for Linux!
#load the packages
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("tcltk", "ICC")
#OR
install.packages("ICC") # Install ape package if you don't have it already
library(ICC) # Load the ape package

##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILE MUST BE THERE! 
#Windows
choose.dir()
#Linux. You will need the tcltk package for this command:
setwd(tk_choose.dir())


##5. FIGURE OUT THE NUMBER OF SAMPLES FOR THE REPETIBILITY AND HOW MANY TIMES MEASURE THE VARIABLES:
n1 = Nest("h", w = 0.1, ICC = 0.9, k = 5, alpha = 0.05)
n1

#"h" = hypothetical test, no preliminary results
# w = the value of the confidence interval, this value must be as low as possible. Lower values more repeats and sample are need!
# ICC = ICC value I want to find, + or - the confidence interval. Greater values,  better Repeatability
#k = number of repeats I want to do per variable. How many times I want to measure the same variable in the same skull.
#alpha = the p-value for the test (95%)
# In the example above, I would like to know how many skulls I should measure for a hypothetic test, knowing that I will measure them 5 times with a power (acuracy) to be able to detect an ICC of 0.9 + or - 0.1 (0.8 - 1.0). In other words, how many specimens I should measure to the test be capable of detecting an ICC ranging from 0.8 to 1.0 (maximum value). Result = 34. I measured more, 40 skulls, 5 times each. No problem if you measure a higher number.

#DO THE MEASUREMENTS TO THE FOLLOWING STEPS:


##6. READ THE TABLE WITH THE MEASUREMENTS:
dados= read.table("Repeat.csv", header = T, sep=",")
head(dados)

# example: Col = variables; row = 5 repeats per individuals 
#Sample	RB	IOC	ZB	MB	GSL	NL	RL	OL	CD
#MJ002	6.65	9.6	22.29	17.34	46.45	17.66	18.39	11.99	13.51
#MJ002	6.59	9.76	22.31	17.34	46.51	17.98	18.44	11.91	13.19
#MJ002	6.44	9.69	22.27	17.29	46.43	17.89	18.35	11.97	13.75
#MJ002	6.68	9.64	22.33	17.33	46.45	18.07	18.45	11.93	13.38
#MJ002	6.6	9.69	22.28	17.41	46.5	17.97	18.39	11.91	13.43
#MJ007	6.38	9.15	20.03	16.76	41.96	15.98	15.63	11.33	13.6
#MJ007	6.51	9.23	19.99	16.74	42.04	15.84	15.35	11.34	13.36
#MJ007	6.5	9.3	20.03	16.69	41.95	15.71	15.63	11.38	13.08
#MJ007	6.67	9.15	20.03	16.8	41.84	16.23	15.33	11.4	13.04
#MJ007	6.1	9.21	20.05	16.79	42.01	15.88	15.23	11.45	13.09


##7. LOAD THE FUNCTION FOR THIS SCRIPT. THE 'R' FILE MUST BE IN THE FOLDER CHOSE ON THE PREVIOUS STEP.
source("repeatability.R")

##8. CALCULATE THE REPEATIBILITY:
repetibilidade(dados$Sample, dados)

#dados$Sample = columns with specimens ID
#dados = table

#The results will be save in the folder under the name: "Repetibilidade_Tabela.csv".

