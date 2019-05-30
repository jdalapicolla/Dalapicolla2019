###############################################################################
######################## SPECIES DISTRIBUTION MODELING ########################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
#################### STEP 06: SELECTING ENVIRONMENTAL LAYERS ##################

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
#A. ENVIRONMENTAL LAYERS EDITED IN THE STEP 04.
#B. UNBIASED OCCURRENCE POINTS IN THE STEP 05.


##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("raster", "rgdal", "pacman","tcltk","vegan", "dismo", "car"))
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("raster", "rgdal", "tcltk", "vegan", "dismo", "car")

#####################################################################################
########################### CORRELATION APPROACH ####################################
#####################################################################################

#####################################################################################
####   FOLLOWING THE METHOD "PCA-based clustering" FROM Dormann et al. (2013)    ####
####   Collinearity a review of methods to deal with it and a simulation study   ####
#### evaluating their performance. ECOGRAPHY. KEY REFERENCES: Booth et al. (1994)####
#####################################################################################

##4. CHOOSE A FOLDER FOR ANALYSIS. THE LAYERS FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())
#Any computer
setwd('C:/Users/Johnny/Downloads/bio_10m_bil') # change the path


##5. CORRELATION AMONG LAYERS:
#load all layers
layers = stack(sapply(list.files(pattern='asc$', recursive = F), raster))
#calculete correlation:
correlation = layerStats(layers, "pearson", na.rm = T)
correlation = correlation$`pearson correlation coefficient`
correlation = round(correlation,2)
#save the result
dir.create("correlation") 
setwd("./correlation") 
write.csv(correlation, file = "correlacao_camadas_R.csv")


##6. PCA OF THE LAYERS
#extraction the values of the layers
raw_values = values(layers)
#remove the NA. If the layers are all from the same source, as WorldClim, the lines will have NA in the same positions in all layers
raw_values = na.omit(raw_values)
raw_values = as.data.frame(raw_values)
head(raw_values)
summary(raw_values) #see if there are NAs
#transform the data: standard because there are different units
values_trans = decostand(raw_values, method="standardize")
summary(values_trans)
#PCA
pca = prcomp (values_trans)
#% PC
std_list = pca$sdev^2
std_list_pct = std_list / sum (std_list) * 100
std_list_pct_sum = round(std_list_pct, 2)
write.csv(std_list_pct_sum, file = "PC_variance.csv")
# layers and PCs
contr= as.data.frame(pca$rotation)
write.csv(contr, file = "variables_PCA.csv")


##7. SUM OF THE PCS TO DISCOVERY THE NUMEBER OF PC THAT EXPLAIN 90% OF THE VARIATION
n.pc = length(std_list_pct_sum)
sum_pca = matrix(NA, n.pc, 2)
for (i in 1:n.pc) {sum_pca[i,] = c(i, sum(std_list_pct_sum[1:i]))}
colnames(sum_pca) = c("PCs","% Variance")
sum_pca
pc = 4 #number of PC with 90.07% of the variance
write.csv(sum_pca, file = "PC_variance_sum.csv", row.names = F)


##8.CHOOSE FOR THE PC WITH AROUND 90% OF VARIANCE (Wang et al. 2016 - Study on selecting sensitive environmental variables in modelling species spatial distribution), THE VARIABLES THAT HAVE VALUES BIGGER >0.32 - 10% OF THE VARIANCE (Tabachnick and Fidell 1989). THEY WILL FORM THE CLUSTERS
tab = abs(contr) #don't matter the direction
lista = list()
for (i in 1:pc)
{
  linhas=tab[tab[, i] > 0.32, ]
  lista[[i]]=row.names(linhas)
}

linhas.resultado=unlist(lista)
linhas.resultado=unique(linhas.resultado)
write.table(linhas.resultado, file = "Variables_10%_PCA.csv", sep = "\n", row.names = F, col.names = F)


##9. PERFORM A GLM WITH OCCURRENCE AND BACKGROUND POINTS TO CHOOSE THE VARIABLES FROM THE CLUSTER THAT CAN EXPLAIN THE OCCURRENCE POINTS, FOLLOWING THE R2 VALUES.
setwd("..")
predictors = stack(linhas.resultado) #you should be in the folder with your variables
projection(predictors) = CRS('+proj=longlat +datum=WGS84')
predictors
pontos = read.csv("ROB_A_UNBIAS.csv") #occurrence points clean and unbiased
#only lat/lon columns
pontos = pontos[,2:3]
head(pontos)
setwd("./correlation")

#GLM
valores.varia = extract(predictors, pontos, na.rm = T)
valores.varia = valores.varia[which(valores.varia[ ,1] !=0), ]
set.seed(0)
backgr = randomPoints(predictors, 500)
absvals = extract(predictors, backgr)
pb = c(rep(1, nrow(valores.varia)), rep(0, nrow(absvals)))
sdmdata = data.frame(cbind(pb, rbind(valores.varia, absvals)))
head(sdmdata)
varia2=colnames(sdmdata)
varia2=varia2[-1]
tabela= as.data.frame(matrix(NA, length(varia2), 2))
colnames(tabela) = c("R2", "p-value")
row.names(tabela) = varia2
loop=0
for (k in varia2)
{
  loop= loop+1
  m1 <- lm(pb ~ sdmdata[,k], data=sdmdata)
  test=summary(m1)
  r=test$adj.r.squared
  p.value= test$coefficients[2,4]
  tabela[loop,]= cbind(r, p.value)
}
tabela = round(tabela,3)
tabela.validos = tabela[tabela$`p-value`<0.05, ]
tabela.validos = tabela.validos[ order(tabela.validos[,1], decreasing = T), ]
tabela.validos
write.csv(tabela.validos, "Selected_Variables_ROB_A.csv")
setwd("..")

#My results:
#BIO_05

##10. VERIFY MANUALLY THE CORRELATION AMONG THE SELECTED VARIABLES IN THE FILE "Selected_Variables". IF THERE ARE NO VARIABLES OR FEW YOU CAN ADD THE VARIABLES FROM CLUSTER "Variables_10%_PCA.csv" TO ANALIZE THE CORRELATION "variaveis_10%_PCA.csv". MORE THAN 0.7 OF CORRELATION, ONE OF THE VARIABLE MUST BE ELIMINATED

#My results analyzing all 10% variables:
#BIO_15
#BIO_05
#BIO_08
#BIO_14
#BIO_18
#BIO_12
#BIO_04

#########################################################################################
################################### VIF APPROACH ########################################
#########################################################################################
#Zuur et al. 2010 for regression with morphological and ecological data, no exactly with SDM or ENM

dir.create("vif") 

predictors = layers
projection(predictors) = CRS('+proj=longlat +datum=WGS84')
pontos = read.csv("ROB_A_UNBIAS.csv")
#only lat/lon columns
pontos = pontos[,2:3]
head(pontos)
setwd("./vif")

#create a table to GLM
valores.varia = extract(predictors, pontos, na.rm = T)
valores.varia = valores.varia[which(valores.varia[ ,1] !=0), ]
set.seed(0)
backgr = randomPoints(predictors, 500)
absvals = extract(predictors, backgr)
pb = c(rep(1, nrow(valores.varia)), rep(0, nrow(absvals)))
sdmdata = data.frame(cbind(pb, rbind(valores.varia, absvals)))
head(sdmdata)

#GLM with all 19 variables
m1 <- lm(pb ~ sdmdata$bio1.asc+sdmdata$bio2.asc+sdmdata$bio3.asc+sdmdata$bio4.asc+sdmdata$bio5.asc+sdmdata$bio6.asc+sdmdata$bio7.asc+sdmdata$bio8.asc+sdmdata$bio9.asc+sdmdata$bio10.asc+sdmdata$bio11.asc+sdmdata$bio12.asc+sdmdata$bio13.asc+sdmdata$bio14.asc+sdmdata$bio15.asc+sdmdata$bio16.asc+sdmdata$bio17.asc+sdmdata$bio18.asc+sdmdata$bio19.asc, data=sdmdata)

#calculate the vif for all model:
vif(m1)

#save the results
write.csv(vif(m1), "VIF_initial_ROB_A.csv")

#if this error appears "Error in vif.default(m1) : there are aliased coefficients in the model". In this context, ''alias'' refers to the variables that are linearly dependent on others (i.e. cause perfect multicollinearity).You should eliminate one of the perfect co-variables. Use the alias() function:
alias(m1) #co-variables are that ones with -1 and 1 in the table. This case bio6  and bio5. I eliminated the bio6. No special reason, you can choose for a biological reason or because it is important in the environmental analyzed in the study.

#GLM with 18 variables
m1 <- lm(pb ~ sdmdata$bio1.asc+sdmdata$bio2.asc+sdmdata$bio3.asc+sdmdata$bio4.asc+sdmdata$bio5.asc+sdmdata$bio7.asc+sdmdata$bio8.asc+sdmdata$bio9.asc+sdmdata$bio10.asc+sdmdata$bio11.asc+sdmdata$bio12.asc+sdmdata$bio13.asc+sdmdata$bio14.asc+sdmdata$bio15.asc+sdmdata$bio16.asc+sdmdata$bio17.asc+sdmdata$bio18.asc+sdmdata$bio19.asc, data=sdmdata)

#calculate the vif for all model:
vif(m1) # I have to eliminate the biggest number, in the case bio10.

write.csv(vif(m1), "VIF_initial_ROB_A.csv") #if there is a error in the first try

##GLM with 17 variables
m1 <- lm(pb ~ sdmdata$bio1.asc+sdmdata$bio2.asc+sdmdata$bio3.asc+sdmdata$bio4.asc+sdmdata$bio5.asc+sdmdata$bio8.asc+sdmdata$bio9.asc+sdmdata$bio11.asc+sdmdata$bio12.asc+sdmdata$bio13.asc+sdmdata$bio14.asc+sdmdata$bio15.asc+sdmdata$bio16.asc+sdmdata$bio17.asc+sdmdata$bio18.asc+sdmdata$bio19.asc, data=sdmdata)

#calculate the vif for all model:
vif(m1) # I have to eliminate the biggest number, in the case bio1.

######################################################################################
#use this one to eliminate one variable per round:
m1 <- lm(pb ~ sdmdata$bio1.asc+sdmdata$bio2.asc+sdmdata$bio3.asc+sdmdata$bio4.asc+sdmdata$bio5.asc+sdmdata$bio6.asc+sdmdata$bio7.asc+sdmdata$bio8.asc+sdmdata$bio9.asc+sdmdata$bio10.asc+sdmdata$bio11.asc+sdmdata$bio12.asc+sdmdata$bio13.asc+sdmdata$bio14.asc+sdmdata$bio15.asc+sdmdata$bio16.asc+sdmdata$bio17.asc+sdmdata$bio18.asc+sdmdata$bio19.asc, data=sdmdata)
#calculate the vif for all model:
vif(m1)
######################################################################################

#Keep eliminating variables after variable, until all variable have vif<2. In my case I selected 5 valiables:
#BIO_04
#BIO_05
#BIO_12
#BIO_15
#BIO_18

#save the results
write.csv(vif(m1), "VIF_final_ROB_A.csv")

setwd("..")

######################################################################################
#comparing results:
#CORRELATION: 7 VARIABLES
#BIO_15
#BIO_05
#BIO_08
#BIO_14
#BIO_18
#BIO_12
#BIO_04

#VIF: 5 VARIABLES
#BIO_04
#BIO_05
#BIO_12
#BIO_15
#BIO_18

#Difference was two variables more in Correlation approach (BIO_08 and BIO_14).
#You should test both set of variables and choose one with the best indexes in the next step, complexity of the models.