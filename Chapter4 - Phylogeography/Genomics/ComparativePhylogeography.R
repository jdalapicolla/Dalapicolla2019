###############################################################################
######################### PHYLOGEOGRAPHIC ANALYSES ############################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2018 ##########################
#################     STEP 02: COMPARATIVE PHYLOGEOGRAPHY     #################
##################  SUMMARY STATISTICS, AMOVA, PCA, PROCRUSTES  ###############

###############################################################################
######################### JERONYMO DALAPICOLLA, 2019 ##########################
###################### LABORATORIO DE MAMIFEROS - LAMA ########################
####################### ESALQ - UNIVERSIDADE DE SAO PAULO #####################
###############################################################################


##1.FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/

##2.INPUTS FOR THIS STEP:
#A. THE FILE ".structure.tsv" CLEANED AFTER POPULATIONS STEP IN STACKS ORGANIZED BY POPULATIONS/GROUPS.
#B. GEOGRAPHICAL COORDENATES FOR INDIVIDUALS IN DECIMAL DEGREES.
#C. SHAPEFILES OR RASTERS FOR MAPS OF THE STUDY AREAS.

##3. INSTALL AND LOAD THE PACKAGES#
#install the packages
#install.packages(c("adegenet", "ade4", "vegan", "raster", "apex", "pegas", "mmod", "poppr", "rgdal", "rgeos", "ggplot2", "reshape2", "ape", "aspace","tcltk"))
#load the packages
library("adegenet")
library("ade4")
library("vegan")
library("raster")
library("apex")
library("pegas")
library("mmod")
library("poppr")
library("rgdal")
library("rgeos")
library("ggplot2")
library("reshape2")
library("ape")
library("aspace")
library("tcltk") #for Linux

#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("adegenet", "ade4", "vegan", "raster", "apex", "pegas", "mmod", "poppr", "rgdal", "rgeos", "ggplot2", "reshape2", "ape", "aspace", "tcltk")

##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())


#######################################################################################
################################# PRE-ANALYSIS ########################################
#######################################################################################

##5A. OPEN THE STRUCTURE FILE (structure.tsv) IN A TEXT EDITOR AND REMOVE THE FIRST LINE WITH THE NAME OF LOCI AND SNPS
##5B. SAVE THE FILE WITH A NEW NAME, WITH EXTENSION ".STR". IF THE FILE IS TOO LARGE YOU CAN USE THE FOLLOWING STEPS##

##5A. Open the file skiping the first to lines of the file and don't read the file as factors:
x = read.table("batch_1.structure.tsv", sep = "\t", skip = 2, as.is = T)
str(x)

#5B. Verify the name of populations and order of individuals, always in the second column:
x$V1
unique(x$V2)

#5C. Replace each population by numbers, not letters is allowed in some softwares:
x$V2[x$V2=="GALVEZ"] = "1"
x$V2[x$V2=="MADRE"] = "2"
x$V2[x$V2=="BOLIVIA"] = "3"
x$V2[x$V2=="AMAZONAS"] = "4"
x$V2[x$V2=="ACRE"] = "5"
x$V2[x$V2=="IQUITOS"] = "6"

#verify
unique(x$V2)

#5D. Save the file as a str file:
write.table(x, file="output_bre.str", sep="\t", col.names = F, row.names = F, quote = F)


#######################################################################################
##################################### ANALYSES ########################################
#######################################################################################


##6.CALCULATE THE NUMBER OF LOCI IN THE STRUCTURE FILE 
#load the file as a table
x = read.table('output_bre.str', header = F)
#the first colunm is the name of samples, and the second the population number. One individual is representing in the structure file by two lines (each line is an allele), so the loci = (colunms number - 2).
n.loc = (ncol(x) - 2)
n.ind = (nrow(x)/2)


##7. INPUT THE STRUCTURE FILE AS GENIND OBJECT
#change the name of file, number of samples "n.ind", number of loci "n.loc".
#col.lab = number the column with the sample names
#col.pop = number the column with the population names
#onerowperind = if each line represents one individual (TRUE) or if the alleles are in differents lines (FALSE)
#NA.char = character in the matrix that represents NA (missing data), in the PLINK and POPULATIONS files is 0, the default is -9. BE CAREFUL!
#ask = if the function should ask for optional informations about the dataset. Interactive.
input = read.structure('output_bre.str', n.ind = 16, n.loc = 5050, col.lab = 1, col.pop=2, onerowperind = F, row.marknames = 0, NA.char = 0, ask = FALSE)

#########################################
####### SUMMARY STATISTICS ##############
#########################################

#8A.REMOVE MISSING DATA FROM GENIND OBJECT:
input2 = missingno(input, type = "loci")
nloci = length(input2$loc.fac)/2 #number of remained loci/SNP

##8B. CALCULUTE THE SUMMARY STATTISTICS
#for global analysis. For Phi_st No missing data allowed!!!!
FST_stat = diff_stats(input2, phi_st = T)
write.csv(data.frame(FST_stat$global, nloci), file = "FST_global_bre.csv")


#plotting the results
#create the data frame
per.locus = melt(FST_stat$per.locus, varnames = c("Locus", "Statistic"))
stats = c("Hs", "Ht", "Gst", "Gprime_st", "D", "D", "Phi_st")
glob = data.frame(Statistic = stats, value = FST_stat$global)
head(per.locus)
head(glob)

#create the theme for ggplot graphs
theme_pca = theme(legend.text = element_text(face = "italic",
                                             colour="black",
                                             family = "Helvetica",
                                             size = rel(1)), 
                  axis.title = element_text(colour="black",
                                            family = "Helvetica",
                                            size = rel(1.2)), 
                  axis.text = element_text(family = "Helvetica",
                                           colour = "black",
                                           size = rel(1)), 
                  axis.line = element_line(size = 1,colour = "black"), 
                  axis.ticks = element_line(colour="black",size = rel(1)),
                  
                  panel.grid.minor = element_blank(), 
                  panel.background = element_rect(fill = "whitesmoke"), 
                  panel.grid.major = element_line(colour="black",size = rel(0.2), linetype = "dotted"),
                  legend.key = element_blank(), 
                  legend.title = element_text(colour = "black",
                                              size = rel(1.5),
                                              family = "Helvetica"), 
                  plot.title = element_text(colour = "black",
                                            face = "bold",
                                            hjust = 0.5, #alingment
                                            size = rel(1.7),
                                            family = "Helvetica"))

#plot
pp = ggplot(per.locus, aes(x = Statistic, y = value)) +
  geom_boxplot() +
  geom_point() +
  geom_point(size = rel(3), color = "red", data = glob) +
  ggtitle("Estimates of Population Differentiation") +
  theme_pca

plot(pp) #verify graph

#save as svg
ggsave("FST_stat_bre.svg", plot = plot(pp), device = "svg", path = NULL, scale = 1, width = 20, height = 15, units = "cm", dpi = 300, limitsize = TRUE)

#Hs and Ht: heterozygosity expected for this population with and without the sub-populations defined in the input data respectively. Global is a mean.
#Gst_Nei = Gst in per locus or Gst_est in global.
#Hedrick's G"st: is Gprime_st both areas.
#Jost's D = D per loci and two global estimates. 'D_het' uses the averages of Hs and Ht across all loci while 'D_mean' takes the harmonic mean of all loci. It will be NA with negative values.
#Because all of these statistics are estimated from estimators of HS and HT, it’s possible to get negative values for each of these differentiation measures. Populations can’t be negatively differentiated, so you should think of these as estimates of a number close to zero (it’s up to you and your reviewers to decide if you report the negative numbers of just zeros).

#by populations. Number meanings are in 5C.
D = pairwise_D(input2, linearized = FALSE)
D = as.matrix(D)
write.csv(D, file = "D_groups_bre.csv")

GST_N = pairwise_Gst_Nei(input2, linearized = FALSE)
GST_N = as.matrix(GST_N)
write.csv(GST_N, file = "GST_N_groups_bre.csv")

GST_H = pairwise_Gst_Hedrick(input2, linearized = FALSE)
GST_H = as.matrix(GST_H)
write.csv(GST_H, file = "GST_H_groups_bre.csv")

#by populations and global to other indexes. See the help for poppr() for acronyms' meanings:
sum_stat = poppr(input2)
write.csv(sum_stat, file = "Sum_Stat_groups_bre.csv")

#8C. Calculate the bootstrap for FST family:
#create 100 matrices
bs = chao_bootstrap(input2, nreps = 100)

#per index:
D_boots = summarise_bootstrap(bs, D_Jost)
write.csv(D_boots$summary.global.het, file = "D_boots_bre_archs.csv")

GST_N_boots = summarise_bootstrap(bs, Gst_Nei)
write.csv(GST_N_boots$summary.global.het, file = "GST_N_boots_bre_archs.csv")

GST_H_boots = summarise_bootstrap(bs, Gst_Hedrick)
write.csv(GST_H_boots$summary.global.het, file = "GST_H_boots_bre_archs.csv")


##If you having problems with quatile function in mmod, try to uninstall the package and install the updated version from github with devtools package:
#remove.packages("mmod")
#devtools::install_github("dwinter/mmod", ref="bs_fix")


#########################################
################ AMOVA ##################
#########################################

#9A.create levels for AMOVA:

populations = c("GALVEZ", "GALVEZ", "GALVEZ", "MADRE", "MADRE", "BOLIVIA", "BOLIVIA", "BOLIVIA", "AMAZONAS", "AMAZONAS", "ACRE", "ACRE", "ACRE","IQUITOS", "IQUITOS", "IQUITOS")

endemism = c("Inambari", "Inambari", "Inambari", "South", "South", "South", "South", "South", "Napo", "Napo", "Inambari", "Inambari", "Inambari", "Napo", "Napo", "Napo")

ecorregion = c("Southwest","Southwest","Southwest","Southwest","Southwest", "Bolivian", "Bolivian","Bolivian", "Ucayali", "Ucayali", "Varzea","Varzea","Varzea","Varzea","Varzea","Varzea")

archs = c("IQUSER", "IQUSER", "IQUSER", "SOUTH", "SOUTH", "SOUTH", "SOUTH", "SOUTH","MARAND","MARAND", "IQUSER", "IQUSER","IQUSER", "IQUSER","IQUSER", "IQUSER")


#9B. Insert in the levels as strata
strata(input2) = data.frame(populations, endemism, ecorregion, archs)

#9C. Calculate the AMOVA

#####################
##### ENDEMISM #####
#####################
amova_noMiss = poppr.amova(input2, hier = ~ endemism/populations, squared = TRUE, within = FALSE ,quiet = FALSE, nperm=1000)

#Phi in this order:
#Phi-samples-total = (Phi_ST)
#Phi-samples-interfluves = (Phi_SC)
#Phi-interfluves-total = (Phi_CT)

#9D. Save results
DF2 = data.frame(amova_noMiss$statphi$Phi)
DF2= rbind(DF2, c("Total"))
phi_name = c("Phi_ST", "Phi_SC", "Phi_CT", "Total")

res_amova = data.frame(amova_noMiss$results, amova_noMiss$componentsofcovariance$Sigma,amova_noMiss$componentsofcovariance$`%`, DF2, phi_name)
write.csv(res_amova, file = "AMOVA_res_bre_endemism.csv", row.names = T)

#9E. Perform permutation process to test for significance
set.seed(1999)
signif = randtest(amova_noMiss, nrepet = 10000)
plot(signif)

write.csv(data.frame(signif$names, signif$obs, signif$expvar, signif$alter,signif$pvalue), file= "AMOVA_res_perm_bre_endemism.csv", row.names = T)

ggsave("AMOVA_res_perm_bre_endemism.svg", plot = plot(signif), device = "svg", path = NULL, scale = 1, width = 10, height = 30, units = "cm", dpi = 300, limitsize = TRUE)

######################
##### ECORREGION #####
######################
amova_noMiss = poppr.amova(input2, hier = ~ ecorregion/populations, squared = TRUE, within = FALSE ,quiet = FALSE, nperm=1000)

#Phi in this order:
#Phi-samples-total = (Phi_ST)
#Phi-samples-interfluves = (Phi_SC)
#Phi-interfluves-total = (Phi_CT)

#9D. Save results
DF2 = data.frame(amova_noMiss$statphi$Phi)
DF2= rbind(DF2, c("Total"))
phi_name = c("Phi_ST", "Phi_SC", "Phi_CT", "Total")

res_amova = data.frame(amova_noMiss$results, amova_noMiss$componentsofcovariance$Sigma,amova_noMiss$componentsofcovariance$`%`, DF2, phi_name)
write.csv(res_amova, file = "AMOVA_res_bre_ecorregion.csv", row.names = T)

#9E. Perform permutation process to test for significance
set.seed(1999)
signif = randtest(amova_noMiss, nrepet = 10000)
plot(signif)

write.csv(data.frame(signif$names, signif$obs, signif$expvar, signif$alter,signif$pvalue), file= "AMOVA_res_perm_bre_ecorregion.csv", row.names = T)

ggsave("AMOVA_res_perm_bre_ecorregion.svg", plot = plot(signif), device = "svg", path = NULL, scale = 1, width = 10, height = 30, units = "cm", dpi = 300, limitsize = TRUE)


######################
######## ARCHS #######
######################
amova_noMiss = poppr.amova(input2, hier = ~ archs/populations, squared = TRUE, within = FALSE ,quiet = FALSE, nperm=1000)

#Phi in this order:
#Phi-samples-total = (Phi_ST)
#Phi-samples-interfluves = (Phi_SC)
#Phi-interfluves-total = (Phi_CT)

#9D. Save results
DF2 = data.frame(amova_noMiss$statphi$Phi)
DF2= rbind(DF2, c("Total"))
phi_name = c("Phi_ST", "Phi_SC", "Phi_CT", "Total")

res_amova = data.frame(amova_noMiss$results, amova_noMiss$componentsofcovariance$Sigma,amova_noMiss$componentsofcovariance$`%`, DF2, phi_name)
write.csv(res_amova, file = "AMOVA_res_bre_archs.csv", row.names = T)

#9E. Perform permutation process to test for significance
set.seed(1999)
signif = randtest(amova_noMiss, nrepet = 10000)
plot(signif)

write.csv(data.frame(signif$names, signif$obs, signif$expvar, signif$alter,signif$pvalue), file= "AMOVA_res_perm_bre_archs.csv", row.names = T)

ggsave("AMOVA_res_perm_bre_archs.svg", plot = plot(signif), device = "svg", path = NULL, scale = 1, width = 10, height = 30, units = "cm", dpi = 300, limitsize = TRUE)

#########################################
################ PCA ####################
#########################################

##10. SCALE THE DATA MATRIX
#PCA can be done based on the covariance matrix as well as the correlation matrix. Scaling the data matrix means that all variables have zero mean and unit variance (also known as "normalizing", "studentisizing", "z-scoring"), makes the two approaches identical. This is because the covariance between two normalized variables *is* the correlation coefficient.
#For individual-based PCA we scale the genotype values across individuals at each locus to have unit variance. For the analysis all variables have mean = 0. 
#center = TRUE means the mean will be 0
#scale = TRUE means the unit of variance will be scaled from 0
#NA.method = "mean", the NA's will be replace by the mean allele frequencies;
#Warning message: "In .local(x, ...) : Some scaling values are null. Corresponding alleles are removed." This is the expected behaviour when loci are non-polymorphic. Why would you like to keep them? Especially in the context of a PCA, where it will creates variables with a variances of zero, and therefore Inf values when scaling?
input_scaled = scaleGen (input2, center = TRUE, scale = TRUE, NA.method = "mean")


##11. PCA
#center = TRUE means centring by the mean, in this case is 0
#scannf	= FALSE, indicates the screeplot should not be displayed
#nf = 3, if scannf = FALSE, an integer must be indicated. It will be the number of kept axes in the PCA
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE, nf = 3)


##12. SAVE THE RESULTS AS .CSV
#coordenates in the PCA by axis
write.csv(pca_input$li, file = "PCA_xis_coord_bre_archs.csv")
#% of PC variation
perc_pca = as.data.frame(pca_input$eig)
soma = sum (perc_pca)
perc_pca[,2] = (perc_pca/soma)*100
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
write.csv(perc_pca, file = "PCA_contribution_pc_eig_bre_archs.csv")

#13. CALCULATE THE % CONTRIBUTION OF EACH PC
PC1 = paste0("PC1 (",round(pca_input$eig[1]/sum(pca_input$eig)*100,2),"%)")
PC2 = paste0("PC2 (",round(pca_input$eig[2]/sum(pca_input$eig)*100,2),"%)")
PC3 = paste0("PC3 (",round(pca_input$eig[3]/sum(pca_input$eig)*100,2),"%)")


##14. SAVE/VIEW THE RESULTS AS .PDF TO EDIT AND PUBLISH
#create a dataframe with the results + grouping factor
df_out = as.data.frame(pca_input$li)
df_out$POP = populations #factor to group
head(df_out)

p = ggplot (df_out,aes(x=Axis1,y=Axis2, shape = POP)) +geom_point(size=rel(5)) +theme_pca +ggtitle("P. brevicauda") +guides(shape=guide_legend("Populations")) +labs(y=PC2, x = PC1, colour = "POP") +scale_shape_manual(labels = c("ACRE", "AMAZONAS", "BOLIVIA", "GALVEZ", "IQUITOS", "MADRE"), values = c(0,1,2,3, 4, 5)) #labels must be in alphabetical order

#verify the graph
plot(p)

#save as svg
ggsave("PCA_bre.svg", plot = plot(p), device = "svg", path = NULL, scale = 1, width = 20, height = 15, units = "cm", dpi = 300, limitsize = TRUE)


#########################################
############# PROCRUSTES ################
#########################################

##15.IMPORT THE LON-LAT MATRIX FOR THE SAMPLES
coord_mat = read.table('coord.txt', header = F, sep = "\t")
head(coord_mat)
tail(coord_mat)
str(coord_mat)
coord_mat2= coord_mat
coord_mat = coord_mat[ ,c(14:15)]

#Should be as matrix without row and cols names 
coord_mat = as.matrix(coord_mat)

##16. CALCULATE THE CORRELATION AND P-VALUE
test_pro_1 = protest (X=coord_mat, Y=pca_input, scale=T, symmetric=T, permutations=10000)
#p-value is "Significance"
#correlation after permutation is "Correlation in a symmetric Procrustes rotation"

#save the result in a csv file:
proc=cbind.data.frame (test_pro_1$ss, test_pro_1$scale, test_pro_1$signif, test_pro_1$permutations)
colnames(proc) = c("Procrustes Sum of Squares (m12 squared)", "Correlation in a symmetric Procrustes rotation", "p-value", "Number of Permutation")
rownames(proc) = c("bre")
write.table(proc, "PROCRUSTES_RES_bre.csv",sep = ",", quote = F, row.names = F)

#use this to get coordinates for the PCA value
test_pro_2 = procrustes (X=coord_mat, Y=pca_input, scale=T, symmetric=F)

#view the arrows
plot (test_pro_2, main=NULL, ar.col='black', xlab='', ylab='')

##17. SAVE/VIEW THE RESULTS AS .PDF TO EDIT AND PUBLISH
##USE A MAP OF YOUR STUDY AREA TO THIS STEP
##CHANGE THE FOLDER IF IT IS NECESSARY

#import raster map (asc, tiff, ESRI)
#map = raster("cinza.asc") #single band raster as asc
map2 = stack("NoSea_CL.tif") #multi band raster as satelite images 
#map3 = raster("relevo.asc") 
#map4 = raster("alt")

#import shapefiles of countries and states
fr.states = shapefile("Brazil_Estados.shp")
fr.countries = shapefile("ne_10m_admin_0_countries.shp")
fr.rivers = shapefile("ne_10m_rivers_lake_centerlines.shp")
fr.archs = shapefile("archs.shp")
fr.stable = shapefile("StableBRE.shp")

#denide an extension for the map
#use x as longitude; y as latitude; same length for a square shape
#redefine the maps and shapefiles sizes according to the study area
#define lat and lon limits, use square shapes to avoid blank spaces like 20 degrees
ext = extent(-80, -60, -20, 0)
axis.x = c(-80, -60)
axis.y = c(-20, 0)

#test the redefinition
#plot(crop(map, ext)) #single band raster as asc
plotRGB(crop(map2, ext)) #multi band raster as satelite images 

#redefine the maps and shapefiles sizes according to the study area
fr.states2 = crop(fr.states, ext)
fr.countries2 = crop(fr.countries, ext)
fr.rivers2 = crop(fr.rivers,ext)
fr.archs2 = crop(fr.archs,ext)
fr.stable2 = crop(fr.stable, ext)

#CREATE THE PDF FILE
pdf("./PROCRUSTES_MAP_bre_stables.pdf", onefile=FALSE)

#add maps in color
plot.new()
par(oma=c(5,5,2,2)) # Add margins
plot.window(xlim= axis.x, ylim= axis.y)
plotRGB(map2, r = 1, g = 2, b = 3, ext= ext, alpha =120, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
axis(1, at=seq(-80,-60, by=5),  cex.axis=1)
axis(2, at=seq(-20, 0, by=5),  cex.axis=1, las=1)
mtext(side=1, text='Longitude', line=2.2, cex=1.15)
mtext(side=2, text='Latitude', line=2.5 , cex=1.15)

#add states and countries maps
plot(fr.rivers2, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lty="solid", lwd=0.05, col="deepskyblue3")
plot(fr.countries2, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=0.3)
plot(fr.states2, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=0.3)
plot(fr.archs2, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lty="dotted", lwd=3)
plot(fr.stable, add=TRUE, axes = F, bg = FALSE, col = rgb(0,0,0,0.3), bty ="n", xaxt='n', ann=FALSE, yaxt='n', lty="blank")


#add maps in grey scale
#plot(map, col=grey(0:8/8), axes = F, bg = FALSE, bty ="n", frame.plot=FALSE, xaxt='n', ann=FALSE, yaxt='n', legend=F, ext= extent(-70, -45, -20, 0)) #option1
#plot(map, col=gray.colors(1000, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL), axes = F, bg = FALSE, bty ="n", frame.plot=FALSE, xaxt='n', ann=FALSE, yaxt='n', legend=F, ext= extent(-70, -45, -20, 0)) #option2
#axis(1, at=seq(-70,-45, by=5),  cex.axis=1); #axes
#axis(2, at=seq(-20, 0, by=5),  cex.axis=1, las=1); #axes

#arrows
arrows(test_pro_2$Y[,1]+test_pro_2$translation[1,1], test_pro_2$Y[,2]+test_pro_2$translation[1,2], test_pro_2$X[,1]+test_pro_2$translation[1,1], test_pro_2$X[,2]+test_pro_2$translation[1,2], length=0.05, code = 0, cex=0.1, lwd=2) #code is the style of arrow, 0 = without head

#points by populations
points(test_pro_2$Y[1:3,1]+test_pro_2$translation[1,1],test_pro_2$Y[1:3,2]+test_pro_2$translation[1,2], col = 'black', bg = 'darkslategrey', cex = 2.5, pch=21)
points(test_pro_2$Y[4:5]+test_pro_2$translation[1,1],test_pro_2$Y[4:5,2]+test_pro_2$translation[1,2], col = 'black', bg ='red', cex = 2.5, pch=21)
points(test_pro_2$Y[6:8]+test_pro_2$translation[1,1],test_pro_2$Y[6:8,2]+test_pro_2$translation[1,2], col = 'black', bg ='blue', cex = 2.5, pch=21)
points(test_pro_2$Y[9:10]+test_pro_2$translation[1,1],test_pro_2$Y[9:10,2]+test_pro_2$translation[1,2], col = 'black', bg ='yellow', cex = 2.5, pch=21)
points(test_pro_2$Y[11:13]+test_pro_2$translation[1,1],test_pro_2$Y[11:13,2]+test_pro_2$translation[1,2], col = 'black', bg ='pink', cex = 2.5, pch=21)
points(test_pro_2$Y[14:16]+test_pro_2$translation[1,1],test_pro_2$Y[14:16,2]+test_pro_2$translation[1,2], col = 'black', bg ='brown', cex = 2.5, pch=21)

#points by locality, one color
#points(test_pro_2$X[,1]+test_pro_2$translation[1,1],test_pro_2$X[,2]+test_pro_2$translation[1,2], col = 'black', cex = 1, pch=17)

#point by locality, multiple color
points(test_pro_2$X[1:3,1]+test_pro_2$translation[1,1],test_pro_2$X[1:3,2]+test_pro_2$translation[1,2], col = 'black', bg = 'darkslategrey', cex = 2.5, pch=24)
points(test_pro_2$X[4:5,1]+test_pro_2$translation[1,1],test_pro_2$X[4:5,2]+test_pro_2$translation[1,2], col = 'black', bg ='red', cex = 2.5, pch=24)
points(test_pro_2$X[6:8,1]+test_pro_2$translation[1,1],test_pro_2$X[6:8,2]+test_pro_2$translation[1,2], col = 'black', bg ='blue', cex = 2.5, pch=24)
points(test_pro_2$X[9:10,1]+test_pro_2$translation[1,1],test_pro_2$X[9:10,2]+test_pro_2$translation[1,2], col = 'black', bg ='yellow', cex = 2.5, pch=24)
points(test_pro_2$X[11:13,1]+test_pro_2$translation[1,1],test_pro_2$X[11:13,2]+test_pro_2$translation[1,2], col = 'black', bg ='pink', cex = 2.5, pch=24)
points(test_pro_2$X[14:16,1]+test_pro_2$translation[1,1],test_pro_2$X[14:16,2]+test_pro_2$translation[1,2], col = 'black', bg ='brown', cex = 2.5, pch=24)

#create a legend box for your graph
legend = c("GALVEZ", "MADRE", "BOLIVIA", "AMAZONAS", "ACRE", "IQUITOS")
col = c('darkslategrey', 'red', 'blue', "yellow", "pink", "brown")
legend(x=-79.5, y=-15, legend = legend, cex = 0.8, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col, pch= 22, pt.cex=1.5)

#scale bar and north arrows 
#scalebar(500, xy=c(-46.5,-16.5), type="bar", divs=4, below = "Km", lonlat = T, lwd=2, adj=c(0,-1.25)) #adj first number is lateral aligment and second is the horizontal aligment
scalebar(250, xy=c(-67,-19.5), type="line", lonlat = T, lwd=3, label = "250 Km", font=2)
#north.arrow(xb=-51.5, yb=-16.25, len=0.2, lab="N", col='black', font.sub=2)

dev.off()


############################################
######## PROCRUSTES - PERMUTATION ##########
############################################

##18. DEFINE SOME VARIABLES 
#number of pops 
pop_counter = 6

#vector with number of individuals per population (in the same order of .str file)
input2$pop #to verify
samp_vec = c(3,2,3,2,3,3)

#vector with pop identifiers. See 5C to names and order:
pop_vec = c("GALVEZ", "MADRE", "BOLIVIA", "AMAZONAS", "ACRE", "IQUITOS")

##19. ANALYSES REMOVING ONE POPULATION BY LOOP
####this block does a Procrustes analysis on the total dataset first, then deletes
#one pop at a time and reruns the procrustes (starting at pop #1). the stats are
#written to a file (testeBC_statsfile.txt), which must already exist in the 
#folder
pro = protest(X=coord_mat,Y=pca_input$li,scale=T,symmetric=T,permutations=10000)

#asin_d(pro$rotation[2,1])
#Since the R default is to compute trigonometric functions on angular 
#measurements stored in radians,function asin_d performs the conversion 
#from degrees, reducing the need to do so a priori, outside the function.

# define other variables
total = pro$scale
rotationallpops = asin_d(pro$rotation[2,1])
t_prime = NULL
t_prime2 = NULL
rotation = NULL
rotation2 = NULL
signif = NULL
difference = NULL
max_diff = 0
count = 1

for(m in 1:pop_counter) {
  ## delete from coordinate matrix those related to each of the eliminated pop
  temp_mat = coord_mat[-(count:(count+samp_vec[m]-1)),]
  ## delete from pca$li matrix the results related to each of the eliminated pop
  temp_pca = pca_input$li[-(count:(count+samp_vec[m]-1)),]
  ## delete from the matrix of mean scaled genotypes those related to each of the eliminated pop
  red = input_scaled[-(count:(count+samp_vec[m]-1)),]
  ## does pca with reduced matrix and protest with the new dataset (without each eliminated pop)
  pca1 = dudi.pca(red, center = T, scale = T, scannf=F,nf = 2)
  pro = protest(X=temp_mat,Y=pca1$li,scale=T,symmetric=T,permutations=10000)
  t_prime2 = c(t_prime2,pro$scale)
  rotation2 = c(rotation2, asin_d(pro$rotation[2,1]))
  signif = c(signif, pro$signif)
  difference = c(difference,(pro$scale-total) )
  
  pro = protest(X=temp_pca,Y=pca1$li,scale=T,symmetric=T,permutations=10000)
  t_prime = c(t_prime,pro$scale)
  rotation = c(rotation, asin_d(pro$rotation[2,1]))
  
  
  count = count + samp_vec[m]
}

##20. SAVE RESULTS AS A FIGURE. 
#Each color represent one population in the same order of .str file
#order difference values to use as ylim
sort(difference)
#makes barplot
setEPS()
postscript("PROCRUSTES_PERM_bre.eps")
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(2,2,1,0));
plot.window(xlim=c(0, 8), ylim=c(-0.3, 0.3));
col=c('darkslategrey', 'red', 'blue', "yellow", "pink", "brown")
barplot (difference, ylim = c(-0.3, 0.3)  ,main = "Relative difference in t0", col=col)
legend (x="topright", ncol=2, pch=21, bty='n', cex= 0.8, pt.cex=2,legend = pop_vec, col = "black", pt.bg = col)
#graph with names and without colors
#barplot (difference, names.arg = pop_vec, cex.names = 0.5)
abline(h=0)
dev.off()

#Makes data frame with all relevant information and saves it in a .csv file
test_data <- data.frame(t0 = total, Angle = rotationallpops,POP = pop_vec, GEO = t_prime2, DIFF = difference, ROT_GEO = rotation2, P = signif, PCA = t_prime, ROT_PCA = rotation)
write.csv(test_data, file="PROCRUSTES_PERM_STAT_bre.csv")