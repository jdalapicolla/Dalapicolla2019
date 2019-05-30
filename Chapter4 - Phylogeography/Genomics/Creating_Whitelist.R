###############################################################################
######################### PHYLOGEOGRAPHIC ANALYSES ############################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
#######################     STEP 01: DATA CLEANING     ########################
#############################  SEGREGATING SITES  #############################

###############################################################################
################# Based on Andrea Thomaz and Sarp Kaya scripts ################
###################### KNOWLES LAB, UNIVERSITY OF MICHIGAN ####################
###############################################################################

#November 17, 2015 - Andrea Thomaz for Stacks version 1.35
#Modified by Sarp Kaya - May 20th, 2017 for Stacks version 1.46 or above
#read .vcf file from stacks output (WITHOUT write_random_SNP flag) for:
#plot the frequency of variable sites per position along all loci
#calculate theta based on number of segregating sites and individuals to create blacklist to delete very variable loci.


##1. FOR OTHER STEPS, PLEASE GO TO: https://github.com/jdalapicolla/

##2.INPUTS FOR THIS STEP:
#A. ".vcf" file from stacks output (WITHOUT write_random_SNP flag)

##3. INSTALL AND LOAD THE PACKAGES#
#install the packages
install.packages(c("plyr", "pegas","tcltk"))
#load the packages
library("plyr")
library("pegas")
library("tcltk") #for Linux

#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("plyr", "pegas","tcltk")

##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())


##5. CHANGE THE SCIENTIFIC NUMBERS
# To close scientific numerics because it creates problem in whitelist, if you want to open use options(scipen = 0)
options(scipen = 999) 


##6. READ VCF
data = read.table('batch_1.vcf', header = FALSE, sep = "\t")
head(data[1:20,1:20])


##7. CONVERT THE VCF IN A DATA FRAME:
#sequence length
seq_len = 140 #MODIFY HERE according the length of the sequence. I sequenced 150pb, with 10 pb of barcodes = 140

#selecting loci ID and position from column $V
#this is for ID column problem in stacks v.1.46  
loci_num = as.numeric(sub('_.*', '', data$V3))
pos2 = as.numeric(sub('.*_', '', data$V3))

#creates dataframe with loci ID, the variable positions and the number of individuals in each loci
new_data = data.frame(loci_ID = loci_num,
                       pos_vcf1 = data[,2],
                       pos = pos2, 
                       #pos_dea = (data[,2] - seq_len*(loci_num-1))-2,
                       ind = rowSums(data[,10:length(data)] != "./.:0:.,."))

head(new_data)
min(new_data$pos)#should always be position 5 (first positions are the adapters)
max(new_data$pos)#should always be position 139
length(unique(new_data$loci_ID))#how many loci do I have
length(new_data$loci_ID)


##8. GRAPH OF VARIABLE SITES ALONG THE LOCI
par(mar = rep(2, 4))
#saving graph with frequency of variable sites along the loci
pdf("./SNPdistr_BRE_ARC_pos140bp.pdf")
hist(new_data$pos, xlim = c(-1,seq_len), breaks = c(seq(-1, seq_len-1, by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites');
abline(2200, 0, col = "red")#helps to find where starts to increase toward the end, last positions have strong increase
abline(v = 127, col = "red")#helps to figure out where to cut off before increase in bad calls
#move the lines around to visualize depending on my case
dev.off()


##9. BASE ON THE GRAPH, CHOOSE HOW MANY POSITION TO DELETE FROM THE END
to_del = 13 #how many sites to delete in the end of the sequence
#13 is based on the 127 I chose for the ab line above
seq_len_cut = seq_len - to_del


##10. CREATE A WHITELIST TO EXCLUDE (to_del) POSITIONS
whitelist = subset(new_data, pos < seq_len_cut)[,c(1,3,4)]
pdf("./SNPdistr_pos_cutto127bp.pdf")
hist(whitelist$pos, xlim = c(0,seq_len_cut), breaks = c(seq(-1, seq_len_cut -1 , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites');
dev.off()


##11. CALCULATING THETA FOR ALL LOCI
var.sites = count(whitelist, "loci_ID")
length(var.sites$loci_ID)
max(var.sites$freq)
theta_calc = merge(unique(whitelist[,-2]), var.sites, by = "loci_ID")
theta_calc$theta = 0
head(theta_calc)
for (i in 1:length(theta_calc$theta)){
  theta_calc[i,4] = theta.s(theta_calc$freq[i], theta_calc$ind[i])/seq_len_cut
}


##12. CALCULATING 95% QUANTILE TO EXCLUDE LOCI EXTREMELY VARAIBLE
quant = quantile(theta_calc$theta, probs=0.95) #set the value to be 
quant
pdf("./theta127bp.pdf")
hist(theta_calc$theta)
abline(v = quant, col="red")
dev.off()


##13. NUMBER OF MUTATIONS IN A LOCUS
#what is the maximum number of mutations in a locus
max(theta_calc$freq) #max theta before 
x = subset(theta_calc, theta < quant)
max(x$freq) #max theta after, make sure is realistic for a 140 bp sequence
#think about what mutation rate the spp might have


##14. CREATING A BLACKLIST OF LOCI TO REMOVE
blacklist = subset(theta_calc, theta > quant)[,1]
#write.table(blacklist, file="blacklist.txt", sep = '\n', row.names = F, col.names = F)


##15. CREATING A WHITELIST FOR RE-RUN POPULATIONS STEP IN STACKS
#removes the blacklist from the whitelist and write off white list
whitelist$blacklist = match(whitelist$loci_ID, blacklist, nomatch = 0)
whitelist_final = subset(whitelist, blacklist == 0)
length(unique(whitelist_final$loci_ID)) #number of unique loci, this is the number I need to get out with "write random loci"
length(whitelist_final$loci_ID) #number of snps
write.table(whitelist_final[,1:2], file="whitelist_BRE.txt", sep = '\t', row.names = F, col.names = F)