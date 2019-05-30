###############################################################################
######################## SPECIES DELIMITATION ANALYSES ########################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################
################################     iBPP    ##################################
########################### STEP 02: MOLECULAR DATA ###########################


###############################################################################
################ Based on Joyce R. Prado and J. P. Huang scripts ##############
###################### KNOWLES LAB, UNIVERSITY OF MICHIGAN ####################
###############################################################################

##1. FOR OTHER STEPS, PLEASE GO TO: https://github.com/jdalapicolla/

##2.INPUTS FOR THIS STEP:
#A. ".loci" FILE FROM IPYRAD STEP 7.
#B. ".csv" FILE WITH ALL SAMPLES PER LINEAGES/PUTATIVE SPECIES YOU NEED TO TEST. SEE THE EXAMPLES FILES.

##3. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())


##4. CHANGE THE SAMPLE NAMES FROM LOCI FILE TO A PATTERN PER LINEAGE/PUTATIVE SPECIES. EACH LINEAGE IN A DIFFERENT .csv FILE. THE NUMBER IN FRONT OF THE SAMPLE NAME WILL REFER TO A DIFFERENT LINEAGE.
#read the .loci file
dat = scan('proe_78.loci', what = 'character', sep = '\n')
#read all lineages in csv files in the folder
file.list = list.files( , pattern = '.csv')

#replace the name of sample to a new pattern. Like "AMNH272698" to ">11_AMNH272698". ">" to indicate the beginning of sample, "11_" to the lineage, clade1, lineage 1, and a "_" to separate the name of sample. Each csv file will be consider a new lineage
for(i in 1:length(file.list)){
  temp.tax = scan(file.list[i], what = 'character', sep ='\n')
  for(iter in 1:length(temp.tax)){
    s.pattern = temp.tax[iter]
    temp.p = paste('>',i+10, '_', sep = '')
    rep.pattern = paste(temp.p, s.pattern, sep = '')
    dat = gsub(s.pattern, rep.pattern, dat)
  }
}

s.pattern = '//'
rep.pattern = '//   '
dat = gsub(s.pattern, rep.pattern, dat)

#save the results
write.table(dat, file = 'HOP_changed.loci', sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)

##5. SPLIT THE SAMPLES BY LOCUS
#You need to load the "name_changed.loci" because the file need to be without quotes or row/col names. DO NOT USE the dat object from step 4 again:
data.loci = scan('HOP_changed.loci', what = 'character', sep = '\n');#used the edited .loci file
head(data.loci)

#create a vector with lineages codes. 
species.vec = c(">11_",">12_", ">13_")

#split the loci
break.lines = grep('//', data.loci)
index.lines = c(0, break.lines)
loci.count = 0 #number of loci counter. will print the number of loci at the end


##6. CREATE THE INPUT FILE FOR iBPP: 
#number of loci
length(break.lines) 


for(ii in 1:50442){#the number of loci in the .loci file. Can be obtained by typing length(break.lines).
  
  num.1 <- index.lines[ii] + 1;
  num.2 <- index.lines[ii + 1] - 1;
  
  temp.text <- data.loci[c(num.1:num.2)];
  
  num.sp <- 0;
  ind.vec <- NULL;
  for(iter in 1:length(species.vec)){
    temp.ind.vec <- grep(species.vec[iter], temp.text);
    ind.vec <- c(ind.vec, temp.ind.vec);
    num.sp <- num.sp + (length(temp.ind.vec) > 0);
    
    if(num.sp == length(species.vec)){
      loci.count <- loci.count + 1;
      
      subset.temp <- temp.text[ind.vec];
      
      string.ed <- subset.temp[1];
      find.p <- '>11_[[:alnum:]]+[[:space:]]+'#for calc length of the loci; data here are all 110#change the regex according to how you name your samples
      sub.pa <- '';
      temp.dna.seq <- sub(find.p,sub.pa,string.ed);
      dna.seq <- strsplit(temp.dna.seq, '');
      loci.length <- length(dna.seq[[1]])
      n.ind <- length(subset.temp);
      first.line <- paste(n.ind, loci.length, sep = ' ');#first line for bpp each locus
      
      new.string <- NULL;
      a11.vec <- grep(species.vec[1], subset.temp);
      count <- 1;
      for(iii in a11.vec){
        s.pattern <- '>11_[[:alnum:]]+';#regex for the 1st species
        bpp.tail <- paste('^', count, sep = '');
        rep.pattern <- paste('	HG', bpp.tail, sep = '');#how you want the sample to be named in bpp file
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a12.vec <- grep(species.vec[2], subset.temp);
      count <- 1;
      for(iii in a12.vec){
        s.pattern <- '>12_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+3, sep = '');#count + 1 because there are 1 19 samples
        rep.pattern <- paste('	PRO', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }  
      a13.vec <- grep(species.vec[3], subset.temp);
      count <- 1;
      for(iii in a13.vec){
        s.pattern <- '>13_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+76, sep = '');#count + 1 because there are 1 19 samples
        rep.pattern <- paste('	TEP', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }   
   
      
      temp.bpp.file <- c(first.line, new.string, '\n');
      write.table(temp.bpp.file, append = TRUE, file='HOP_ibpp.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);#change the file name here
      
    }
  }
  
}
cat('a total of', loci.count, 'loci');#let you know how many loci are there saved in the bpp file