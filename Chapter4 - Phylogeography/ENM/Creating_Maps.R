###############################################################################
############################# CREATE MAPS IN R ################################
################################## TUTORIAL ###################################
###############################################################################

######################### JERONYMO DALAPICOLLA, 2019 ##########################

###############################################################################
######################### JERONYMO DALAPICOLLA, 2019 ##########################
###################### LABORATORIO DE MAMIFEROS - LAMA ########################
####################### ESALQ - UNIVERSIDADE DE SAO PAULO #####################
###############################################################################

##1.FOR MORE SCRIPTS, PLEASE GO TO: https://github.com/jdalapicolla/


##2.INPUTS FOR THIS STEP:
#A. MAPS (RASTER OR SHAPEFILES).
#B. LONG-LAT COORDENATES FOR ALL YOUR SAMPLES.


##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("raster", "rgdal", "rgeos", "GISTools", "pacman"))
#load the packages
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("raster", "tcltk", "rgdal", "rgeos", "GISTools")

##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILES MUST BE THERE! 
#Windows
choose.dir()
#Linux
setwd(tk_choose.dir())
#Any computer
setwd('C:/Users/Johnny/Downloads/bio_10m_bil') # change the path


##5. IMPORT THE MAPS (asc, tiff, ESRI)
#map = raster("cinza.asc") #single band raster as asc
map2 = stack("NoSea_CL.tif") #multi band raster as satelite images 
#map3 = raster("relevo.asc")
#map4 = raster("alt")

#import shapefiles of countries and states
fr.states = shapefile("Brazil_Estados.shp")
fr.countries = shapefile("ne_10m_admin_0_countries.shp")


##6. DEFINE AN EXTENSION FOR THE MAP
#use x as longitude; y as latitude; same length for a square shape
#redefine the maps and shapefiles sizes according to the study area
#define lat and lon limits, use square shapes to avoid blank spaces like 15 degrees of latitude and 15 degrees of longitude
ext = extent(-75, -60, -15, 0) #define extension
axis.x = c(-75, -60) #define axis
axis.y = c(-15, 0) #define axis

#Verify the area, modify long and lat above if necessary
#plot(crop(map, ext)) #single band raster as asc
plotRGB(crop(map2, ext)) #multi band raster as satelite images 

#redefine and shapefiles sizes according to the study area
fr.states2 = crop(fr.states, ext)
fr.countries2 = crop(fr.countries, ext)


##7. SAVE THE RESULTS IN A PDF FILE. IT IS EASIER TO EDIT IN A CORELDRAWN, ILLUSTRATOR OR INKSCAPE IF IT IS NECESSARY
pdf("./MAP_PROECHIMYS_CERRADO.pdf", onefile = FALSE) #change the name of pdf #onefile is to avoid a blank page be printed before your map. Sometimes the pdf file has 2 pages, 1 blank and other with the map. #For testing you can avoid this command, and after you sure about the map you can run the commands again with this pdf command line.  


##8. ADD MAPS IN COLOR
plot.new() #open a new area to draw
par(oma=c(5,5,2,2)) # Add margins
plot.window(xlim= axis.x, ylim= axis.y) #add the axis to the maps
plotRGB(map2, r = 1, g = 2, b = 3, ext= ext, alpha =120, stretch="lin", bty="n", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness; # yaxs is the position of the grid, "r" is 1 map unit deslocated from the map; # r, g, b indicate the number of band is read, green and blue.
axis(3, at=seq(-75,-60, by=5),  cex.axis=1) #add grids
axis(2, at=seq(-15, 0, by=5),  cex.axis=1, las=1) #las is the orientation of the text/numbers
#mtext(side=1, text='Longitude', line=2.2, cex=1.15) #is the label, not necessary to regular maps, only procrustes.
#mtext(side=2, text='Latitude', line=2.5 , cex=1.15)

#add states and countries maps
plot(fr.countries2, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n')
plot(fr.states2, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n')


##9. ADD MAPS IN GRAY SCALE
plot(map, col=grey(0:8/8), axes = F, bg = FALSE, bty ="n", frame.plot=FALSE, xaxt='n', ann=FALSE, yaxt='n', legend=F, ext= ext) #option1

plot(map, col=gray.colors(1000, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL), axes = F, bg = FALSE, bty ="n", frame.plot=FALSE, xaxt='n', ann=FALSE, yaxt='n', legend=F, ext= ext) #option2; 1000 is the number of shades of gray, and 0.3 and 0.9 the range the where the number of colors will be create. 0 is white, 1 is black.
axis(3, at=seq(-70,-45, by=5),  cex.axis=1); #grids
axis(2, at=seq(-20, 0, by=5),  cex.axis=1, las=1);


##10. IMPORT THE LON-LAT MATRIX FOR THE SAMPLES
coord_mat = read.table('coord.txt', header = F, sep = "\t") #my table with 19 colunms, separeted by tab and WITHOUT HEADER
#verify the table
head(coord_mat)
tail(coord_mat)
str(coord_mat)

coord_mat2= coord_mat #save a copy if you need the raw table in the future

coord_mat = coord_mat[coord_mat[,9] =="P. longicaudatus", ] #select the lines of the target species, in the 9th column.
coord_mat = coord_mat[ ,c(18:19)] #select the long and lat colunms, in my case, long is 18th and lat the 19th.
coord_mat = as.matrix(coord_mat) #transform in a matrix.


##11. ADD POINTS TO THE MAP
#points by localities, each locality is a different color.
points(coord_mat[1:2,1], coord_mat[1:2,2], col = 'black', bg = 'darkslategrey', cex = 1.5, pch=21) #long is x the first
points(coord_mat[3:9,1], coord_mat[3:9,2], col = 'black', bg ='pink', cex = 1.5, pch=24)
#add how many points you want to, changing the color (bg), size (cex), and symbols (pch)

#points by locality, same color to all
points(coord_mat[ ,1], coord_mat[ ,2], col = 'black', cex = 1, pch=17)


##12. ADD A LEGEND BOX FOR YOUR MAP
legend = c("Headwaters Juruá River", "Lower Central Juruá  River", "Madre de Dios and Beni rivers", "Mouth Juruá River", "Upper Central Juruá River")
col = c('darkslategrey', 'pink', 'blue', 'red', 'yellow') #same order than the names
legend(x=-66.3, y=-12, legend = legend, cex = 0.8, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col, pch= c(22, 21, 24, 17, 23), pt.cex=1.5) #pt.cex is tha size of the symbols; xy is the position of the left corner; pch is in the same order than names.


#13. ADD THE SCALE BAR AND NORTH ARROW 
#scalebar(500, xy=c(-46.5,-16.5), type="bar", divs=4, below = "Km", lonlat = T, lwd=2, adj=c(0,-1.25)) #adj first number is lateral aligment and second is the horizontal aligment
scalebar(250, xy=c(-73,-14.5), type="line", lonlat = T, lwd=3, label = "250 Km", font=2)
north.arrow(xb=-51.5, yb=-16.25, len=0.2, lab="N", col='black', font.sub=2) # I prefer add a north arrow by myself on Inkscape

#I don't know how to put a zoom retangle in the map, so I prepare this map a part and add the zoom retangle manualy using Inkscape.


dev.off() #windows will be closed or pdf will be created.
