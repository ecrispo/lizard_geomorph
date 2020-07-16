install.packages("devtools")
devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
install.packages("Morpho")
library(geomorph)
#Make sure to set the working directory - where specimen photos are saved

filelist <- c("imageLake5.jpg", "imageLake48.jpg", "imagePond5.jpg", "imagePond6.jpg")

#Works better to digitize in R than in RStudio
digitize2d(filelist,
  nlandmarks=15,
  scale = 10,
  tpsfile = "newFile.tps" ,
  MultScale = FALSE,
  verbose = TRUE
)

#You now have a tps file called "newFile" that contains the data you just collected.
#See page 17 of Quick Guide to Geomorph:
myData<-readland.tps("newFile.tps", specID = "ID", readcurves = FALSE, warnmsg = TRUE)

#Generalized Procrustes Analysis
#Starts on page 29 of Quick Guide to Geomorph 
#Starts on page 63 of Zelditch
myGPA<-gpagen(myData, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
       max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
summary(myGPA)
myGPA$coords
myGPA$Csize
plotAllSpecimens(myGPA$coords,mean=T)
#Try mean=T option to plot the consensus landmark positions
?plotAllSpecimens

#Data Analysis page 43 of Quick Guide to Geomorph
myGPA$site=c("lake","lake","pond","pond")
is.factor(myGPA$site)
myGPA$site<-as.factor(myGPA$site)
myGPA$site

gdf1 <- geomorph.data.frame(coords = myGPA$coords, site = myGPA$site,
                           logcs = log(myGPA$Csize))
gdf1
procD.lm(coords ~ site + logcs, data = gdf1, iter = 999,
         RRPP = FALSE, print.progress = T)
res <- procD.lm(coords ~ logcs,data=gdf1)
res

#EXAMPLE
data("plethodon")
plethodon$site
is.factor(plethodon$site)
data(plethodon) # example dataset
Y.gpa <- gpagen(plethodon$land, print.progress = FALSE) # GPA-alignment
pleth.gdf <- geomorph.data.frame(shape = Y.gpa$coords,
                           site = plethodon$site, species = plethodon$species) # make geomorph data frame
# permutation option 1: randomize raw values
procD.lm(shape ~ species * site, data = pleth.gdf, iter = 999,
         RRPP = FALSE, print.progress = FALSE)



