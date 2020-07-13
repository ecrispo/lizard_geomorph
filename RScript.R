install.packages("devtools")
devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
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
#There is no partial or relative warps in geomorph
myGPA<-gpagen(myData, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
       max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
summary(myGPA)
myGPA$coords
myGPA$Csize
plotAllSpecimens(myGPA$coords,mean=F,links=fish.gpa$links)
#Try mean=T option to plot the consensus landmark positions
?plotAllSpecimens

#Data Analysis page 43 of Quick Guide to Geomorph
#Relative warps on page 228 of Schlager (Morpho package)
install.packages("Morpho")
library(Morpho)
