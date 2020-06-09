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
#For some reason I got missing data (NA) for some landmarks.

#You now have a tps file called "newFile" that contains the data you just collected.
#See page 17 of Quick Guide to Geomorph:
readland.tps("newFile2.tps", specID = "ID", readcurves = FALSE, warnmsg = TRUE)
#For some reason I am getting warning message that there is no scale for some images (SCALE=)
#Deleted missing landmark data and changed number of landmarks, and no it is working okay
#Can estimate missing landmarks - see page 26 of Quick Guide to Geomorph

#Generalized Procrustes Analysis
#Starts on page 29 of Quick Guide to Geomorph 

