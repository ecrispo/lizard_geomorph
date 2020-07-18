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
myGPA<-gpagen(myData, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
       max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
summary(myGPA)
myGPA$coords
myGPA$Csize
plotAllSpecimens(myGPA$coords,mean=T)
#Try mean=T option to plot the consensus landmark positions
?plotAllSpecimens



#Data Analysis page 43 of Quick Guide to Geomorph
#First create a new variable that includes the site
myGPA$site=c("lake","lake","pond","pond")
is.factor(myGPA$site)
myGPA$site<-as.factor(myGPA$site)
myGPA$site

#Create a data frame
gdf1 <- geomorph.data.frame(coords = myGPA$coords, site = myGPA$site,
                           logcs = log(myGPA$Csize))
gdf1

#Perform a linear model, MANCOVA
#Response variable is the landmark coordinates 
#Explanatory variable is the site
#Covariate is the log of the centroid size
model1<-procD.lm(coords ~ site + logcs, data = gdf1, iter = 999)
anova(model1)

#View coordinates in a different format
y<-two.d.array(myGPA$coords)
y
model2<-procD.lm(y~gdf1$site+gdf1$logcs,iter=999)
anova(model2)

#Use model 3 as it is the simplest:
model3<-procD.lm(myGPA$coords~myGPA$site+log(myGPA$Csize),iter=999)
anova(model3)

#Examine effects of centroid size alone (ALLOMETRY).
res <- procD.lm(coords ~ logcs,data=gdf1)
summary(res)

#Discuss type I/II/II SS and what permutation of residuals means
#Discuss why log the centroid size (more scatter as size increases)
hist(myGPA$Csize,xlab="Centroid Size",ylab="Number of fish",main="")
hist(log(myGPA$Csize),xlab="Log transformed centroid size",ylab="Number of fish",main="")
#Discuss coordinates versus partial warps (see chapter 5 of Zelditch)


plotTangentSpace(myGPA$coords,groups=myGPA$site,legend=T,label=T)


