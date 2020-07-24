#Install the needed packages
install.packages("devtools")
install.packages("Rtools")
install.packages("RRPP")
devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
#The next two packages are optional
install.packages("lattice")
install.packages("Morpho")
#Activate the geomorph package and check which version you're using
#Note that the Undergraduate Guide was constructed using geomorph version 3.2.1
library(geomorph)
packageVersion('geomorph')
#View the vignettes, i.e. worked out geomorph examples
?vignettes
vignette("geomorph.PCA")
vignette("geomorph.assistance")
vignette("geomorph.functions")
vignette("geomorph.digitize3D") 
#Make sure to set the working directory - where specimen photos are saved
#Create a list of your specimen images
filelist <- c("imageLake5.jpg", "imageLake48.jpg", "imagePond5.jpg", "imagePond6.jpg")
#Place landmarks on your specimen images (this is called digitizing)
#Works better to digitize in R than in RStudio
digitize2d(filelist,
  nlandmarks=15,
  scale = 10,
  tpsfile = "newFile.tps" ,
  MultScale = FALSE,
  verbose = TRUE
)
#You now have a tps file called "newFile" that contains the data you just collected.
#Read the data in R or RStudio
myData<-readland.tps("newFile.tps", specID = "ID", readcurves = FALSE, warnmsg = TRUE)
#Generalized Procrustes Analysis is next
#Starts on page 63 of Zelditch et al Second Edition
myGPA<-gpagen(myData, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
       max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
summary(myGPA)
myGPA$coords
myGPA$Csize
#Plot new standardized landmark coordinates after GPA
plotAllSpecimens(myGPA$coords,mean=T)
#Try mean=T option to plot the consensus landmark positions or mean=F to exclude them
?plotAllSpecimens
#Data Analysis next
#First create a new variable that includes the site and convert it to a factor
myGPA$site=c("lake","lake","pond","pond")
is.factor(myGPA$site)
myGPA$site<-as.factor(myGPA$site)
myGPA$site
#Perform a linear model, MANCOVA
#Response variable is the landmark coordinates 
#Explanatory variable is the site
#Covariate is the log of the centroid size
#You have a few options for data frames and analyses, all yielding same results:
#OPTION 1
gdf1 <- geomorph.data.frame(coords = myGPA$coords, site = myGPA$site,
                           logcs = log(myGPA$Csize))
gdf1
model1<-procD.lm(coords ~ site + logcs, data = gdf1, iter = 999)
anova(model1)
#OPTION 2
y<-two.d.array(myGPA$coords)
y
model2<-procD.lm(y~gdf1$site+gdf1$logcs,iter=999)
anova(model2)
#OPTION 3
#Use option 3 as it is the simplest:
model3<-procD.lm(myGPA$coords~myGPA$site+log(myGPA$Csize),iter=999)
anova(model3)
#Examine effects of centroid size alone (ALLOMETRY).
modelAllo <- procD.lm(coords ~ logcs,data=gdf1)
summary(modelAllo)
#What if we want to test for an interaction between site and size?
modelInteraction<-procD.lm(myGPA$coords~myGPA$site*log(myGPA$Csize),iter=999)
anova(modelInteraction)
#Understand what type I/II/II SS and what permutation of residuals means
#Understand why log the centroid size (more scatter as size increases)
hist(myGPA$Csize,xlab="Centroid Size",ylab="Number of fish",main="")
hist(log(myGPA$Csize),xlab="Log transformed centroid size",ylab="Number of fish",main="")
#What we really want is for data to be normally distributed within each group
#and for variances to be similar.
#Can use lattice package to compare groups:
library(lattice)
histogram(~myGPA$Csize|myGPA$site)
#Let's do a PCA and plot the coordinates
#Calculate the PC and determine the proportion of the variation explained by each
PCA <- gm.prcomp(myGPA$coords)
summary(PCA)
plot(PCA, pch = 22, bg = c(rep("red", 2), rep("blue", 2)), cex = 1.5)
#After using the following code, you need to click on figure where you want legend to appear
legend(locator(1),levels(myGPA$site),pch=15, cex=0.8, col=c("red","blue"))
#Can also use the following code if you prefer this figure format.
#It displays differences in shape along each of the two axes. These are called deformation grids
plotTangentSpace(myGPA$coords,groups=myGPA$site,legend=T,label=T)
plotTangentSpace(myGPA$coords, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, label = TRUE,
                 groups = myGPA$site, legend = TRUE)
#If we want to only visualize the shape changes represented by PC scores:
res <- plotTangentSpace(myGPA$coords, groups = factor(paste(myGPA$site)))
ref<-mshape(myGPA$coords)
#The first shows us what the minimum PC1 score represents in terms of shape change:
plotRefToTarget(M1=ref, M2=res$pc.shapes$PC1min, method="TPS")
#The second shows us what the maximum PC1 score represents in terms of shape change:
plotRefToTarget(M1=ref, M2=res$pc.shapes$PC1max, method="TPS")
#You see that PC1 represents a shift in the tail from a downward tilt to an upward tilt.
#It is most likely due to error in placing the fish in a straight plane before taking photos.
#And the next two lines of code do the same for PC2:
plotRefToTarget(M1=ref, M2=res$pc.shapes$PC2min, method="TPS")
plotRefToTarget(M1=ref, M2=res$pc.shapes$PC2max, method="TPS")
#You see that PC2 represents a shortening and fattening of the fish.
#Get get a list of PC scores for each specimen in Morpho package
library(Morpho)
RW<-relWarps(myData, scale = TRUE, CSinit = TRUE, alpha = 0, 
         orp = TRUE,  noalign = FALSE)
RW
#Using code above, the Var values for expected (exVar) and cumulative (cumVar) 
#variance should be same values you got using geomorph.
#The bescores are individual PC scores for each specimen.
#Plot them below, as above:
plot(RW$bescores,pch = 22, bg = c(rep("red", 2), rep("blue", 2)), cex = 1.5)

