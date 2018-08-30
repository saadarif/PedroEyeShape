#A script to to sliding, gpa, PCA and extraction of PCA or Shape data
#using geomorph package
#26/04/2018 SA
#This updated and some of the unnecessary ANOVA stuff has been removed.


#Need to install the geomorph package, this only needs to be done ONCE
#install.packages("geomorph")

#update to the latest stable release of geomorph
devtools::install_github("geomorphR/geomorph",ref="Stable")

#You may have to install the devtools package first if you get an error above (with install.packages("devtools"))

#load the package to use it
library("geomorph")


#---------------------------------------------------------------------------------------------------------------------
#1. Preparing files

#Most of this only needs to be done once and files can be reused unless landmarks are changed

#Set the working directory to where the landmark files (.tps) are stored
#this is the directory on my computer
setwd("/media/data_disk/PROJECTS/Saad/Pedro_morph/PedroEyeShape/data/")

#Read in the appropriate landmark file
#these files are in the tps format
#The specID option tells the function where the sample names are located
#I joined the 4 individual files into a single one called "All.tps"
head <- readland.tps("Z4T10All.tps", specID = "ID")

#First we need to define the semilandmarks and this only needs to be done once
#unless we change landmarks later on
#define.sliders(head, nsliders=30)
#this saves a file with the identity of the semi-lanmarks in the current
#working directory as "curveslide.csv", which we will need to read in to do the sliding
curves <- as.matrix(read.csv("../curveslide.csv", header=T))
#Similarly we can define a links file for help visualization
#define.links(head)
#I have saved this file as a text file and can be read in and used later
links <- as.matrix(read.table("../links.txt"))

#Finally we need to define pairs of landmarks for object symmetry
# I just do this manually below
pairs <- cbind(c(34, 33, 1, 8, 7, 6, 5, 4, 3, 2, 14, 13, 12, 11, 10, 9, 40, 41, 42, 35, 30, 29 ), 
               c(38, 39,15,22,21,20,19,18,17,16, 28, 27, 26, 25, 24 ,23, 43, 45, 44, 37, 32, 31))

#finally we need to extract indvidual names so we can use these later
x <- dimnames(head)[[3]]
ind <- sapply(strsplit(x, '\t'), function(x){x[1]})
#lets also create a vector for pops
pops <- as.factor(sapply(strsplit(x, '_'), function(x){paste0(x[1], "_", x[2])}))
#It's probably best to check all these additional files to make sure they are correct
#as we keep resuing them
#---------------------------------------------------------------------------------------------------------------------
#2. Sliding landmarks and doing the Procrustes fit

#We are now read to slide our landmarks
#We need the information on the semi-landmarks here
head.gpa <- gpagen(head, curves= curves)

#Now the sliding of the landmarks is done and we can do the actual fit with object symmetry
#we need to specify the individuals and provide the pairing of landmarks
head.sym <- bilat.symmetry(A = head.gpa , ind=ind, object.sym = TRUE, 
                             land.pairs=pairs, RRPP = TRUE, iter = 499)

#Do the PCA
headPCA <- plotTangentSpace(A=head.sym$symm.shape, axis1=1, axis2=2,groups=pops, labels=ind) 

#This plot does suggest that PC1 might just be an error although it is difficult to look at deformation
#grids as they are upside down NEED to FIX

#PC 2 explains ~20% of the variation and seems to distinguish between strains
#extract PC2 score fo use as a phenotype column  for R/QTL etc.
PC2Scores <-  headPCA$pc.scores[,2]

#Generate deformation grids
#Calculate consensus shape
consensus <- mshape(head.sym$symm.shape)
#find indvidual with max PC1 score and plot is magnified by 2x to exaggerate deformation
par(mfrow=(c(1,2)))
x <-which.max(headPCA$pc.scores[,1])
plotRefToTarget(consensus,head.sym$symm.shape[,,x], method="vector", mag =2, links=links)
#do the same for other end of PC1
x <-which.min(headPCA$pc.scores[,1])
plotRefToTarget(consensus,head.sym$symm.shape[,,x], method="vector", mag =2, links=links)
#You can save the plots as pdf and add them whichever way is suitable

#Do the same for PC2
x <-which.max(headPCA$pc.scores[,1])
plotRefToTarget(consensus,head.sym$symm.shape[,,x], method="TPS", mag =2, links=links)
#do the same for other end of PC1
x <-which.min(headPCA$pc.scores[,1])
plotRefToTarget(consensus,head.sym$symm.shape[,,x], method="TPS", mag =2, links=links)

#Do the same for PC3
x <-which.max(headPCA$pc.scores[,2])
plotRefToTarget(consensus,head.sym$symm.shape[,,x], method="TPS", mag =2, links=links)
#do the same for other end of PC1
x <-which.min(headPCA$pc.scores[,2])
plotRefToTarget(consensus,head.sym$symm.shape[,,x], method="TPS", mag =2, links=links)

#------------------------------------------------------------------------
#Better PCA plots
#colour scheme for pops
col.pops <- rainbow(length(levels(pops)))
#You can also assign the names by hand using HEX codes or names of colours. 

#This vector needs to have dimension labels
names(col.pops) <- levels(pops)

#Using match() we can generate a vector of length(n) assigning a colour to each specimen
col.pops <- col.pops[match(pops, names(col.pops))]

#make axes labels using the output of head pca above
xlab <- paste("Principal Component 1 ", "(", round(headPCA$pc.summary$importance[2,1]*100, 1), "%)", sep="")
ylab <- paste("Principal Component 2 ", "(", round(headPCA$pc.summary$importance[2,2]*100, 1), "%)", sep="")

#make the plot
dev.off()
plot(headPCA$pc.scores[,2], headPCA$pc.scores[,3], pch=21, cex=1, bg=col.pops, xlab=xlab, ylab=ylab, asp=T)
legend("topleft", legend= unique(pops), pch=19,  col=unique(col.pops))
#You can change the axes above and save the plot as a pdf


#-------------------------------------------------------
#Regressing size on shape and testing for pop, sex differences

#I saved the size data as a .csv file

eyesize <- read.csv("Z4T10Parentals_EyeArea.csv", header=F)

#First we need to reorder the indiviuals for eye size in the same was as for shape
eyesize$V1=as.character(eyesize$V1)
eyesize$sqrtSize=sqrt(eyesize$V2)
#First we need to get rid of any potential whitespace in the individual vector we made earlier
ind<- trimws(ind)
#No reorder the indiviuals in the area file to the same order as the shape file
eyesize <-eyesize[order(match(eyesize[,1],ind)),]
#double check this 
View(cbind(eyesize$V1,ind))

#We need to create a relevant dataframe as above
#Note you may want to consider the square root of eye area if you want to use it as size
gdf4 <- geomorph.data.frame(shape = headPCA$pc.scores[,2:dim(headPCA$pc.scores)[2]], size = as.numeric(eyesize$V2),
                            pop = pops_only , sex = sex) # geomorph data frame

head.size.model<-procD.lm(shape ~ size*sex*pop ,  data=gdf4, iter = 9999, RRPP = F)
alloplot<-plot(head.size.model, type = "regression", pch = 21, predictor=gdf4$size, bg=col.pops, reg.type="Pred")
summary(head.size.model)

#references of min and max of regression score (instead of eyesize)
#We could plot deformation grids (from the consensus) of the smallest versus the biggest individualspar(mfrow=(c(1,2)))
x <-which.max(alloplot$RegScore)
plotRefToTarget(consensus,head.sym$symm.shape[,,x], method="TPS", mag =2, links=links)
#do the same for other end of PC1
y <-which.min(alloplot$RegScore)
plotRefToTarget(consensus,head.sym$symm.shape[,,y], method="TPS", mag =2, links=links)
#deforemation grid from min to max
plotRefToTarget(head.sym$symm.shape[,,y],head.sym$symm.shape[,,x], method="TPS", mag =2, links=links)
plotRefToTarget(head.sym$symm.shape[,,x],head.sym$symm.shape[,,y], method="TPS", mag =2, links=links)

#Alternative to above - should work now
head.size.model2<-procD.allometry(shape ~ size, ~sex*pop , logsz=F,  data=gdf4, iter = 9999, RRPP = F)
plot(head.size.model2, method = "Pred", gp.labels=T, warpgrids = F, shapes=F)
summary(head.size.model2)

#Another   test of homogeneity of slopes, including pairwise slopes comparisons - only slightly different from above
head.size.model3 <- advanced.procD.lm(shape ~ size+sex+pop, f2=~size*sex*pop, data=gdf4, groups=~sex*pop, slope=~size, iter=9999)
summary(head.size.model3)
