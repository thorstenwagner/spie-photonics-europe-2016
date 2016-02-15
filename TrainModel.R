#The MIT License (MIT)
#
#Copyright (c) 2016 Thorsten Wagner (wagner@biomedical-imaging.de)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

############################
###  LOAD/INSTALL LIBS   ###
############################
if (!require("randomForest",character.only = TRUE))
{
      install.packages("randomForest",dep=TRUE,repos='http://mirrors.softliste.de/cran/')
}

if (!require("caret",character.only = TRUE))
{
      install.packages("caret",dep=TRUE,repos='http://mirrors.softliste.de/cran/')
}
if (!require("scales",character.only = TRUE))
{
      install.packages("scales",dep=TRUE,repos='http://mirrors.softliste.de/cran/')
}
lvl.to.num <- function(a){
	return(as.numeric(as.character(a)));
}

addFeatureVector <- function(feature.frame,feature.name,feature.selection,vector){
	eval.parent(substitute(vector <- c(vector,lvl.to.num(feature.frame[[feature.name]][feature.selection]))))
}

############################
###      LOAD DATA       ###
############################
reload=TRUE;
if(reload){
	#Load free tracks
	load("tracks_free.RData");
	positions.free <- help[[1]];
	features.free <- help[[2]];
	gc();
	print("Free loaded");
	#Load confined tracks
	load("tracks_confined.RData");
	positions.confined <- help[[1]];
	features.confined <- help[[2]];
	gc();
	print("Confined loaded");
	#Load anomalous tracks
	load("tracks_anomalous.RData"); 
	positions.anomalous <- help[[1]];
	features.anomalous <- help[[2]];
	gc();
	print("Anomalous loaded");
	#Load active tracks
	load("tracks_active.RData"); 
	positions.active <- help[[1]];
	features.active <- help[[2]];
	gc();
	print("Active loaded");
}


############################
# PREPARE DATA  STRUCTURE  #
############################
tracklength.short=61;
tracklength.middle=181;
tracklength.long=541;

#Indexes for free tracks
free.isShortTrack <- features.free$lengths==tracklength.short;
free.isMiddleLongTrack <- features.free$lengths==tracklength.middle;
free.isLongTrack <- features.free$lengths==tracklength.long;
free.noDrift <-features.free$driftVelo==0;
free.lowDrift<-features.free$driftVelo==0.27;
free.middleDrift<-features.free$driftVelo==0.80;
free.highDrift<-features.free$driftVelo==2.40;

#Indexes for confined tracks
confined.isShortTrack <- features.confined$lengths==tracklength.short;
confined.isMiddleLongTrack <- features.confined$lengths==tracklength.middle;
confined.isLongTrack <- features.confined$lengths==tracklength.long;
confined.noDrift <-features.confined$driftVelo==0;
confined.noDrift <-features.free$driftVelo==0;
confined.lowDrift<-features.confined$driftVelo==0.27;
confined.middleDrift<-features.confined$driftVelo==0.80;
confined.highDrift<-features.confined$driftVelo==2.40;

#Indexes for anomalous tracks
anomalous.isShortTrack <- features.anomalous$lengths==tracklength.short & features.anomalous$exVolFracs==0.95;
anomalous.isMiddleLongTrack <- features.anomalous$lengths==tracklength.middle & features.anomalous$exVolFracs==0.94;
anomalous.isLongTrack <- features.anomalous$lengths==tracklength.long & features.anomalous$exVolFracs==0.66;
anomalous.noDrift <-features.anomalous$driftVelo==0;

#Indexes for active tracks
active.isShortTrack <- features.active$lengths==tracklength.short;
active.isMiddleLongTrack <- features.active$lengths==tracklength.middle;
active.isLongTrack <- features.active$lengths==tracklength.long;
active.lowDrift<-features.active$driftVelo==0.27;
active.middleDrift<-features.active$driftVelo==0.80;
active.highDrift<-features.active$driftVelo==2.40;


#Combined indexes 
free.NoDrift.Long <- free.noDrift & free.isMiddleLongTrack;
confined.NoDrift.Long <- confined.noDrift & confined.isMiddleLongTrack;
anomalous.NoDrift.Long <- anomalous.noDrift & anomalous.isMiddleLongTrack;

free.middleDrift.Long <- free.middleDrift & free.isMiddleLongTrack;
confined.middleDrift.Long <- confined.middleDrift & confined.isMiddleLongTrack;
active.middleDrift.Long <- active.middleDrift & active.isMiddleLongTrack;



#Combine data for visualization

types <- factor(c(	rep("Free",500),
				  	rep("SUB.",1000),
					rep("SUPER",1500)
					))


###Elongation
elong <- lvl.to.num (features.free$elong[free.NoDrift.Long]);
addFeatureVector(features.confined,"elong",confined.NoDrift.Long,elong);
addFeatureVector(features.anomalous,"elong",anomalous.NoDrift.Long,elong);
addFeatureVector(features.active,"elong",active.middleDrift.Long,elong);
addFeatureVector(features.free,"elong",free.middleDrift.Long,elong);
addFeatureVector(features.confined,"elong",confined.middleDrift.Long,elong);

###Fractal dimension
fd <- lvl.to.num (features.free$fd[free.NoDrift.Long]);
addFeatureVector(features.confined,"fd",confined.NoDrift.Long,fd);
addFeatureVector(features.anomalous,"fd",anomalous.NoDrift.Long,fd);
addFeatureVector(features.active,"fd",active.middleDrift.Long,fd);
addFeatureVector(features.free,"fd",free.middleDrift.Long,fd);
addFeatureVector(features.confined,"fd",confined.middleDrift.Long,fd);

###Long time / short time diffusion coefficient ratio
LtStRatio <- lvl.to.num (features.free$LtStRatio[free.NoDrift.Long]);
addFeatureVector(features.confined,"LtStRatio",confined.NoDrift.Long,LtStRatio);
addFeatureVector(features.anomalous,"LtStRatio",anomalous.NoDrift.Long,LtStRatio);
addFeatureVector(features.active,"LtStRatio",active.middleDrift.Long,LtStRatio);
addFeatureVector(features.free,"LtStRatio",free.middleDrift.Long,LtStRatio);
addFeatureVector(features.confined,"LtStRatio",confined.middleDrift.Long,LtStRatio);

combinedData <- data.frame(ELONG=elong,FD=fd,LTST.RATIO=LtStRatio);

############################
###   SHOW SCATTERPLOT   ###
############################

free.range <- 1:500;
sub.range <- 501:1500;
super.range <- 1501:3000;
whichModes <- c(free.range,sub.range,super.range)
gc.ramp <- hue_pal()(3);
#col=types[whichModes],
pairs(combinedData[whichModes,], 
		main="Normal vs. sup. vs. super",oma=c(10,3,3,3),
		labels=c("EL","FD","LTST-DC-RATIO"),
		col = c(gc.ramp[1], gc.ramp[2], gc.ramp[3])[types[whichModes]]) 
par(xpd = TRUE) 
legend("bottom", fill = unique(c(gc.ramp[1], gc.ramp[2], gc.ramp[3])[types[whichModes]]), legend = c("norm. diffusion","subdiffusion","superdiffusion"),horiz=TRUE) 

############################
### CLASSIFICATION       ###
############################

#Shuffle features and select 90% for training and 10% for testing
set.seed(1234);
rand.shuffle <- order(runif(length(combinedData[[1]])));
combinedData.shuffled <- combinedData[rand.shuffle,];
types.shuffled <- types[rand.shuffle];
numberTrainingData <- as.integer(length(combinedData[[1]])*0.9);
numberTestData <- length(combinedData[[1]])-numberTrainingData;

features.training <- combinedData.shuffled[1:numberTrainingData,];
types.training <- types.shuffled[1:numberTrainingData];
features.test <-  combinedData.shuffled[(numberTrainingData+1):(numberTrainingData+numberTestData),];
types.test <- types.shuffled[(numberTrainingData+1):(numberTrainingData+numberTestData)];

if(FALSE){
	#CutOff: If all classes have the probability, prefer normal diffusion
	noc = nlevels(types); # number of classes
	coFree = 1/(noc+1)
	coRest = (1-coFree)/(noc-1);
	cutoff <- c(coFree,rep(coRest,noc-1)); 
	features.model <- randomForest(features.training,types.training,mtry=1,cutoff=cutoff)
} else {
	load("randomForestModel.RData");
}


# Print some performance measures
print(features.model);
features.predict <- predict(features.model,features.test);
t <- table(features.predict,types.test);
print(t);
print(prop.table(t));
print(prop.table(table(features.predict==types.test)));
print(confusionMatrix(features.predict,types.test));

if(FALSE){
	save(features.model,file="randomForestModel.RData")
}






