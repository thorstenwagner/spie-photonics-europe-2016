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

library("randomForest");

lvl.to.num <- function(a){
	return(as.numeric(as.character(a)));
}

makePrediction <- function(features,model,title){
	###Elongation
	elong <- lvl.to.num (features$elong);

	###Fractal dimension
	fd <- lvl.to.num (features$fd);


	###Long time / short time diffusion coefficient ratio
	LtStRatio <- lvl.to.num (features$LtStRatio);


	combinedData <- data.frame(ELONG=elong,FD=fd,LTST.RATIO=LtStRatio);

	features.predict <- predict(model,combinedData);
	
	print(title);
	print(table(features.predict))
}


load("randomForestModel.RData");

##########################################################
###                   PREDICTION                       ###
##########################################################

load("tracks_nta_many_free.RData");
features <- help[[2]];
makePrediction(features,features.model,"### PREDICTION OF FREE DIFFUSION 100nm POLYSTYRENE PARTICLES (NTA) ###");

load("tracks_cyto_confined.RData");
features <- help[[2]];
makePrediction(features,features.model,"### PREDICTION OF CONFINED DIFFUSION 50nm GOLD PARTICLES (DARKFIELD) ###");

load("tracks_cyto_active.RData");
features <- help[[2]];
makePrediction(features,features.model,"### PREDICTION OF ACTIVE DIFFUSION 50nm GOLD PARTICLES (DARKFIELD) ###");

load("tracks_laser_confined.RData");
features <- help[[2]];
makePrediction(features,features.model,"### PREDICTION OF CONFINED DIFFUSION 50nm GOLD PARTICLES (LASER-SCANNING) ###");

load("tracks_laser_active.RData");
features <- help[[2]];
makePrediction(features,features.model,"### PREDICTION OF ACTIVE DIFFUSION 50nm GOLD PARTICLES (LASER-SCANNING) ###");
