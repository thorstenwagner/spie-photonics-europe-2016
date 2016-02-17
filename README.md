# Data for my SPIE Photonics Europe 2016 publication

This the data and scripts to reproduce the results my publication at the SPIE Photonics West 2016:

*T. Wagner, A. Kroll, M. Wiemann, H.-G. Lipinski, "Classification of nanoparticle diffusion processes in vital cells by a multifeature random forests approach: Application to simulated data, darkfield and confocal laser scanning microscopy", SPIE Photonics Europe, Brussels, Belgium, 2016*

##Trajectory data
#### Simulated trajectory data
* tracks_active.RData: Simulated active trajectories
* tracks_anomalous.RData: Simulated anomalous trajectories
* tracks_confined.RData: Simulated confined trajectories
 
#### Real world trajectory data
* tracks_cyto_active.RData:  Active diffusion trajectories of 50 nm gold in V79 fibroblasts recorded by darkfield microscopy
* tracks_cyto_confined.RData: Confined diffusion trajectories of 50 nm gold in V79 fibroblasts recorded by darkfield microscopy
* tracks_laser_active.RData: Active diffusion trajectories of 50 nm gold in V79 fibroblasts recorded by confocal laser scanning microscopy
* tracks_laser_confined.RData: Confined diffusion trajectories of 50 nm gold in V79 fibroblasts recorded by confocal laser scanning microscopy
* tracks_nta_many_free.RData: Free diffusion trajectories of 100 nm polystyrene particles in V79 fibroblasts recorded by confocal laser scanning microscopy

##Objects
* randomForestModel.RData: RObject with trained random forest model

##Scripts
* TrainModel.R: Script to train the model and evaluate it with simulated trajectories.
* ClassifiyTracks.R: Script which classifies real world trajectories
