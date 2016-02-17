package de.biomedical_imaging;

import java.util.ArrayList;

import javax.vecmath.Point3d;

import de.biomedical_imaging.traJ.ExportTools;
import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.TrajectoryUtil;
import de.biomedical_imaging.traJ.features.AspectRatioFeature;
import de.biomedical_imaging.traJ.features.ElongationFeature;
import de.biomedical_imaging.traJ.features.FractalDimensionFeature;
import de.biomedical_imaging.traJ.features.MeanSquaredDisplacmentCurvature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;
import de.biomedical_imaging.traJ.features.ShortTimeLongTimeDiffusioncoefficentRatio;
import de.biomedical_imaging.traJ.features.ShortTimeLongTimeSCDFFeature;
import de.biomedical_imaging.traJ.features.SplineCurveDynamicsFeature;
import de.biomedical_imaging.traJ.features.SplineCurveSpatialFeature;
import de.biomedical_imaging.traJ.features.StandardDeviationDirectionFeature;
import de.biomedical_imaging.traJ.simulation.AbstractSimulator;
import de.biomedical_imaging.traJ.simulation.ActiveTransportSimulator;
import de.biomedical_imaging.traJ.simulation.AnomalousDiffusionScene;
import de.biomedical_imaging.traJ.simulation.AnomalousDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.CentralRandomNumberGenerator;
import de.biomedical_imaging.traJ.simulation.ConfinedDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.FreeDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.ImmobileSphereObstacle;
import de.biomedical_imaging.traJ.simulation.SimulationUtil;

public class GenerateSetOfTracks2 {
	/*
	 *  This script generates 500 tracks per diffusion mode, length, drift type
	 *  
	 *  There are 3 different diffusion modes:
	 *   - Free diffusion
	 *   - Annomalous diffusion
	 *   - Confined diffusion (Radius is set to 5µm)
	 *   - Active transport
	 *   
	 *  There are four different drift types:
	 *   - No Drift (0µm/s)
	 *   - Constant drift (low speed, 0.27µm/s)
	 *   - Constant drift (middle speed 0.8µm/s, ermittelt aus Partikelbewegung in Capping3Color video)
	 *   - Constant drift (high speed , 2.4µm/s)
	 *   
	 *  There three different lengths
	 *   - Short (2 seconds)
	 *   - Middle (6 secondos)
	 *   - Long (18 seconds)
	 *  
	 *  The for anomalous diffusion, a scene of obstacles is generated randomly in dependence of the simulated tracklength.
	 *  The scene is generated in such a way, that the probabilty of interaction is roughly the same for all tracklengths.
	 *  The diffusion coefficient is set to 9.02*10^-14 [m^2/s]. This corresponds to a 50 nm particle at 22C° in water
	 */
	enum SIM_TYPE {
		FREE,
	    ANOMALOUS,
	    CONFINED,
	    ACTIVE
	}
	private static CentralRandomNumberGenerator r;
	
	public static void main(String[] args) {
		
		r  = CentralRandomNumberGenerator.getInstance();
		r.setSeed(22);
		//General Parameters
		double diffusioncoefficient = 9.02*Math.pow(10,-14); //[m^2/s];
		double timelag = 1.0/30; //s
		int dimension = 2;
		int[] tracklengths = {2,6,18}; //s;
		int numberOfTracks = 500;
		
		//Active transport / drift
		double[] driftspeed = {0, 0.27*Math.pow(10, -6),0.8*Math.pow(10, -6),2.4*Math.pow(10, -6)}; // m/s
		double angleVelocity = Math.PI/4.0; //rad/s
		
		//Confined diffusion parameters
		double radius_confined = 1*Math.pow(10, -6); // m;
		double probNonInteraction = 0.1;
		
		//Anomalous Diffusion parameters
		double[] excludedVolumeTresholds  = new double[3]; //0.5;
		excludedVolumeTresholds[0] = 0.2;
		excludedVolumeTresholds[1] = 0.4;
		excludedVolumeTresholds[2] = 0.6;
		double diameterMean = 1*Math.pow(10, -6);
		double diameterSD = 0.5/3 *Math.pow(10, -6);
		ArrayList<AnomalousDiffusionScene> scenes = new ArrayList<AnomalousDiffusionScene>();
		double[] size = {10.0*Math.pow(10, -6),10.0*Math.pow(10, -6)}; // 10 x 10 µm
		double probNonInteractionAnom = 0.05;
		for (double tlength : tracklengths) {
			
			AnomalousDiffusionScene scene = new AnomalousDiffusionScene(size, 2);
			double radiusNonInteract = Math.sqrt(-1*Math.log(1-probNonInteraction)*(4*diffusioncoefficient*tlength));
			System.out.println("Radius: " + radiusNonInteract);
			double p = scene.estimateProbNonInteraction(radiusNonInteract);
			while(p>probNonInteractionAnom){
				double radius = (diameterMean + r.nextGaussian()*diameterSD)/2;
				double[] pos = {r.nextDouble()*size[0],r.nextDouble()*size[1]};
				ImmobileSphereObstacle obstacle = new ImmobileSphereObstacle(pos, radius,2);
				scene.addObstacle(obstacle);
				p = scene.estimateProbNonInteraction(radiusNonInteract);

			}
			scenes.add(scene);
			
			System.out.println("Scene exVolThresh: " + scene.estimateExcludedVolumeFraction() + " Tracklength: " + tlength);
		}
		System.out.println("Scenes generated");
		
		
		SIM_TYPE[] types = SIM_TYPE.values();
		ArrayList<Trajectory> trajectorys_free = new ArrayList<Trajectory>();
		ArrayList<Trajectory> trajectorys_confined = new ArrayList<Trajectory>();
		ArrayList<Trajectory> trajectorys_active = new ArrayList<Trajectory>();
		ArrayList<Trajectory> currentList = null;
		System.out.println("Generation of Free, Confined and Active Tracks");
		for (int tracklength : tracklengths) {
			int numberOfSteps = (int)(tracklength * 1/timelag);
			for (double drift : driftspeed) {
				for (SIM_TYPE type : types) {
					for(int i = 0 ; i < numberOfTracks; i++){
						if(type != SIM_TYPE.ANOMALOUS){
							AbstractSimulator sim = null;
							String typestring = "";
							typestring += type.toString();
							Trajectory t = null;
							
							switch (type) {
							case FREE:
								sim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, dimension, numberOfSteps);
								typestring += ",D_"+String.format("%.4f", diffusioncoefficient*Math.pow(10, 12)).replace(",", ".")+
										",dt_"+String.format("%.2f", timelag).replace(",", ".");
								currentList = trajectorys_free;
								break;
							case CONFINED:
								radius_confined = Math.sqrt(-1*Math.log(1-probNonInteraction)*(4*diffusioncoefficient*tracklength));
								sim = new ConfinedDiffusionSimulator(diffusioncoefficient, timelag, radius_confined, dimension, numberOfSteps);
								currentList = trajectorys_confined;
								typestring += ",D_"+String.format("%.4f", diffusioncoefficient*Math.pow(10, 12)).replace(",", ".")+
										",dt_"+String.format("%.3f", timelag).replace(",", ".")+"_r_"+String.format("%.1f", radius_confined*Math.pow(10, 6)).replace(",", ".");
								break;
							case ACTIVE:
								sim = new ActiveTransportSimulator(drift, angleVelocity, timelag, dimension, numberOfSteps);
								currentList = trajectorys_active;
								typestring += ",DriftVelocity_"+String.format("%.3f", drift*Math.pow(10, 6)).replace(",", ".")+
										",DriftAngleVelocity_"+String.format("%.3f", angleVelocity).replace(",", ".");
								break;					
							default:
								break;
							}
							if( (type==SIM_TYPE.ACTIVE && drift==0)==false){
								t = sim.generateTrajectory();
				
								if(drift!=0 && type != SIM_TYPE.ACTIVE){
									ActiveTransportSimulator ats = new ActiveTransportSimulator(drift, angleVelocity, timelag, dimension, numberOfSteps);
									Trajectory active = ats.generateTrajectory();
									t = TrajectoryUtil.combineTrajectory(t, active);
									typestring += ",DriftVelocity_"+String.format("%.3f", drift*Math.pow(10, 6)).replace(",", ".")+
											",DriftAngleVelocity_"+String.format("%.3f", angleVelocity).replace(",", ".");
								}
								t.setType(typestring);
								currentList.add(t);
							}
						}

					}
				}
			}
		}
		System.out.println("Start anomalous track generation");
		ArrayList<Trajectory> trajectorys_anomalous = new ArrayList<Trajectory>();

		for(int j = 0; j < tracklengths.length; j++) {
			int tracklength = tracklengths[j];
			
			int numberOfSteps = (int)(tracklength * 1/timelag);
			AnomalousDiffusionScene scene = scenes.get(j);
	
					for(int i = 0 ; i < numberOfTracks; i++){
						String typestring = "";
						AnomalousDiffusionSimulator sim = new AnomalousDiffusionSimulator(diffusioncoefficient, timelag, dimension, numberOfSteps, scene);
						sim.setStartPoint(getAnomalousStartPosition(scene,diameterMean,diameterSD)); 
						Trajectory t = sim.generateTrajectory();
						t = sim.generateTrajectory();
						typestring = SIM_TYPE.ANOMALOUS.toString();
						typestring += ",D_"+String.format("%.4f", diffusioncoefficient*Math.pow(10, 12)).replace(",", ".")+
								",dt_"+String.format("%.3f", timelag).replace(",", ".")+
								",exVolFrac_"+ String.format("%.2f", scene.estimateExcludedVolumeFraction()).replace(",", ".");

						t.setType(typestring);
						trajectorys_anomalous.add(t);
						System.out.println("Anom Track generated");
					
					}
				
			
		}
		
		
		//Add features
		
		int numberOfSegmentsSplineFit = 7;
		int numberOfPointsForShortTimeLongTimeRatio = 3;
		
		
		//Free
		for (Trajectory t : trajectorys_free) {
			int maxTimeLagPowerLaw1 = t.size()/20;
			int timelagForDirectionDeviationLong = t.size()/20; 
			t.addFeature(new AspectRatioFeature(t));
			t.addFeature(new ElongationFeature(t));
			t.addFeature(new FractalDimensionFeature(t));
			t.addFeature(new MeanSquaredDisplacmentCurvature(t));
			t.addFeature(new PowerLawFeature(t, 1, maxTimeLagPowerLaw1));
			t.addFeature(new StandardDeviationDirectionFeature(t, timelagForDirectionDeviationLong));
			t.addFeature(new SplineCurveDynamicsFeature(t, numberOfSegmentsSplineFit, 1));
			t.addFeature(new SplineCurveSpatialFeature(t, numberOfSegmentsSplineFit));
			t.addFeature(new ShortTimeLongTimeDiffusioncoefficentRatio(t, numberOfPointsForShortTimeLongTimeRatio));
			t.addFeature(new ShortTimeLongTimeSCDFFeature(t, numberOfPointsForShortTimeLongTimeRatio));
		
		}
		//Confined
		for (Trajectory t : trajectorys_confined) {
			int maxTimeLagPowerLaw1 = t.size()/20;
			int timelagForDirectionDeviationLong = t.size()/20; 
			t.addFeature(new AspectRatioFeature(t));
			t.addFeature(new ElongationFeature(t));
			t.addFeature(new FractalDimensionFeature(t));
			t.addFeature(new MeanSquaredDisplacmentCurvature(t));
			t.addFeature(new PowerLawFeature(t, 1, maxTimeLagPowerLaw1));
			t.addFeature(new StandardDeviationDirectionFeature(t, timelagForDirectionDeviationLong));
			t.addFeature(new SplineCurveDynamicsFeature(t, numberOfSegmentsSplineFit, 1));
			t.addFeature(new SplineCurveSpatialFeature(t, numberOfSegmentsSplineFit));
			t.addFeature(new ShortTimeLongTimeDiffusioncoefficentRatio(t, numberOfPointsForShortTimeLongTimeRatio));
			t.addFeature(new ShortTimeLongTimeSCDFFeature(t, numberOfPointsForShortTimeLongTimeRatio));
		}
		//Anomalous
		for (Trajectory t : trajectorys_anomalous) {
			int maxTimeLagPowerLaw1 = t.size()/20;
			int timelagForDirectionDeviationLong = t.size()/20; 
			t.addFeature(new AspectRatioFeature(t));
			t.addFeature(new ElongationFeature(t));
			t.addFeature(new FractalDimensionFeature(t));
			t.addFeature(new MeanSquaredDisplacmentCurvature(t));
			t.addFeature(new PowerLawFeature(t, 1, maxTimeLagPowerLaw1));
			t.addFeature(new StandardDeviationDirectionFeature(t, timelagForDirectionDeviationLong));
			t.addFeature(new SplineCurveDynamicsFeature(t, numberOfSegmentsSplineFit, 1));
			t.addFeature(new SplineCurveSpatialFeature(t, numberOfSegmentsSplineFit));
			t.addFeature(new ShortTimeLongTimeDiffusioncoefficentRatio(t, numberOfPointsForShortTimeLongTimeRatio));
			t.addFeature(new ShortTimeLongTimeSCDFFeature(t, numberOfPointsForShortTimeLongTimeRatio));
		}
		//Active
		for (Trajectory t : trajectorys_active) {
			int maxTimeLagPowerLaw1 = t.size()/20;
			int timelagForDirectionDeviationLong = t.size()/20; 
			t.addFeature(new AspectRatioFeature(t));
			t.addFeature(new ElongationFeature(t));
			t.addFeature(new FractalDimensionFeature(t));
			t.addFeature(new MeanSquaredDisplacmentCurvature(t));
			t.addFeature(new PowerLawFeature(t, 1, maxTimeLagPowerLaw1));
			t.addFeature(new StandardDeviationDirectionFeature(t, timelagForDirectionDeviationLong));
			t.addFeature(new SplineCurveDynamicsFeature(t, numberOfSegmentsSplineFit, 1));
			t.addFeature(new SplineCurveSpatialFeature(t, numberOfSegmentsSplineFit));
			t.addFeature(new ShortTimeLongTimeDiffusioncoefficentRatio(t, numberOfPointsForShortTimeLongTimeRatio));
			t.addFeature(new ShortTimeLongTimeSCDFFeature(t, numberOfPointsForShortTimeLongTimeRatio));
		}
		
		
		System.out.println("Tracks generated");
		ExportTools.exportTrajectories(trajectorys_free, "/home/thorsten/tracks_free.xml");
		ExportTools.exportTrajectories(trajectorys_confined, "/home/thorsten/tracks_confined.xml");
		ExportTools.exportTrajectories(trajectorys_anomalous, "/home/thorsten/tracks_anomalous.xml");
		ExportTools.exportTrajectories(trajectorys_active, "/home/thorsten/tracks_active.xml");
		System.out.println("Tracks exported");
		
		for(int i = 0; i < scenes.size(); i++){
			ExportTools.exportAnomalousSceneAsXML(scenes.get(i), "/home/thorsten/obstacles_tl_"+tracklengths[i]+".xml");
		}

	}
	
	public static Point3d getAnomalousStartPosition(AnomalousDiffusionScene scene,double diameterMean, double diameterSD){
		int sdfact = 3;
		double lowerBoundX = diameterMean + sdfact*diameterSD;
		double lowerBoundY = diameterMean + sdfact*diameterSD;
		double upperBoundX = scene.getSize()[0] - (diameterMean + sdfact*diameterSD);
		double upperBoundY = scene.getSize()[1] - (diameterMean + sdfact*diameterSD);
		double x = (upperBoundX-lowerBoundX)*r.nextDouble() + lowerBoundX;
		double y = (upperBoundY-lowerBoundY)*r.nextDouble() + lowerBoundY;
		
		while(scene.checkCollision(new double[]{x,y} )){
			x = (upperBoundX-lowerBoundX)*r.nextDouble() + lowerBoundX;
			y = (upperBoundY-lowerBoundY)*r.nextDouble() + lowerBoundY;
		}
		
		return new Point3d(x, y, 0);
	}
	
	public static Trajectory addDrift(Trajectory t, double drift, double angleVelocity, double timelag, int dimension, int numberOfSteps){
		ActiveTransportSimulator ats = new ActiveTransportSimulator(drift, angleVelocity, timelag, dimension, numberOfSteps);
		Trajectory active = ats.generateTrajectory();
		t = TrajectoryUtil.combineTrajectory(t, active);
		return t;
		
	}

}
