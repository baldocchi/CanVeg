# CanVeg
CANVEG is a one-dimensional, multi-layer biosphere-atmosphere gas exchange model to compute water vapor, CO2 and isoprene flux densities. The model consists of coupled micrometeorological and eco-physiological modules.  The micrometeorological modules compute leaf and soil energy exchange, turbulent diffusion, scalar concentration profiles and radiative transfer through the canopy.  Environmental variables, computed with the micrometeorological module, in turn, drive the physiological modules that compute leaf photosynthesis, stomatal conductance, transpiration and leaf, bole and soil/root respiration, as well as isoprene emission rates. The micrometeorological turbulence diffusion is based on Lagrangian theory that produces a dispersion matrix.

Canveg.m is the main program. It calls specific parameter subroutines for the canopy being studied. Here we call input met files, set up latitude and longitude, set canopy, vegetaion, and soil conditions, etc

Details about the model are in the file Canveg_documentation

The architectural flow chart of the code runs as follows.  

Subroutines	

[prm]=parameter_alfalfa();	Reads all the model parameters as a structure

[DIJ]=DispCanveg_v2a(prm);  	%Compute the Dispersion matrix. It only changes with LAI, so it can be run once, and output can be saved and read for future runs

inmet=csvread(‘AlfMetBouldinInput.csv');
%	Input meteorological conditions; check code for inputs and their order. Day, hour, solar radiation, air temperature, humidity, wind speed, friction velocity, pressure, CO2 concentration are main inputs. %Also can add aerodynamic height, inferred LAI, and water table.

[soil]=fSetSoilAlfalfa(met,prm); 
	%Set and input soil parameters

[quantum,nir,ir,Qin,rnet,Sun,Shade]=fZeroArrays(prm);
	Zero and pre-allocate Arrays for faster execution
 
sunang= fSunAngle(prm.time_zone,prm.lat_deg,
prm.long_deg, met.day, met.hhour);
	
% Compute sun angles
[sunrad]=fDiffuse_Direct_Radiation (met.rglobal,met.parin,met.P_kPa, sunang.sine_beta);
	%Computes fraction of direct and diffuse light to apply the model on sun and shaded leaf fractions

[leafang]=LeafAngle(sunang, prm);   
	Computes the G function, the direction cosine between leaf normal and the sun for prescribed leaf angle distributions

[prof]=initial_profile_Matrix(met,prm); 
	Pre-allocates arrays for profiles of scalars and fluxes
[prof.wind]=fUZ_Matrix(met,prm);
	Computes empirical wind speed profile in canopy to estimate boundary layer conductances
[quantum]=fRadTranCanopy_Matrix(sunang,leafang,
quantum,waveband,prm);
	Computes Radiative transfer for specific wave bands, eg PAR, NIR. If you have spectral information you look at narrow wave bands, SIF, etc
Start Iteration Loop {…	
[ir]= fIR_RadTranCanopy_Matrix(leafang,ir,
quantum,met,Sun,Shade, prm); 
	Computes IR fluxes based on information on leaf temperature
[Qin]=fQin_Matrix(quantum,nir,ir,prm);
	Computes incoming longwave and shortwave energy absorbed by each layer
[Sun,Shade]=fEnergy_Carbon_Fluxes_Matrix(Sun, Shade, Qin, quantum, met, prof, prm);
	Computes carbon, water and energy fluxes for each layer for sun and shaded leaf fractions. This master sub routine calls routines for leaf boundary layer resistance, leaf energy balance and photosynthesis.  Leaf -air temperature is used to assess whether the boundary layer is affected by convection or free convection.
[soil]=fSoilEnergyBalanceMatrix(quantum, nir,ir, met, prof, prm, soil,j);
	Compute Soil Energy balance Fluxes 
[prof.Tair_K]=fConcMatrix(prof.H,soilflux, prof.delz, Dij,met,met.T_air_K, prm, fact.heatcoef);
	Compute sources and sink strengths and use them and Dij to compute scalar profiles for air temperature, humidity and CO2.
End Iteration Loop	
	
	Sum Layers and compute Canopy fluxes weighted for sun/shade fractions
	Plot and visualize
	
	


