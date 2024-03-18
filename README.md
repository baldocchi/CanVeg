# CanVeg
CANVEG is a one-dimensional, multi-layer biosphere-atmosphere gas exchange model to compute water vapor, CO2 and isoprene flux densities. The model consists of coupled micrometeorological and eco-physiological modules.  The micrometeorological modules compute leaf and soil energy exchange, turbulent diffusion, scalar concentration profiles and radiative transfer through the canopy.  Environmental variables, computed with the micrometeorological module, in turn, drive the physiological modules that compute leaf photosynthesis, stomatal conductance, transpiration and leaf, bole and soil/root respiration, as well as isoprene emission rates. The micrometeorological turbulence diffusion is based on Lagrangian theory that produces a dispersion matrix.

Canveg.m is the main program. It calls specific parameter subroutines for the canopy being studied. Here we call input met files, set up latitude and longitude, set canopy, vegetaion, and soil conditions, etc

Details about the model are in the file Canveg_documentation

The architectural flow chart of the code runs as follows.  

%Subroutines	

%Reads all the model parameters as a structure. You will need to edit the parameter file for the canopy
% you are working with or create a new one. It will include information on latitude and longitude, plant structure
% and function, soil conditions, parameters for the stomatal conductance and photosynthesis model, etc, etc. One stop
% shopping for key model information relevant to the canopy of interest.

[prm]=parameter_alfalfa();

%Compute the Dispersion matrix. It only changes with LAI and canopy  height, so it can be run 'once' for that conditions of structure and function
% Output can be saved and read for future runs

[DIJ]=DispCanveg_v2a(prm);  

%  Input meteorological conditions; check code for inputs and their order. 
%  Day, hour, solar radiation, air temperature, humidity, wind speed, friction velocity, pressure, ...
% CO2 concentration are main inputs. %Also can add aerodynamic height, inferred LAI, and water table.

inmet=csvread(‘AlfMetBouldinInput.csv');

%Set and input soil parameters

[soil]=fSetSoilAlfalfa(met,prm); 
	
%Zero and pre-allocate Arrays for faster execution

[quantum,nir,ir,Qin,rnet,Sun,Shade]=fZeroArrays(prm);
	
 % Compute sun angles
 
sunang= fSunAngle(prm.time_zone,prm.lat_deg, prm.long_deg, met.day, met.hhour);
	
%Computes fraction of direct and diffuse light to apply the model on sun and shaded leaf fractions

[sunrad]=fDiffuse_Direct_Radiation (met.rglobal,met.parin,met.P_kPa, sunang.sine_beta);

%Computes the G function, the direction cosine between leaf normal and the sun for prescribed leaf angle distributions

[leafang]=LeafAngle(sunang, prm);   
	
%Pre-allocates arrays for profiles of scalars and fluxes

[prof]=initial_profile_Matrix(met,prm); 

%Computes empirical wind speed profile in canopy to estimate boundary layer conductances	

[prof.wind]=fUZ_Matrix(met,prm);
	
%Computes Radiative transfer for specific wave bands, eg PAR, NIR. If you have spectral information you look at narrow wave bands, SIF, etc

[quantum]=fRadTranCanopy_Matrix(sunang,leafang, quantum,waveband,prm);

[nir]=fRadTranCanopy_Matrix(sunang,leafang, nir,waveband,prm);
	
%Start Iteration Loop {…

%Computes IR fluxes based on information on leaf temperature

[ir]= fIR_RadTranCanopy_Matrix(leafang,ir,quantum,met,Sun,Shade, prm); 

 %Computes incoming longwave and shortwave energy absorbed by each layer
 
[Qin]=fQin_Matrix(quantum,nir,ir,prm);

%Computes carbon, water and energy fluxes for each layer for sun and shaded leaf fractions. 
%This master sub routine calls routines for leaf boundary layer resistance, leaf energy balance 
%and photosynthesis.  Leaf -air temperature is used to assess whether the boundary layer is affected by convection or free convection.	

[Sun,Shade]=fEnergy_Carbon_Fluxes_Matrix(Sun, Shade, Qin, quantum, met, prof, prm);
	
%Compute Soil Energy balance Fluxes 

[soil]=fSoilEnergyBalanceMatrix(quantum, nir,ir, met, prof, prm, soil,j);

%Compute sources and sink strengths and use them and Dij to compute scalar profiles for air temperature, humidity and CO2.

[prof.Tair_K]=fConcMatrix(prof.H,soilflux, prof.delz, Dij,met,met.T_air_K, prm, fact.heatcoef);
	
%End Iteration Loop	
	
%Sum Layers and compute Canopy fluxes weighted for sun/shade fractions
	
%Plot and visualize
	
	


