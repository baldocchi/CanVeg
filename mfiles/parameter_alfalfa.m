function [prm]=parameter_alfalfa()

%  to do.. make LAI, dff and Markov vary with height, like deciduous forest
% update gsmin from Harley and Baldocchi

% Beland found error in Kc and Ko, they were flipped, fixed


    prm.version='Aug 3 2023';
    
   
  % prm.filename='d:\CanAlfalfa\AlfTwitchellMetInput.csv'; 
  %  prm.fluxdata='d:\CanAlfalfa\AlfMetTwitchellFluxes.csv';
    
  % prm.filename='D:\Canalfalfa\AlfBouldinMetInput.csv'; 
  %prm.filename='D:\CanAlfalfa\AlfBouldinMetInputL_4.csv'; 
  %prm.filename='d:\CanAlfalfa\AlfBouldinMetInputL_3to4.csv'; 
  %prm.filename='d:\CanAlfalfa\AlfBouldinMetInputL_2to3.csv'; 
  
  
  
  % 2019 d141-157, LAI ~ 4, base case
   prm.filename='D:\Canalfalfa\AlfMetBouldinInput.csv';
   prm.fluxdata='d:\CanAlfalfa\AlfMetBouldinFluxes.csv';
   
  
   
 
   % prm.fluxdata='d:\CanAlfalfa\AlfMetBouldinFluxesL_4.csv';
   % prm.fluxdata='d:\CanAlfalfa\AlfMetBouldinFluxesL_3to4.csv';
   % prm.fluxdata='d:\CanAlfalfa\AlfMetBouldinFluxesL_2to3.csv';
   
   %   prm.Dij='D:\Canalfalfa\Dij_Alfalfa.csv';
     
      prm.Dij='D:\Canalfalfa\Dij_AlfalfaL5.csv';   % D70 to 105 2020
      
      %prm.Dij='D:\Canalfalfa\Dij_AlfalfaL4.csv';
      %prm.Dij='D:\Canalfalfa\Dij_AlfalfaL3.csv';
      %prm.Dij='D:\Canalfalfa\Dij_AlfalfaL2_3.csv'
      %prm.Dij='D:\Canalfalfa\Dij_AlfalfaL2_3.csv';
      %prm.Dij ='/Users/baldocchi/DataFolder/Canalfalfa/Dij.csv';
   %   prm.Dij= 'D:\Canalfalfa\Dij_AlfalfaTwitchell.csv';  % run with case LAI , ht 0.8 m
   
   
        prm.time_zone=-8;     % time zone
        
    % location, alfalfa, Twitchell Island
 %    prm.lat_deg= 38.1074;   % latitude  
 %    prm.long_deg= -121.6469; % longitude  
    
    % Bouldin Island
    prm.lat_deg= 38.0991538;
    prm.long_deg= -121.49933;
    
    prm.title={'CanVeg, Alfalfa'};
    prm.stomata='Amphi';     % amphistomatous = 2; hypostomatous = 1
    prm.hypo_amphi=2;  % hypostomatous, 1, amphistomatous, 2
    
    % canopy structure
    prm.Veg='Alfalfa';
    
   % prm.LAI = 3.6;  % leaf area index  May day 141 to 157, 2019, Bouldin Alfalfa from LI2200
   % prm.veg_ht=0.78;  % vegetation canopy height, aerodynamic height, Bouldin, 2018
    
   prm.LAI = 5.0;  % leaf area index day 70 to 105, 2020, Bouldin Alfalfa from LI2200
   prm.veg_ht=0.80;  % vegetation canopy height, aerodynamic height, Bouldin, 2018
    
    
    
%     prm.LAI = 4;  % leaf area index >4  Bouldin Alfalfa from LI2200
%    prm.veg_ht=1.00;  % vegetation canopy height, aerodynamic height, Bouldin, 2018

    
    
    
 %    prm.LAI =2.5;   % sort data 2 to 3 lai  
  %   prm.veg_ht=0.60; % aerodynamic height for LAI range
    
%    prm.LAI =2;    % lear area index Twitchell, July 1 2015
 %    prm.veg_ht=0.65;  % Twitchell
   
    prm.dff= 0.1; % leaf area of each layer
    prm.nlayers=prm.LAI/prm.dff;  % number of canopy layers
   
   
    prm.dff=ones(prm.nlayers,1)*prm.dff;
    
   prm.sumlai=prm.LAI-cumsum(prm.dff);
   
   prm.sumlai(prm.sumlai <0)=0;
   
   
    
    prm.meas_ht= 5;   % measurement height
    prm.dht=0.6*prm.veg_ht; %zero plane displacement height, m
    prm.z0 =0.1*prm.veg_ht;  % aerodynamic roughness height, m
    prm.markov=0.95;   % markov clumping factor, 0 to 1, from LI 2200
   
    prm.markov=ones(prm.nlayers,1)*prm.markov;
    
    prm.dff_clmp=prm.dff./prm.markov;
   
       
    prm.jtot=prm.nlayers;    % number of layers
    prm.jktot=prm.jtot+1;    % number of layers plus 1
    prm.jtot3=prm.nlayers*3;  % 3 times layers for Dij
    
    prm.dht_canopy=prm.veg_ht/prm.jtot;  % height increment of each layer
    prm.ht_atmos=prm.meas_ht-prm.veg_ht;
    n_atmos_layers=50;
    prm.dht_atmos=prm.ht_atmos/n_atmos_layers;
    prm.nlayers_atmos=prm.jtot + floor(prm.ht_atmos/prm.dht_atmos);

    prm.hrs=48;   % half hour periods per day

% calculate height of layers in the canopy and atmosphere

for i=1:prm.jtot
prm.zht(i,1)=i*prm.dht_canopy;
prm.delz(i,1)=prm.dht_canopy;
end

j=0;
for i=prm.jtot+1:prm.nlayers_atmos+1
    j=j+1;
    prm.zht(i,1)=j* prm.dht_atmos + prm.veg_ht;
    prm.delz(i,1)=prm.dht_atmos;
end

     % divide by height of the layers in the canopy
     
  %   prm.adens= (prm.dff(1:prm.nlayers) ./ prm.delz(1:prm.nlayers))';   % lai density
    
     prm.adens=(prm.dff(1:prm.nlayers) ./ prm.dht_canopy)';
    
%      optical properties PAR wave band
   
    
%     Mobasheri, Mohammad Reza, and Sayyed Bagher Fatemi. 2013. 
%     'Leaf Equivalent Water Thickness assessment using reflectance at 
%     optimum wavelengths', Theoretical and Experimental Plant Physiology, 25: 196-202.


    prm.par_reflect = .05;    % 
    prm.par_trans = .05;
    prm.par_absorbed = (1. - prm.par_reflect - prm.par_trans);
      
    prm.par_soil_refl = .05;    % 
    


%         optical properties NIR wave band

    prm.nir_reflect = 0.60;   % 0.4 .based on field and Ocean Optics...wt with Planck Law. High leaf N and high reflected NIR
    prm.nir_trans = 0.20;     % 0.4 ...UCD presentation shows NIR transmission is about the same as reflectance; 80% NIR reflectance 700 to 1100 nm
    
    prm.nir_soil_refl = 0.10; % Ocean Optics spectra May 2, 2019 after cutting
    prm.nir_absorbed = (1. - prm.nir_reflect - prm.nir_trans);
    
% van der Tol, Christiaan, et al. 2016. 'A model and measurement comparison of diurnal
% cycles of sun-induced chlorophyll fluorescence of crops', Remote Sensing of Environment, 186: 663-77.

% corroborating information from 

%     Mobasheri, Mohammad Reza, and Sayyed Bagher Fatemi. 2013. 'Leaf Equivalent Water Thickness assessment
%     using reflectance at optimum wavelengths', Theoretical and Experimental Plant Physiology, 25: 196-202.
 
    
    %Mohamed, E. S., A. M. Saleh, A. B. Belal, and Abd Allah Gad. 2018. 
    %'Application of near-infrared reflectance for quantitative assessment of soil properties', 
    % The Egyptian Journal of Remote Sensing and Space Science, 21: 1-14.

  
    % Leaf angle distributions
    
    % options include spherical, planophile, erectophile, uniform, plagiophile,extremophile...
    
    prm.leafangle='spherical';
          
    % IR fluxes
    
    prm.sigma = 5.670367E-08 ;       % Stefan Boltzmann constant W m-2 K-4
    
%     Bo-Hui Tang & Zhao-Liang Li (2019) 
%     Evaluation of an ASTER emissivity product with field spectral radiance measurements for natural surfaces,
%     International Journal of Remote Sensing, 
%     40:5-6, 1709-1723, DOI: 10.1080/01431161.2018.1455243
    
    prm.ep = .985;                    % emissivity of leaves
    prm.epm1=0.015;                   % 1- ep  
    prm.epsoil = .98;                % Emissivity of soil  
    prm.epsigma=5.5566e-8;           % ep*sigma 
    prm.epm1=0.02;                   % 1- ep  
    prm.epsigma2 = 11.1132e-8;       % 2*ep*sigma
    prm.epsigma4 = 22.2264e-8;       %  4.0 * ep * sigma
    prm.epsigma6 = 33.3396e-8;       % 6.0 * ep * sigma 
    prm.epsigma8 = 44.448e-8;        % 8.0 * ep * sigma 
    prm.epsigma12= 66.6792e-8;       % 12.0 * ep * sigma
    
            
    prm.ir_reflect= 1-prm.ep ;        % tranmittance is zero, 1 minus emissivity 
    prm.ir_trans= 0 ;
    prm.ir_soil_refl= 1-prm.epsoil;
    prm.ir_absorbed=prm.ep ;     
        
          
        % Universal gas constant  

    prm.rugc = 8.314;              % J mole-1 K-1 
    prm.rgc1000 = 8314;            % gas constant times 1000.
    prm.Cp = 1005;                 % specific heat of air, J kg-1 K-1
    
% parameters for the Farquhar et al photosynthesis model

%          Consts for Photosynthesis model and kinetic equations.
%          for Vcmax and Jmax.  Taken from Harley and Baldocchi (1995, PCE)
 
       % carboxylation rate at 25 C temperature, umol m-2 s-1 
       % these variables were computed from A-Ci measurements and the
       % Sharkey Excel spreadsheet tool
    
    prm.vcopt =  171;  %carboxylation rate at 25 C temperature, umol m-2 s-1; from field measurements Rey Sanchez alfalfa
    prm.jmopt =  259;  %electron transport rate at 25 C temperature, umol m-2 s-1, field measuremetns Jmax = 1.64 Vcmax
    prm.rd25 = 2.68;  % dark respiration at 25 C, rd25= 0.34 umol m-2 s-1, field measurements 
    
    prm.hkin = 200000.0;    % enthalpy term, J mol-1
    prm.skin = 710.0;       % entropy term, J K-1 mol-1
    prm.ejm = 55000.0;      % activation energy for electron transport, J mol-1
    prm.evc = 55000.0;      % activation energy for carboxylation, J mol-1
    

%       Enzyme constants & partial pressure of O2 and CO2
%       Michaelis-Menten K values. From survey of literature.
    

        prm.ko25 = 274.6;   % kinetic coef for O2 at 25 C, microbars  
        prm.kc25 = 419.8;   % kinetic coef for CO2 at 25C,  millibars 
        prm.o2 = 210.0;     % oxygen concentration  mmol mol-1  
        
        
        prm.ekc = 80500.0;     % Activation energy for K of CO2; J mol-1  
        prm.eko = 14500.0;     % Activation energy for K of O2, J mol-1
        prm.erd = 38000.0;     % activation energy for dark respiration, eg Q10=2 
        prm.ektau = -29000.0;  % J mol-1 (Jordan and Ogren, 1984)
        prm.tk_25 = 298.16;    % absolute temperature at 25 C
       
        
        % Ps code was blowing up at high leaf temperatures, so reduced opt
        % and working perfectly
        
        prm.toptvc = 303.0;    % optimum temperature for maximum carboxylation, 311 K
        prm.toptjm = 303.0;    % optimum temperature for maximum electron transport,311 K
            
        prm.kball = 8.17;% % Ball-Berry stomatal coefficient for stomatal conductance, data from Camilo Rey Sanchez bouldin Alfalfa
       
        prm.bprime=0.05;   % intercept leads to too big LE0.14;    % mol m-2 s-1 h2o..Camilo Rey Sanchez..seems high
        
        prm.bprime16 = prm.bprime/1.6;  % intercept for CO2, bprime16 = bprime / 1.6;
        
        
        % simple values for priors and first iteration
        
        prm.rsm = 145.0;  % Minimum stomatal resistance, s m-1.
        prm.brs=60.0;     % curvature coeffient for light response
       

        prm.qalpha = .22;       %  leaf quantum yield, electrons 
        prm.qalpha2 = 0.0484;   % qalpha squared, qalpha2 = pow(qalpha, 2.0);
        
        prm.lleaf = 0.04;       % leaf length, m, alfalfa, across the trifoliate	
        

      
               % // Diffusivity values for 273 K and 1013 mb (STP) using values from Massman (1998) Atmos Environment
                %// These values are for diffusion in air.  When used these values must be adjusted for
                %// temperature and pressure
                
                %// nu, Molecular viscosity 
                

        prm.nuvisc = 13.27;   % // mm2 s-1
        prm.nnu = 0.00001327;  %// m2 s-1

        %// Diffusivity of CO2 

        prm.dc = 13.81;        % // mm2 s-1
        prm.ddc = 0.00001381;  % // m2 s-1

        %//   Diffusivity of heat

        prm.dh = 18.69;         %// mm2 s-1
        prm.ddh = 0.00001869;   %// m2 s-1
            

        %//  Diffusivity of water vapor 

        prm.dv = 21.78;        % // mm2 s-1
        prm.ddv = 0.00002178;   %// m2 s-1 
 

        %// Diffusivity of ozone  

        prm.do3=14.44;          % mm2 s-1
        prm.ddo3 = 0.00001444; % m2 s-1
    
        prm.betfac=1.5;    %the boundary layer sheltering factor from Grace
        
        %Constants for leaf boundary layers
       
        prm.lfddh=prm.lleaf/prm.ddh;

        %Prandtl Number

        prm.pr = prm.nuvisc / prm.dh;
        prm.pr33 = prm.pr^0.33;

           
        prm.lfddv=prm.lleaf/prm.ddv;

        %SCHMIDT NUMBER FOR VAPOR  

        prm.sc = prm.nuvisc / prm.dv;
        prm.sc33 = prm.sc^0.33;

        %  SCHMIDT NUMBER FOR CO2  

        prm.scc = prm.nuvisc / prm.dc;
        prm.scc33 = prm.scc^0.33;


        %Grasshof Number  

        prm.grasshof=9.8*prm.lleaf^3/prm.nnu^2;
        
        prm.Mair = 28.97;              
	    prm.dLdT = -2370.;   
        
        
        prm.extinct=2; % extinction coefficient wind in canopy

           
      % Dispersion Matrix Lagrangian model
      
      prm.npart=500000;   % number of random walk particles, use about 10,000 for testing, up to 1M for smoother profiles

end