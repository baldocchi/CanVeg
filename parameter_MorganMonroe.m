function [prm]=parameter_MorganMonroe()

      
    prm.version='Feb 28 2023';
    
    % debug, Beland found Kc and Ko flipped, fixed
    
   % data from Fluxnet 2015 for Morgan Monroe, IN
   
   % data for summer 2010, days 150 to 200
    
   prm.filename = 'C:\Users\baldocchi\Documents\MATLAB\CanVeg\MorganMonroeInput.csv';
    
   prm.fluxdata='C:\Users\baldocchi\Documents\MATLAB\CanVeg\MorganMonroeFlux.csv';


    % location,Morgan Monroe, IN
    
    prm.time_zone=-6;     % time zone
    prm.lat_deg= 	39.3232 ;
    prm.long_deg=  -86.4131;
    
    prm.title={'CanVeg, DeciduousForest'};
    prm.stomata='Hypo';     % amphistomatous = 2; hypostomatous = 1
    prm.hypo_amphi=1;  % hypostomatous, 1, amphistomatous, 2
    
    % canopy structure
    prm.Veg='DeciduousForest';
    
     prm.hrs=24;  % half hour periods per day
     
     
     laidata='Lidar';
     laidata='Spline';
     
     
     switch laidata
         case 'Lidar'
             
             lai=readtable('C:\Users\baldocchi\Documents\MATLAB\CanVeg\profile Morgan.xlsx');
        
    prm.LAI=sum(lai.LAI);   % LIDAR LAI for each layer
    
    
    % prm.sumlai=prm.LAI:-prm.dff:0;  % cumulative LAI from top to bottom   
    
    prm.sumlai=prm.LAI-cumsum(lai.LAI);
     
    prm.veg_ht=35.3;  % vegetation canopy height, LIDAR
 
  %  prm.dff= 0.1; % leaf area of each layer
    prm.dff=lai.LAI;    % each layer has a different LAI increment
  
   % prm.nlayers=int32(prm.LAI/prm.dff);  % number of canopy layers
    prm.nlayers=length(lai.LAI);   % number of layers from LIDAR..top level was 0
        
    prm.meas_ht= 46;   % measurement height
    prm.dht=0.6*prm.veg_ht; %zero plane displacement height, m
    prm.z0 =0.1*prm.veg_ht;  % aerodynamic roughness height, m
    %prm.markov= 0.8 Beland;   % markov clumping factor, 0 to 1, from LI 2200
    prm.markov=lai.clumpingFactor;  % clumping factor for each height
    
       
    prm.jtot=int32(prm.nlayers);    % number of layers..had to int them. int8 did not work for 145 layers
    prm.jktot=prm.jtot+1;    % number of layers plus 1
    prm.jtot3=int32(prm.nlayers*3);  % 3 times layers for Dij
    
    prm.dht_canopy=double(prm.veg_ht)/double(prm.jtot);  % height increment of each layer
    prm.ht_atmos=prm.meas_ht-prm.veg_ht;
    n_atmos_layers=30;   % number of atmospheric layers
    prm.dht_atmos=prm.ht_atmos/n_atmos_layers;   % increment of atmospheric layer from top of canopy to tower
    
    
    %prm.nlayers_atmos=prm.jtot + floor(prm.ht_atmos/prm.dht_atmos);    % number of canopy and atmospheric layers
     prm.nlayers_atmos=int32(prm.nlayers + n_atmos_layers);
     
     prm.zht=zeros(prm.nlayers_atmos+1,1);

% calculate height of layers in the canopy and atmosphere

% for i=1:prm.jtot
% prm.zht(i,1)=double(i)*prm.dht_canopy;
% prm.delz(i,1)=prm.dht_canopy;
% end

j=0;
for i=prm.jtot+1:prm.nlayers_atmos+1
    j=j+1;
    prm.zht(i,1)= j .* prm.dht_atmos + prm.veg_ht;
    prm.delz(i,1)=prm.dht_atmos;
end

prm.zht(1:prm.jtot,1)=lai.height(1:prm.jtot,1);
prm.delz(1:prm.jtot,1)=prm.zht(2:prm.jtot+1)-prm.zht(1:prm.jtot);


   prm.adens= (prm.dff(1:prm.nlayers) ./ prm.delz(1:prm.nlayers));   % lai density
   prm.adens=prm.adens';
    

             
         case 'Spline'
             
       % Spline fit of Martin Beland's  profiles of LAI and omega from Morgan Monroe from LiDAR
      
    
    lai=readtable('C:\Users\baldocchi\Documents\MATLAB\CanVeg\profile Morgan4thOrderSmooth.xlsx');
    
    
    prm.veg_ht=35;
     
     % spline fit LAI density, 1 m increments....but then integral does not
     % equal LAI....need to weight by number of layers in original profile
     % with number of layers in spline fit profile.  I am assuming 1 m
     % increments of canopy height, as the vegetation is about 35 m.  the
     % original profile has 114 layers that were about 0.3 m in difference

    laispline = spline(lai.height,(114/35)*lai.LAIsmooth,(1:35))';
    
    prm.nlayers=length(laispline);
    prm.nlevels=prm.nlayers+1;
    prm.jtot=prm.nlayers;
    prm.jktot=prm.nlevels;
   
     prm.LAILidar=sum(lai.LAI);   % LIDAR LAI for each layer
     prm.LAIsmooth=sum(lai.LAIsmooth);
     prm.LAIspline=sum(laispline);
     
     prm.LAI=prm.LAIspline;
          
     prm.sumlai=prm.LAI-cumsum(laispline);
     
    
    % prm.sumlai=prm.LAI:-prm.dff:0;  % cumulative LAI from top to bottom   
    
    %prm.sumlai=prm.LAI-cumsum(lai.LAI);
     
   % prm.veg_ht=35.3;  % vegetation canopy height, LIDAR
 
  %  prm.dff= 0.1; % leaf area of each layer
   % prm.laiz=lai.LAI;    % each layer has a different LAI increment
  
  
        
    prm.meas_ht= 46;   % measurement height
    prm.dht=0.6*prm.veg_ht; %zero plane displacement height, m
    prm.z0 =0.1*prm.veg_ht;  % aerodynamic roughness height, m
    %prm.markov= 0.8 Beland;   % markov clumping factor, 0 to 1, from LI 2200
    %prm.markov=lai.clumpingFactor;  % clumping factor for each height
    
    %% spline fit Markov

    markov = spline(lai.height,lai.clumpingFactor,(1:35));
    prm.markov=markov';
    
       
    prm.jtot=int32(prm.nlayers);    % number of layers..had to int them. int8 did not work for 145 layers
    prm.jktot=prm.jtot+1;    % number of layers plus 1
    prm.jtot3=int32(prm.nlayers*3);  % 3 times layers for Dij
    
    prm.dht_canopy=double(prm.veg_ht)/double(prm.jtot);  % height increment of each layer
    prm.ht_atmos=prm.meas_ht-prm.veg_ht;
    n_atmos_layers=30;   % number of atmospheric layers
    prm.dht_atmos=prm.ht_atmos/n_atmos_layers;   % increment of atmospheric layer from top of canopy to tower
    
    
    %prm.nlayers_atmos=prm.jtot + floor(prm.ht_atmos/prm.dht_atmos);    % number of canopy and atmospheric layers
     prm.nlayers_atmos=int32(prm.nlayers + n_atmos_layers);
     
     prm.zht=zeros(prm.nlayers_atmos+1,1);

% calculate height of layers in the canopy and atmosphere

% for i=1:prm.jtot
% prm.zht(i,1)=double(i)*prm.dht_canopy;
% prm.delz(i,1)=prm.dht_canopy;
% end
  %prm.zht(1:prm.jtot,1)=lai.height(1:prm.jtot,1);     % lidar profile
     prm.zht(1:prm.jtot,1)=1:prm.jtot;                    % spline fit

     prm.delz(1:prm.jtot,1)=ones(prm.jtot,1);
     prm.dff(1:prm.jtot,1)=laispline;
j=0;
for i=prm.jtot+1:prm.nlayers_atmos+1
    
    j=j+1;
    prm.zht(i,1)= j .* prm.dht_atmos + prm.veg_ht;
    prm.delz(i,1)=prm.dht_atmos;
end

         prm.adens= (prm.dff(1:prm.nlayers) ./ prm.delz(1:prm.nlayers));   % lai density
         prm.adens=prm.adens';
         otherwise             
            
             
     end
     
     
     
    
 
    
    






%      optical properties PAR wave band
   
    
%     Mobasheri, Mohammad Reza, and Sayyed Bagher Fatemi. 2013. 
%     'Leaf Equivalent Water Thickness assessment using reflectance at 
%     optimum wavelengths', Theoretical and Experimental Plant Physiology, 25: 196-202.


    prm.par_reflect = .05;    % 
    prm.par_trans = .05;
    prm.par_absorbed = (1. - prm.par_reflect - prm.par_trans);
      
    prm.par_soil_refl = .05;    % 
    


%         optical properties NIR wave band

    prm.nir_reflect = 0.35;   
    prm.nir_trans = 0.35;     
    
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
    prm.ep = .98;                    % emissivity of leaves
    prm.epm1=0.02;                   % 1- ep  
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
    
    prm.vcopt = 34;  % Riccuito 60;  %carboxylation rate at 25 C temperature, umol m-2 s-1; from Wilson
    prm.jmopt = 54;  %120;  %electron transport rate at 25 C temperature, umol m-2 s-1, field measuremetns Jmax = 1.6 Vcmax
    prm.rd25 = 0.34;  % dark respiration at 25 C, rd25= 0.34 umol m-2 s-1, field measurements Harley ddb
    
    prm.hkin = 200000.0;    % enthalpy term, J mol-1
    prm.skin = 710.0;       % entropy term, J K-1 mol-1
    prm.ejm = 55000.0;      % activation energy for electron transport, J mol-1
    prm.evc = 55000.0;      % activation energy for carboxylation, J mol-1
    

%         Enzyme constants & partial pressure of O2 and CO2
%         Michaelis-Menten K values. From survey of literature.

        prm.ko25 = 274.6;   % kinetic coef for O2 at 25 C, microbars  
        prm.kc25 = 419.8;   % kinetic coef for CO2 at 25C,  millibars 
        prm.o2 = 210.0;     % oxygen concentration  mmol mol-1  
        
        
        prm.ekc = 80500.0;     % Activation energy for K of CO2; J mol-1  
        prm.eko = 14500.0;     % Activation energy for K of O2, J mol-1
        prm.erd = 38000.0;     % activation energy for dark respiration, eg Q10=2 
        prm.ektau = -29000.0;  % J mol-1 (Jordan and Ogren, 1984)
        prm.tk_25 = 298.16;    % absolute temperature at 25 C
        prm.toptvc = 311.0;    % optimum temperature for maximum carboxylation
        prm.toptjm = 311.0;    % optimum temperature for maximum electron transport
            
        prm.kball = 9.5;		% Ball-Berry stomatal coefficient for stomatal conductance, 
               
                                % this value is close to values in Miner et al review and plots I see in papers by Franks et al
                                % Xu annd Baldocchi and Medlyn et al
                                
                                % Franks, P. J., et al. (2017). "Stomatal Function across Temporal and Spatial Scales: 
                                % Deep-Time Trends, Land-Atmosphere Coupling and Global Models." Plant Physiology 174(2): 583., 
                                
                                % Medlyn, B. E. et al 2011. Reconciling the optimal and empirical approaches to 
                                % modelling stomatal conductance. Global Change Biology 17:2134-2144.Medlyn et al
                                % the smaller yields better fluxes day and night
                                
        prm.bprime=0.0175;      % Harley and Baldocchi, 1995 PCe found
                                % intercept was 17.5 mmol H2o m-2 s-1 for
                                % oak...0.0175 mol m-2 s-1/39 mol m-3...m/s
                                % bprime, intercept of Ball-Berry model, mol m-2 s-1 
        
        prm.bprime16 = prm.bprime/1.6;  % intercept for CO2, bprime16 = bprime / 1.6;
        
        
        % simple values for priors and first iteration
        
        prm.rsm = 145.0;  % Minimum stomatal resistance, s m-1.
        prm.brs=60.0;     % curvature coeffient for light response
       

        prm.qalpha = .22;       %  leaf quantum yield, electrons 
        prm.qalpha2 = 0.0484;   % qalpha squared, qalpha2 = pow(qalpha, 2.0);
        
        prm.lleaf = 0.1;       % leaf length, m, oak leaf
        

      
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
        
        prm.extinct=1.5;

           
      % Dispersion Matrix Lagrangian model
      
       prm.Dij='C:\Users\baldocchi\Documents\MATLAB\CanVeg\Dij_DeciduousForestn35.csv';
      
      prm.npart=100000;   % number of random walk particles, use about 10,000 for testing, up to 1M for smoother profiles

end