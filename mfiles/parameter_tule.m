function [prm]=parameter()

    
    
    prm.version='Feb 28 2023';
    
    % debug Beland found Kc and Ko were flipped, fixed
    
  %  prm.filename='/Users/baldocchi/DataFolder/CanTule/WP_input_smoke_2018.csv';
     prm.filename='D:\CanTule\TuleInMet.csv';
    
     prm.fluxdata='d:\Cantule\Tuleflux.csv';
   % prm.fluxdata='/users/baldocchi/datafolder/Cantule/Tuleflux.csv'
    
    % location, West Pond, Twitchell Island
    
    prm.time_zone=-8;     % time zone
    prm.lat_deg= 38.1074;   % latitude
    prm.long_deg= -121.6469; % longitude  
    
    prm.title={'Canveg, tule, West Pond'};
    
    prm.stomata= 'Hypo';     % amphistomatous = 2; hypostomatous = 1
    prm.hypo_amphi=1;  % hypostomatous, 1, amphistomatous, 2
    
    % canopy structure
    prm.Veg='Tule';
    prm.LAI = 3;  % leaf area index
    prm.dff= 0.1; % leaf area of each layer
    prm.nlayers=prm.LAI/prm.dff;  % number of canopy layers
    prm.dff=ones(prm.nlayers,1)*prm.dff;
    prm.sumlai=prm.LAI-cumsum(prm.dff);
    prm.sumlai(prm.sumlai <0)=0;
    
    
    
    
    prm.veg_ht=3;     % vegetation canopy height
    prm.meas_ht= 5;   % measurement height
    prm.dht=0.6*prm.veg_ht; %zero plane displacement height, m
    prm.z0 =0.1*prm.veg_ht;  % aerodynamic roughness height, m
    prm.markov=0.9;   % markov clumping factor
    prm.markov=ones(prm.nlayers,1)*prm.markov;
    
    prm.jtot=prm.nlayers;    % number of layers
    prm.jktot=prm.jtot+1;    % number of layers plus 1
    prm.jtot3=prm.nlayers*3;  % 3 times layers for Dij
    
    prm.dht_canopy=prm.veg_ht/prm.jtot;  % height increment of each layer
    prm.ht_atmos=prm.meas_ht-prm.veg_ht;
    n_atmos_layers=50;
    prm.dht_atmos=prm.ht_atmos/n_atmos_layers;
    prm.nlayers_atmos=prm.jtot + floor(prm.ht_atmos/prm.dht_atmos);

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

      prm.adens= (prm.dff(1:prm.nlayers) ./ prm.delz(1:prm.nlayers))';   % lai density
    
    
    
    %      optical properties PAR wave band
    %      after Norman (1979) and NASA/NIST report
    
    % for tule assuming the 'soil' is the water level

    prm.par_reflect = .10;
    prm.par_trans = .07;
    prm.par_soil_refl = .05;
    prm.par_absorbed = (1. - prm.par_reflect - prm.par_trans);


%         optical properties NIR wave band
%         after Norman (1979) and NASA report

    prm.nir_reflect = 0.45;
    prm.nir_trans = 0.05;      % they seem think and little with transmit
    prm.nir_soil_refl = 0.10;
    prm.nir_absorbed = (1. - prm.nir_reflect - prm.nir_trans);
    
    prm.leafangle='erectophile';
   
    % IR fluxes
    
    prm.sigma = 5.670367E-08 ;           % Stefan Boltzmann constant W m-2 K-4
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
        
        
        % leaf Ps model
        
        % Universal gas constant  

    prm.rugc = 8.314;              % J mole-1 K-1 
    prm.rgc1000 = 8314;            % gas constant times 1000.
    prm.Cp = 1005;                 % specific heat of air


%          Consts for Photosynthesis model and kinetic equations.
%          for Vcmax and Jmax.  Taken from Harley and Baldocchi (1995, PCE)
 
       % carboxylation rate at optimal temperature, umol m-2 s-1 
    
       % Amax = 27.4 μmol m−2 s−1 Knapp and Yavitt, 1995 Aquatic Botany
    prm.vcopt = 128.8;  %81, carboxylation rate at optimal temperature, umol m-2 s-1; from lit, Vcmax = 4.7 Amax, Wilson
    prm.jmopt = 211;  %133, electron transport rate at optimal temperature, umol m-2 s-1 Jmax = 1.64 Vcmax
    prm.rd25 = .22;  % dark respiration at 25 C, rd25= 0.34 umol m-2 s-1 
    
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
            
        prm.kball = 9.5;		%  Ball-Berry stomatal coefficient for stomatal conductance
        prm.bprime = .0175;	    %  intercept of Ball-Berry model, mol m-2 s-1 
        prm.bprime16 = 0.0109375;  % intercept for CO2, bprime16 = bprime / 1.6;
 
        prm.rsm = 145.0;  % Minimum stomatal resistance, s m-1.
        prm.brs=60.0;     % curvature coeffient for light response
        prm.stomata= 'Hypo';     % amphistomatous = 2; hypostomatous = 1
        prm.hypo_amphi=1;  % hypostomatous, 1, amphistomatous, 2

        prm.qalpha = .22; %   leaf quantum yield, electrons 
        prm.qalpha2 = 0.0484;   % qalpha squared, qalpha2 = pow(qalpha, 2.0);
        
        prm.lleaf = .02;       %// leaf length, m, width

      
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

        %PRANDTL NUMBER  

        %Prandtl Number

        prm.pr = prm.nuvisc / prm.dh;
        prm.pr33 = prm.pr^0.33;

        % DIFFUSIVITY OF WATER VAPOR, m2 s-1  
       
        prm.lfddv=prm.lleaf/prm.ddv;

        %SCHMIDT NUMBER FOR VAPOR  

        prm.sc = prm.nuvisc / prm.dv;
        prm.sc33 = prm.sc^0.33;

        %  SCHMIDT NUMBER FOR CO2  

        prm.scc = prm.nuvisc / prm.dc;
        prm.scc33 = prm.scc^0.33;


        %Grasshof Number  

        prm.grasshof=9.8*prm.lleaf^3/prm.nnu^2;
        
        prm.Mair = 29.;              
	    prm.dLdT = -2370.;   

       % add soil properties
      % or call soil set up
      
      % Dispersion Matrix Lagrangian model
      
      prm.npart=1000000;   % number of random walk particles, use about 10,000 for testing, up to 1M for smoother profiles

      prm.extinct=3; % extinction coefficient wind in canopy

end