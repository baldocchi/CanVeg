function [soil] = fSetSoilAlfalfa(met,prm, soil)

% 
%   Routines, algorithms and parameters for soil moisture were from 
% 
%   Campbell, G.S. 1985. Soil physics with basic. Elsevier
% 
%   updated to algorithms in Campbell and Norman and derived from Campbell et al 1994 Soil Science

%   Also confirming based on new Soil Physics in Python by Marco Bittelli

%    http://marcobittelli.com/?page_id=767


%    Nov 29, 2020
		
        

%        Soil water content 

         soil.water_content_15cm =met.soilmoisture;    % //  measured at 10 cmwater content of soil m3 m-3  


 %       // Water content of litter. Values ranged between 0.02 and 0.126

        soil.water_content_litter = .0;   %// assumed constant but needs to vary


  %       // soil content
  
  
  
  %%   fraction porosity + mineral + organic = 1
  
  %   airborne fraction = porosity - volumetric water content
  
        soil.bulkdensity=1.06;     % g cm-3   Data from Tyler Anthony
        soil.bulkdensity_kg_m3 = soil.bulkdensity *100*100*100/1000;
        
        soil.pore_fraction = 1-soil.bulkdensity/2.65;    %// from alfalfa, 1 minus ratio bulk density 1.00 g cm-3/2.65 g cm-3, density of solids
        

        soil.clay_fraction = .3;      %//  Clay fraction   
		soil.peat_fraction = 0.08;    %//  SOM = a C; C = 4.7%, a = 1.72  Kuno Kasak, 2019 Bouldin Alfalfa
       
        soil.mineral_fraction= 1-soil.pore_fraction -soil.peat_fraction ; % // from bulk density asssuming density of solids is 2.65
        
        soil.air_fraction = soil.pore_fraction - met.soilmoisture;

		soil.Cp_water= 4180;  % // J kg-1 K-1, heat capacity
		soil.Cp_air =  1065;
		soil.Cp_org = 1920; 
		soil.Cp_mineral = 870;

		soil.K_mineral= 2.5;  % // W m-1 K-1, thermal conductivity
		soil.K_org= 0.8;
		soil.K_water= 0.25;
     
        soil.dt = 20.;        %  // Time step in seconds
 
        soil.mtime = floor(1800 /soil.dt);   %// time steps per half hour
        
        soil.n_soil=10;  %number of soil layers
        soil.n_soil_1=soil.n_soil+1;  %number of soil levels
        soil.n_soil_2=soil.n_soil+2;  %number of soil levels
        
        
        % pre allocate space

         soil.Cp_soil=zeros(prm.nn,soil.n_soil);
         soil.K_soil=zeros(prm.nn,soil.n_soil);
         soil.cp_soil=zeros(prm.nn,soil.n_soil);
         soil.k_conductivity_soil=zeros(prm.nn,soil.n_soil);
         
         soil.evap=zeros(prm.nn,1);
         soil.heat=zeros(prm.nn,1);
         soil.lout=zeros(prm.nn,1);
         soil.rnet=zeros(prm.nn,1);
         soil.gsoil=zeros(prm.nn,1);

           %     //  Assign soil layers and initial temperatures
           %     //  and compute layer heat capacities and conductivities

           
         soil.T_soil = ones(prm.nn,soil.n_soil_1) .* met.Tsoil;
         soil.T_soil_old=soil.T_soil;
         
        % Compute soils depths, from 0 to base
        % we define n layers and n+1 levels, including top and bottom
        % boundaries
        
        % define Geometric Depths in Soil 
        
        soil.z_soil=zeros(soil.n_soil_1,1);
        soil.dz=zeros(soil.n_soil_1,1);
        soil.d2z=zeros(soil.n_soil,1);
      
        
        depth=0.15;   % depth of Tsoil lower boundary condition, n+1
                       
        soil.z_soil(1)=0;
        
        % exponential change in soil depth
        
        %nfact= [1 2 4 8 16 32 64 128 256 512]   
              
        nfact=2.^(1:soil.n_soil);
        
        dz=depth/sum(nfact);
                       
        soil.dz = dz .* nfact;  % array of delta z
        
        
        for i=1:soil.n_soil
        soil.z_soil(i)=soil.dz(i)+soil.dz(i);
        end
        
       
        soil.z_soil(2:soil.n_soil_1)=soil.z_soil(1:soil.n_soil);
         soil.z_soil(1)=0;
        
        for i=2:soil.n_soil
        soil.d2z(i) = soil.z_soil(i+1)-soil.z_soil(i-1);
        end
         soil.d2z(1)= soil.dz(1);
        
       
       aream2=1;  % meter squared area
       
       % soil.vol(1)=0;
        soil.vol = aream2 * soil.dz;   % soil volume
        
      
         
         % initialize soil temperature to the deep base temperature at 15 cm
         
         soil.T_soil=ones(prm.nn,soil.n_soil_2) .*met.Tsoil +273.15;
         
         % initialize upper boundary temperature as air temperature
         % in later iterations this is reset to the radiative surface
         % temperature
         
         soil.T_soil_up_boundary=met.T_air_K;
         soil.sfc_temperature=soil.T_soil_up_boundary;
         
          soil.sfc_temperature_old= soil.sfc_temperature;
         
         soil.bulk_density=ones(prm.nn,soil.n_soil) .* 0.83;   % /soil  bulk density for the alfalfa, g cm-3
         
        
         % Call Subroutine to compute Heat Capacity and Thermal Conductivity of the 
         % Soil, weighted for mineral, water, air, organic matter
         
         [soil]=HeatCapacityThermalConductivityMatrix(met,soil);   
         
        
%                 // Heat capacity and conductivity.

         % rho Cp Volume/ 2 dt

         soil.cp_soil = soil.Cp_soil .* soil.d2z' ./ (2* soil.dt);  % transpose dz to make 2d matrix
         
         % K/dz
 	     soil.k_conductivity_soil = soil.K_soil ./ soil.dz;
         
         soil.k_conductivity_soil_bound=soil.k_conductivity_soil(:,1);

         
         %soil.k_conductivity_soil(:,1)=0;
         soil.k_conductivity_soil(:,12)=0;
     
      					
         soil.lout=prm.epsigma.*met.T_air_K.^4; %initialization
         soil.llout=soil.lout;
         
         soil.evap=zeros(prm.nn,1); % initialization
        
          [soil.resistance_h2o]=  Soil_Sfc_Res(soil.water_content_15cm);            
          %soil.resistance_h2o=2000;
          
          % lower boundary
          soil.T_soil_low_bound=met.Tsoil+273.15;
          
          soil.T_soil(:,soil.n_soil_1)= soil.T_soil_low_bound;
          
    

end




