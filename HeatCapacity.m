function [soil] =HeatCapacity(met, soil)

% Heat Capacity

% from Bitelli Soil Physics in Python

% bulkDensity=0.83;
% waterContent=0.2;
% 
% Cp= (2.4e6 * bulkDensity / 2650 + 4.18e6 * waterContent);   
    
% from Campbell Soil Physics in Basic



%         met.soilmoisture=0.2;
%         met.P_Pa=101300;  % pressure Pa
%         met.air_density_mole = 0.0414;
%         met.air_density=1.2;
%         met.dest= 189; % Pa C-1 at 25 C
%         
%         dest=met.dest;
%         
        
%       // soil content

%         soil.clay_fraction = .3;      %//  Clay fraction   
% 		soil.peat_fraction = 0.129;    %//  SOM = a C; C = 7.5%, a = 1.72
% 		soil.pore_fraction = 0.687;    %// from alfalfa, 1 minus ratio bulk density 0.83 g cm-3/2.65 g cm-3, density of solids
% 		soil.mineral_fraction= 0.558; % // from bulk density asssuming density of solids is 2.65
% 
% 		soil.air_fraction = soil.pore_fraction - met.soilmoisture;
% 
% 		soil.Cp_water= 4180;  % // J kg-1 K-1, heat capacity
% 		soil.Cp_air =  1065;
% 		soil.Cp_org = 1920; 
% 		soil.Cp_mineral = 870;
% 
% 		soil.K_mineral= 2.5;  % // W m-1 K-1, thermal conductivity
% 		soil.K_org= 0.8;
% 		soil.K_water= 0.25;
        
        fw=1./(1+power((met.soilmoisture/0.15),-4));  %// terms for Stefan flow as water evaporates in the pores
        
        soil.K_air= 0.024 + 44100.*2.42e-5.*fw.*met.air_density_mole.* met.dest./met.P_Pa;

        k_fluid=soil.K_air + fw .*(soil.K_water-soil.K_air);
% 
		wt_air=2./(3.*(1+.2.*(soil.K_air./k_fluid-1))) + 1./(3.*(1+(1-2*.2).*(soil.K_air./k_fluid -1)));
		wt_water=2./(3.*(1+.2.*(soil.K_water./k_fluid-1))) + 1./(3.*(1+(1-2.*.2).*(soil.K_water./k_fluid -1)));
		wt_mineral=2./(3*(1+.2.*(soil.K_mineral./k_fluid-1))) + 1 ./(3.*(1+(1-2*.2).*(soil.K_mineral./k_fluid -1)));
		wt_org=2./(3.*(1+.2.*(soil.K_org./k_fluid-1))) + 1./(3.*(1+(1-2*.2).*(soil.K_org./k_fluid -1)));

		Cp_soil_num= ( met.air_density .* soil.Cp_air .* soil.air_fraction + 1000.000 .* soil.Cp_water .* met.soilmoisture +... 
			1300.000 .* soil.Cp_org .* soil.peat_fraction + 2650.000 .* soil.Cp_mineral .* soil.mineral_fraction);
    
         soil.Cp_soil=Cp_soil_num./( met.air_density .*  soil.air_fraction + 1000.000 .*  met.soilmoisture + ...
			1300.000  .* soil.peat_fraction + 2650.000 .* soil.mineral_fraction);
    
    
end

