function [soil] = ThermalConductivity(soil, met)

% Sept 27, 2020
% merging Campbell Soil Physics in Basic with
% Bitelli Soil Physics in Python to Matlab


%soil thermal conductivity according to Campbell et al. Soil Sci. 158:307-313


%     ga = 0.088 ;               %deVries shape factor; assume same for all mineral soils
%     thermalConductivitysolid = 2.5;  %average thermal conductivity of soil minerals [W/mC]
%     atmPressure = met.P_kPa;                    %[kPa] 
%  
%     q = 7.25 * clay + 2.52;              %regression from soils in Campbell et al. 1994
%     xwo = 0.33 * clay + 0.078;           %regression from soils in Campbell et al. 1994
    
     %       // soil content

%         soil.clay_fraction = .3;      %//  Clay fraction   
% 		soil.peat_fraction = 0.129;    %//  SOM = a C; C = 7.5%, a = 1.72
% 		soil.pore_fraction = 0.687;    %// from alfalfa, 1 minus ratio bulk density 0.83 g cm-3/2.65 g cm-3, density of solids
% 		soil.mineral_fraction= 0.558; % // from bulk density asssuming density of solids is 2.65
    
    
        % thermal conductivity code from Campbell and Norman

   
        fw=1./(1+power((met.soilmoisture/0.15),-4));  %// terms for Stefan flow as water evaporates in the pores
        soil.K_air= 0.024 + 44100.*2.42e-5.*fw.*met.air_density_mole.* met.dest./met.P_Pa; 
     
 		k_fluid=soil.K_air + fw .*(soil.K_water-soil.K_air);
 
 		wt_air=2./(3.*(1+.2.*(soil.K_air./k_fluid-1))) + 1./(3.*(1+(1-2*.2).*(soil.K_air./k_fluid -1)));
      
 		wt_water=2./(3.*(1+.2.*(soil.K_water./k_fluid-1))) + 1./(3.*(1+(1-2.*.2).*(soil.K_water./k_fluid -1)));
       
 		wt_org=2./(3.*(1+.2.*(soil.K_org./k_fluid-1))) + 1./(3.*(1+(1-2*.2).*(soil.K_org./k_fluid -1)));  
        
        wt_mineral=2./(3*(1+.2.*(soil.K_mineral./k_fluid-1))) + 1 ./(3.*(1+(1-2*.2).*(soil.K_mineral./k_fluid -1)));
                           
         
        % compute thermal conductivty weighting by mineral, water, air and
        % organic fractions
    
                       
%          	% numerator
        K_soil_num=soil.mineral_fraction .* wt_mineral.*soil.K_mineral + soil.air_fraction .* wt_air.*soil.K_air +...
			soil.water_content_15cm .* wt_water.*soil.K_water + soil.peat_fraction .* wt_org.*soil.K_mineral;
% 
		soil.K_soil=K_soil_num./(soil.mineral_fraction .* wt_mineral + soil.air_fraction .* wt_air +...
			soil.water_content_15cm .* wt_water + soil.peat_fraction .* wt_org);
                       
                    
    
end

