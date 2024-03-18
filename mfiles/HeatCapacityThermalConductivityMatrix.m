function [soil] = HeatCapacityThermalConductivityMatrix(met,soil)

% Nov 17, 2020


% compute heat capacity of the soil and thermal conductivity

% merging Campbell Soil Physics in Basic with
% Bitelli Soil Physics in Python to Matlab


%soil thermal conductivity according to Campbell et al. Soil Sci. 158:307-313

% 
% 
% 		// adopt new equations from Campbell and Norman and Campbell 1994, after the basic book was written
%         
%           //      C1 = .65 - .78 * soil.bulk_density[I] + .6 * soil.bulk_density[I] * soil.bulk_density[I];
%           //      C2 = 1.06 * soil.bulk_density[I];  // corrected according to Campbell notes
%            //     C3= 1. + 2.6 / sqrt(soil.mineral_fraction);     
%           //      C4 = .03 + .1 * soil.bulk_density[I] * soil.bulk_density[I];
% 
% 			// soil conductivity needs debugging ? or are units for bulk of soil kg m-3
% 			//soil.k_conductivity_soil[I] = (C1 + C2 *
% 			soil.water_content_15cm - (C1 - C4) * ...
%                 exp(-pow((C3 * soil.water_content_15cm), 4.))) / (soil.z_soil[I + 1] - soil.z_soil[I]);
         
  

    
        % thermal conductivity code from Campbell and Norman

   
        fw=1./(1+power((met.soilmoisture/0.15),-4));  %// terms for Stefan flow as water evaporates in the pores
        soil.K_air= 0.024 + 44100.*2.42e-5.*fw.*met.air_density_mole.* met.dest./met.P_Pa; 
     
 		k_fluid=soil.K_air + fw .*(soil.K_water-soil.K_air);
 
 		wt_air=2./(3.*(1+.2.*(soil.K_air./k_fluid-1))) + 1./(3.*(1+(1-2*.2).*(soil.K_air./k_fluid -1)));
      
 		wt_water=2./(3.*(1+.2.*(soil.K_water./k_fluid-1))) + 1./(3.*(1+(1-2.*.2).*(soil.K_water./k_fluid -1)));
       
 		wt_org=2./(3.*(1+.2.*(soil.K_org./k_fluid-1))) + 1./(3.*(1+(1-2*.2).*(soil.K_org./k_fluid -1)));  
        
        wt_mineral=2./(3*(1+.2.*(soil.K_mineral./k_fluid-1))) + 1 ./(3.*(1+(1-2*.2).*(soil.K_mineral./k_fluid -1)));
                           
         
        % compute heat capacity and thermal conductivty weighting by mineral, water, air and
        % organic fractions
        
        % rho_s Cs dT/dt = d(k dT/dz)/dz
        
               
        
       
        % soil density rho (kg m-3) *  Heat Capacity, J kg-1 K-1 -> J m-3
        % K-1
        
        % factors of 1000 is kg H2O m-3, or density of mineral soil (2650) and peat (1300) in
        % terms of kg m-3
        
        soil.Cp_soil= (met.air_density .*  soil.air_fraction .* soil.Cp_air +  1000.000 .* soil.Cp_water .* met.soilmoisture +... 
			1300.000 .* soil.Cp_org .* soil.peat_fraction + 2650.000 .* soil.Cp_mineral .* soil.mineral_fraction);
    
%         % from Python Code and Simpler Campel function
        % soil.Cp_soil =  2.4e6*soil.bulkdensity + 4.18e6 .*met.soilmoisture; 
        
    
    
          % thermal conductivity, W m-1  K-1
          
        K_soil_num=soil.mineral_fraction .* wt_mineral.*soil.K_mineral + soil.air_fraction .* wt_air.*soil.K_air +...
			soil.water_content_15cm .* wt_water.*soil.K_water + soil.peat_fraction .* wt_org.*soil.K_mineral;
 
 		soil.K_soil=K_soil_num./(soil.mineral_fraction .* wt_mineral + soil.air_fraction .* wt_air +...
 			soil.water_content_15cm .* wt_water + soil.peat_fraction .* wt_org);
                       
                            
    
end

