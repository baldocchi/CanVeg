function [soil] = fSetSoil(met)

% 
%   Routines, algorithms and parameters for soil moisture were from 
% 
%   Campbell, G.S. 1985. Soil physics with basic. Elsevier
% 
%   updated to algorithms in Campbell and Norman and derived from Campbell et al 1994 Soil Science
% 
% 
%   Need to adjust for clay and organic fractions.. Need to adjust heat capacity and conductivity for peat

		

%        Soil water content 

        soil.water_content_15cm = met.soilmoisture;    % //  measured at 10 cmwater content of soil m3 m-3  


 %       // Water content of litter. Values ranged between 0.02 and 0.126

        soil.water_content_litter = .0;   %// assumed constant but needs to vary


  %       // soil content

        soil.clay_fraction = .3;      %//  Clay fraction   
		soil.peat_fraction = 0.129;    %//  SOM = a C; C = 7.5%, a = 1.72
		soil.pore_fraction = 0.687;    %// from alfalfa, 1 minus ratio bulk density 0.83 g cm-3/2.65 g cm-3, density of solids
		soil.mineral_fraction= 0.558; % // from bulk density asssuming density of solids is 2.65

		soil.air_fraction = soil.pore_fraction - soil.water_content_15cm;

		soil.Cp_water= 4180;  % // J kg-1 K-1, heat capacity
		soil.Cp_air =  1065;
		soil.Cp_org = 1920; 
		soil.Cp_mineral = 870;

		soil.K_mineral= 2.5;  % // W m-1 K-1, thermal conductivity
		soil.K_org= 0.8;
		soil.K_water= 0.25;
		

		% thermal conductivity code from Campbell and Norman

		fw=1./(1+power((soil.water_content_15cm/0.15),-4));  %// terms for Stefan flow as water evaporates in the pores
		
		%desdt=DESDT(input.ta+273.15);
        
        
		soil.K_air= 0.024 + 44100*2.42e-5*fw*met.air_density_mole*met.desdt/met.press_Pa;

		k_fluid=soil.K_air + fw *(soil.K_water-soil.K_air);

		wt_air=2/(3*(1+.2*(soil.K_air/k_fluid-1))) + 1/(3*(1+(1-2*.2)*(soil.K_air/k_fluid -1)));
		wt_water=2/(3*(1+.2*(soil.K_water/k_fluid-1))) + 1/(3*(1+(1-2*.2)*(soil.K_water/k_fluid -1)));
		wt_mineral=2/(3*(1+.2*(soil.K_mineral/k_fluid-1))) + 1/(3*(1+(1-2*.2)*(soil.K_mineral/k_fluid -1)));
		wt_org=2/(3*(1+.2*(soil.K_org/k_fluid-1))) + 1/(3*(1+(1-2*.2)*(soil.K_org/k_fluid -1)));

		Cp_soil_num= ( met.air_density * soil.Cp_air * soil.air_fraction + 1000.000 * soil.Cp_water * soil.water_content_15cm +... 
			1300.000 * soil.Cp_org * soil.peat_fraction + 2650.000 * soil.Cp_mineral * soil.mineral_fraction);

		soil.Cp_soil=Cp_soil_num/( met.air_density *  soil.air_fraction + 1000.000 *  soil.water_content_15cm + ...
			1300.000  * soil.peat_fraction + 2650.000 * soil.mineral_fraction);
		
		soil.K_soil_num=soil.mineral_fraction * wt_mineral*K_mineral + soil.air_fraction * wt_air*K_air +...
			soil.water_content_15cm * wt_water*K_water + soil.peat_fraction * wt_org*K_mineral;

		soil.K_soil=K_soil_num/(soil.mineral_fraction * wt_mineral + soil.air_fraction * wt_air +...
			soil.water_content_15cm * wt_water + soil.peat_fraction * wt_org);
       
        soil.dt = 20.;        %  // Time step in seconds
 
        soil.mtime = floor(3600 /soil.dt);   %// time steps per hour


           %     //  Assign soil layers and initial temperatures
           %     //  and compute layer heat capacities and conductivities

  
        for I=0:n_soil
        
        IP1=I+1;    
        IM1=I-1;
        soil.z_soil(IP1) = soil.z_soil(I) + .005 * 1.5^IM1;
        soil.T_soil(I) = met.Tsoil;

            %   // assign bulk densities for litter and soil

        if (soil.z_soil(IP1) < z_litter)
        soil.bulk_density(IP1) = .074 ;  % // litter  
        else
        soil.bulk_density(IP1) = 0.83;   % // soil  bulk density for the alfalfa, g cm-3
        end
        
        end % }  //  next I  

        for I=1:10
        
%                 // Heat capacity and conductivity.
% 
%        		          
% 			// assuming in the numeric code it is bulk density times Cp as values in code have
% 			// Campbell and Norman have rhos Cs = rhoa Cpa + .. correctly in the book, pg 119.
% 
% 			// use weighted Cp check units kg m-3 times J kg-1 K-1

        soil.cp_soil(I) = soil.Cp_soil * (soil.z_soil(I + 1) - soil.z_soil(I - 1)) / (2. * soil.dt);


% 		// adopt new equations from Campbell and Norman and Campbell 1994, after the basic book was written
%         
%           //      C1 = .65 - .78 * soil.bulk_density[I] + .6 * soil.bulk_density[I] * soil.bulk_density[I];
%           //      C2 = 1.06 * soil.bulk_density[I];  // corrected according to Campbell notes
%            //     C3= 1. + 2.6 / sqrt(soil.mineral_fraction);     
%           //      C4 = .03 + .1 * soil.bulk_density[I] * soil.bulk_density[I];
% 
% 			// soil conductivity needs debuting ?? or are units for bulk of soil kg m-3
% 			//soil.k_conductivity_soil[I] = (C1 + C2 * soil.water_content_15cm - (C1 - C4) * exp(-pow((C3 * soil.water_content_15cm), 4.))) / (soil.z_soil[I + 1] - soil.z_soil[I]);
%         
	     soil.k_conductivity_soil(I) = soil.K_soil / (soil.z_soil(I + 1) - soil.z_soil(I));
     
						
        end %// NEXT I 
        

end

