function [Rsoil] = Soil_Sfc_Res(wg)
             

%            // Camillo and Gurney model for soil resistance 
%            // Rsoil= 4104 (ws-wg)-805, ws=.395, wg=0 
%            // ws= 0.395
% 
% 	       // wg is at 10 cm, use a linear interpolation to the surface, top cm, mean between 0 and 2 cm

	       wg0 = 1 .* wg/10;

%        //    y=4104.* (0.395-wg0)-805.;
% 
% 	   // model of Kondo et al 1990, JAM for a 2 cm thick soil

	   Rsoil = 3e10 .* power((0.395-wg0), 16.6);

	
    end

