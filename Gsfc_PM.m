function [ GsPM ] = Gsfc_PM(met,Can)  
% Canopy conductance Penman Monteith 
%  12/26/2020 
% D. Baldocchi

% input data
% data.LE_gf,data.RNET, data.G, vpd, data.TA, data.ustar, data.ubar );

% need to create input.Avail

Cp = 1005;
Press=met.P_Pa; % Pa

Tk=met.T_air_K;



psych=Cp * Press ./(fLambda(Tk)*0.622); % psychrometric constant

s=fdESdT(Tk);  % slope of saturation vapor pressure-temperature

Ram=met.wind ./(met.ustar .* met.ustar);  % boundary layer resistance for momentum

Sc=15.1/24.9;  % Schmidt number
Pr=15.1/22.2;  % Prandtl number

Rb= 2./(.4*met.ustar) * (Sc/Pr)^.6666; % quasi laminar boundary layer resistance for vapor

RaH=Ram+Rb;
GaH=1 ./RaH;  % boundary layer conductance for heat, m/s

vpd=met.vpd_Pa;  % vapor pressure deficit, Pa



% canopy surface conductance, m/s, inverted Penman Monteith Eq
% GsPM=(psych .*GaH .* input.LE_gf) ./(s .*(input.RNET-input.G)+input.rhoa .*Cp .*vpd .* GaH-input.LE_gf .*(s+psych));

GsPM=(psych .*GaH .* Can.LE) ./(s .*(Can.Avail)+met.air_density .*Cp .*vpd .* GaH-Can.LE .*(s+psych));

end



  function [fT]=DESDT(T)

% Slope of the saturation vapor pressure-temperature curve
% 	
% Thermodynamic relationship from HESS    (Pa K-1)
% f(Tk)
param.Rstar = 8.3144; 

ess=100. * exp(54.8781919 - 6790.4985 ./ T - 5.02808 .* log(T));

fT = ess .* LAMBDA(T) * 18. ./ (param.Rstar .* T .* T * 1000.);
  end 
 

function [fT]=ES(T) 
% 	Saturation vapor pressure,  Pa, f(TC)

% fT = 100. * exp(54.8781919 - 6790.4985 ./ T - 5.02808 .* log(T));

fT=613.75*exp(17.502*T./(240.97+T));

end


function [fT]=LAMBDA(T)

% Latent heat of vaporization, J kg-1
% f(Tk)

fT = 3149000 - 2370 * T;
end

