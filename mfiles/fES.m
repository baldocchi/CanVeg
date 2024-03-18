
function [fT]=fES(T) 
% 	Saturation vapor pressure,  Pa, T Kelvin

%fT = 100. * exp(54.8781919 - 6790.4985 ./ T - 5.02808 .* log(T));

Tc=T-273.15;
fT=613.75 * exp(17.502 .*Tc./(240.97 + Tc));



