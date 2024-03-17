function [fT]=fdESdT(T)

% Slope of the saturation vapor pressure-temperature curve
% 	
% Thermodynamic relationship from HESS    (Pa K-1)
%param.Rstar = 8.3144; 

fT = fES(T) .* fLambda(T) * 18. ./ (8.3144 .* T .* T * 1000.);