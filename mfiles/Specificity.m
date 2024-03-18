function [tau]=Specificity(T)
% partioning coefficients for CO2 and O2

% Roupsard et al 1996 Ann Sci Forest 53: 243-256

%Henry Law coefficience mol l-1 bar-1

Kh_co2.a0=78.5e-3;
Kh_co2.a1=-2.89e-3;
Kh_co2.a2=54.7e-6;
Kh_co2.a3=-0.417e-6;

Kh_o2.a0=2.1e-3;
Kh_o2.a1=-57.1e-6;
Kh_o2.a2=1.024e-6;
Kh_o2.a3=-7.503e-9;

T2=T.*T;
T3 = T2 .* T;

Kh_co2.T=Kh_co2.a0 + Kh_co2.a1 * T + Kh_co2.a2 * T2 + Kh_co2.a3 *T3;

Kh_o2.T=Kh_o2.a0 + Kh_o2.a1 * T + Kh_o2.a2 * T2 + Kh_o2.a3 *T3;

% Specificity


% A = Vc - 0.5 Vo - Rd

% A = Vc (1- 0.5 Vo/Vc) -Rd

S=102;  % specificity at 25 C from Rouspard et al

tau=S*Kh_co2.T ./Kh_o2.T;

% Vc/Vo = C/O S KHco2/KHO2

% phi = Vo/Vc
end


