function [ y ] = TL( Z, input, turb )
%
%   Detailed explanation goes here


% /*
%     ****************************************************************
%        This function gives the variation of T(L) with height
%     ****************************************************************
% 
% 
%         Adopt scaling values of Massman and Weil that adjust Tl profiles for different
%         canopy strcutures
% 
%   u* Tl/h = A1 * (z-d)/h;  z > h
% 
%   u* Tl/h = A1 (1-d/h) gamma_3/ (sigmaw(z)/u*); z <= h
% 
%   A1 = A2 (sigmaw(z)/u* gamma_3)^1/2 (1-d/h)^-1/2
% 
%   gamma_3 = sigmaw(h)/u*
% 
% */

%                 // factor of 2 comes from (1 - d/h)^-1/2; (1-0.75)^-1/2  
                


            A1=0.6*sqrt(SIGMA(Z,input, turb)/input.ustar*turb.sigma_h)* 2.;

		
        if (Z <= input.HH)
                

        
%         // The factor 0.25 = 1 - d/h = 1 - 0.75
        
                y = A1* input.HH*0.25*turb.sigma_h/SIGMA(Z,input, turb);
        
                else
                
                
% 				 // u* Tl/h = A1 * (z-d)/h;  z > h	

                y = A1 * (Z-input.DD)/input.ustar;
        
        end
       




end

