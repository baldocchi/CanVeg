function [ y ] = TL( Z, input, turb )
%
%   Detailed explanation goes here


% /*  1/9/2021
%     ****************************************************************
%        This function gives the variation of T(L) with height
%     ****************************************************************
% 
% 
%         Adopt scaling values of Raupach that adjust Tl profiles for different
%         canopy strcutures
% 
%    TL = k(z-d)/(1.25^2 u*, z> h
%    TL constant below h

% in future double check Brunet BLM 2020 and Wohlfarhrt
% 
% */

%                             


      		
        if (Z >= 1.5* input.HH)
                
        
        y=0.4*(Z-input.DD)/(1.56 * input.ustar);
       
            %    y = A1* input.HH*0.25*turb.sigma_h/SIGMA(Z,input, turb);
        
        else
                
                    % TL ~ 0.25 h/u*
                
 				 y=0.25*input.HH/ input.ustar;

              
        
        end
       




end

