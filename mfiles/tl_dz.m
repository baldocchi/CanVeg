function [ y ] = tl_dz( Z, input)
%
%   1/9/2021


% /*
%     ****************************************************************
%        This function gives the variation of T(L) with height
%     ****************************************************************
                
y=zeros(1,length(Z));

                         
        
        y=0.4*(Z-input.DD)/(1.56 * input.ustar);
       
               
                    % TL ~ 0.25 h/u*
       y(Z <= 1.5*input.HH)=0.25*input.HH/ input.ustar;

              
       

end

