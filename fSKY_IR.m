 function [y]=fSKY_IR(T,sunrad,prm)
       
%  Infrared radiation from sky, W m-2, using algorithm from Norman
       
       % alternatively
       % emissivity of the atmosphere from Prata, 1996 QJRMS
       % this function performs matrix calc for the array of met inputs
       %[ep_air]=EMISSIVITY_ATMOS(met.eair_hPa, met.Tair_K); 
       
       % 8/2/2023
       
%        double SKY_IR (double T)
%         {
%        
%         // Infrared radiation from sky, W m-2, using algorithm from Norman
%         
%         double y;
%  Idso and Jackson, 1969 JGR
 
%         y = sigma * pow(T,4.) * ((1. - .261 * exp(-.000777 * pow((273.16 - T), 2.))) * solar.ratrad + 1 - solar.ratrad);
%         
%         return y;
%         }
%        
       % Linair= ep_air .* 5.67e-8 .* met.Tair_K .^4;
       
       
        %sunrad.ratrad=0.8;
        
		product= ( 1. - .261 * exp(-.000777 * (273.16 - T).^ 2. ).* sunrad.ratrad + 1. - sunrad.ratrad);
      
    	y = product .* prm.sigma .* T.^4.; 

end

