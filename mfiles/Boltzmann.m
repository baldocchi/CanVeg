%//////////////////////////////////////////////
 
        %clear all;
        param.rugc = 8.314;              % J mole-1 K-1 
        rgc1000 = 8314;            % gas constant times 1000.


%          Consts for Photosynthesis model and kinetic equations.
%          for Vcmax and Jmax.  Taken from Harley and Baldocchi (1995, PCE)
 

        param.hkin = 200000.0;    % 200000 enthalpy term, J mol-1
        param.skin = 710.0;       % entropy term, J K-1 mol-1
        ejm = 55000.0;      % activation energy for electron transport, J mol-1
        evc = 55000.0;      % activation energy for carboxylation, J mol-1
    

%         Enzyme constants & partial pressure of O2 and CO2
%         Michaelis-Menten K values. From survey of literature.


         kc25 = 274.6;   % kinetic coef for CO2 at 25 C, microbars  
        
         ko25 = 419.8;   % kinetic coef for O2 at 25C,  millibars 
         o2 = 210.0;     % oxygen concentration  mmol mol-1  
         
        ekc = 80500.0;     % Activation energy for K of CO2; J mol-1  
        eko = 14500.0;     % Activation energy for K of O2, J mol-1
        erd = 38000.0;     % activation energy for dark respiration, eg Q10=2 
        ektau = -29000.0;  % J mol-1 (Jordan and Ogren, 1984)
        tk_25 = 298.16;    % absolute temperature at 25 C
        toptvc = 303;%311;    % optimum temperature for maximum carboxylation
        toptjm = 303% 311.0;    % optimum temperature for maximum electron transport
        
         vcopt = 127;   % carboxylation rate at optimal temperature, umol m-2 s-1 
         jmopt = 277;   % electron transport rate at optimal temperature, umol m-2 s-1 
         rd25 = .5;     % dark respiration at 25 C, rd25= 0.34 umol m-2 s-1 
      

%         rt = param.rugc * tlk;               %  product of universal gas constant and abs temperature
% 
%         tprime25 = tlk - tk_25;       % temperature difference

       
       
%      KC and KO are solely a function of the Arrhenius Eq.

for i=275:320
    
    tlki=i;
      rt = param.rugc * tlki;               %  product of universal gas constant and abs temperature

        tprime25 = tlki - tk_25;   
         ttemp = exp((param.skin * tlki - param.hkin) / rt) + 1.0; %  denominator term

        jmax(i-274) = fTBOLTZ(jmopt, ejm, toptjm, tlki);
        vcmax(i-274) = fTBOLTZ(vcopt, evc, toptvc, tlki);
        
        
end

        ratio=vcmax/fTBOLTZ(vcopt, evc, toptvc, 298)

       figure(1)
       clf;
       plot(275:320,ratio,'LineWidth',3)
       
      
       
       figure(2)
       clf;
       plot(275:320,jmax,'LineWidth',3)

