function [soil] = SoilRespiration(Ac,Tsoil,soilmoisture,veght,Rd,prm)
% 12/4/2020

% Use statistical model derived from ACE analysis on alfalfa

VegType=prm.Veg;

 switch VegType
        case 'Alfalfa'

% gpp, Tsoil, height, soil moisture; Phi/Resp;    % double check if ACE
% alfalfa was with Tk or Tc
 
x=[Ac+Rd,Tsoil,soilmoisture,veght];


b1=[-1.51316616076344;0.673139978230031;-59.2947930385706;-3.33294857624960];
b2=[1.38034307225825;3.73823712105636;59.9066980644239;1.36701108005293];
b3=[2.14475255616910;19.9136298988773;0.000987230986585085;0.0569682453563841]; 


phi =b1' + b2' .* x./(b3'+x);

phisum=sum(phi,2);

%   Detailed explanation goes here

b0_0=4.846927939437475;
b0_1=2.22166756601498;
b0_2=-0.0281417120818586;

soil.Respiration =b0_0+b0_1.*phisum + b0_2.*phisum .* phisum - Rd;   % Reco-Rd = Rsoil


     case 'Tule'

% if Tule

watertable=soilmoisture;  % the model was based on water table in meter and measured in cm

TC=Tsoil-273.15;  % ACE for Tules was calculated with Centigrade

x=[Ac+Rd,TC,watertable];   % for tule soilmoisture is water table, Tsoil is water temperature


b1=[1.62630740309143;-0.0972080624769812;-3.34991935356261];
b2=[0.667332480459099;11306589.2650554;0.00123789249208206];
b3=[1.02973143355825;101496007.349346;1.13804193403370];


phi =b1' + b2' .* x./(b3'+x);

% phiest(i,:)=b1(i)+b2(i) .*x(i,:) ./(b3(i)+ x(i,:));

phisum=sum(phi,2);


b0_0=2.90604140359958;
b0_1=1.54051093524537;
b0_2= 0.329925286533385;

%soil.Respiration

soil.Respiration =b0_0+b0_1.*phisum + b0_2.*phisum .* phisum - Rd;   % Reco-Rd = Rsoil

     case 'Savanna'

         % Tang and Baldocchi, June

         % open space
         soil.RespOpen=-0.12+0.029 .* (Tsoil-273.15);

         % under tree
          soil.RespTree=5.33+0.040 .* (Tsoil-273.15);

          %weighted

          soil.Respiration=0.4*soil.RespOpen + 0.6*soil.RespTree;

     case 'DeciduousForest'


     otherwise

 end
end

