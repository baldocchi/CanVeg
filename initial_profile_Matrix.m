function [prof]=initial_profile_Matrix(met,prm)

% 3/31/2021
% initialize the profiles of Tair, eair, CO2 in the atmosphere and canopy

% make the layer heights in the canopy equal to number of layers/height
% above the canopy, assign the increments 0.1 m
% this is different from CanOak where I kept dz constant and changed dLAI
% for each layer. Here I keep dLAI constant and change dz for each layer. I
% think this is better for computing rad transfer and source sinks in the
% canopy

%     prm.dht_canopy=prm.veg_ht/prm.jtot;
%     prm.ht_atmos=prm.meas_ht-prm.veg_ht;
%     n_atmos_layers=50;
%     prm.dht_atmos=prm.ht_atmos/n_atmos_layers;
%     prm.nlayers_atmos=prm.jtot + floor(prm.ht_atmos/prm.dht_atmos);
% 
% 
% 
%     % consider look up tables for functional calls
% 
%     for i=1:prm.jtot
%     prm.zht(i)=i*prm.dht_canopy;
%     end
% 
%     j=0;
%     for i=prm.jtot+1:prm.nlayers_atmos
%         j=j+1;
%         prm.zht(i)=j* prm.dht_atmos + prm.veg_ht;
%     end



prof.nlayers=prm.nlayers_atmos; %prm.jtot + ht_atmos/dht_atmos; % number of atmospheric layers

prof.co2=ones(prm.nn,prof.nlayers);
prof.Tair_K=ones(prm.nn,prof.nlayers);
prof.Told_K=ones(prm.nn,prof.nlayers);
prof.eair_Pa=ones(prm.nn,prof.nlayers);
prof.eair_old_Pa=ones(prm.nn,prof.nlayers);
prof.wind=zeros(prm.nn,prm.jtot);

prof.co2=prof.co2 .* met.CO2;
prof.Tair_K=prof.Tair_K .* met.T_air_K;
prof.Told_K=prof.Told_K .* met.T_air_K;
prof.eair_Pa=prof.eair_Pa .* met.eair_Pa;
prof.eair_old_Pa=prof.eair_Pa;


prof.Tsfc=zeros(prm.nn,prm.jtot);
prof.H=zeros(prm.nn,prm.jtot);
prof.LE=zeros(prm.nn,prm.jtot);
prof.Ps=zeros(prm.nn,prm.jtot);

% height of canopy layers

prof.zht=prm.zht;
prof.delz=prm.delz;




end