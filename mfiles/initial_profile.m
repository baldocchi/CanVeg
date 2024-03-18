function [prof]=initial_profile(met,prm)

% initialize the profiles of Tair, eair, CO2 in the atmosphere and canopy

% make the layer heights in the canopy equal to number of layers/height
% above the canopy, assign the increments 0.1 m
% this is different from CanOak where I kept dz constant and changed dLAI
% for each layer. Here I keep dLAI constant and change dz for each layer. I
% think this is better for computing rad transfer and source sinks in the
% canopy


ht_atmos=prm.meas_ht-prm.veg_ht;
dht_canopy=prm.veg_ht/prm.jtot;

dht_atmos=0.1;
prof.nlayers=prm.jtot + ht_atmos/dht_atmos; % number of atmospheric layers


prof.zht=zeros(prof.nlayers,1);
prof.co2=zeros(prof.nlayers,1);
prof.Tair_K=zeros(prof.nlayers,1);
prof.Told_K=zeros(prof.nlayers,1);
prof.eair_Pa=zeros(prof.nlayers,1);

prof.co2(:,1)=met.CO2(1,1);
prof.Tair_K(:,1)=met.T_air_K(1,1);
prof.Told_K(:,1)=met.T_air_K(1,1);
prof.eair_Pa(:,1)=met.eair_Pa(1,1);

% height of canopy layers

for i=1:prm.jtot
prof.zht(i)=i*dht_canopy;
end

% height of atmospheric layers

j=0;
for i=prm.jtot+1:prof.nlayers
    j=j+1;
   prof.zht(i)=j* dht_atmos + prm.veg_ht;
end

end