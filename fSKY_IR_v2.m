function [RlongIn]=fSKY_IR_v2(met,sunrad,prm)
% New subroutine for longwave incoming energy
% 8/3/2023


% Choi, Minha, Jennifer M. Jacobs, and William P. Kustas. 2008.
% 'Assessment of clear and cloudy sky parameterizations for daily downwelling 
% longwave radiation over different land surfaces in Florida, USA'
% , Geophysical Research Letters, 35.


% Brunt Eq for clear sky

met.ea_mb=met.eair_Pa/100;

Rldc=(0.605 + 0.048 .* met.ea_mb .^0.5)*5.67e-8 .* met.T_air_K .^4;

% fraction of clouds

c= 1-met.rglobal./sunrad.Rgpotential;

% c is undefined at night so fix it based on the c value of clouds for that
% day

c(c<0)=0;
c(c==inf)=NaN;

cc=reshape(c,48,prm.ndays);

ccday=nanmean(cc,1);

for j=1:prm.ndays
tf=isnan(cc(:,j));
cc(tf,j)=ccday(j);
c=reshape(cc,48*prm.ndays,1);

end

% long wave with clear and clouds
RlongIn=Rldc .*(1-c) + c .* 5.67e-8 .* met.T_air_K .^4;
end

