
function [y]=psi(zoverl)

 y(zoverl==0)=0;
 
 y=4.7*zoverl(zoverl>0);

 x=(1-15*zoverl(zoverl<0))^.25;
 y=2*log((1+x)/2)+log((1+x*x)/2)-2*atan(x)+pi/2;
    
end


