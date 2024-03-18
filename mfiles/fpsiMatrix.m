
function [y]=fpsiMatrix(zoverl)

 y=zeros(length(zoverl),1);
 
y(zoverl==0)=0;
 
 y(zoverl>0)=4.7 .*zoverl(zoverl>0);

 x=(1-15*zoverl(zoverl<0)) .^ 0.25;
 y(zoverl<0)=2 .*log((1+x)/2)+log((1+x .* x)/2)-2 .*atan(x)+pi/2;
    
end


