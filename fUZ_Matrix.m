function [wnd]=fUZ_Matrix(met,prm)

% 3/15/2021

 LL=(prm.meas_ht -prm.dht)./met.zL;  % Monin Obukhov Length scale

% U(Z) inside the canopy 
  
  % use u(z)=u(h)exp(alpha(z/h-1)) following Cionco
  
  % need to adjust for z/L
  
  % Psi =f(met.zL)
  % log wind profile         
  %velocity(i)=(fvelocity/k)*(log(abs((heightarray(i)-d)/roughness))-psi(heightarray(i)-d)/LL);



  if LL == inf
  UH=met.wind .* log((prm.veg_ht-prm.dht)./prm.z0)./log((prm.meas_ht -prm.dht)/prm.z0);
  else
 
  tst= fpsiMatrix((prm.meas_ht -prm.dht)./LL);
  tst1=fpsiMatrix((prm.veg_ht -prm.dht)./LL);
      
  UH=met.wind .* (log((prm.veg_ht-prm.dht)./prm.z0)-tst1)./(log((prm.meas_ht -prm.dht)/prm.z0) - tst);
  end
  
   
  wndexp= exp(prm.extinct *(prm.zht(1:prm.jtot)./prm.veg_ht -1));
  wnd=UH .* wndexp';

end

  
       


     