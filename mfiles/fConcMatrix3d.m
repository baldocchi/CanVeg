        
function [ccnc]=fConcMatrix(Source,soilflux,delz, Dij,met, cref,prm, factor)
 % Subroutine to compute scalar concentrations from source
 % estimates and the Lagrangian dispersion matrix

 
 % Dec 9, 2020
 % this sub has been modified to include new variable soil energy and 
 % mas fluxes
 
 ustar=met.ustar;

% met.zL=0;  %Monin Obukhov, compute later, zero for now

 
       % // Compute concentration profiles from Dispersion matrix
       
        ustar_ref = 1;
       
        ustfact = ustar_ref./ustar;   %   ustar_ref /ustar;  
 
 % loop through ustfact for each met run
 
 cc=zeros(prm.jtot,1);
 disperzl=zeros(prm.nlayers_atmos,prm.jtot);
 soilbnd=zeros(prm.nn,1);
 cncc=zeros(prm.nn,prm.nlayers_atmos);
 
 disper=zeros(prm.nn,prm.nlayers_atmos,prm.jtot);
 
                  
 sumcc=zeros(prm.nn,prm.nlayers_atmos); 
 

 
% loop k otherwise I need a 3D disperision matrix Dij(k,i,j) if I am to perform matrix multiplications
% something to think about in the future fConMatrix3d...

% initialize matrices

 
 
%  CC is the differential concentration (Ci-Cref)
% 
% 
%         Ci-Cref = SUM (Dij S DELZ), units mg m-3 or mole m-3
% 
%          S = dfluxdz/DELZ
% 
%         note delz values cancel
% 
%         scale dispersion matrix according to friction velocity
 
                   disper(:,:,:)=ustfact .* Dij;
 
  
                  disperzl(met.zL < 0) = disper(met.zL < 0) .* (0.973*-0.7182)./(met.zL(met.zL < 0) -0.7182);
                  
                  disperzl(met.zL >=0)=disper(met.zL >=0) * (-0.31 * met.zL(met.zL >=0) + 1.00);
                  
 
                    sumcc=sum(delz'.* disperzl.*Source);  
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  for k=1:length(ustfact)
%
%         CC is the differential concentration (Ci-Cref)
% 
% 
%         Ci-Cref = SUM (Dij S DELZ), units mg m-3 or mole m-3
% 
%          S = dfluxdz/DELZ
% 
%         note delz values cancel
% 
%         scale dispersion matrix according to friction velocity
 
                   disper=ustfact(k) .* Dij;
 
%   updated Dispersion matrix alfalfa
% 
                  if (met.zL(k) < 0)
                  disperzl = disper * (0.973*-0.7182)/(met.zL(k) -0.7182);
                  else
                  disperzl=disper * (-0.31 * met.zL(k) + 1.00);
                  end
             
      
       
for i=1:prm.nlayers_atmos

     

% /*
%         CC is the differential concentration (Ci-Cref)
% 
% 
%         Ci-Cref = SUM (Dij S DELZ), units mg m-3 or mole m-3
% 
%          S = dfluxdz/DELZ
% 
%         note delz values cancel
% 
%         scale dispersion matrix according to friction velocity
% */

  %      disper = ustfact(k) .* Dij;              %      // units s/m

       % // scale dispersion matrix according to Z/L

      
       %  sumccprod= delz' .* disperzl .* Source(k,:);  %  needed to transpose delz for matrix
     
       sumcc(i)=sum(delz(1:prm.jtot)'.* disperzl(i,1:prm.jtot).*Source(k,1:prm.jtot));  
      
     
      %  // scale dispersion matrix according to Z/L


        dispersoil = ustfact(k) .* Dij(i,1);
        
       % disperzlsoil=dispersoil;
                
%                       
                 if(met.zL < 0)
                 disperzlsoil = dispersoil .* (0.97*-0.7182)./(met.zL -0.7182);
                 else
                 disperzlsoil=dispersoil .* (-0.31 .* met.zL + 1.00);
                 end
                
       % // add soil flux to the lowest boundary condition

        soilbnd(k)=soilflux(k) .*disperzlsoil(k) ./factor(k);  

        cc(i)=sumcc(i) ./factor(k)+soilbnd(k);
        

end % // next i z level

% // factor to adjust Dij with alternative u* values
                                             
 %       //  Compute scalar profile below reference
 
       
        for i=1:prm.nlayers_atmos 
        ccnc(k,i) = cc(i) + cref(k) - cc(prm.nlayers_atmos);
        end


 end %//next k, half-hour met data
 


end


