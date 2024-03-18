        
function [ccnc]=Conc(Source, Dij,ustar, cref,prm)
 % Subroutine to compute scalar concentrations from source
 % estimates and the Lagrangian dispersion matrix

 cc=zeros(prm.jtot);
 cncc=zeros(prm.jtot3);

 zL=0;  %Monin Obukhov, compute later, zero for now
 soilflux=0;  % zero for now, later add in soil fluxes
 
       % // Compute concentration profiles from Dispersion matrix
       
        ustar_ref = 1;
       
 ustfact = ustar_ref/ustar;   %   ustar_ref /ustar;  
 
 
       
for i=1:prm.jtot3

        sumcc=zeros(prm.jtot3,1);

    for j=1:prm.jtot

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

        disper = ustfact .* Dij;              %      // units s/m

       % // scale dispersion matrix according to Z/L

              
       
               % // updated Dispersion matrix (Dec, 2002)

                if(zL < 0)
                disperzl = disper * (0.679* zL -0.5455)/(zL -0.5462);
                else
                disperzl=disper;
                end

                if (i < 50)
                delz=0.06; % for 0 to 50; 0.1 51-70
                else
                delz=0.1;
                end
        
        sumcc(i) = sumcc(i-1) + delz * disperzl * Source(j);

        
     end % // next j 


      %  // scale dispersion matrix according to Z/L


        disper = ustfact .* Dij(:,1);
                
                      
                if(met.zl < 0)
                disperzl = disper * (0.679* met.zl -0.5455)/(met.zl -0.5462);
                else
                disperzl=disper;
                end

                
       % // add soil flux to the lowest boundary condition

        soilbnd=soilflux*disperzl/factor;

        cc(i)=sumcc(i)/factor+soilbnd;

end % // next i

% // factor to adjust Dij with alternative u* values
                                              %  // Note that disperion matrix was computed using u* = 0.405
 %       //  Compute scalar profile below reference
 
       
        for i=1:jtot3 
        cncc(i) = cc(i) + cref - cc(izref);
        end



end




        
  

