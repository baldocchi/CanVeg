function [Qin]=fQin_Matrix(quantum,nir,ir,prm)

% Nov 28, 2020



%          Available energy on leaves for evaporation.
%          Values are average of top and bottom levels of a layer.
%          The index refers to the layer.  So layer 3 is the average
%          of fluxes at level 3 and 4.  Level 1 is soil and level
%          j+1 is the top of the canopy. layer jtot is the top layer
%
%          Sunlit and shaded Absorbed Energy on the Top and Bottom of Leaves
% 
%          rad.sun_normal = rad.normal * rad.absorbed;
%                 
%          rad.sh_abs= (rad.dn_flux + rad.up_flux)* rad.absorbed ;
%          
%          rad.sun_abs= rad.sun_normal + rad.shade;


% rethink the balance of in minus out rather than in times absorbed..
          
     vis.sun_abs=quantum.sun_abs/4.6;     % convert umol m-2 s-1 PPFD to W m-2 in visible
     vis.sh_abs=quantum.sh_abs/4.6;
     
     % Qin = Rg(1-rho) + ep Lin
     
     % IR shade is on the upper and lower level of leaves
     
     %ir.shade=ir.shade*prm.ep;    % IR in minus IR reflected (IR - (1-ep)IR)
     
     % noticed in code ir.shade was multiplied by prm.ep twice, eg in
     % fIR_RadTranCanopy_MatrixV2 and here. removed prm.ep in the IR
     % subroutine
     
     Qin.sun_abs = nir.sun_abs + vis.sun_abs + ir.shade*prm.ep;   
     Qin.shade_abs = nir.sh_abs + vis.sh_abs + ir.shade*prm.ep;
      
%      
%       figure(4)
%       %clf
%       
%       plot(mean(Qin.sun_abs,1),1:prm.jtot,'LineWidth', 2)
%       hold on
%       plot(mean(Qin.shade_abs,1),1:prm.jtot,'LineWidth', 2)
%       
%       
%       legend('sun','shade')
%       xlabel('Radiation Flux Density')
%       ylabel('Canopy Depth')
%       title('Qin flux density, W m-2');
%       hold on;
%      
  
end
         
 


