function [ir]=fIR_RadTranCanopy_Matrix_Bonan(leafang,ir,rad,met,Sun,Shade, prm)
    % Dennis Baldocchi
    % Nov 4, 2020
    
        

% ------------------------------------------------------------
%                      IR_RadTranCanopy
% 
%     This subroutine computes the flux density of diffuse
%     radiation in the infrared waveband.  The Markov model is usedr
%     to compute the probability of beam penetration through clumped foliage.
%     
%       The algorithms of Norman (1979) are used.
%       
%         Norman, J.M. 1979. Modeling the complete crop canopy.
%         Modification of the Aerial Environment of Crops.
%         B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.

%       updated using Numerics in Bonan Climate Change and Terrestrial Ecosystem
%       Modeling
% -----------------------------------------------------------------
 

% may have some errors in layering. Double checking with Bonan
tic

        SUP=zeros(prm.nn,prm.jktot);
        SDN=zeros(prm.nn,prm.jktot);
        temp=zeros(prm.nn,prm.jktot);
        temp1=ones(prm.nn,prm.jktot);
        
        % set upper boundary condition, in_dn was initialized as ones
        % then assume IR at all layers is equal to the input from above for
        % first iteration
        
         ir.ir_dn = temp1 .* ir.in;   % level jktot
         ir.ir_up = temp1 .* ir.in;   % level jktot
         
%         compute IR radiative source flux as a function of
%         leaf temperature weighted according to
%         sunlit and shaded fractions
%         
%         source=ep*sigma*(laisun*tksun^4 + laish*tksh^4)
        
         ir.IR_source_sun = rad.prob_beam(1:prm.jtot)  .* power(Sun.Tsfc,4.);
          ir.IR_source_shade = rad.prob_shade(1:prm.jtot)  .* power(Shade.Tsfc,4.);
          ir.IR_source = prm.epsigma * (ir.IR_source_sun + ir.IR_source_shade);
          
           SUP(:,prm.jktot:-1:2) = ir.IR_source(:,prm.jtot:-1:1) .* (1. - leafang.integ_exp_diff);
 
%        Intercepted IR that is radiated downward

         SDN(:,prm.jtot:-1:1) = ir.IR_source(:,prm.jtot:-1:1) .* (1. - leafang.integ_exp_diff);

 
          
% Numerical Algorithms in Bonan from Norman

%         Lin,i = Lin,i+1 * (tau,diff, i+1) + Lout,i * ((1-tau,diff,i+1)*
%         emiss + IRsource,i+1 *(1-tau,diff,i+1)

%         Lout, i+1= Lout,i*(tua,diff,i+1)+ Lin,i+1((1-tau,diff,i+1)*emiss) +
%         IRsource,i+1 *(1-tau,diff,i+1)

%         SUP=SDN = IRsource,i+1 *(1-tau,diff,i+1)

% note SUP and SDN for i+1 is a function of source, T^4, i+1

          %SUP(:,prm.jktot:-1:2) = ir.IR_source(:,prm.jtot:-1:1) .* (1. - leafang.integ_exp_diff);
        
 
%        Intercepted IR that is radiated downward

        
         
         reflc_lay_IR = (1 - leafang.integ_exp_diff) * prm.epm1;
         
         
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%continue debug and revising code
      
      
   for k=1:4 % iterate as you need to know ups and downs to compute ups and downs
    
%         Downward IR radiation, sum of that from upper layer that is transmitted
%         and the downward source generated in the upper layer.


          for j = prm.jtot:-1:2
               jj=j+1;
       
%          Lin,i = Lin, i+1 * (tau,diff, i+1) + Lout, i * ((1-tau,diff,i+1)* emiss
%          + IRsource,i+1 *(1-tau,diff,i+1)     
           
       

          ir.ir_dn(:,j) = ir.ir_dn(:,jj) .* (leafang.integ_exp_diff)  + ...
          ir.ir_up(:,j) * reflc_lay_IR...
          + SDN(:,jj);
         
           end
           
        
 
        % need to solve for soil temperature. For debugging assume air
        % temperature or a constant lower boundary condition
        
         soil.sfc_temperature=met.Tsoil;
         emiss_IR_soil = prm.epsigma .* power((soil.sfc_temperature + 273.16),4.);
         
         
         % call fSoilEnergyBalance
         %soil.T_Kelvin;
         % emiss_IR_soil=soil.lout 
 
         SUP(:,1) = ir.ir_dn(:,2) .* (1. - prm.epsoil);
         
         ir.ir_up(:,1) = emiss_IR_soil + SUP(:,1);

           for j=2:prm.jktot
               jj=j-1;
               
               
%         Lout, i+1= Lout,i*(tau,diff,i+1)+ Lin,i+1((1-tau,diff,i+1)*emiss) +
%         IRsource,i+1 *(1-tau,diff,i+1)   

        
               
          ir.ir_up(:,j) = ir.ir_up(:,jj) .* leafang.integ_exp_diff...
          +   ir.ir_dn(:,j) .* reflc_lay_IR...
          + SUP(:,j);
           end
        
   end
         
%           for K = 1:4
%           for j=prm.jtot:-1:1
%               jj=j+1;
%             reflc_lay_IR = (1 - leafang.integ_exp_diff) * prm.epm1;
%         
%           ir.ir_dn(:,j) = leafang.integ_exp_diff * ir.ir_dn(:,jj) +...
%           ir.ir_up(:,j) * reflc_lay_IR + SDN(:,j);
%           end  % next j
%           
%           SUP(:,1) = ir.ir_dn(:,1) * (1 - prm.epsoil);
%           ir.ir_up(:,1) = emiss_IR_soil + SUP(:,1);
%           
%           reflc_lay_IR = (1 - leafang.integ_exp_diff) * (prm.epm1);
%           
%           for j=2:prm.jktot
%               jj=j-1;
%           ir.ir_up(:,j) = reflc_lay_IR * ir.ir_dn(:,j) +...
%               ir.ir_up(:,jj) * leafang.integ_exp_diff + SUP(:,j);
%           end % next j
%                   
%           end  %  next K   
%         

      
        
      figure(6)
      clf
      
      plot(mean(ir.ir_up,1),1:prm.jktot,'LineWidth', 1)
      hold on
      plot(mean(ir.ir_dn,1),1:prm.jktot,'LineWidth', 1)
      
      
      legend('up','down')
      xlabel('Radiation Flux Density')
      ylabel('Canopy Depth')
      title('IR flux density, W m-2');
      toc
end

      
      
      
      
      
      
      
      

