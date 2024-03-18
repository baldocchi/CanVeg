function [ir]=fIR_RadTranCanopy_Matrix(leafang,ir,rad,soil,Sun,Shade, prm)
    % Dennis Baldocchi
    % Jan 27, 2021
    
        

% ------------------------------------------------------------
%                      IR_RadTranCanopy
% 
%     This subroutine computes the flux density of diffuse
%     radiation in the infrared waveband.  The Markov model is used
%     to compute the probability of beam penetration through clumped foliage.
%     
%       The algorithms of Norman (1979) are used.
%       
%         Norman, J.M. 1979. Modeling the complete crop canopy.
%         Modification of the Aerial Environment of Crops.
%         B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.
% -----------------------------------------------------------------
% 

% Need to add info on soil surface temperature and IR out later

% debug
% was having problem where IR.Shade of layer one was only 1/2 of the total
% flux. fixed it and the IR profile looks better, too. Had a kink before at
% the lower layer

tic

        SUP=zeros(prm.nn,prm.jktot);
        SDN=zeros(prm.nn,prm.jktot);
        temp=zeros(prm.nn,prm.jktot);
        temp1=ones(prm.nn,prm.jktot);
        
        % set upper boundary condition, in_dn was initialized as ones
        
        ir.ir_dn = temp1 .* ir.in;
%         compute IR radiative source flux as a function of
%         leaf temperature weighted according to
%         sunlit and shaded fractions
%         
%         source=ep*sigma*(laisun*tksun^4 + laish*tksh^4)

                  
          ir.IR_source_sun = rad.prob_beam(:,1:prm.jtot)  .* power(Sun.Tsfc,4.);
          ir.IR_source_shade = rad.prob_shade(:,1:prm.jtot)  .* power(Shade.Tsfc,4.);
          ir.IR_source = prm.epsigma * (ir.IR_source_sun + ir.IR_source_shade);
          
          
%         Intercepted IR that is radiated up

%            for j = prm.jtot:-1:1
%                jj=j+1;

          SUP(:,prm.jktot:-1:2) = ir.IR_source(:,prm.jtot:-1:1) .* (1. - leafang.integ_exp_diff(prm.jtot:-1:1)');
 
%        Intercepted IR that is radiated downward

         SDN(:,prm.jtot:-1:1) = ir.IR_source(:,prm.jtot:-1:1) .* (1. - leafang.integ_exp_diff(prm.jtot:-1:1)');
         
      
          % end
    
%         Downward IR radiation, sum of that from upper layer that is transmitted
%         and the downward source generated in the upper layer.
% 
%          REMEMBER LEVEL JJ IS AFFECTED BY temperature OF LAYER
%          ABOVE WHICH IS JJ + 1

%           %  see if I can decrement arrays

          for j = prm.jtot:-1:1
               jj=j+1;
         ir.ir_dn(:,j) = leafang.integ_exp_diff(j) .* ir.ir_dn(:,jj) + SDN(:,j);
          end
           
        
        % emiss_IR_soil = prm.epsigma .* power((soil.sfc_temperature + 273.16),4.);
        
        emiss_IR_soil = prm.epsigma .* power(soil.sfc_temperature,4.);
 
         SUP(:,1) = ir.ir_dn(:,1) .* (1. - prm.epsoil);
         
         ir.ir_up(:,1) = emiss_IR_soil + SUP(:,1);

           for j=2:prm.jktot
               jj=j-1;
          ir.ir_up(:,j) = leafang.integ_exp_diff(jj) .* ir.ir_up(:,jj) + SUP(:,j);
           end
        
           reflc_lay_IR = (1 - leafang.integ_exp_diff) .* (prm.epm1);% ep-1
         
          for K = 1:4
          for j=prm.jtot:-1:1
              jj=j+1;
                   
          ir.ir_dn(:,j) = leafang.integ_exp_diff(j) * ir.ir_dn(:,jj) +...
          ir.ir_up(:,j) .* reflc_lay_IR(j)' + SDN(:,j);
          end  % next j
          
          SUP(:,1) = ir.ir_dn(:,1) * (1 - prm.epsoil);
          ir.ir_up(:,1) = emiss_IR_soil + SUP(:,1);
          
      
          
          for j=2:prm.jktot
              jj=j-1;
          ir.ir_up(:,j) = reflc_lay_IR(jj)' .* ir.ir_dn(:,j) +...
              ir.ir_up(:,jj) .* leafang.integ_exp_diff(jj) + SUP(:,j);
          end % next j
                  
          end  %  next K   
        

      

% IR flux onto top and bottom of layers
for i=prm.jktot:-1:2
ii=i-1;
%ir.shade(:,i)=ir.ir_dn(:,i) + ir.ir_up(:,ii); 
ir.shade(:,ii)=ir.ir_dn(:,i) + ir.ir_up(:,ii); 


end


%ir.shade=ir.shade*prm.ep;    % IR in minus IR reflected (IR - (1-ep)IR)
        
      figure(6)
          
      plot(nanmean(ir.ir_up(:,1:prm.jtot),1),prm.sumlai,'+-','LineWidth', 1)
      ax=gca;
      set(ax, 'ydir','reverse');
      hold on
      plot(nanmean(ir.ir_dn(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 1)
      legend('up','down')
      xlabel('Radiation Flux Density')
      ylabel('Canopy Depth')
      title('IR flux density, W m-2');
      hold on;
      toc
end

      
      
      
      
      
      
      
      

