function [ir]=fIR_RadTranCanopy_MatrixV2(leafang,ir,rad,soil,Sun,Shade, prm)
    % Dennis Baldocchi
    % Feb 23, 2021
    
        

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

%     And adaption algorithms from Bonan, Climate Change and Terrestrial Ecosystem
%     Modeling of Norman's model
% -----------------------------------------------------------------
%   assume IR transmission by leaf is zero and reflectance is 1- emissivity


%fixed problem



tic

        ir.ir_up=zeros(prm.nn,prm.jktot);
        ir.ir_dn=zeros(prm.nn,prm.jktot);
        temp=zeros(prm.nn,prm.jktot);
        temp1=ones(prm.nn,prm.jktot);
        
        % set upper boundary condition, in_dn was initialized as ones
        
          if ir.flag ==1 
        
          ir.ir_dn = temp1 .* ir.in;
          ir.ir_up = temp1 .* ir.in;
          else
              
              ir.ir_dn(:,prm.jktot)=ir.in;
              ir.ir_up(:,prm.jktot)=ir.in;
          end
          
%         compute IR radiative source flux as a function of
%         leaf temperature weighted according to
%         sunlit and shaded fractions
%         
%         source=ep*sigma*(laisun*tksun^4 + laish*tksh^4)

           
           ir.IR_source_sun = rad.prob_beam(:,1:prm.jtot)  .* power(Sun.Tsfc,4.);
           ir.IR_source_shade = rad.prob_shade(:,1:prm.jtot)  .* power(Shade.Tsfc,4.);
           ir.IR_source = prm.epsigma * (ir.IR_source_sun + ir.IR_source_shade);
          
         

          % compute downward IR from the Top down         
          
          % following Bonan and consider forward and backward scattering of
          % IR, scat = (1-ep)/2; normally we assume reflectance is 1 -ep
          % and transmission is zero, but that only allows backward
          % scattering and my IR fluxes are not so great leaving the canopy
          
          scat=(1-prm.ep)/2;
          
          
          
            for i=prm.jtot:-1:1   
             ii=i+1;
   
%             ir.ir_dn(:,i) = ir.ir_dn(:,ii) .* leafang.integ_exp_diff(i)+ ...
%                         ir.ir_up(:,i) .*  (1-leafang.integ_exp_diff(i)) * (1-prm.ep) + ...
%                         ir.IR_source(:,i) .*(1-leafang.integ_exp_diff(i));

              ir.ir_dn(:,i) = ir.ir_dn(:,ii) .* (leafang.integ_exp_diff(i) + (1-leafang.integ_exp_diff(i)) *scat) +...
                        ir.ir_up(:,i) .*  (1-leafang.integ_exp_diff(i)) * scat + ...
                        ir.IR_source(:,i) .*(1-leafang.integ_exp_diff(i));

            end
            
            % Compute upward IR from the bottom up with the lower boundary
            % condition based on soil temperature
            
            ir.ir_up(:,1) =(1-prm.epsoil)* ir.ir_dn(:,1) + prm.epsoil .* prm.sigma .* power(soil.sfc_temperature,4);
            
             
            for i =1:prm.jtot
                ii=i+1;
                    
%              ir.ir_up(:,ii) = ir.ir_up(:,i) .* leafang.integ_exp_diff(i)+  ...
%                         ir.ir_dn(:,ii) .*  (1-leafang.integ_exp_diff(i)) * (1-prm.ep) + ...
%                       ir.IR_source(:,i) .*(1-leafang.integ_exp_diff(i)) ;      

            ir.ir_up(:,ii) = ir.ir_up(:,i) .* (leafang.integ_exp_diff(i)+ (1-leafang.integ_exp_diff(i)) *scat) + ...
                        ir.ir_dn(:,ii) .*  (1-leafang.integ_exp_diff(i)) * scat + ...
                      ir.IR_source(:,i) .*(1-leafang.integ_exp_diff(i)) ; 
                    
            end
            
         
          ir.shade=zeros(prm.nn,prm.jtot);
          
          % IR shade on top + bottom of leaves
          
         % ir.shade=(ir.ir_up(:,1:prm.jtot)+ir.ir_dn(:,2:prm.jktot));
         
          ir.shade=(ir.ir_up(:,1:prm.jtot)+ir.ir_dn(:,1:prm.jtot));  % test to see if Rnet sum = Rnet top
           
          % C code
          %ir_shade = solar.ir_dn[JJ] + solar.ir_up[JJ]; 
          
          
          % Bonan shows IR Balance for ground area
          
          ir.balance = (1-leafang.integ_exp_diff)' .*((ir.ir_up(:,1:prm.jtot)+ir.ir_dn(:,2:prm.jktot)) * prm.ep -...
              2.*ir.IR_source);
        
      figure(666)
          
      plot(nanmean(ir.ir_up(:,1:prm.jtot),1),prm.sumlai,'+-','LineWidth', 1)
      ax=gca;
      set(ax, 'ydir','reverse');
      hold on
      plot(nanmean(ir.ir_dn(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 1)
      ax=gca;
      set(ax, 'ydir','reverse');
      legend('up','down')
      xlabel('Radiation Flux Density')
      ylabel('Canopy Depth')
      title('IR flux density, W m-2');
      hold on;
      toc
end

      
      
      
      
      
      
      
      

