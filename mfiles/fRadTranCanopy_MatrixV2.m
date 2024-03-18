function [rad]=fRadTranCanopy_MatrixV2(sunang,leafang,rad,waveband,prm)

% Dennis Baldocchi
    % July 28, 2023
    
      
    % converting radiative transfer code of Norman, from C in Canveg to
    % Matlab
    
    % using Equations from Bonan to get the right layering for scattering
 

% ------------------------------------------------------------
%                      RadTranCanopy  V2
% 
%     This subroutine computes the flux density of direct and diffuse
%     radiation in the near infrared waveband.  The Markov model is used
%     to compute the probability of beam penetration through clumped foliage.
%     
%       The algorithms of Norman (1979) are used.
%       
%         Norman, J.M. 1979. Modeling the complete crop canopy.
%         Modification of the Aerial Environment of Crops.
%         B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.
% -----------------------------------------------------------------
% 
% I removed some of the do loops and took advantage of Matlabs matrix
% algebra

% this version is more general and one can assign either par or nir
% coeffients for the calculations.  And if one has spectral data for narrow
% bands one can run the model for that, too

% It also considers different leaf angle distributions and clumping

% information on sun angles is passed by sunang.*

% information on the incoming radiation for the PAR or NIR bands is passed
% through rad.*
    
% select the waveband for these computations of radiative transfer through the canopy. 
% We need to assign different transmission and reflectance coefficients
% for each waveband


% Jiang Peishi is working with my code and found some errors and typos
    %waveband='par';
    
    %waveband='nir';
    
     switch waveband
%         
         case 'par'
%             
         rad.reflect= prm.par_reflect;  % reflectance of leaf
         rad.trans=prm.par_trans;       % transmittances of leaf
         rad.soil_refl=prm.par_soil_refl; % soil reflectances
         rad.absorbed=prm.par_absorbed;   % fraction absorbed
         niter=5;                         % number of iterations  
         
         case 'nir'
             
         rad.reflect= prm.nir_reflect;
         rad.trans= prm.nir_trans ;
         rad.soil_refl=prm.nir_soil_refl ;
         rad.absorbed=prm.nir_absorbed ;
         niter=50;                       % number of interations
            
         otherwise
             
         rad.reflect= 0.25 ;
         rad.trans= 0.25 ;
         rad.soil_refl= 0.25;
         rad.absorbed= 0.5 ;
         rad.beam= 0.5*rglobal;
         rad.diffuse= 0.5* rglobal;
             
     end
     
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
     % initialize the rad fluxes in case it is night 
     
       rad.dn=zeros(prm.nn,prm.jktot);  
       rad.up=zeros(prm.nn,prm.jktot);
       rad.sun=zeros(prm.nn,prm.jktot);
       rad.shade=zeros(prm.nn,prm.jktot);

      
       beam=ones(prm.nn,prm.jktot);
            
    
      
 % incident radiation above the canopy
 
 
   %   compute the amount of diffuse light that is transmitted through a
     %   layer. diffuse radiation through gaps and absorbed radiation
     %   transmitted through leaves

        
       % rad.incoming = rad.inbeam+ rad.indiffuse;
        
        % will need to loop for this..lets run matrix and avoid
        % run code only for non zero light levels
        %if(rad.incoming(indat) >= 1. || sunang.sine_beta(indat) >=0)
        %gotgo NIRNIGHT;
        
        % fraction of beam radiation above the canopy at layer jktot
        %fraction_beam = rad.beam(indat) / (rad.beam(indat) + rad.diffuse(indat));
       
        rad.fraction_beam = rad.inbeam ./ rad.incoming;
        
        % at night if rad.incoming is zero then set fraction_beam to zero
        TF=isnan(rad.fraction_beam);
        rad.fraction_beam(TF)=0;
                
        beam(:,prm.jktot) = beam(:,prm.jktot) .* rad.fraction_beam;   % compute beam from initial ones
        
                  
       
       Tbeam=ones(prm.nn,prm.jktot);
       rad.P0=ones(prm.nn,prm.jktot);
       rad.dn = zeros(prm.nn,prm.jktot);
       rad.up=zeros(prm.nn,prm.jktot);
       
       % Compute Beam Radiation transfer and Probablity of Beam
       
     
        exp_direct = exp(-(prm.dff .* prm.markov)' .* (leafang.Gfunc ./ sunang.sine_beta));
        exp_direct(:,prm.jktot)=1;
        exp_direct(sunang.sine_beta<=0,:)=0;
        
        rad.P0=cumprod(exp_direct,2,'reverse');
        
        Tbeam=rad.P0;
        
        rad.prob_beam(:,prm.jtot)=ones;
         
             
       
        % probability of beam transfer through canopy with clumped
        % vegetation
        
        % Check Lemeur and Blad, Gutschick..
        
        % Pbeam= - cos(theta)/G dP0/dL = markov exp(-markov G L/cose(theta))
        
        %rad.prob_beam=prm.markov' .* Tbeam(:,1:prm.jtot);
        
        % rad.prob_beam=prm.markov' .* exp(-prm.markov' .* prm.sumlai' .* (leafang.Gfunc ./ sunang.sine_beta));
        
      
        rad.prob_beam=rad.P0;   
        
        rad.prob_beam=prm.markov' .* rad.P0(:,1:prm.jtot);  
        
        rad.prob_beam(sunang.sine_beta <=0,:)=0;
        
        rad.prob_beam(rad.prob_beam>1)=1;
        
        rad.prob_shade=1-rad.prob_beam;
       
        rad.prob_beam(sunang.sine_beta<0,:) = 0;
        rad.prob_shade(sunang.sine_beta<0,:)=1;    
                    
        
       % compute diffuse and complementary radiation 
       
       % from Norman 1979, as recoded by Bonan 2019
        
       rad.dn(:,prm.jktot)=(1-rad.fraction_beam).* rad.incoming;   % downward diffuse radiation at the top of the canopy
       
       for i=prm.jtot:-1:1    
             ii=i+1;

            rad.dn(:,i) = rad.dn(:,ii) .* (leafang.integ_exp_diff(i)+ (1-leafang.integ_exp_diff(i))* rad.trans) + ...
                        rad.up(:,i) .*  (1-leafang.integ_exp_diff(i)) * rad.reflect + ...
                        rad.incoming .* rad.fraction_beam .* (Tbeam(:,ii) .* (1- exp_direct(:,i)) * rad.trans);

                    
             rad.up(:,ii) = rad.up(:,i) .* (leafang.integ_exp_diff(i)+ (1-leafang.integ_exp_diff(i))* rad.trans) + ...
                        rad.dn(:,ii) .*  (1-leafang.integ_exp_diff(i)) * rad.reflect + ...
                        rad.incoming .* rad.fraction_beam .* (Tbeam(:,ii) .* (1- exp_direct(:,i)) * rad.reflect);       
                    
        end
       
       
          rad.up(:,1) = (rad.dn(:,1) + rad.incoming .* rad.fraction_beam .* Tbeam(:,1)) * rad.soil_refl;
          
           
          % iterate
          
          for k=1:niter   % #n..is more needed for NIR?
              
              
        
       for i=prm.jtot:-1:1
             ii=i+1;
%           rad.dn(:,i) = rad.dn(:,ii) .* (leafang.integ_exp_diff(i)+ (1-leafang.integ_exp_diff(i))* rad.trans) + ...
%                         rad.up(:,i) .*  (1-leafang.integ_exp_diff(i)) * rad.reflect + ...
%                         rad.incoming .* rad.fraction_beam .* (Tbeam(:,ii) .* (1- exp_direct(:,i)) * rad.reflect);

% had typo rad.reflect instead of rad.trans

           rad.dn(:,i) = rad.dn(:,ii) .* (leafang.integ_exp_diff(i)+ (1-leafang.integ_exp_diff(i))* rad.trans) + ...
                         rad.up(:,i) .*  (1-leafang.integ_exp_diff(i)) * rad.reflect + ...
                         rad.incoming .* rad.fraction_beam .* (Tbeam(:,ii) .* (1- exp_direct(:,i)) * rad.trans);

                   
             rad.up(:,ii) = rad.up(:,i) .* (leafang.integ_exp_diff(i)+ (1-leafang.integ_exp_diff(i))* rad.trans) + ...
                        rad.dn(:,ii) .*  (1-leafang.integ_exp_diff(i)) * rad.reflect + ...
                        rad.incoming .* rad.fraction_beam .* (Tbeam(:,i) .* (1- exp_direct(:,i)) * rad.reflect);  
       end
       
       
          %rad.up(:,1) = (rad.dn(:,1) + Tbeam(:,1)) * rad.soil_refl;
          
          rad.up(:,1) = (rad.dn(:,1) + rad.incoming .* Tbeam(:,1)) * rad.soil_refl;
         
              
          end  % next k
          
          
          % this new version is in Flux Density so may not need multiply by
          % rad.incoming
          
        
          

%     // Compute radiation flux densities

               
     %   upward diffuse radiation flux density, on the horizontal
     
     rad.up_flux = rad.up;
     
     %   downward beam radiation flux density, incident on the horizontal

      rad.beam_flux = Tbeam .* rad.incoming .* rad.fraction_beam;
      
      
      % downward diffuse radiation flux density on the horizontal
      
       rad.dn_flux = rad.dn ;
       
       %  total downward radiation, incident on the horizontal
       
       rad.total = rad.beam_flux + rad.dn_flux;
     
       
       
       % compute radiation absorbed on the direct and diffuse leaves of
       % each layer
     
         % normal radiation on sunlit leaves of given angle
     
         % should be radiation at top of the canopy
         
         rad.normal = rad.beam_flux(:,prm.jktot) .* leafang.Gfunc ./ sunang.sine_beta;
         
         
         
         % Gutschick
         
         % rad.normal = rad.beam_flux(:,prm.jktot) .* leafang.Gfunc;
         
         rad.normal(rad.normal <~ 0)=0;
         
        
         % amount of radiation absorbed on the sun and shade leaves
                       
         rad.sun_normal_abs = rad.normal * rad.absorbed;
                
         %rad.sh_abs= (rad.dn_flux(:,2:prm.jktot) + rad.up_flux(:,1:prm.nlayers))* rad.absorbed ;
         
         
         rad.sh_abs= (rad.dn_flux(:,1:prm.nlayers) + rad.up_flux(:,1:prm.nlayers))* rad.absorbed ;
         
               
         % C code
         
         %solar.quantum_sh[JJ] = (solar.par_down[JJ] + solar.par_up[JJ])*solar.par_absorbed;     /* umol m-2 s-1  */
         
         % Norman 1979
         % rad.sun_abs= absorption * (Rdiffuse in(j) + Rdiffuse out(j)) + absorption * R(N) *fraction beam * G/cos(zenith)
         % rad.sun_abs= rad.sun_normal_abs(:,1:prm.nlayers) + rad.sh_abs;
        
         rad.sun_abs= rad.sun_normal_abs + rad.sh_abs;
         
             
                
         % beam on soil
         
        %  rad.sun_abs(:,1)=(rad.beam_flux(:,1)+ rad.shade(:,1))  ;
        % rad.sun_abs(:,1)=(rad.beam_flux(:,1)+ rad.shade(:,1)) * (1-rad.soil_refl) ;
        
                
               
        
        
      if (waveband == 'par')
         figure(111)
         
          clf
          H=axes;
          get(H);
      end
      
      if (waveband == 'nir')
          figure(112)
          clf
      end
      
%       plot(mean(rad.sun_normal_abs,1),1:prm.jktot,'LineWidth', 2)
%       hold on
     
      plot(mean(rad.up_flux(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2);
      ax=gca;
      set(ax, 'ydir','reverse');
      
      hold on
      plot(mean(rad.dn_flux(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2)
      hold on
      plot(mean(rad.beam_flux(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2)
      hold on
      plot(mean(rad.total(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2)
      hold on
      plot(mean(rad.sun_abs(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2)
      hold on
      plot(mean(rad.sh_abs(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2)
     
     
      
      legend('up','down','beam','total','sun abs','shade abs')
      xlabel('Radiation Flux Density')
      ylabel('Canopy Cumulative LAI')
      title(waveband)

      
        
  
 end
        
