function [rad]=fRadTranCanopy_Matrix(sunang,leafang,rad,waveband,prm)

% Dennis Baldocchi
    % Jan 8, 2021
    
      
    % converting radiative transfer code of Norman, from C in Canveg to
    % Matlab
% 
% /*
% ------------------------------------------------------------
%                      RadTranCanopy
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

% I removed the vertical variation in leaf area density and leaf angles with height that
% was in CANVEG. So this is a simpler set of algorithms
% 

% discovered the algorithm for ADUM in C was different from  Norman
% (1979) and in my J Applied Ecol paper.
% But it is like the code in versions of CUPID I have and found on the internet

% changed increments of SUP and SDN from J to JJP1, according to CUPID code.
% This is now giving me some complementary scattering of NIR in the upper
% canopy which I was not getting with previous code, but saw in the J
% Applied Ecology paper.

% continue to investigate simple scattering terms in Norman 1980 with
% multilevel model

% site geometric information, structural information on the canopy and optical properties of the
% vegeation are passed through information in prm.*

% information on sun angles is passed by sunang.*

% information on the incoming radiation for the PAR or NIR bands is passed
% through rad.*
    
% select the waveband for these computations of radiative transfer through the canopy. 
% We need to assign different transmission and reflectance coefficients
% for each waveband



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
         
         case 'nir'
             
         rad.reflect= prm.nir_reflect;
         rad.trans= prm.nir_trans ;
         rad.soil_refl=prm.nir_soil_refl ;
         rad.absorbed=prm.nir_absorbed ;
            
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

      
      TBEAM=ones(prm.nn,prm.jktot);
      beam=ones(prm.nn,prm.jktot);
      
%       reflectance_layer=zeros(1,prm.jktot);
%       transmission_layer=zeros(1,prm.jktot);

      %%%% do I need 2D matrix 
      ADUM=zeros(1,prm.jktot);
      
       SUP=zeros(prm.nn,prm.jktot);
       SDN=zeros(prm.nn,prm.jktot);
       
    
      
 % incident radiation above the canopy
 
        
       % rad.incoming = rad.inbeam+ rad.indiffuse;
        
        % will need to loop for this..lets run matrix and avoid
        % run code only for non zero light levels
        %if(rad.incoming(indat) >= 1. || sunang.sine_beta(indat) >=0)
        %gotgo NIRNIGHT;
        
        % fraction of beam radiation above the canopy at layer jktot
        %fraction_beam = rad.beam(indat) / (rad.beam(indat) + rad.diffuse(indat));
       
        fraction_beam = rad.inbeam ./ rad.incoming;
        
        % at night if rad.incoming is zero then set fraction_beam to zero
        TF=isnan(fraction_beam);
        fraction_beam(TF)=0;
                
        beam(:,prm.jktot) = beam(:,prm.jktot) .* fraction_beam;   % compute beam from initial ones
        
     %   TBEAM(:,prm.jktot) = TBEAM(:,prm.jktot) .* fraction_beam; % compute TBEAM from initial ones
        
       
        % SDN is beam radiation transmitted downward between two layers
        % bottom boundary condition is zero
        
       

%        Compute probability of penetration for direct and
%        diffuse radiation for each layer in the canopy
% 
%        Level 1 is the soil surface and level jktot is the
%        top of the canopy.  layer 1 is the layer above
%        the soil and layer jtot is the top layer.

%        compute the probability of diffuse radiation penetration through the
%        hemisphere.  this computation is not affected by penubra
%        since we are dealing only with diffuse radiation from a sky
%        sector.

%      Integrated probability of diffuse sky radiation penetration
%      for each layer

%      then compute the amount of diffuse radiation reflected by each layer

        reflectance_layer = (1. - leafang.integ_exp_diff) * rad.reflect;
        
       
     %   compute the amount of diffuse light that is transmitted through a
     %   layer. diffuse radiation through gaps and absorbed radiation
     %   transmitted through leaves

         transmission_layer = (1. - leafang.integ_exp_diff) * rad.trans + leafang.integ_exp_diff;
        
                 
    %     compute the probability of beam penetration through gaps of a
    %     layer
    
       exp_direct = exp(-(prm.dff .* prm.markov)' .* (leafang.Gfunc ./ sunang.sine_beta));
       
       exp_direct(sunang.sine_beta<=0,:)=0;
       
        
      %  prm.sumlai=prm.LAI:-prm.dff:0;  % cumulative LAI from top to bottom   
       
       
           

        %// Probability of beam penetration.
        
        
%       exp_direct = exp(-prm.dff * prm.markov* leafang.Gfunc(indat) / sunang.sine_beta(indat));

        % cumulative prob of beam pentration with cumulative LAI
        PEN2 = exp(-(prm.sumlai .*prm.markov)' .* (leafang.Gfunc ./ sunang.sine_beta));
        
        PEN2(sunang.sine_beta<0,:)=0;

      
        rad.prob_beam = prm.markov' .*PEN2;  % sinB/G dP0/dL
        
        
        % prob_beam at level  should be one
        rad.prob_beam(:,prm.jktot)=1;

       % // probability of shade

       
        QU = 1.0 - rad.prob_beam;

        QU(QU > 1)=1;
        QU(QU < 0)=0;


        %// probability of umbra

        rad.prob_shade = QU;
        
         rad.prob_beam(sunang.sine_beta<0,:) = 0;
         rad.prob_shade(sunang.sine_beta<0,:)=1;    
        
        
%        Beam transmission through each layer of dLAI

       for j=prm.jtot:-1:1
       beam(:,j) = beam(:,j+1) .* exp_direct(:,j);
       end

   %    beam(:,prm.jtot:-1:1) = beam(:,prm.jktot:-1:2) .* exp_direct;
       
             
        TBEAM = beam;
       
        
        % first scattering of beam intercepted, reflected and transmitted
        % in up and down directions through Layer of leaves

%         SUP(:,prm.jktot:-1:2) = (TBEAM(:,prm.jktot:-1:2) - TBEAM(:,prm.jtot:-1:1)) * rad.reflect;
%         SDN(:,prm.jtot:-1:1) = (TBEAM(:,prm.jktot:-1:2) - TBEAM(:,prm.jtot:-1:1)) * rad.trans;
%         
        
        for j=prm.jtot:-1:1
    
             dTbeam=(TBEAM(:,j+1) - TBEAM(:,j));
             SUP(:,j+1) =  dTbeam .* rad.reflect;
             SDN(:,j) = dTbeam .* rad.trans;
         end  %  // next J


%      initiate scattering using the technique of NORMAN (1979).
%      scattering is computed using an iterative technique.
%      Initially compute up and down diffuse fields without any direct 
%      complementary radiation. Confirming the code by consulting a version
%      of CUPID  http://soils.wisc.edu/facstaff/wayne/cupid/radiat.html
%      subroutine Radiat (curadia.f)

         SUP(:,1) = TBEAM(:,1) * rad.soil_refl;  % bottom boundary for beam reflected up at soil
         rad.dn(:,prm.jktot) = 1. - fraction_beam; % fraction of radiation that is diffuse
         ADUM(1) = rad.soil_refl;  % bottom boundary condition for ADUM
         
         % first step to compute diffuse fields and scattering is assume
         % beam fraction is zero
         
         % an array of transmission_layer and reflectance_layer was
         % computed.  they did not change much so for expedience I am
         % computing the mean for this loop

          
         
          TLAY2 = transmission_layer .* transmission_layer;
          
          for i=2:prm.jktot
              ii=i-1;
          ADUM(i) = (ADUM(ii) .* TLAY2(ii)') ./ (1. - (ADUM(ii) .* reflectance_layer(ii)')) + reflectance_layer(ii)';             
          end
         
         for i=prm.jtot:-1:1
             ii=i+1;
         rad.dn(:,i) = rad.dn(:,ii) .* transmission_layer(i)' ./ (1. - ADUM(ii) * reflectance_layer(i)') + SDN(:,i);
         rad.up(:,ii) = ADUM(ii) .* rad.dn(:,ii) + SUP(:,ii);
         end

        
       rad.up(:,1) = rad.soil_refl * rad.dn(:,1) + SUP(:,1);
       
% Iterate
       for k=1:5
           
 %       lower boundary: upward radiation from soil

       
  
         
%     Iterative calculation of upward diffuse and downward beam +
%     diffuse radiation to compute scattering

           
        for j=2:prm.jktot   %i=prm.jtot:-1:1
            %ii=i+1;
            jj=prm.jktot-j+1;
            jjp1=jj+1;
            
%        rad.dn(:,i)= transmission_layer * rad.dn(:,ii) + rad.up(:,i) * reflectance_layer + SDN(:,i);
        rad.dn(:,jj)= transmission_layer(jj)' * rad.dn(:,jjp1) + rad.up(:,jj) * reflectance_layer(jj)' + SDN(:,jj);
        end  % next j
      
       
       
%       upward radiation at soil is reflected beam and downward diffuse  

       rad.up(:,1) = (rad.dn(:,1) + TBEAM(:,1)) * rad.soil_refl;

       % iterate
       for i=2:prm.jktot
           ii=i-1;
      rad.up(:,i)= reflectance_layer(ii)' * rad.dn(:,i) + rad.up(:,ii) * transmission_layer(ii)' + SUP(:,i);
       end
       
       end    % next k
       
  

%     // Compute radiation flux densities


     %rad.total(indat,:) = rad.inbeam(indat) + rad.indiffuse(indat);
     llai = prm.LAI;

          
     %   upward diffuse radiation flux density, on the horizontal
     
     rad.up_flux = rad.up .* rad.incoming;
     
     %   downward beam radiation flux density, incident on the horizontal


    %%%%% beam is not the same as prob_beam...? double check
    
    
    %%%%%getting too low of beam flux in computations
     
      rad.beam_flux = beam .* rad.incoming;
      
       
      % downward diffuse radiation flux density on the horizontal
      
       rad.dn_flux = rad.dn .* rad.incoming;
       
       %  total downward radiation, incident on the horizontal
       
       rad.total = rad.beam_flux + rad.dn_flux;
     
     
     %// normal radiation on sunlit leaves
     
     % should be radiation at top of the canopy
         
         rad.normal = rad.beam_flux(:,prm.jktot) .* leafang.Gfunc ./ sunang.sine_beta;
         rad.normal(rad.normal <~ 0)=0;
         
         
         
        
         % amount of radiation absorbed on the sun and shade leaves
         
         % rad.sun_normal = rad.normal .* rad.absorbed;
         
         rad.sun_normal_abs = rad.normal * rad.absorbed;
                
         rad.sh_abs= (rad.dn_flux(:,2:prm.jktot) + rad.up_flux(:,1:prm.nlayers))* rad.absorbed;
     
     
        % rad.sun_abs= rad.sun_normal_abs(:,1:prm.nlayers) + rad.sh_abs;
        
         rad.sun_abs= rad.sun_normal_abs + rad.sh_abs;
         
         % beam on soil
         
         rad.sun_abs(:,1)=rad.beam_flux(:,1)+ rad.shade(:,1) ;
        
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
      
      legend('up','down','beam','total')
      xlabel('Radiation Flux Density')
      ylabel('Canopy Cumulative LAI')
      title(waveband)

        
  
 end
        
