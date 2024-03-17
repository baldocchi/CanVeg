function [soil]=fSoilEnergyBalanceMatrix(quantum, nir, ir,met, prof,prm,soil,indx)

                
%                  The soil energy balance model of Campbell has been adapted to
%                  compute soil energy fluxes and temperature profiles at the soil
%                  surface.  The model has been converted from BASIC to C.  We
%                  also use an analytical version of the soil surface energy
%                  balance to solve for LE, H and G.
% 
%                  Combine surface energy balance calculations with soil heat 
%                  transfer model to calculate soil conductive and convective heat 
%                  transfer and evaporation rates.  Here, only the deep temperature 
%                  is needed and G, Hs and LEs can be derived from air temperature 
%                  and energy inputs.
% 
%                 Soil evaporation models by Kondo, Mafouf et al. and 
%                 Dammond and Simmonds are used. Dammond and Simmonds for example 
%                 have a convective adjustment to the resistance to heat transfer.
%                 Our research in Oregon and Canada have shown that this consideration
%                 is extremely important to compute G and Rn_soil correctly.


% mar 2, 2021

% was finding I was not getting energy balance closure, while
% MainSoilPhysics was.  Fixed a few things like quadratic equation for Tsfc
% was missing
        
              


        soil.n_soil_1=soil.n_soil+1;
         soil.n_soil_2=soil.n_soil+2;
         
         a_soil=zeros(prm.nn,soil.n_soil_2);
            b_soil=zeros(prm.nn,soil.n_soil_2);
               c_soil=zeros(prm.nn,soil.n_soil_2);
                  d_soil=zeros(prm.nn,soil.n_soil_2);
                        soil.storage=zeros(prm.nn,1);

       % // kv_soil is the water vapor transfer coef for the soil

        soil.water_content_sfc=met.soilmoisture;  % // at the soil surface

          

             %   //  Compute soilevap as a function of energy balance at the soil
              %  //  surface. Net incoming short and longwave energy

       
           %     // radiation balance at soil in PAR band, W m-2

             soil.par = (quantum.beam_flux(:,1) + quantum.dn_flux(:,1) - quantum.up_flux(:,1))./4.6 ;
%       
%         
%               %  // radiation balance at soil in NIR band, W m-2
% 
                 soil.nir = nir.beam_flux(:,1) + nir.dn_flux(:,1) - nir.up_flux(:,1);
% 
%                % // net incoming solar radiation balance at soil  and incoming terrestrial, W m-2
%         
                soil.Qin = soil.par + soil.nir + ir.ir_dn(:,1)*prm.epsoil;

% check if energy balance of soil is for Qin or Rnet and if this Rnet is Qin              
                
     %   // set air temperature over soil with lowest air layer, filtered

                soil.T_air = prof.Tair_K(:,2);
                
                % Initial value need to calc later
                if indx==1
                soil.sfc_temperature=soil.T_air;
                end




   %   // Compute Rh_soil and rv_soil from wind log profile for lowest layer


         % // wind speed one layer above soil

  
      u_soil=prof.wind(:,2);
   
    %   Rh_soil = 32.6 / u_soil; // bound. layer rest for heat above soil
    %   Rv_soil = 31.7 / u_soil; // bound. layer rest for vapor above soil
        


%// Stability factor from Daamen and Simmonds 1996 WRR
% older parameterization was blowing up. Went back to original paper and
% coded it

           
         soil.z0=0.02;  % one millimeter

         Ram_soil=log(prm.zht(2)/soil.z0).^2 ./(0.4*0.4 .*u_soil);
         
         %// Stability factor from Daamen and Simmonds 

     
         stabdel=5.*9.8*(prm.zht(2))*(soil.sfc_temperature-soil.T_air)./((soil.T_air).*u_soil.*u_soil);
      %   stabdel=0;     


        
        facstab(stabdel > 0)=power(1.+stabdel(stabdel > 0),-0.75);
        
        facstab(stabdel <= 0)=power(1.+stabdel(stabdel <= 0),-2.);
        


         
        facstab(facstab < .01)=.01;
        

     
        facstab(facstab > 1)=1;
       
        
         Rh_soil=Ram_soil .*facstab';
        
       % Rh_soil = 98.*facstab ./ u_soil;
      
         
        Rh_soil(Rh_soil > 1000.)=1000.;
        
              
        Rh_soil(Rh_soil < 5.)=5.;
        

  
        

        Rv_soil=Rh_soil;
        
        soil.Rv_soil=Rv_soil;
        soil.Rh_soil=Rh_soil;
        

           %     //  kcsoil is the convective transfer coeff for the soil. (W m-2 K-1)
           
       

        kcsoil = (prm.Cp .* met.air_density) ./ Rh_soil;

           %     // soil surface conductance to water vapor transfer

         kv_soil = 1. ./ (Rv_soil + soil.resistance_h2o);
                                
   
 
   %   Compute products of absolute air temperature...or Tsoil..check
   %   derivation  It is air temperature from the linearization
   
        tk1= prof.Tair_K(:,1);
   
        tk2 = tk1 .*tk1 ;
        tk3 = tk2 .* tk1;
        tk4 = tk3 .* tk1;

     %   // Longwave emission at air temperature, W m-2

%         soil.lout = prm.epsoil .*prm.sigma .* tk4;
%         soil.llout=soil.lout;
        
          %   // Slope of the vapor pressure-temperature curve, Pa/C
      %  //  evaluate as function of Tk

        dest = fdESdT(tk1);

       

      %  // Second derivative of the vapor pressure-temperature curve, Pa/C
      %  // Evaluate as function of Tk


         d2est = fd2ESdT(tk1);
         
         fact.latent=fLambda(tk1);
        
         est=fES(tk1);
         
         %vpdsoil=met.vpd_Pa;
         vpdsoil=est-prof.eair_Pa(:,1);
         
       % soil.heat = kcsoil .* (soil.T_soil(:,1)- soil.T_air_K);     % place holder later put Prof.Tair(1)

      % call finite difference routine to solve Fourier Heat Transfer equation for soil
      
     
               
      [soil]=FiniteDifferenceMatrix(soil,prm);
 

  % compute storage
    tmparray=zeros(prm.nn,1);
    tmparray(2:prm.nn,1) =(soil.T_soil(2:prm.nn,1) - soil.T_soil_old(1:prm.nn-1,1));
    
    
   %soil.gsoil is computed in FiniteDifferenceMatrix
    
    soil.storage = soil.cp_soil(:,1) .* tmparray;
    soil.gsoil = soil.gsoil + soil.storage;
         
         
      %          // coefficients for latent heat flux density

        lecoef = met.air_density .* .622 .* fact.latent .* kv_soil ./ (met.P_kPa*1000);
  
      
      %  // The quadratic coefficients for the solution to
                
       %         //   a LE^2 + b LE +c =0

        repeat = kcsoil + 4. .* prm.epsoil.* prm.sigma .* tk3;

        acoeff = lecoef .* d2est ./ (2. .* repeat);
        acoef = acoeff;
        
        bcoef = -(repeat) - lecoef .* dest + acoeff .* (-2.*soil.Qin +2 .* soil.llout+2.*soil.gsoil);
        
        ccoef = (repeat) .* lecoef .* vpdsoil + lecoef .* dest .* (soil.Qin - soil.llout-soil.gsoil) + ...
                             acoeff .* ((soil.Qin .* soil.Qin) + soil.llout .* soil.llout+soil.gsoil.*soil.gsoil - ...
                                 2 .*soil.Qin .* soil.llout-2 .*soil.Qin.*soil.gsoil+2.*soil.gsoil .*soil.llout);


     product = bcoef .* bcoef - 4 .* acoef .* ccoef;
                             
     
    
     
     
%// LE1 = (-BCOEF + (BCOEF ^ 2 - 4 * ACOEF * CCOEF) ^ .5) / (2 * ACOEF) 

       le1= (-bcoef + power(product,.5)) ./ (2 .* acoef);

%// LE2 = (-BCOEF - (BCOEF * BCOEF - 4 * acoef * CCOEF) ^ .5) / (2 * acoef) 


       
        le2= (-bcoef - power(product,.5)) ./ (2 .* acoef);
       
        %        // latent energy flux density over soil, W m-2
        
        soil.evap=real(le2);
        
        
        
         % // solve for Ts using quadratic solution  
      
        att = 6 .* prm.epsoil .*prm.sigma .* tk2 + d2est .* lecoef ./ 2;

        btt = 4 .* prm.epsoil .* prm.sigma .* tk3 + kcsoil + lecoef .* dest;

        ctt = -soil.Qin + soil.llout+soil.gsoil + lecoef .* vpdsoil;


         %       // IF (BTLF * BTLF - 4 * ATLF * CTLF) >= 0 THEN /

        product = btt .* btt - 4. .* att .* ctt;

       
        %        // T_sfc_K = TAA + 
         %       //     (-BTLF + SQR(BTLF * BTLF - 4 * ATLF * CTLF)) / (2 * ATLF) 

        if (product >= 0.)
        soil.sfc_temperature = met.T_air_K + (-btt + sqrt(product)) ./ (2. .* att);
        else
        soil.sfc_temperature=met.T_air_K;
        end
      
        
   % // Soil surface temperature, K

                soil.T_Kelvin=soil.sfc_temperature;
                
                soil.T_soil_up_boundary=soil.sfc_temperature;
                
                soil.del_Tk =soil.sfc_temperature-soil.T_air;
        
      %          // IR emissive flux density from soil, W m-2

           soil.lout_sfc = prm.epsoil .* prm.sigma .*power(soil.sfc_temperature,4);
           
           
         %%%%%%%%%%%%%%%%%%  
           
        
        % dT = (Q -LE - Gsoil -  ep sigma Ta^4)/( rho Cp gh + 4 ep sigma Ta^3)
        
         soil.del_Tk= (soil.Qin - soil.evap - soil.gsoil - soil.llout) ./repeat;

         soil.sfc_temperature= soil.T_air + soil.del_Tk;
         soil.T_Kelvin=soil.sfc_temperature;
         soil.lout_sfc = prm.epsoil .* prm.sigma .*power(soil.sfc_temperature,4);
         
         
         %      // Sensible heat flux density over soil, W m-2


             soil.heat=soil.del_Tk .* kcsoil;
             
             
           soil.rnet=soil.Qin-soil.lout_sfc;  
      

% soil energy balance test
%  figure(99)
%  clf
%  plot(soil.rnet,(soil.evap+soil.heat+soil.gsoil),'.')


 
end