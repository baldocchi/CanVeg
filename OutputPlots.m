function [] = OutputPlots(prm,Can,soil,met,prof,Veg)

% 8/2, 2023

% list of plots to evaluate the model

 ttime=met.day +met.hhour/24;
 
    % compare computations with CanVeg C version
    
    switch prm.Veg
        case 'Tule'
    
    data=readtable(prm.fluxdata);
   
     
    LE=data.Var3;
      H=data.Var4;
       GPP=data.Var5;
       Rnet=data.Var6;
        albedo=data.Var7;
           Fco2=data.Var8;
             Gsoil=data.Var9;
    
    
   
        case 'Alfalfa'
   
            % output fluxes from CanalfalfaInput.m
            
     data=readtable(prm.fluxdata);
     
     
  
% Bouldin_Island_alfalfa.m
% spectral corrections applied to prm.fluxdata

% outflx=data.DOY(use);
% outflx=horzcat(outflx,data.time(use));
% outflx=horzcat(outflx,data.LE_gf(use));
% outflx=horzcat(outflx,data.H_gf(use));
% outflx=horzcat(outflx,data.gpp_ANNnight(use));
% outflx=horzcat(outflx,Metdata.RNET(use));
% outflx=horzcat(outflx,albedo(use));
% outflx=horzcat(outflx,data.wc_gf(use));
% outflx=horzcat(outflx,Gsoil_crt(use));


% output1=DOY(use);
% output1=horzcat(output1,hhmm(use));
% output1=horzcat(output1,alf.RNET(use));
% output1=horzcat(output1,alf.LE_gf(use));
% output1=horzcat(output1,alf.H_gf(use));
% output1=horzcat(output1,alf.wc_gf(use));
% output1=horzcat(output1,-alf.gpp_Reichstein(use));
% output1=horzcat(output1,albedo(use));
% output1=horzcat(output1,Gsoil(use));
% output1=horzcat(output1,Trad(use));
% output1=horzcat(output1,NIRv(use));
% output1=horzcat(output1,LongIn(use));



%strg={'DOY','hhmm','Rnet','LE','H','Fco2','GPP','albedo','Gsoil','Trad','NIRv','LongIn'};


% output fluxes
strg={'DOY','hhmm','LE','H','GPP','Rnet','albedo','Fco2','Gsoil','Trad','NIRv', 'LongIn'};
      
      LE=data.Var4;
      H=data.Var5;
      GPP=data.Var7;  % make positive sign
      Rnet=data.Var3;
      albedo=data.Var8;
      Fco2=data.Var6;   
      Gsoil=data.Var9;
      Trad=data.Var10;
      NIRv=data.Var11;
      LongIn=data.Var12;
      
% corrections applied to prm.fluxdata
         LE=LE/.8494;
         H=H/1.0186;
         Fco2=Fco2/0.792;   % spectral correction is close to one
       
       
    
          

     
        case 'DeciduousForest'
            
             data=readtable(prm.fluxdata);

      
     LE=data.Var3;
      H=data.Var4;
       GPP=data.Var5;
       Rnet=data.Var6;
        albedo=data.Var7;
           Fco2=data.Var8;
             Gsoil=data.Var9;
             
             LE(LE==-9999)=NaN;
             H(H==-9999)=NaN;
             GPP(GPP==-9999)=NaN;
             Rnet(Rnet==-9999)=NaN;
             Fco2(Fco2==-9999)=NaN;
             Gsoil(Gsoil==-9999)=NaN;


        case 'Savanna'
            
             data=readtable(prm.fluxdata);

%             strg={'DOY','hhmm','Rnet','LE','H','Fco2','Gsoil','GPP'};


        Rnet=data.Var3;
        LE=data.Var4;
        H=data.Var5;
        Fco2=data.Var6;
        Gsoil=data.Var7;
        GPP=data.Var8;
             
             
             LE(LE==-9999)=NaN;
             H(H==-9999)=NaN;
             Rnet(Rnet==-9999)=NaN;
             Fco2(Fco2==-9999)=NaN;
        
        otherwise
     
    end
    
     % plot LAI profile
    fnum=0;
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(prm.dff,prm.zht(1:prm.nlayers),'-.','LineWidth',2)
    ylabel('height')
    xlabel('LAI per layer')
    
    
     fnum=fnum+1;
    figure(fnum)
    clf
    plot(prm.adens,prm.zht(1:prm.nlayers),'-.','LineWidth',2)
    ylabel('height')
    xlabel('LAI density per layer, m2 m-3')
    
        
    fnum=fnum+1;
    figure(fnum)
    clf
    EBClosure=Rnet - Gsoil- LE - H;
    
    testEB=find(abs(EBClosure) < 50);
    
    histogram(EBClosure,30)
    xlabel('Energy Balance Closure Rnet -Gsoil - LE  - H) measured')
    
    nanmedian( EBClosure)
    
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(Rnet,Gsoil+LE + H,'.')
    xlabel('Rnet measured, W m-2')
    ylabel('H + LE + Gsoil, measured, W m-2')
    
    
    

 [ GsPM ] = Gsfc_PM(met,Can); 
      
      
      
      fnum=fnum+1;
      figure(fnum)
      clf
      plot(GsPM,Veg.gs,'.','MarkerSize',10)
      xlim([0 0.025])
      ylim([0 0.025])
      xlabel('Canopy Conductance, m/s')
      ylabel('Canopy Stomatal Conductance, m/s')
      title(prm.title)
    
      %xlim([0 0.015])
      
    % plot mean diurnal Ps
    
      
  
     Divpd=reshape(Veg.vpd,prm.hrs,prm.ndays);
    
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(Divpd,2))  
    xlabel('hour')
    ylabel ('vpd Pa')
    title(prm.title);
    
  %  GsPM(GsPM > 0.025 | GsPM < 0)=NaN;
        
    Digs=reshape(Veg.gs,prm.hrs,prm.ndays);
    DigsPM=reshape(GsPM,prm.hrs,prm.ndays);
    
   
     
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(Digs,2))  
    hold on
     plot(1:prm.hrs,nanmean(DigsPM,2))  
    xlabel('hour')
    ylabel ('gs m/s')
    title(prm.title);
    ylim([0 0.4])
    legend('gs','gsPM')
    
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(1./Digs,2))  
    xlabel('hour')
    ylabel ('rs m/s')
    title(prm.title);
    
 %    Can.Rnet=[Can.Rnet;NaN];
     Dirnet=reshape(Can.Rnet,prm.hrs,prm.ndays);
     DiRnet=reshape(Rnet,prm.hrs,prm.ndays);
   
     
     Dirtsfc=reshape(Veg.Tsfc,prm.hrs,prm.ndays);
     Dirtair=reshape(met.T_air_K,prm.hrs,prm.ndays);
     
  fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(Dirtsfc,2)) 
    hold on
    plot(1:prm.hrs,nanmean(Dirtair,2)) 
    xlabel('hour')
    ylabel ('T C')
    title(prm.title);
    legend('Tsfc','Tair')
    
    
    fnum=fnum+1;
    figure(fnum)
    clf
    histogram(Dirtsfc- Dirtair,50)
    xlabel('Tsfc - Tair')
    title(prm.title);
    
   
      Bowen=H./LE;
      
      Bowen(abs(Bowen) > 2.5)=NaN;
      
      LEcorr =(Rnet -Gsoil) ./ (1+Bowen);
      
      LEcorr(LEcorr > 700 | LEcorr < -100)=NaN;
      
     fnum=fnum+1;
    figure(fnum)
      clf
      histogram(LEcorr,30)
      xlim([-100 700])
      xlabel('LE Measured, Corrected with Bowen Ratio, W m-2')
      title(prm.title)
    
    % with spectral corrected LE 
     fnum=fnum+1;
    figure(fnum)
      clf
      plot(LE,Can.LE, '.')
      hold on
          
     plot(LEcorr,Can.LE, '.')
      
      xlabel(' LE measured')  
      ylabel(' LE calculated')  
      title(prm.title);
      xlim([0 max(Can.LE)])
      ylim([0  max(Can.LE) ])
      legend('Eddy Cov','Bowen Ratio')
      
      
      % with spectral corrected LE 
     fnum=fnum+1;
    figure(fnum)
      clf
      plot(LE,Can.LE, '.')
      
      
      xlabel(' LE measured, spectral correction')  
      ylabel(' LE calculated')  
      title(prm.title);
      xlim([0 max(Can.LE)])
      ylim([0  max(Can.LE) ])
      legend('Eddy Cov')
      
      
      
      
      % check spectral attenutation for LE
      
      %%%% what if we use better energy balance closure,  testEB
      
     fnum=fnum+1;
    figure(fnum)
      clf
      plot(LE(testEB),Can.LE(testEB), '.')
      
     xlabel(' LE measured,  conditional on abs(EBClosure) < 50);')  
      ylabel(' LE calculated')  
      title(prm.title);
      xlim([0 max(Can.LE)])
      ylim([0  max(Can.LE) ])
      
      
      
    DiLE=reshape(Can.LE,prm.hrs,prm.ndays);
    DiLEec=reshape(LE,prm.hrs,prm.ndays);
    DiLEBowen=reshape(LEcorr,prm.hrs,prm.ndays);
    
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(DiLE,2),'LineWidth',3)  
    hold on
    plot(1:prm.hrs,nanmean(DiLEec,2),'LineWidth',3)  
    % hold on
    % plot(1:prm.hrs,nanmean(DiLEBowen,2),'LineWidth',3)  
    xlabel('hour')
    ylabel ('LE W m-2')
    title(prm.title);
    legend('Computed','Measured')
        
      
    
    DiH=reshape(Can.H,prm.hrs,prm.ndays);
    DiHec=reshape(H,prm.hrs,prm.ndays);
    
   
    
   fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(DiH,2))  
    hold on
    plot(1:prm.hrs,nanmean(DiHec,2))  
    xlabel('hour')
    ylabel ('H W m-2')
    title(prm.title);
    legend('Computed','Measured')
        
    fnum=fnum+1;
    figure(fnum)
    clf
    histogram(H-Can.H,25)
    xlabel('H measured - calculated, W m-2')
    
    nanmean(H-Can.H)
    
    
    fnum=fnum+1;
    figure(fnum)
    clf
    histogram(LE-Can.LE,25)
    xlabel('LE measured - calculated, W m-2')
      
     fnum=fnum+1;
    figure(fnum)
      clf
      plot(abs(GPP(Veg.Gpp<60)),Veg.Gpp(Veg.Gpp<60), '.')
      xlabel(' Ac measured \mu mol m^{-2} s^{-1}')  
      ylabel(' Ac calculated')  
       title(prm.title);
       
    % xlim([0 50])
    % ylim([0 50])
       
       
    fnum=fnum+1;
    figure(fnum)
      clf
      plot(H,Can.H, '.')
%       hold on
%       Hcorr =(Rnet-Gsoil)./(1 + 1/Bowen);
%       plot(Hcorr,Can.H, '.')
      xlabel(' H measured')  
      ylabel(' H calculated')  
      title(prm.title);
      legend('Eddy Cov')
      
      
      ccc=sum(prof.Rnet .* prof.delz(1:prm.jtot)',2);
      
      fnum=fnum+1;
    figure(fnum)
      clf
      plot(Rnet,Can.Rnet_calc, '.')
      hold on
      plot(Rnet,Can.Rnet, '.')

      xlabel(' Rnet measured')  
      ylabel(' Rnet calculated, top of canopy')  
      legend('Rnet top of canopy','Rnet Integ')
       title(prm.title);
              
       xlim([-100 900])
       
      fnum=fnum+1;
    figure(fnum)
      clf
      plot(Can.Rnet_calc,Can.Rnet, '.')
   
      xlabel(' Rnet calculated, top of canopy')  
      ylabel(' Can.Rnet')  
       title(prm.title);
              
       xlim([-100 900])
       
       
       fnum=fnum+1;
    figure(fnum)
        clf
      histogram(Rnet-Can.Rnet, 25)
       xlabel('Rnet measured - calculated, W m-2')
       
      fnum=fnum+1;
    figure(fnum)
      clf
      plot(nanmean(prof.Tair_K,1),prof.zht(1:prm.nlayers_atmos));
      xlabel('Tair, K')
      ylabel('ht m')
      
      fnum=fnum+1;
    figure(fnum)
      clf
       plot(nanmean(prof.H,1),prof.zht(1:prm.jtot));
       xlabel('dH/dz')
         ylabel('ht m')
      
      
      
     fnum=fnum+1;
    figure(fnum)
      clf
      plot(nanmean(prof.eair_Pa,1),prof.zht(1:prm.nlayers_atmos));
      xlabel('eair, Pa')
      ylabel('ht m')
      
      fnum=fnum+1;
    figure(fnum)
      clf
       plot(nanmean(prof.LE,1),prof.zht(1:prm.jtot));
       xlabel('dLE/dz')
         ylabel('ht m')
      
      
      
     fnum=fnum+1;
    figure(fnum)
      clf
      plot(nanmean(prof.co2,1),prof.zht(1:prm.nlayers_atmos));
      xlabel('CO2, ppm')
      ylabel('ht m')
      
      
     fnum=fnum+1;
    figure(fnum)
      clf
      plot(nanmean(prof.Ps,1),prof.zht(1:prm.jtot));
      xlabel('Ps, \mu mol m-3 s-1')
      ylabel('ht m')
      
      
      
      
     fnum=fnum+1;
    figure(fnum)
      clf
      plot(nanmean(prof.Rnet,1),prof.zht(1:prm.jtot));
      xlabel('Rnet, W m-3 s-1')
      ylabel('ht m')
      
      
      
      
     DiGPP=reshape(abs(GPP),prm.hrs,prm.ndays); 
     DiPs=reshape(Veg.Gpp,prm.hrs,prm.ndays);  
    
   
    
    fnum=fnum+1;
    figure(fnum)
     clf
    plot(1:prm.hrs,nanmean(DiPs,2),'LineWidth',3)  
    hold on
    plot(1:prm.hrs,nanmean(DiGPP,2),'LineWidth',3)  
    
    xlabel('hour')
    ylabel ('Canopy Photosynthesis, umol m-2 s-1')
    %title('CanVeg, Tule, West Pond')
    title(prm.title); 
     legend('Computed','Measured')

    fnum=fnum+1;
    figure(fnum)
    clf
    plot(Fco2,Can.NEE,'.','MarkerSize',12)
    xlabel('NEE measured')
    ylabel('NEE computed')
   
   
    
      
    
     DiLE=reshape(Can.LE,prm.hrs,prm.ndays);
     
     DiET=reshape(LE,prm.hrs,prm.ndays); 
    
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(DiLE,2))  
    hold on
    plot(1:prm.hrs,nanmean(DiET,2))  
    xlabel('hour')
    ylabel ('LE W m-2')
    legend('Computed','Measured')
    title(prm.title);
 
    
      VegType=prm.Veg;
     
     switch VegType
         case 'Tule'
      DiGec=reshape(met.Gwater,prm.hrs,prm.ndays);
         case 'Alfalfa'
         DiGec=reshape(Gsoil,prm.hrs,prm.ndays);
         case 'DeciduousForest'
              DiGec=reshape(Gsoil,prm.hrs,prm.ndays);
         
         otherwise
    end 
    
   
       
      
   fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(DiRnet,2),'LineWidth', 1.25);
    hold on
    plot(1:prm.hrs,nanmean(DiLEec,2),'LineWidth', 1.25)  
    hold on
    plot(1:prm.hrs,nanmean(DiHec,2),'LineWidth', 1.25)  
     hold on
  %  plot(1:prm.hrs,nanmean(DiGec,2),'LineWidth', 1.25) 
    xlabel('hour')
    ylabel ('Energy Flux Density, W m-2')
    title(prm.title);
    legend('Rnmeas','LEmeas','Hmeas','Gsoil/watermeas')
     
    
    
     
    
    % fnum=fnum+1;
    % figure(fnum)
    % clf
    % 
    %   Dialbedo=reshape(Can.albedo_calc,prm.hrs,prm.ndays); 
    %   DiAlbedo=reshape(albedo,prm.hrs,prm.ndays); %measured
    % 
    %   DiAlbedo(DiAlbedo <0 | DiAlbedo > 0.6)=NaN;
    % 
%          plot(12:36,nanmean(Dialbedo(12:36,:),2),'LineWidth', 1.25);
%          hold on
%            plot(12:36,nanmean(DiAlbedo(12:36,:),2),'LineWidth', 1.25);

  % plot(8:16,nanmean(Dialbedo(8:16,:),2),'LineWidth', 1.25);
  %        hold on
  %          plot(8:16,nanmean(DiAlbedo(8:16,:),2),'LineWidth', 1.25);
  % 
  %   ylabel('albedo')
  %   legend('Computed','Measured')
           
           
            % Canopy Photosynthesis vs reflected NIR  
           
           fnum=fnum+1;
           figure(fnum)
           clf
           plot(Can.NIR_refl,Veg.Gpp,'.')
                 
           
           xlabel('NIR reflected by Vegetation W m^{-2}')
           ylabel('Canopy Assimilation \mu mol m^{-2} s^{-1}')
           title(prm.title);
           
           
      
                      
            % Canopy conductance vs reflected NIR  
           
           fnum=fnum+1;
          figure(fnum)
           clf
           plot(Can.NIR_refl,Veg.gs,'.')
           xlabel('NIR reflected by Vegetation, W m^{-2}')
           ylabel('Canopy Conductance m s^{-1}')
          
           ylim([0 .03])
         % Flux weighted Tair vs Trad
                  
         Taero =nanmean((prof.H(:,1:prm.nlayers) .* prof.Tair_K(:,1:prm.nlayers))./(prof.H(:,1:prm.nlayers)),2);
         
         
          nanmean( nanmean( Veg.Tsfc-Can.Trad))
         
        fnum=fnum+1;
    figure(fnum)
         clf
         plot(Can.Trad,Taero,'.')
         
       xlabel('Tsfc, K')
       ylabel('Taero, H weighted Tair')
       title(prm.title)
       
       
     nanmean( nanmean( Veg.Tsfc-Taero))
     
     
     
     fnum=fnum+1;
    figure(fnum)
         clf
         plot(Can.Trad,Can.Taero,'.')
         
       xlabel('Tsfc, K')
       ylabel('Taero, Dij * H weighted Tair')
       title(prm.title)
     
     
     
 
 %%% plot soil temperature profiles, with std/sqrt(nn)
 
 sqrtnn=sqrt(prm.nn);
 
 
 
 
 switch VegType
    case 'Alfalfa'
        
      fnum=fnum+1;
    figure(fnum)
       clf
       plot(Gsoil,soil.gsoil,'.')
 
 xlabel('Gsoil, measured W m-2')
 ylabel('Gsoil, computed, W m-2')
 title (prm.title)
      
        
        fnum=fnum+1;
    figure(fnum)
 clf
 histogram(Gsoil-soil.gsoil,25)
 xlabel('Gsoil measured- Gsoil calculated, W m-2')
 nanmean(Gsoil-soil.gsoil)
 
 fnum=fnum+1;
    figure(fnum)
 clf
 xmeas=[nanmean(met.Tsoil_2);nanmean(met.Tsoil_4);nanmean(met.Tsoil_8);nanmean(met.Tsoil)];
 errxmeas=[nanstd(met.Tsoil_2);nanstd(met.Tsoil_4);nanstd(met.Tsoil_8);nanstd(met.Tsoil)];
 errxmeas =errxmeas/sqrtnn;
 ymeas=[.02,.04,.08,.16];
 errorbar(xmeas,ymeas,errxmeas,'horizontal','LineWidth',2)
 %plot(xmeas,ymeas,'.-');
 
 hold on
 errxcalc=nanstd(soil.T_soil(:,1:11))./sqrtnn;
 xcalc=nanmean(soil.T_soil(:,1:11))-273.15;
 errorbar(xcalc,soil.z_soil,errxcalc,'horizontal','LineWidth',2)
 
 ax=gca;
 set(ax, 'ydir','reverse');
 
  %plot(nanmean(soil.T_soil(:,1:11))-273.15,soil.z_soil,'+-');
  legend('measured','calculated');
  ylabel('soil depth, m');
  xlabel('Soil Temperature, C');
  
  
    
     canopyht=reshape(met.zcanopy,prm.hrs,prm.ndays);
     daymat=reshape(met.day,prm.hrs,prm.ndays);
     dday=nanmean(daymat);
     
     fnum=fnum+1;
    figure(fnum)
     clf
     plot( dday, nanmean(canopyht(10:38,:),1),'.','MarkerSize',10)
     ylabel('canopy height, m')
     title('aerodynamic canopy height')
     xlabel('day of year')
     
     
     canopy_lai=reshape(met.LAI,prm.hrs,prm.ndays);
     
     fnum=fnum+1;
    figure(fnum)
     clf
     plot( dday, nanmean(canopy_lai(10:38,:),1),'.','MarkerSize',10)
     ylabel('canopy LAI, m')
     title('Bouldin Alfalfa')
     xlabel('day of year')
     
     
     % 
  
 
  case 'Tule'
     
      use=(met.WaterTable > 16);
         
       fnum=fnum+1;
    figure(fnum)
       clf
       plot(met.Gwater(use),soil.gsoil(use),'.')
 
 xlabel('Gsoil, measured as storage in water column W m-2')
 ylabel('Gsoil, computed, W m-2')
 title (prm.title)
      
  fnum=fnum+1;
    figure(fnum)
 clf
 histogram(met.Gwater(use)-soil.gsoil(use),25)
 xlabel('Gsoil measured- Gsoil calculated, W m-2')
 nanmean(met.Gwater-soil.gsoil)
 title (prm.title)
      
     fnum=fnum+1;
    figure(fnum)
 clf
 xmeas=[nanmean(met.Twater_4cm);nanmean(met.Twater_8cm);nanmean(met.Twater)];
 errxmeas=[nanstd(met.Twater_4cm);nanstd(met.Twater_8cm);nanstd(met.Twater)];
 errxmeas =errxmeas/sqrtnn;
 ymeas=[.04,.08,.16];
 errorbar(xmeas,ymeas,errxmeas,'horizontal','LineWidth',2)
 %plot(xmeas,ymeas,'.-');
 
 hold on
 errxcalc=nanstd(soil.T_soil(:,1:11))./sqrtnn;
 xcalc=nanmean(soil.T_soil(:,1:11))-273.15;
 errorbar(xcalc,soil.z_soil,errxcalc,'horizontal','LineWidth',2)
 
  %plot(nanmean(soil.T_soil(:,1:11))-273.15,soil.z_soil,'+-');
  legend('measured','calculated');
  ylabel('soil depth, m');
  xlabel('Soil Temperature, C');
  title (prm.title)
  
  fnum=fnum+1;
    figure(fnum)
 
  plot(ttime,met.WaterTable,'.')
  ylabel('Water Table, cm')
  xlabel('Day-hour')
  title (prm.title)
  
  
    
     canopyht=reshape(met.zcanopy,prm.hrs,prm.ndays);
     daymat=reshape(met.day,prm.hrs,prm.ndays);
     dday=nanmean(daymat);
     
     fnum=fnum+1;
    figure(fnum)
     clf
     plot( dday, nanmean(canopyht(10:38,:),1))
     ylabel('canopy height, m')
     title('aerodynamic canopy height')
     xlabel('day of year')
  
     case 'DeciduousForest'
         
           fnum=fnum+1;
    figure(fnum)
       clf
       plot(Gsoil,soil.gsoil,'.')
 
 xlabel('Gsoil, measured W m-2')
 ylabel('Gsoil, computed, W m-2')
 title (prm.title)
      
        
        fnum=fnum+1;
    figure(fnum)
 clf
 histogram(Gsoil-soil.gsoil,25)
 xlabel('Gsoil measured- Gsoil calculated, W m-2')
 nanmean(Gsoil-soil.gsoil)
 
        case 'Savanna'
         
           fnum=fnum+1;
    figure(fnum)
       clf
       plot(Gsoil,soil.gsoil,'.')
 
 xlabel('Gsoil, measured W m-2')
 ylabel('Gsoil, computed, W m-2')
 title (prm.title)
      
        
        fnum=fnum+1;
    figure(fnum)
 clf
 histogram(Gsoil-soil.gsoil,25)
 xlabel('Gsoil measured- Gsoil calculated, W m-2')
 nanmean(Gsoil-soil.gsoil)
 
      
      
      otherwise   
end
 
 
           % test energy balance closure
           
           x=Can.Rnet ;
           y=Can.H + Can.LE + soil.gsoil;
           
           
        fnum=fnum+1;
    figure(fnum)
         clf
         plot(x,y,'.')
         xlabel('Rnet, calculated, W m-2')
         ylabel('H + LE + Gsoil, Canveg')
         
         
             
    
% plot soil T profile and measurements


     DiLEsoil=reshape(soil.evap,prm.hrs,prm.ndays);
        DiHsoil=reshape(soil.heat,prm.hrs,prm.ndays);
           DiGsoil=reshape(soil.gsoil,prm.hrs,prm.ndays);
              DiRnsoil=reshape(soil.rnet,prm.hrs,prm.ndays);
              
              fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(Dirnet,2),'LineWidth', 1.25);
    hold on
    plot(1:prm.hrs,nanmean(DiLE,2),'LineWidth', 1.25)  
    hold on
    plot(1:prm.hrs,nanmean(DiH,2),'LineWidth', 1.25)  
     hold on
    plot(1:prm.hrs,nanmean(DiGsoil,2),'LineWidth', 1.25) 
    xlabel('hour')
    ylabel ('Energy Flux Density, W m-2')
    title(prm.title);
    legend('Rn','LE','H','Gsoil')

              
              
               fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(Dirnet,2),'LineWidth', 1.25);
    hold on
    plot(1:prm.hrs,nanmean(DiLE,2),'LineWidth', 1.25)  
    hold on
    plot(1:prm.hrs,nanmean(DiH,2),'LineWidth', 1.25)  
    hold on
    plot(1:prm.hrs,nanmean(DiGsoil,2),'LineWidth', 1.25)  
    xlabel('hour')
    ylabel ('Energy Flux Density, W m-2')
    title(prm.title);
    legend('Rn','LE','H','Gsoil')
      
              
              
              
                fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(DiRnsoil,2),'LineWidth', 1.25);
    hold on
    plot(1:prm.hrs,nanmean(DiLEsoil,2),'LineWidth', 1.25)  
    hold on
    plot(1:prm.hrs,nanmean(DiHsoil,2),'LineWidth', 1.25)  
    hold on
    plot(1:prm.hrs,nanmean(DiGsoil,2),'LineWidth', 1.25)  
    xlabel('hour')
    ylabel ('Energy Flux Density, W m-2')
    %title('CanVeg, Tule, West Pond')
    title(prm.title);
    legend('Rnsoil','LEsoil','Hsoil','Gsoil')
    
    
   fnum=fnum+1;
    figure(fnum)
    clf
    plot(nanmean(soil.T_soil(:,1:soil.n_soil_1)),-soil.z_soil(1:soil.n_soil_1),'.-', 'LineWidth',2)
    ylabel('soil depth, m')
    xlabel('Soil Temperature, K')
   
    
 
    
     DiFco2=reshape(Fco2,prm.hrs,prm.ndays); 
     DiNEE=reshape(Can.NEE,prm.hrs,prm.ndays); 
     
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(1:prm.hrs,nanmean(DiFco2,2))  
    hold on
    plot(1:prm.hrs,nanmean(DiNEE,2))  
    legend('NEE measured','NEE calculated')
    ylabel('NEE \mu mol m-2 s-1')
    
    
    % histogram of upper leaf temperature. look at clumping effect
    fnum=fnum+1;
    figure(fnum)
    clf
    
    n1=prm.nlayers- floor(1/ prm.dff);
    histogram(prof.Tsfc(:,n1:prm.nlayers));
    xlabel('Tsfc, upper 1 m2/m2 lai layers, K');
    nanmean(nanmean(prof.Tsfc(:,n1:prm.nlayers)))
    mean(kurtosis(prof.Tsfc(:,n1:prm.nlayers)))
    
    
    
    % equilibrium evaporation
    
    LEeq= fdESdT(met.T_air_K)./(fdESdT(met.T_air_K)+ 66.5) .*(Rnet-Gsoil);
    LEeqmat=reshape(LEeq,prm.hrs,prm.ndays);
    
    LE_LEeq=nansum(DiLEec)./nansum(LEeqmat);
    
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(LE_LEeq,'.-','MarkerSize',10)
    
    hold on
    
    LEeqCan= fdESdT(met.T_air_K)./(fdESdT(met.T_air_K)+ 66.5) .*(Can.Rnet-Can.Gsoil);
    LEeqCanmat=reshape(LEeqCan,prm.hrs,prm.ndays);
    
    LE_LEeqCan=nansum(DiLE)./nansum(LEeqCanmat);
    
    plot(LE_LEeqCan,'.-','MarkerSize',10)
      
    title ('Bouldin alfalfa')
    xlabel('days')
    ylabel('Priestley Taylor Coefficient, based on Daily Sums')
    legend('Measured','Computed')
    
    nanmedian(LE_LEeq)
     nanmedian(LE_LEeqCan)
       
    % water use efficiency
    
    fnum=fnum+1;
    figure(fnum)
    clf
    plot(LE,abs(GPP),'.')
    hold on
    plot(Can.LE,Veg.Gpp,'+r')
    xlabel('LE W m-2')
    ylabel('GPP \mu mole m-2 s-1')
    title('water use efficiency, alfalfa')
    legend('measured','Canveg')
    
    
    
        
     Trad=met.Trad;
     
    fnum=fnum+1;
    figure(fnum)
     clf
   plot(Trad,Can.Trad,'.')
   xlabel('T_{rad} measured K')
   ylabel('T_{rad} computed K')
     title (prm.title)
     
     
    fnum=fnum+1;
    figure(fnum)
     clf
   histogram(Trad-Can.Trad,25)
   xlabel('T_{rad} measured - T_{rad} computed K')
   title (prm.title)
   
   median(Trad-Can.Trad)
      
   
   
     Taero=Can.Taero;
       Trad=met.Trad;
       
      %save('d:\CanAlfalfa\TaeroTrad_L5.mat','Taero','Trad');
      %save('d:\CanAlfalfa\TaeroTrad_L4.mat','Taero','Trad');
    
      % Does Rnet decay in the canopy with exponential function, eg Beer's
      % Law
    
%      x=reshape(quantum.P0,48,17,51);
%      xx=nanmean(nanmean(x(20:32,:,1:50),2));
%      
%    
%      y=reshape(prof.Rnet,48,17,50);
%      yy=nanmean(nanmean(y(14:36,:,:),2));
%     
%      yyy=reshape(yy,1,50);
%      xxx=reshape(xx,1,50);
%           
%      figure(58)
%      clf
%      plot(xxx,yyy,'.')
%      xlabel('P_0')
%      ylabel('Rnet, W m-2')
%      title('CanVeg, Alfalfa, 0700-1700')
     
    
     
%      x=reshape(quantum.P0,48,17,51);
%      xx=x(:,:,1:50);
%    
%      use=find(xx>0);
%      
%      y=reshape(prof.Rnet,48,17,50);
%      
%      figure(588)
%      clf
%      plot(xx(use),y(use),'.')
%      xlabel('P_0')
%      ylabel('Rnet, W m-2')
%       title('CanVeg, Alfalfa, P0>0')
   
    
     
     
       %  profile of Tsun and Tshade
    
%       fnum= fnum+1;
%      figure(fnum)
%      clf
%      plot(nanmean(Sun.Tsfc(:,1:prm.jtot),1),prm.sumlai','+-','LineWidth', 1)
%      ax=gca;
%       set(ax, 'ydir','reverse');
%      hold on
%      plot(nanmean(Shade.Tsfc(:,1:prm.jtot),1),prm.sumlai','+-','LineWidth', 1)
%      ax=gca;
%       set(ax, 'ydir','reverse');
%        xlabel('Surface Temperature')
%       ylabel('Canopy Depth')
%       legend('sun','shade');
      
      
      
      % water use efficiency test
      
      
      IWUE=met.llambda .* Veg.Gpp ./Can.LE .* (met.vpd_Pa./met.P_Pa);
      IWUE_veg=met.llambda .* Veg.Gpp ./Veg.LE .* (met.vpd_Pa./met.P_Pa);
      
      iWUE=Veg.Gpp ./ Veg.gs;
      
      
      fnum= fnum+1;
      figure(fnum)
      clf
      plot(IWUE, iWUE, '.')
      xlabel('IWUE')
      ylabel('iWUE')
      title(VegType)
      ylim([0 2500])
      xlim([0 2500])
      
      
      
      fnum=fnum+1 ;
      figure(fnum)
      clf
      plot(IWUE, IWUE_veg, '.')
      xlabel('IWUE_{Can}')
      ylabel('IWUE_{veg}')
      title(VegType)
      ylim([0 2500])
      xlim([0 2500])
      
              
      
      [ GsPM ] = Gsfc_PM(met,Can) ; 
      
      
      iWUE_Gsfc=Veg.Gpp ./ GsPM;
        
      fnum=fnum+1 ;
      figure(fnum)
      clf
      plot(IWUE, iWUE_Gsfc, '.')
      xlabel('IWUE')
      ylabel('iWUE_Gsfc')
      title(VegType)
      ylim([0 2500])
      xlim([0 2500])
      
      
       fnum=fnum+1 ;
      figure(fnum)
      clf
      plot(iWUE, iWUE_Gsfc, '.')
      xlabel('iWUE_{Gstom}')
      ylabel('iWUE_{Gsfc}')
      title(VegType)
      ylim([0 2500])
      xlim([0 2500])
     
       fnum=fnum+1 ;
      figure(fnum)
      clf
      plot(met.LongIn,LongIn,'.')
      ylabel('IRin, measured, W m-2')
      xlabel('IRin, calculated, W m-2')
      
      
       fnum=fnum+1 ;
      figure(fnum)
      clf
      plot(nanmean(reshape(LongIn,48,17),2),'.');
      
      
     fnum=fnum+1 ;
     figure(fnum)
     clf
     histogram( met.LongIn-LongIn)
     xlabel('IRin calc - IRin measured W m-2')
  
end

