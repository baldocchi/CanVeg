function [boundary_layer_res] = fBoundary_Resistance_Matrix(prof, met,TLF, prm)

% 10/20/2020  

% 	BOUNDARY_RESISTANCE
% 
% 	This subroutine computes the leaf boundary layer
% 	resistances for heat, vapor and CO2 (s/m).
% 
% 	Flat plate theory is used, as discussed in Schuepp (1993) and
% 	Grace and Wilson (1981).
% 
% 	We consider the effects of turbulent boundary layers and sheltering.
% 	Schuepp's review shows a beta factor multiplier is necessary for SH in
% 	flows with high turbulence.  The concepts and theories used have been
% 	validated with our work on HNO3 transfer to the forest.
% 
% 
% 	Schuepp. 1993 New Phytologist 125: 477-507
% 
% 
% 	Diffusivities have been corrected using the temperature/Pressure algorithm in Massman (1998)


      boundary_layer_res.heat=zeros(prm.nn,prm.jtot);
      boundary_layer_res.vapor=zeros(prm.nn,prm.jtot);
      boundary_layer_res.co2=zeros(prm.nn,prm.jtot);
      
      Sh_heat=zeros(prm.nn,prm.jtot);
      Sh_vapor=zeros(prm.nn,prm.jtot);
      Sh_CO2=zeros(prm.nn,prm.jtot);

      Tref = TLF; 
   
    deltlf= Tref- prof.Tair_K(:,1:prm.jtot);  % make sure K

    graf = prm.grasshof .* deltlf ./ prof.Tair_K(:,1:prm.jtot); 
    
    graf(deltlf< 0)=0;
    
        
    nnu_T_P = prm.nnu .* (101.3 ./ met.P_kPa).*(prof.Tair_K(:,1:prm.jtot) ./ 273.16).^1.81;
    
   % tst=UZ(indat,zzz(1:prm.jktot),met,prm);
   
   % compute profile of UZ
    
    Re = prm.lleaf .* prof.wind(:,1:prm.jtot) ./ nnu_T_P;

    Re5 = Re .^0.5;
	Re8 = Re .^0.8;

%	if (Re > 14000.)
	
% 		turbulent boundary layer
 
% 		SH = .036 * Re8 * pr33*betfact;
% 		SHV = .036 * Re8 * sc33*betfact;
% 		SHCO2 = .036 * Re8 * scc33*betfact;

        Res_factor = 0.036.*Re8.*prm.betfac;

		Sh_heat(Re > 14000.) = Res_factor(Re > 14000.) * prm.pr33;
		Sh_vapor(Re > 14000.) = Res_factor(Re > 14000.)* prm.sc33;
		Sh_CO2(Re > 14000.) = Res_factor(Re > 14000.) * prm.scc33;
 %   end
	
	% laminar
	
  %  if (Re <= 14000)

		Res_factor_lam = 0.66*Re5*prm.betfac;
		
% 		laminar sublayer
 
% 		SH = .66 * Re5 * pr33*betfact;
% 		SHV = .66 * Re5 * sc33*betfact;
% 		SHCO2 = .66 * Re5 * scc33*betfact;
	

		Sh_heat(Re <= 14000) = Res_factor_lam(Re <= 14000) * prm.pr33;
		Sh_vapor(Re <= 14000) = Res_factor_lam(Re <= 14000) * prm.sc33;
		Sh_CO2(Re <= 14000) = Res_factor_lam(Re <= 14000) * prm.scc33;
  %  end
   
	%If there is free convection

     tsst = graf ./ (Re .* Re);
    
%     if (tsst > 1.)
% 	
% 
% 		%     Compute Grashof number for free convection
% 
% 		if (graf < 100000.)
% 			GR25 = .5 * graf  .^ 0.25;
%         end
%         
% 		if (graf >= 100000.)
% 			GR25 = .13 .* graf.^ 0.33;
%         end
% 
%         Sh_heat = prm.pr33 .* GR25;
% 		Sh_vapor = prm.sc33 .* GR25;
% 		Sh_CO2 = prm.scc33 .* GR25;
% 
%     end
    
     tstgr=(tsst > 1 & graf < 100000);
     GR25 = .5 * graf  .^ 0.25;
     Sh_heat(tstgr) = prm.pr33 .* GR25(tstgr);
     Sh_vapor(tstgr) = prm.sc33 .* GR25(tstgr);
	 Sh_CO2(tstgr) = prm.scc33 .* GR25(tstgr);
     
     tstgra=(tsst > 1 & graf >= 100000);
     GR33 = .13 .* graf.^ 0.33;
     Sh_heat(tstgra) = prm.pr33 .* GR33(tstgra);
     Sh_vapor(tstgra) = prm.sc33 .* GR33(tstgra);
	 Sh_CO2(tstgra) = prm.scc33 .* GR33(tstgra);
     
    %// lfddx=lleaf/ddx


	%// Correct diffusivities for temperature and pressure

	ddh_T_P = prm.ddh .* (101.3 ./ met.P_kPa).*power((prof.Tair_K(:,1:prm.jtot) / 273.16), 1.81);
	ddv_T_P = prm.ddv .* (101.3 ./ met.P_kPa).*power((prof.Tair_K(:,1:prm.jtot) / 273.16), 1.81);
	ddc_T_P = prm.ddc .* (101.3 ./ met.P_kPa).*power((prof.Tair_K(:,1:prm.jtot) / 273.16), 1.81);

	boundary_layer_res.heat = (prm.lleaf ./ (ddh_T_P .* Sh_heat)); 
	boundary_layer_res.vapor = (prm.lleaf ./ (ddv_T_P .* Sh_vapor));
	boundary_layer_res.co2 = (prm.lleaf ./ (ddc_T_P .* Sh_CO2));



end  % end of function

	

	

