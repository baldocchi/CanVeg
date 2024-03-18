
function [energy]=fEnergy_Leaf_1D_Matrix(boundary_layer_res,Qin,met,prof,gs,prm)

%  [ans]=fEnergy_Leaf_Hypo_Matrix(boundary_layer_res,Qin.sun_abs,met,prof,Sun.gs,prm);
   
% Energy Balance subroutine
% leaf version, one side
% Matlab CanVeg

% D Baldocchi
% Nov 27,2020

% debug..CanVeg has factor of 2 larger coef for T..check and check.. I am
% explicity computing top and bottom energy balance, so I may need to drop
% factors of 2 for H and Lout, as I am getting too negative H and too cool
% Tsfc, reducing Rnet

% 	Initialize Parameters
	
% 	Incoming short and longwave radiation

   Q_In=Qin;  %net incoming short wave and incoming Longwave
   
   Rs=1 ./gs;  % stomatal resistance 
   
   % think about computing ke for top and bottom and if hypo or amphi..to
   % do
   
   ke=1 ./(boundary_layer_res.vapor + Rs);
   %ke=2 ./(boundary_layer_res.vapor + 2* Rs);  % amphistomatous leaf
   %ke=(gs .*boundary_layer_res.vapor) ./ (gs + boundary_layer_res.vapor);
   
   met.P_Pa=1000 .* met.P_kPa;
   
   tk2=met.T_air_K .* met.T_air_K;
   tk3=tk2 .* met.T_air_K;
   tk4=tk3 .* met.T_air_K;
   
   % Using derivation for 2 sided leaves
    
	%energy.Lout = prm.ep * prm.sigma * met.T_air_K .^4.0;
       
	llout = prm.epsigma .* tk4;
   
	
% 	Compute gas and model coefficient

    [es]=fES(met.T_air_K);    
    met.vpd_Pa=es-met.eair_Pa;
    
    var.air_density = met.P_kPa .* prm.Mair ./ (prm.rugc .* met.T_air_K); 

	[dest]=fdESdT(met.T_air_K);
	[d2est]=fd2ESdT(met.T_air_K);
    [llambda]=fLambda(met.T_air_K);
    
    
    % need to double check derivation of LE
   
	% lecoef=var.air_density*.622.*llambda./(met.P_Pa .*(boundary_layer_res.vapor + Rs));
    
     lecoef=var.air_density*.622.*llambda .*ke./ met.P_Pa ;
	
     hcoef=var.air_density.*prm.Cp ./boundary_layer_res.heat;
     hcoef2 = 2 * hcoef;
     
     
     qrad=Q_In;
	
     %  4.0 * ep * sigma
     repeat = hcoef + prm.epsigma4 .* tk3;
     
     Acoef=lecoef.*d2est./(2.*repeat);
     Acoef=Acoef/4;
    
     Bcoef=-(repeat) - lecoef .* dest ./ 2. + Acoef .* (-qrad / 2. + llout);

     
      Ccoef=Acoef .* (Q_In .* Q_In -2*Q_In .*llout + llout .*llout) +...
      lecoef.*(met.vpd_Pa.*repeat+dest.*(Q_In-llout));
     
      Ccoef = repeat .* lecoef .* met.vpd_Pa + lecoef .* dest .* (qrad / 2. - llout)...
      + Acoef .* ((qrad .* qrad) / 4. + llout .* llout - qrad .* llout);
  
   
  
     
    
    %  solve for LE
% 
%  a LE^2 + bLE + c = 0

le1 = (-Bcoef + (Bcoef .*Bcoef - 4. * Acoef .* Ccoef).^.5) ./ (2. * Acoef);
le2 = (-Bcoef - (Bcoef .* Bcoef - 4. * Acoef .* Ccoef).^.5) ./ (2. * Acoef);


% if getting complex numbers..need real root

energy.LE=real(le2);

% Derivation for one sided
% 
%  aT^2 + bT + c =0

aT = 12. * prm.epsigma .* met.T_air_K .* met.T_air_K + d2est .* lecoef / 2;  
bT = 8. * prm.epsigma .*tk3 + hcoef + dest .* lecoef;
cT =  -qrad +2* llout + lecoef .* met.vpd_Pa ;



% quadratic equation

%del_Tk=(-bT + (bT .* bT - 4. * aT .* cT).^.5) ./ (2. * aT);

% need alternative solution in how I derive Ts-Ta
del_Tk=(-bT - (bT .* bT - 4. * aT .* cT).^.5) ./ (2. * aT);

met.Tsfc_K= met.T_air_K + real(del_Tk);


% one sided leaf

%  H is sensible heat flux density

        energy.H = hcoef .* real(del_Tk);


% Lout is longwave emitted energy

        energy.Lout =prm.epsigma* met.Tsfc_K .^4.;
        energy.Qin=Q_In;
        energy.Rnet=Q_In-energy.Lout;
        energy.Tsfc=met.Tsfc_K;
        energy.esTair=es;
        [energy.esTsfc]=fES(met.Tsfc_K);
end

	

%   void ENERGY_BALANCE (double qrad, double *tsfckpt, double taa, double rhovva, double rvsfc,
%             double stomsfc, double *lept, double *H_leafpt, double *lout_leafpt)
% {
%     
%     
% /* 
%         ENERGY BALANCE COMPUTATION
% 
%         A revised version of the quadratic solution to the leaf energy balance relationship is used.
% 
%         Paw U, KT. 1987. J. Thermal Biology. 3: 227-233
% 
% 
%          H is sensible heat flux density on the basis of both sides of a leaf
%          J m-2 s-1 (W m-2).  Note KC includes a factor of 2 here for heat flux
%          because it occurs from both sides of a leaf.
% */
%        
%              
%         double est, ea, tkta, le2;
%         double tk2, tk3, tk4;
%         double dest, d2est;
%         double lecoef, hcoef, hcoef2, repeat, acoeff, acoef;
%         double bcoef, ccoef, product;
%         double atlf, btlf, ctlf,vpd_leaf,llout;
%         double ke;
% 
%              
%              
%         tkta=taa;   // taa is already in Kelvin
% 		        
%         est = ES(tkta);  //  es(T)  Pa 
%       
% 
%         // ea  = RHOA * TAA * 1000 / 2.165
% 
%          ea = 1000 * rhovva * tkta /2.170;   // vapor pressure above leaf, Pa rhov is kg m-3
% 
% 
% 
%         // Vapor pressure deficit, Pa
% 
% 
%         vpd_leaf = est - ea;
% 
%         if (vpd_leaf < 0.)
%         vpd_leaf = 0;
% 
% 
%         // Slope of the vapor pressure-temperature curve, Pa/C
%         // evaluate as function of Tk
% 
% 
%         dest = DESDT(tkta);
% 
% 
%         // Second derivative of the vapor pressure-temperature curve, Pa/C
%         // Evaluate as function of Tk
% 
% 
%          d2est = DES2DT(tkta);
% 
% 
%         // Compute products of air temperature, K
%         
%         tk2 = tkta * tkta;
%         tk3 = tk2 * tkta;
%         tk4 = tk3 * tkta;
% 
%        
% 
%         // Longwave emission at air temperature, W m-2
% 
% 
%         llout = epsigma * tk4;
% 
% /*
% 
%         Coefficient for latent heat flux
% 
%         Oaks evaporate from only one side. They are hypostomatous. 
%         Cuticle resistance is included in STOM.
% 
% */
% 
% 
%         ke = 1./ (rvsfc + stomsfc);  // hypostomatous
% 
% 		// ke = 2/ (rvsfc + stomsfc);  // amphistomatous
% 
%         lecoef = met.air_density * .622 * fact.latent * ke / met.press_Pa;
% 
% 
%         // Coefficients for sensible heat flux
% 
% 
%         hcoef = met.air_density*cp/bound_layer_res.heat;
%         hcoef2 = 2 * hcoef;
% 
%         
%        // The quadratic coefficients for the a LE^2 + b LE +c =0
% 
% 
%         repeat = hcoef + epsigma4 * tk3;
% 
%         acoeff = lecoef * d2est / (2. * repeat);
%         acoef = acoeff / 4.;
% 
%         bcoef = -(repeat) - lecoef * dest / 2. + acoeff * (-qrad / 2. + llout);
% 
%         ccoef = repeat * lecoef * vpd_leaf + lecoef * dest * (qrad / 2. - llout) + acoeff * ((qrad * qrad) / 4. + llout * llout - qrad * llout);
% 
% 
%         // LE1 = (-BCOEF + (BCOEF ^ 2 - 4 * ACOEF * CCOEF) ^ .5) / (2 * ACOEF) 
% 
%         product = bcoef * bcoef - 4. * acoef * ccoef;
% 
%         // LE2 = (-BCOEF - (BCOEF * BCOEF - 4 * acoef * CCOEF) ^ .5) / (2. * acoef) 
% 
%         
%         le2= (-bcoef - sqrt(product)) / (2. * acoef);
%   
%         *lept=le2;  // need to pass pointer out of subroutine
% 
% 
%         // solve for Ts using quadratic solution  
% 
%       
%         // coefficients to the quadratic solution
%         
%         atlf = epsigma12 * tk2 + d2est * lecoef / 2.;
% 
%         btlf = epsigma8 * tk3 + hcoef2 + lecoef * dest;
% 
%         ctlf = -qrad + 2 * llout + lecoef * vpd_leaf;
% 
% 
%         // IF (BTLF * BTLF - 4 * ATLF * CTLF) >= 0 THEN 
% 
%         product = btlf * btlf - 4 * atlf * ctlf;
% 
%        
%         // T_sfc_K = TAA + (-BTLF + SQR(BTLF * BTLF - 4 * ATLF * CTLF)) / (2 * ATLF) 
% 
%         if (product >= 0)
%         *tsfckpt = tkta + (-btlf + sqrt(product)) / (2 * atlf);
%         else
%         *tsfckpt=tkta;
% 
% 
%        if(*tsfckpt < -230 || *tsfckpt > 325)
%        *tsfckpt=tkta;
% 
%         // long wave emission of energy
% 
%         *lout_leafpt =epsigma2*pow(*tsfckpt,4);
% 
%         // H is sensible heat flux 
% 
%         *H_leafpt = hcoef2 * (*tsfckpt- tkta);
% 
% 
% return;
% }
% 
% 
% 
%   
