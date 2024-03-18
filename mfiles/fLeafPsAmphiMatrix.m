 function [ps]=fLeafPsAmphiMatrix(Iphoton,cca, Tlk, rb_co2, P_kPa, eair_Pa,prm)
 

%         8-31-2021, Leaf_photosynthesis model, Amphistomatous leaves   

% converting to matrix


%         DENNIS BALDOCCHI
%         Ecosystem Science Division
%         Department of Environmental Science, Policy and Management
%         & Berkeley Atmospheric Science Center
%         345 Hilgard Hall
%         University of California, Berkeley
%         Berkeley, CA 94720-3110
%         baldocchi@berkeley.edu
%         510-642-2874
       
%         
%          This program solves a cubic equation to calculate
%          leaf photosynthesis.  This cubic expression is derived from solving
%          five simultaneous equations for A, PG, cs, CI and GS.
%          Stomatal conductance is computed with the Ball-Berry model.
%          The cubic derivation assumes that b', the intercept of the Ball-Berry
%          stomatal conductance model, is non-zero.
% 
%           Gs = k A rh/cs + b'
% 
%           We also found that the solution for A can be obtained by a quadratic equation
%           when Gs is constant or b' is zero.
% 
% 
%             The derivation is published in:
% 
%             Baldocchi, D.D. 1994. An analytical solution for coupled leaf photosynthesis
%             and stomatal conductance models. Tree Physiology 14: 1069-1079.
% 
% 
% -----------------------------------------------------------------------
% 
%           A Biochemical Model of C3 Photosynthesis
% 
%             After Farquhar, von Caemmerer and Berry (1980) Planta.
%             149: 78-90.
% 
%         The original program was modified to incorporate functions and parameters
%         derived from gas exchange experiments of Harley, who paramertized Vc and J in
%         terms of optimal temperature, rather than some reference temperature, eg 25C.
% 
%         Program calculates leaf photosynthesis from biochemical parameters
% 
%         rd25 - Dark respiration at 25 degrees C (umol m-2 s-1)
%         tlk - leaf temperature, Kelvin
%         jmax - optimal rate of electron transport
%         vcopt - maximum rate of RuBP Carboxylase/oxygenase
%         iphoton - incident photosynthetically active photon flux (mmols m-2 s-1)
%    
%             note: Harley parameterized the model on the basis of incident PAR
% 
%         gs - stomatal conductance (mol m-2 s-1), typically 0.01-0.20
%         pstat-station pressure, bars
%         aphoto - net photosynthesis  (umol m-2 s-1)
%         ps - gross photosynthesis (umol m-2 s-1)
%         aps - net photosynthesis (mg m-2 s-1)
%         aphoto (umol m-2 s-1)
% 
% --------------------------------------------------
% 
%         iphoton is radiation incident on leaves
% 
%         The temperature dependency of the kinetic properties of
%         RUBISCO are compensated for using the Arrhenius and
%         Boltzmann equations.  From biochemistry, one observes that
%         at moderate temperatures enzyme kinetic rates increase
%         with temperature.  At extreme temperatures enzyme
%         denaturization occurs and rates must decrease.
% 
%         Arrhenius Eq.
% 
%         f(T)=f(tk_25) exp(tk -298)eact/(298 R tk)), where eact is the
%         activation energy.
% 
%         Boltzmann distribution
% 
%         F(T)=tboltz)
% 
% 
%         Define terms for calculation of gross photosynthesis, PG
% 
%         PG is a function of the minimum of RuBP saturated rate of
%         carboxylation, Wc, and the RuBP limited rate of carboxylation, Wj.
%         Wj is limiting when light is low and electron transport, which
%         re-generates RuBP, is limiting.  Wc is limiting when plenty of RuBP is
%         available compared to the CO2 that is needed for carboxylation.
% 
%         Both equations take the form:
% 
%         Vc - 0.5 Vo = PG-photorespiration= (a CI-a d)/(e CI + b)
% 
%         PG-photorespiration=min[Wj,Wc] (1-gamma/Ci)
% 
%         Wc=Vcmax Ci/(Ci + Kc(1+O2/Ko))
% 
%         Wj=J Ci/(4 Ci + 8 gamma)
% 
%         Ps kinetic coefficients from Harley at WBW.
% 
%         Gamma is the CO2 compensation point
% 
%      Information on the leaf photosynthetic parameters can be found in:
% 
%      Harley, P.C. and Baldocchi, 1995.Scaling carbon dioxide and water vapor exchange
%      from leaf to canopy in a deciduous forest:leaf level parameterization.
%      Plant, Cell and Environment. 18: 1146-1156.
% 
%      Wilson, K.B., D.D. Baldocchi and P.J. Hanson. 2000. Spatial and seasonal variability of
%      photosynthesis parameters and their relationship to leaf nitrogen in a deciduous forest.
%      Tree Physiology. 20, 565-587.
% 
% 
%      Tests of the model are reported in:
% 
%      Baldocchi, D.D. 1997. Measuring and modeling carbon dioxide and water vapor
%      exchange over a temperate broad-leaved forest during the 1995 summer drought.
%      Plant, Cell and Environment. 20: 1108-1122
% 
%      Baldocchi, D.D. and P.C. Harley. 1995. Scaling carbon dioxide and water vapor
%      exchange from leaf to canopy in a deciduous forest: model testing and application.
%      Plant, Cell and Environment. 18: 1157-1173.
% 
%      Baldocchi, D.D and T.P. Meyers. 1998. On using eco-physiological, micrometeorological
%      and biogeochemical theory to evaluate carbon dioxide, water vapor and gaseous deposition 
%      fluxes over vegetation. Agricultural and Forest Meteorology 90: 1-26.
%      
%      Baldocchi, D.D. Fuentes, J.D., Bowling, D.R, Turnipseed, A.A. Monson, R.K. 1999. Scaling
%      isoprene fluxes from leaves to canopies: test cases over a boreal aspen and a mixed species temperate
%      forest. J. Applied Meteorology. 38, 885-898.
% 
%      Baldocchi, D.D. and K.B.Wilson. 2001. Modeling CO2 and water vapor exchange of a
%      temperate broadleaved forest across hourly to decadal time scales. Ecological Modeling 
%           142: 155-184


%            Tlk=zeros(1,prm.jktot);
           
%          Arrhenius constants
%          Eact for Michaelis-Menten const. for KC, KO and dark respiration
%          These values are from Harley


%         ekc, Activation energy for K of CO2; J mol-1  
%         eko, Activation energy for K of O2, J mol-1
%         erd, Activation energy for dark respiration, eg Q10=2 
%         ektau, J mol-1 (Jordan and Ogren, 1984)
%         tk_25, absolute temperature at 25 C
%         toptvc, optimum temperature for maximum carboxylation
%         toptjm, optimum temperature for maximum electron transport
%             
%         kball, Ball-Berry stomatal coefficient for stomatal conductance
%         for water vapor.  For CO2 we need to divide by the ratios of the
%         diffusivities, 1.6

%         bprime, intercept of Ball-Berry model, mol m-2 s-1 
%         bprime16, intercept for CO2, bprime16 = bprime / 1.6;
%         bprimes applies to when Anet is zero, not when Agross is zero
%  
%         rsm, Minimum stomatal resistance, s m-1.
%         brs, curvature coeffient for light response
% 
%         qalpha, leaf quantum yield, electrons 
%         qalpha2, qalpha squared, qalpha2 = pow(qalpha, 2.0);

       
% debug
% found errors and inconsistancies in computing gross and net
% photosynthesis and application to stomatal conductance

% cleaned code and complete matrix modifications and removed for loops
        
% rank roots #1,#2 and #3 according to the minimum, intermediate and maximum


minroot=ones(prm.nn,prm.jtot) * 1e10;
midroot=zeros(prm.nn,prm.jtot);
maxroot=-ones(prm.nn,prm.jtot)* 1e10;
aphoto=zeros(prm.nn,prm.jtot);
Aps=zeros(prm.nn,prm.jtot);
root1=zeros(prm.nn,prm.jtot);
root2=zeros(prm.nn,prm.jtot);
root3=zeros(prm.nn,prm.jtot);
ps_1=zeros(prm.nn,prm.jtot);
delta_1=zeros(prm.nn,prm.jtot);

PkPa=ones(prm.nn,prm.jtot).*P_kPa;
           
        TlC=Tlk-273.15;
       
        rt = prm.rugc * Tlk;               %  product of universal gas constant and abs temperature

        tprime25 = Tlk - prm.tk_25;       % temperature difference

        ttemp = exp((prm.skin .* Tlk - prm.hkin) ./ rt) + 1.0; %  denominator term
        
    
%      KC and KO are solely a function of the Arrhenius Eq.


        [kct] = fTEMP_FUNC(prm.kc25, prm.ekc, tprime25, 298, Tlk);
        [ko] = fTEMP_FUNC(prm.ko25, prm.eko, tprime25, 298,Tlk);
        [tau]=Specificity(TlC);
        
        
        % fix the Ko and O2 values with same units
        
        ko25_Pa= prm.ko25* 100;  % Pa
        o2_Pa= prm.o2 * 101.3 ;  % Pa
        
        % bc = kct * (1.0 + o2 / ko);
        bc= kct * (1. + o2_Pa/ko25_Pa);

       

%         gammac is the CO2 compensation point due to photorespiration, umol mol-1
%         Recalculate gammac with the new temperature dependent KO and KC
%         coefficients..C at Vc = 0.5 Vo
%         gammac = O2/(2 tau)
%         O2 has units of kPa so multiplying 0.5 by 1000 gives units of Pa


          gammac = 500.0 * prm.o2 ./ tau;   % umol/mol
  
     
%         Scale rd with vcmax and apply temperature
%         correction for dark respiration

        rdzref= prm.vcopt * 0.004657;

        [rd] = fTEMP_FUNC(rdzref, prm.erd, tprime25, 298, Tlk);

        %       reduce respiration by 40% in light according to Amthor 
        
        

        rd(Iphoton > 10)=rd(Iphoton > 10) .* 0.4;
        rd(Iphoton <=10) = rd(Iphoton <=10);   
        
            
        % temperature correction to dark respiration
       
    
%        Apply temperature correction to JMAX and vcmax

%       Vcmax and Jmax were computed using Sharkey tool for 25 C
    

        ratio_jmax = fTBOLTZ(prm.jmopt, prm.ejm, prm.toptjm, Tlk)./...
            fTBOLTZ(prm.jmopt, prm.ejm, prm.toptjm, 298);
            
        jmax=ratio_jmax * prm.jmopt;
                
        ratio_vcmax = fTBOLTZ(prm.vcopt, prm.evc, prm.toptvc,Tlk) ./ ...
            fTBOLTZ(prm.vcopt, prm.evc, prm.toptvc,298);
        
        vcmax=prm.vcopt * ratio_vcmax;


%         T leaf boundary layer resistance
        
%         gb_mole leaf boundary layer conductance for CO2 exchange,
%         mol m-2 s-1         
%  
%         RB has units of s/m, convert to mol-1 m2 s1 to be
%         consistant with R.

      
% 

        rb_mole = rb_co2 .* Tlk .* 101.3* 0.022624 ./(273.15 .* P_kPa);

        gb_mole = 1. ./ rb_mole;

        dd = gammac;
        b8_dd = 8 * dd;


% ***************************************
% 
%         aphoto = PG - rd, net photosynthesis is the difference
%         between gross photosynthesis and dark respiration. Note
%         photorespiration is already factored into PG.
% 
% ****************************************

%         coefficients for Ball-Berry stomatal conductance model
% 
%         Gs = k A rh/cs + b'
% 
%         rh is relative humidity, which comes from a coupled
%         leaf energy balance model
% */

        % [rh_leaf]  = SFC_VPD(tlk, leleaf);
        
        rh_leaf=eair_Pa./fES(Tlk); % need to transpose matrix

        k_rh = rh_leaf .* prm.kball;  %// combine product of rh and K ball-berry


%         Gs from Ball-Berry is for water vapor.  It must be divided
%         by the ratio of the molecular diffusivities to be valid
%         for A and CO2 exchange

        k_rh_co2 = k_rh / 1.6;      % adjust the coefficient for the diffusion of CO2 rather than H2O

        gb_k_rh = gb_mole .* k_rh_co2;

        ci_guess = cca * .7;   % // initial guess of internal CO2 to estimate Wc and Wj

      %        // cubic coefficients that are only dependent on CO2 levels
      
      % set switch with prm.stomata
      
     
     
   
    switch(prm.stomata)
        case 'Hypo'
 
%  % Hypostomatous
          alpha_ps = 1.0 + (prm.bprime16 ./ gb_mole) - k_rh; % transpose gbmole
          bbeta = cca .* (gb_k_rh - 2.0 * prm.bprime16 - gb_mole);
          gamma = cca .* cca .* gb_mole * prm.bprime16;
          theta_ps = gb_k_rh - prm.bprime16;
          
 
        case ('Amphi')
% Amphistomatous       
         alpha_ps = 1.0 + (prm.bprime16 ./ (gb_mole)) - k_rh_co2;
         bbeta = cca .* (2*gb_k_rh - 3.0 * prm.bprime16 - 2*gb_mole);
         gamma = cca .* cca .* gb_mole * prm.bprime16 * 4;
         theta_ps = 2* gb_k_rh - 2 * prm.bprime16;
        otherwise
    end
%                   
%      Test for the minimum of Wc and Wj.  Both have the form:
% 
%         W = (a ci - ad)/(e ci + b)
% 
%         after the minimum is chosen set a, b, e and d for the cubic solution.
% 
%         estimate of J according to Farquhar and von Cammerer (1981)
% 
% 
%         J photon from Harley

%        looking at Collatz these W functions should be in terms of partial
%        pressure, units of Pa
% */

              
        
        j_photon = prm.qalpha .* Iphoton ./...
            sqrt(1. +(prm.qalpha2 * Iphoton .* Iphoton ./ (jmax .* jmax)));
        
       
    
        wj = j_photon .* (ci_guess - dd) ./ (4. .* ci_guess + b8_dd);
        wc = vcmax .* (ci_guess - dd) ./ (ci_guess + bc);

      
        B_ps=zeros(prm.nn,prm.jtot);
        a_ps=zeros(prm.nn,prm.jtot);
        E_ps=zeros(prm.nn,prm.jtot);
        psguess=zeros(prm.nn,prm.jtot);
        
        psguess(wj < wc)=wj(wj < wc);
        B_ps(wj < wc) = b8_dd(wj < wc);
        a_ps(wj < wc) = j_photon(wj < wc);
        E_ps(wj < wc) = 4.0;
        
        psguess(wj>=wc)=wc(wj>=wc);
        B_ps(wj>=wc) = bc(wj>=wc);
        a_ps(wj>=wc) = vcmax(wj>=wc);
        E_ps(wj>=wc) = 1.0;
        
       
%         /*
%         cubic solution:
% 
%          A^3 + p A^2 + q A + r = 0
%         */

% need to transpose some of the  matrices to keep 1 by nlayers

        denom = E_ps .* alpha_ps;

        Pcube = (E_ps .* bbeta + B_ps .* theta_ps - a_ps .* alpha_ps + E_ps .* rd .* alpha_ps);
        Pcube = Pcube./denom;

        Qcube = (E_ps .* gamma + (B_ps .* gamma ./ cca) - a_ps .* bbeta + a_ps .* dd .* theta_ps + E_ps .* rd .* bbeta + rd .* B_ps .* theta_ps);
        Qcube = Qcube./denom;

        Rcube = -a_ps .* gamma + a_ps .* dd .* gamma ./ cca + E_ps .* rd .* gamma + rd .* B_ps .* gamma ./ cca;
        Rcube = Rcube./denom;

       
        
        
        
%         // Use Cubic solution from Numerical Recipes from Press


        P2 = Pcube .* Pcube;
        P3 = P2 .* Pcube;
        Q = (P2 - 3.0 .* Qcube) ./ 9.0;
        R = (2.0 .* P3 - 9.0 .* Pcube .* Qcube + 27.0 .* Rcube) / 54.0;


% /*
%         Test = Q ^ 3 - R ^ 2
%         if test >= O then all roots are real
% */

                    rr=R.*R;
                    qqq=Q.*Q.*Q;
                    
                   tstroots=qqq-rr;       
        
%                 // real roots 

       
                        arg_U = R ./ power(qqq,0.5);

                        ang_L = acos(arg_U);
                        
                        sqrtQ=sqrt(Q);
        
                        root1 = -2.0 .* sqrtQ(tstroots > 0) .* cos(ang_L(tstroots > 0) ./ 3.0) - Pcube(tstroots > 0) ./ 3.0;  
                        root2 = -2.0 .* sqrtQ(tstroots > 0) .* cos((ang_L(tstroots > 0) + pi.*2) ./ 3.0) - Pcube(tstroots > 0) ./ 3.0;
                        root3 = -2.0 .* sqrtQ(tstroots > 0) .* cos((ang_L(tstroots > 0) - pi.*2) ./ 3.0) - Pcube(tstroots > 0) ./ 3.0;
                        
                        
                        % use Matlab calls to find mins, maxs and mids
                        root=[root1,root2,root3];
                        
                        minroot=min(root,[],2);
                        maxroot=max(root,[],2);
                        tst=(root ~=maxroot & root ~= minroot);
                        midroot=sum(tst .*root,2);
                          

% use logical loops


       % // find out where roots plop down relative to the x-y axis
       
            
        
             tst=(minroot >0 & midroot > 0 & maxroot >0);
             %aphoto(tst)=minroot(tst);
              Aps(tst)=minroot(tst);
       
       
             tst=(minroot  < 0 & midroot < 0 & maxroot > 0);
            %aphoto(tst)=maxroot(tst);  
            Aps(tst)=maxroot(tst);    
         

             tst=(minroot  < 0 & midroot >  0 & maxroot >0);
%              aphoto(tst)=midroot(tst);  
               Aps(tst)=midroot(tst);  
             
               

   
        % think I found error and this should be Aps, not aphoto
        % aphoto = Aps-rd
        %  Aps=reshape(root3,prm.nn,prm.jtot);
            
           aphoto = Aps-rd;
            

% /*
%          Here A = x - p / 3, allowing the cubic expression to be expressed
%          as: x^3 + ax + b = 0
% */

     
% 
% /*
%         Ttest for sucrose limitation of photosynthesis, as suggested by
%         Collatz.  Js=Vmax/2
% */
       
      


        j_sucrose = vcmax / 2. -rd;
        
        tst=(j_sucrose < aphoto);
        aphoto(tst)=j_sucrose(tst);
        
       % surface CO2

       cs = cca - aphoto ./ gb_mole;
       rh_leaf(rh_leaf>1)=1;

       
       % this is the stomatal conductance for CO2
       % Ball Berry is computed based on Anet, as parameterized in the
       % field
              
        gs_leaf_mole = (prm.kball .* rh_leaf .* aphoto ./ cs)/1.6 + prm.bprime16;

        
        gs_co2 = gs_leaf_mole;


%         stomatal conductance is mol m-2 s-1
%         convert back to resistance (s/m) for energy balance routine and
%         water


        gs_m_s = 1.6 .*gs_leaf_mole .* Tlk .* 101.3 .* 0.022624 ./(273.15 * P_kPa);
        
        ps.rstom=1./gs_m_s;  
        
       

%         // to compute ci, Gs must be in terms for CO2 transfer


        ci = cs - aphoto ./ gs_co2;
        
         % recompute wj and wc with ci
        
        wj = j_photon .* (ci - dd) ./ (4. * ci + b8_dd);

        wc = vcmax .* (ci - dd) ./ (ci + bc);
        
  % Collatz uses a quadratic model to compute a dummy variable wp to allow
  % for the transition between wj and wc, when there is colimitation.  this
  % is important because if one looks at the light response curves of the
  % current code one see jumps in A at certain Par values
  
  % theta wp^2 - wp (wj + wc) + wj wc = 0
  
  % a x^2 + b x + c = 0
  
  % x = [-b +/- sqrt(b^2 - 4 a c)]/2a

       a=0.98;
       b= -(wj +wc);
       c=wj.*wc;
       
       wp1=(-b + sqrt(b.*b - 4.*a.*c))./(2*a);
       wp2=(-b - sqrt(b.*b - 4.*a.*c))./(2*a);
       
       wp=min(wp1,wp2);

% beta A^2 - A (Jp+Js) + JpJs = 0

       aa = 0.95;
       bb= -(wp + j_sucrose);
       cc = wp .* j_sucrose;
       
       
       Aps1=(-bb + sqrt(bb.*bb - 4*aa.*cc))./(2*aa);
       Aps2=(-bb - sqrt(bb.*bb - 4*aa.*cc))./(2*aa);
       
       Aps=min(Aps1,Aps2);
       
       tst=(Aps < aphoto & Aps > 0);
       
       aphoto(tst)=Aps(tst)-rd(tst);
       
        gs_leaf_mole(tst) = (prm.kball .* rh_leaf(tst) .* aphoto(tst) ./ cs(tst))/1.6 + prm.bprime16;


%         // convert Gs from vapor to CO2 diffusion coefficient


        gs_co2(tst) = gs_leaf_mole(tst);

% /*
%         stomatal conductance is mol m-2 s-1
%         convert back to resistance (s/m) for energy balance routine and
%         water
% */

       % gs_m_s(tst) = 1.6 .*gs_leaf_mole(tst) .* Tlk(tst) .* 101.3 .* 0.022624 ./(273.15 * P_kPa);
        gs_m_s(tst) = 1.6 .*gs_leaf_mole(tst) .* Tlk(tst) .* 101.3 .* 0.022624 ./(273.15 .* PkPa(tst));
        
        ps.rstom(tst)=1./gs_m_s(tst);  % is rstom effectively the hypostomatous leaf resistance on both sides?
        
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
%        % respiration     
%        
       tst=(wj <=rd) | (wc <= rd);   %???
       %tst=(wj <=0) | (wc <= 0);
       
%       
% % %         a quadratic solution of A is derived if gs=ax, but a cubic form occurs
% % %         if gs =ax + b.  Use quadratic case when A is less than zero because gs will be
% % %         negative, which is nonsense
%       
           ps_1(tst) = cca(tst) .* gb_mole(tst) .* gs_co2(tst);
           delta_1(tst) = gs_co2(tst) + gb_mole(tst);
           denom(tst) = gb_mole(tst) .* gs_co2(tst);
%           
           Aquad1(tst) = delta_1(tst) .* E_ps(tst);
           Bquad1(tst) = -ps_1(tst) .* E_ps(tst) - a_ps(tst) .* delta_1(tst)...
               + E_ps(tst) .* rd(tst) .* delta_1(tst) - B_ps(tst) .* denom(tst);
%           
           Cquad1(tst) = a_ps(tst) .* ps_1(tst) - a_ps(tst) .* dd(tst) .* denom(tst)...
               - E_ps(tst) .* rd(tst) .* ps_1(tst) - rd(tst) .* B_ps(tst) .* denom(tst);
%  
           product(tst)=Bquad1(tst) .* Bquad1(tst) - 4.0 .* Aquad1(tst) .* Cquad1(tst);
%           
         % tst1=(product(tst) >= 0);
          sqrprod= sqrt(product(tst));
          aphoto(tst) = (-Bquad1(tst) - sqrprod) ./ (2.0 .* Aquad1(tst));
%          
%          
          %gs_leaf_mole(j,i) = prm.bprime16;
          
          gs_leaf_mole(tst) = prm.bprime16;
          %gs_leaf_mole(tst) = (prm.kball .* rh_leaf(tst) .* aphoto(tst) ./ cs(tst))/1.6 + prm.bprime16;
          
          gs_leaf_mole(gs_leaf_mole<0)= prm.bprime16;  % eliminate few conditions with negative conductance
          
          gs_co2(tst) = gs_leaf_mole(tst) ;
% % 
% % %  /*
% % %        stomatal conductance is mol m-2 s-1
% % %         convert back to resistance (s/m) for energy balance routine and
% % %         water vapor
% % % */
% % 
          gs_m_s(tst) = 1.6 .* gs_leaf_mole(tst) .* Tlk(tst) .* 101.3 .*0.022642 ./(PkPa(tst)*273.15) ;
%          
%          
          cs(tst) = cca(tst) - aphoto(tst) ./ gb_mole(tst);
          ci(tst) = cs(tst) - aphoto(tst) ./ gs_co2(tst);

 %%%%%%%%%%%%%%%%%%%
 
 %output
         
         ps.aphoto=real(aphoto);
         ps.ci=ci;
         ps.gs_co2=real(gs_co2);     % this is the effective hypostomatous resistance 
         ps.gs_m_s=real(gs_m_s);     % conductance for water vapor
         ps.wj=wj;
         ps.wc=wc;
         ps.wp=wp;
         ps.jsucrose=j_sucrose;
         ps.Ag=Aps;
         ps.x1 = root1;
         ps.x2 = root2;
         ps.x3 = root3;
         ps.p = Pcube;
         ps.q = Qcube;
         ps.r = Rcube;
         ps.rd = rd;
 end
       

 
   