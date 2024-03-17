 function [ps]=fLeafPsHypoMatrix(Iphoton,cca, Tlk, rb_co2, P_kPa, eair_Pa,prm)
 

%         8-31-2021, Leaf_photosynthesis model, 

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
%         PG-photorespiration= (a CI-a d)/(e CI + b)
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
%            
%          ps.aphoto=zeros(1,prm.jktot);
%          ps.ci=zeros(1,prm.jktot);
%          ps.gs_co2=zeros(1,prm.jktot);
%          ps.wj=zeros(1,prm.jktot);
%          ps.wc=zeros(1,prm.jktot);
%          ps.wp=zeros(1,prm.jktot);
%          ps.jsucrose=zeros(1,prm.jktot);
%          ps.Ag=zeros(1,prm.jktot);
%          ps.x1 = zeros(1,prm.jktot);
%          ps.x2 = zeros(1,prm.jktot);
%          ps.x3 = zeros(1,prm.jktot);
%          ps.p = zeros(1,prm.jktot);
%          ps.q = zeros(1,prm.jktot);
%          ps.r = zeros(1,prm.jktot);
%          ps.rd = zeros(1,prm.jktot);

%          gs_leaf_mole=zeros(prm.nn,prm.jtot);

         
           
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
%         bprime, intercept of Ball-Berry model, mol m-2 s-1 
%         bprime16, intercept for CO2, bprime16 = bprime / 1.6;
%  
%         rsm, Minimum stomatal resistance, s m-1.
%         brs, curvature coeffient for light response
% 
%         qalpha, leaf quantum yield, electrons 
%         qalpha2, qalpha squared, qalpha2 = pow(qalpha, 2.0);

       
        
% rank roots #1,#2 and #3 according to the minimum, intermediate and maximum


minroot=ones(prm.nn,prm.jtot) * 1e10;
midroot=zeros(prm.nn,prm.jtot);
maxroot=-ones(prm.nn,prm.jtot)* 1e10;
aphoto=zeros(prm.nn,prm.jtot);
root1=zeros(prm.nn,prm.jtot);
root2=zeros(prm.nn,prm.jtot);
root3=zeros(prm.nn,prm.jtot);



        
             
       

           
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
    

        jmax = fTBOLTZ(prm.jmopt, prm.ejm, prm.toptjm, Tlk);
        vcmax = fTBOLTZ(prm.vcopt, prm.evc, prm.toptvc,Tlk);


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
%         APHOTO = PG - rd, net photosynthesis is the difference
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
%         for A

        k_rh = k_rh / 1.6;      % adjust the coefficient for the diffusion of CO2 rather than H2O

        gb_k_rh = gb_mole .* k_rh;

        ci_guess = cca * .7;   % // initial guess of internal CO2 to estimate Wc and Wj

      %        // cubic coefficients that are only dependent on CO2 levels

 % Hypostomatous
        alpha_ps = 1.0 + (prm.bprime16 ./ gb_mole) - k_rh; % transpose gbmole
        bbeta = cca .* (gb_k_rh - 2.0 * prm.bprime16 - gb_mole);
        gamma = cca .* cca .* gb_mole * prm.bprime16;
        theta_ps = gb_k_rh - prm.bprime16;
%         
 % Amphistomatous       
%          alpha_ps = 1.0 + (prm.bprime16 ./ (gb_mole)) - k_rh;
%          bbeta = cca .* (2*gb_k_rh - 3.0 * prm.bprime16 - 2*gb_mole);
%          gamma = cca .* cca .* gb_mole * prm.bprime16 * 4;
%          theta_ps = 2* gb_k_rh - 2 * prm.bprime16;
        
                  
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
        
        
          

%         if wj or wc are less than rd then A would probably be less than zero.  This would yield a
%         negative stomatal conductance.  In this case, assume gs equals the cuticular value. This
%         assumptions yields a quadratic rather than cubic solution for A


        
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

       
        
        
        
%         // Use solution from Numerical Recipes from Press


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
        
                        root1 = -2.0 .* sqrtQ(tstroots > 0) .* cos(ang_L(tstroots > 0) / 3.0) - Pcube(tstroots > 0) / 3.0;  
                        root2 = -2.0 .* sqrtQ(tstroots > 0) .* cos((ang_L(tstroots > 0) + pi*2) / 3.0) - Pcube(tstroots > 0) / 3.0;
                        root3 = -2.0 .* sqrtQ(tstroots > 0) .* cos((ang_L(tstroots > 0) - pi*2) / 3.0) - Pcube(tstroots > 0) / 3.0;
                        
                        
                        
                        tst= (root1 <= root2 & root1 <=root3);
                        minroot(tst)=root1(tst);
                         
                        tst= (root1 <= root2 & root1 <=root3 & root2 <= root3);
                        midroot(tst)=root2(tst);
                        maxroot(tst)=root3(tst);
                        
                        tst= (root1 <= root2 & root1 <=root3 & root2 > root3);
                        midroot(tst)=root3(tst);
                        maxroot(tst)=root2(tst);

                        
                        tst= (root2 <= root1 & root2 <=root3);
                        minroot(tst)=root2(tst);
                         
                        tst= (root2 <= root1 & root2 <=root3 & root1 <= root3);
                        midroot(tst)=root1(tst);
                        maxroot(tst)=root3(tst);
                        
                        tst= (root2 <= root1 & root2 <=root3 & root1 > root3);
                        midroot(tst)=root3(tst);
                        maxroot(tst)=root1(tst);
                        
                        
                        tst= (root3 <= root1 & root3 <=root2);
                        minroot(tst)=root3(tst);
                         
                        tst= (root3 <= root1 & root3 <=root2 & root1 < root2);
                        midroot(tst)=root1(tst);
                        maxroot(tst)=root2(tst);
                        
                        tst= (root3 <= root1 & root3 <=root2 & root1 >= root2);
                        midroot(tst)=root2(tst);
                        maxroot(tst)=root1(tst);
                        

% use logical loops

% debug..  seems I am picking a 'bad' root..for night I was getting maxroot
% at 520 or so rather than -0.35
% check roots and notes from Eva F in canoak
%         for j=1:prm.nn
%         for i=1:prm.jtot
%                         if(root1(j,i) <= root2(j,i) && root1(j,i) <= root3(j,i))
%                         minroot(j,i)=root1(j,i);
%                                 if (root2(j,i) <= root3(j,i))
%                                 midroot(j,i)=root2(j,i);
%                                 maxroot(j,i)=root3(j,i);
%                                 else
%                                 midroot(j,i)=root3(j,i);
%                                 maxroot(j,i)=root2(j,i);
%                                 end
%                         end
% 
% 
%                         if(root2(j,i) <= root1(j,i) && root2(j,i) <= root3(j,i))
%                         minroot(j,i)=root2(j,i);
%                                 if (root1(j,i) <= root3(j,i))
%                                 midroot(j,i)=root1(j,i);
%                                 maxroot(j,i)=root3(j,i);
%                                 else
%                                 midroot(j,i)=root3(j,i);
%                                 maxroot(j,i)=root1(j,i);
%                                 end
%                         end
%                         
%                         if(root3(j,i) <= root1(j,i) && root3(j,i) <= root2(j,i))
%                         minroot(j,i)=root3(j,i);
%                                 if (root1(j,i) < root2(j,i))
%                                 midroot(j,i)=root1(j,i);
%                                 maxroot(j,i)=root2(j,i);
%                                 else
%                                 midroot(j,i)=root2(j,i);
%                                 maxroot(j,i)=root1(j,i);
%                                 end
% 
%                         end  %// end of the if loop for real roots
%                  
%                   

       % // find out where roots plop down relative to the x-y axis
            
        
             tst=(minroot >0 & midroot > 0 & maxroot >0);
             aphoto(tst)=minroot(tst);
       
%             if (minroot(j,i) > 0 && midroot(j,i) > 0 && maxroot(j,i) > 0)
%             aphoto(j,i)=minroot(j,i);
%             end
        
             tst=(minroot  < 0 & midroot < 0 & maxroot > 0);
             aphoto(tst)=maxroot(tst);        

%             if (minroot(j,i) < 0 && midroot(j,i) < 0 && maxroot(j,i) > 0)
%             aphoto(j,i)=maxroot(j,i);
%             end

             tst=(minroot  < 0 & midroot >  0 & maxroot >0);
             aphoto(tst)=midroot(tst);  
             
%             if (minroot(j,i) < 0 && midroot(j,i) > 0 && maxroot(j,i) > 0)
%             aphoto(j,i)=midroot(j,i);
%             end

%         end % end of i loop
%         end % end of j loop
            

% /*
%          Here A = x - p / 3, allowing the cubic expression to be expressed
%          as: x^3 + ax + b = 0
% */

%   aphoto=root3;        
% 
% /*
%         also test for sucrose limitation of photosynthesis, as suggested by
%         Collatz.  Js=Vmax/2
% */
       

      


j_sucrose = vcmax / 2. -rd;
        
        tst=(j_sucrose < aphoto);
        aphoto(tst)=j_sucrose(tst);
        
       % /*
%         Stomatal conductance for water vapor
% 
% 
%         forest are hypostomatous.
%         Hence we don't divide the total resistance
%         by 2 since transfer is going on only one side of a leaf.

      cs = cca - aphoto ./ gb_mole;
       rh_leaf(rh_leaf>1)=1;

       gs_leaf_mole = (prm.kball .* rh_leaf .* aphoto ./ cs) + prm.bprime;


%         // convert Gs from vapor to CO2 diffusion coefficient


        gs_co2 = gs_leaf_mole / 1.6;

% /*
%         stomatal conductance is mol m-2 s-1
%         convert back to resistance (s/m) for energy balance routine
% */

        gs_m_s = gs_leaf_mole .* Tlk * 101.3* 0.022624 ./(273.15 * P_kPa);
        
        ps.rstom=1./gs_m_s;  % is rstom effectively the hypostomatous leaf resistance on both sides?
        
     
        
       

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
       
%        for j=1:prm.nn
%        for i=1:prm.jtot
%        if(Aps(j,i) < aphoto(j,i) && Aps(j,i) > 0)
%            aphoto(j,i)=Aps(j,i) - rd(j,i);
%        end
%        end  % next i
%        end %next j
%        

             
       

      
                
       % respiration     
       
       % aphoto < 0
       
%        tst=(aphoto < 0);
%        
%         gs_leaf_mole(aphoto < 0) = prm.bprime;
%         gs_co2(aphoto < 0) = gs_leaf_mole(aphoto < 0) / 1.6;
%        
%          ps_1(aphoto < 0) = cca(j) * gb_mole(aphoto < 0) * gs_co2(aphoto < 0);
%          delta_1(aphoto < 0) = gs_co2(j,i) + gb_mole(aphoto < 0);
%          denom(aphoto < 0) = gb_mole(aphoto < 0) * gs_co2(aphoto < 0);

       
       %%% may need logical loop like above
       for j=1:prm.nn
       for i=1:prm.jtot
        
        if (wj(j,i) <=rd(j,i)) || (wc(j,i) <= rd(j,i))
     % if(aphoto(j,i) < 0)
            
        gs_leaf_mole(j,i) = prm.bprime;
        gs_co2(j,i) = gs_leaf_mole(j,i) / 1.6;

%  /*
%        stomatal conductance is mol m-2 s-1
%         convert back to resistance (s/m) for energy balance routine
% */

        gs_m_s(j,i) = gs_leaf_mole(j,i) .* Tlk(j,i) * 101.3*0.022642 ./(P_kPa(j)*273.15) ;
        
        

%         a quadratic solution of A is derived if gs=ax, but a cubic form occurs
%         if gs =ax + b.  Use quadratic case when A is less than zero because gs will be
%         negative, which is nonsense

            
         ps_1(j,i) = cca(j) * gb_mole(j,i) * gs_co2(j,i);
         delta_1(j,i) = gs_co2(j,i) + gb_mole(j,i);
         denom(j,i) = gb_mole(j,i) * gs_co2(j,i);

         Aquad1(j,i) = delta_1(j,i) * E_ps(j,i);
         Bquad1(j,i) = -ps_1(j,i) * E_ps(j,i) - a_ps(j,i) * delta_1(j,i)...
             + E_ps(j,i)* rd(j,i) * delta_1(j,i) - B_ps(j,i) * denom(j,i);
         
         Cquad1(j,i) = a_ps(j,i) * ps_1(j,i) - a_ps(j,i) * dd(j,i) * denom(j,i)...
             - E_ps(j,i) * rd(j,i) * ps_1(j,i) - rd(j,i) * B_ps(j,i) * denom(j,i);

         product(j,i)=Bquad1(j,i) * Bquad1(j,i) - 4.0 * Aquad1(j,i) * Cquad1(j,i);

       
        if (product(j,i) >= 0)
        sqrprod= sqrt(product(j,i));
        aphoto(j,i) = (-Bquad1(j,i) - sqrprod) / (2.0 * Aquad1(j,i));
        end 
        

       



% /*
%          Tests suggest that APHOTO2 is the correct photosynthetic root when
%          light is zero because root 2, not root 1 yields the dark respiration
%          value rd.
% */

        cs(j,i) = cca(j) - aphoto(j,i) / gb_mole(j,i);
        ci(j,i) = cs(j,i) - aphoto(j,i) / gs_co2(j,i);
       
         end  % if statement
     
     
       end   % for i loop
       end % end j loop
        
       
       

      
      
        
       
       
       % debug cs -> 0 gs -> inf in second iteration
       
         
        ps.aphoto=real(aphoto);
         ps.ci=ci;
         ps.gs_co2=gs_co2;     % this is the effective hypostomatous resistance ?
         ps.gs_m_s=gs_m_s;
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
       

 
 
%  
%   void PHOTOSYNTHESIS(double Iphoton,double *rstompt, double zzz,double cca,double tlk,
%         double *leleaf, double *A_mgpt, double *resppt, double *cipnt, 
%                 double *wjpnt, double *wcpnt)
%         {
% 
% /*
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
%         gs - stomatal conductance (mols m-2 s-1), typically 0.01-0.20
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
%             Boltzmann distribution
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
%         PG-photorespiration= (a CI-a d)/(e CI + b)
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
% 
%         Jan 14, 1999 Updated the cubic solutions for photosynthesis.  There are
%         times when the restriction that R^2 < Q^3 is violated.  I therefore need
%         alternative algorithms to solve for the correct root.
%  
% ===============================================================

% 
%         double rr,qqq, minroot, maxroot, midroot;
%         
%         rt = rugc * tlk;                // product of universal gas constant and abs temperature
% 
%         tprime25 = tlk - tk_25;       // temperature difference
% 
%         ttemp = exp((skin * tlk - hkin) / rt) + 1.0;  // denominator term
% 
% 		// initialize min and max roots
% 
% 		minroot= 1e10;
% 		maxroot=-1e10;
% 		midroot=0;
% 		root1=0;
% 		root2=0;
% 		root3=0;
% 		aphoto=0;
% 
% 
%         // KC and KO are solely a function of the Arrhenius Eq.
% 
% 
%         kct = TEMP_FUNC(kc25, ekc, tprime25, tk_25, tlk);
%         ko = TEMP_FUNC(ko25, eko, tprime25, tk_25,tlk);
%         tau = TEMP_FUNC(tau25, ektau, tprime25, tk_25,tlk);
% 
%         bc = kct * (1.0 + o2 / ko);
% 
%         if(Iphoton < 1)
%         Iphoton = 0;
% 
% /*
%         gammac is the CO2 compensation point due to photorespiration, umol mol-1
%         Recalculate gammac with the new temperature dependent KO and KC
%         coefficients
% 
%         gammac = .5 * O2*1000/TAU 
% */
% 
%         gammac = 500.0 * o2 / tau;
% 
% /*
%         temperature corrections for Jmax and Vcmax
% 
%         Scale jmopt and VCOPT with a surrogate for leaf nitrogen
%         specific leaf weight (Gutschick and Weigel).
% 
%         normalized leaf wt is 1 at top of canopy and is 0.35
%         at forest floor.  Leaf weight scales linearly with height
%         and so does jmopt and vcmax
%         zoverh=0.65/HT=zh65
% 
% */
%        
%      
%         
%          
% 
%         // growing season, full Ps capacity  (note newer data by Wilson et al shows more
%         // dynamics
% 
%         
%         jmaxz = jmopt ;
%         vcmaxz = vcopt ;
% 
% 		//vcmaxz = leaf.Vcmax;
% 		
%        
% 
%         
% 		
% 
%         
% /*
%         Scale rd with height via vcmax and apply temperature
%         correction for dark respiration
% */
% 
%         rdz=vcmaxz * 0.004657;
% 
% 
%         // reduce respiration by 40% in light according to Amthor 
% 
% 
%         if(Iphoton > 10)
%         rdz *= 0.4;      
% 
%         rd = TEMP_FUNC(rdz, erd, tprime25, tk_25, tlk);
% 
%     
%         // Apply temperature correction to JMAX and vcmax
%     
% 
%         jmax = TBOLTZ(jmaxz, ejm, toptjm, tlk);
%         vcmax = TBOLTZ(vcmaxz, evc, toptvc, tlk);
% 
% /*
%         Compute the leaf boundary layer resistance
%         
%         gb_mole leaf boundary layer conductance for CO2 exchange,
%         mol m-2 s-1         
%  
%         RB has units of s/m, convert to mol-1 m2 s1 to be
%         consistant with R.
% 
%         rb_mole = RBCO2 * .0224 * 1.01 * tlk / (met.pstat * 273.16)
% */
% 
%         rb_mole = bound_layer_res.co2 * tlk * (met.pstat273);
% 
%         gb_mole = 1. / rb_mole;
% 
%         dd = gammac;
%         b8_dd = 8 * dd;
% 
% 
% /***************************************
% 
%         APHOTO = PG - rd, net photosynthesis is the difference
%         between gross photosynthesis and dark respiration. Note
%         photorespiration is already factored into PG.
% 
% ****************************************
% 
%         coefficients for Ball-Berry stomatal conductance model
% 
%         Gs = k A rh/cs + b'
% 
%         rh is relative humidity, which comes from a coupled
%         leaf energy balance model
% */
% 
%         rh_leaf  = SFC_VPD(tlk, zzz, leleaf);
% 
%         k_rh = rh_leaf * sfc_res.kballstr;  // combine product of rh and K ball-berry
% 
% /*
%         Gs from Ball-Berry is for water vapor.  It must be divided
%         by the ratio of the molecular diffusivities to be valid
%         for A
% */
%         k_rh = k_rh / 1.6;      // adjust the coefficient for the diffusion of CO2 rather than H2O
% 
%         gb_k_rh = gb_mole * k_rh;
% 
%         ci_guess = cca * .7;    // initial guess of internal CO2 to estimate Wc and Wj
% 
% 
%        // cubic coefficients that are only dependent on CO2 levels
% 
% 
%         alpha_ps = 1.0 + (bprime16 / gb_mole) - k_rh;
%         bbeta = cca * (gb_k_rh - 2.0 * bprime16 - gb_mole);
%         gamma = cca * cca * gb_mole * bprime16;
%         theta_ps = gb_k_rh - bprime16;
% 
% /*
%         Test for the minimum of Wc and Wj.  Both have the form:
% 
%         W = (a ci - ad)/(e ci + b)
% 
%         after the minimum is chosen set a, b, e and d for the cubic solution.
% 
%         estimate of J according to Farquhar and von Cammerer (1981)
% 
% 
%         J photon from Harley
% */
% 
%               
%         if (jmax > 0)
%         j_photon = qalpha * Iphoton / sqrt(1. +(qalpha2 * Iphoton * Iphoton / (jmax * jmax)));
%         else
%         j_photon = 0;
% 
% 
%         wj = j_photon * (ci_guess - dd) / (4. * ci_guess + b8_dd);
% 
% 
%         wc = vcmax * (ci_guess - dd) / (ci_guess + bc);
% 
% 
% 
%        
% 
% 
% 
%         if(wj < wc)
%         {
% 
%         // for Harley and Farquhar type model for Wj
% 
%         psguess=wj;
%         
%         B_ps = b8_dd;
%         a_ps = j_photon;
%         E_ps = 4.0;
%         }
%         else
%         {
%         psguess=wc;
%           
%         B_ps = bc;
%         a_ps = vcmax;
%         E_ps = 1.0;
%         }
% 
% /*
%         if wj or wc are less than rd then A would probably be less than zero.  This would yield a
%         negative stomatal conductance.  In this case, assume gs equals the cuticular value. This
%         assumptions yields a quadratic rather than cubic solution for A
% 
% */
% 
% 
% 		 // frost and end of leaf photosynthesis and respiration 
% 
%      //   if(time_var.days > time_var.leafdrop)
%      //   {
%      //     wj = 0;
%      //     j_photon = 0;
%      //     wc=0;
%      //     rd = 0;
%      //   }  
% 
%         if (wj <= rd)
%         goto quad;
% 
%         if (wc <= rd)
%         goto quad;
% 
%         /*
%         cubic solution:
% 
%          A^3 + p A^2 + q A + r = 0
%         */
% 
%         denom = E_ps * alpha_ps;
% 
%         Pcube = (E_ps * bbeta + B_ps * theta_ps - a_ps * alpha_ps + E_ps * rd * alpha_ps);
%         Pcube /= denom;
% 
%         Qcube = (E_ps * gamma + (B_ps * gamma / cca) - a_ps * bbeta + a_ps * dd * theta_ps + E_ps * rd * bbeta + rd * B_ps * theta_ps);
%         Qcube /= denom;
% 
%         Rcube = (-a_ps * gamma + a_ps * dd * (gamma / cca) + E_ps * rd * gamma + rd * B_ps * gamma / cca);
%         Rcube /= denom;
% 
% 
%         // Use solution from Numerical Recipes from Press
% 
% 
%         P2 = Pcube * Pcube;
%         P3 = P2 * Pcube;
%         Q = (P2 - 3.0 * Qcube) / 9.0;
%         R = (2.0 * P3 - 9.0 * Pcube * Qcube + 27.0 * Rcube) / 54.0;
% 
% 
% /*
%         Test = Q ^ 3 - R ^ 2
%         if test >= O then all roots are real
% */
% 
%         rr=R*R;
%         qqq=Q*Q*Q;
%         
%                 // real roots 
% 
%         
%                         arg_U = R / sqrt(qqq);
% 
%                         ang_L = acos(arg_U);
%         
%                         root1 = -2.0 * sqrt(Q) * cos(ang_L / 3.0) - Pcube / 3.0;  
%                         root2 = -2.0 * sqrt(Q) * cos((ang_L + PI2) / 3.0) - Pcube / 3.0;
%                         root3 = -2.0 * sqrt(Q) * cos((ang_L -PI2) / 3.0) - Pcube / 3.0;
% 
%                 // rank roots #1,#2 and #3 according to the minimum, intermediate and maximum
%                 // value
% 
% 
%                         if(root1 <= root2 && root1 <= root3)
%                         {minroot=root1;
%                                 if (root2 <= root3)
%                                 {midroot=root2;
%                                 maxroot=root3;}
%                                  else
%                                  {midroot=root3;
%                                  maxroot=root2;}
%                         }
% 
% 
%                         if(root2 <= root1 && root2 <= root3)
%                         {minroot=root2;
%                                 if (root1 <= root3)
%                                 {midroot=root1;
%                                 maxroot=root3;}
%                                  else
%                                  {midroot=root3;
%                                  maxroot=root1;}
%                         }
% 
% 
%                         if(root3 <= root1 && root3 <= root2)
%                                 {minroot=root3;
%                                 if (root1 < root2)
%                                 {midroot=root1;
%                                 maxroot=root2;}
%                                  else
%                                  {midroot=root2;
%                                 maxroot=root1;}
% 
%                 }  // end of the loop for real roots
% 
% 
%         // find out where roots plop down relative to the x-y axis
%             
%         
%         if (minroot > 0 && midroot > 0 && maxroot > 0)
%         aphoto=minroot;
%         
%         
%         if (minroot < 0 && midroot < 0 && maxroot > 0)
%         aphoto=maxroot;
% 
%    
%                 if (minroot < 0 && midroot > 0 && maxroot > 0)
%         aphoto=midroot;
% 
% /*
%          Here A = x - p / 3, allowing the cubic expression to be expressed
%          as: x^3 + ax + b = 0
% */
% 
%                 // aphoto=root3;  // back to original assumption      
% 
% /*
%         also test for sucrose limitation of photosynthesis, as suggested by
%         Collatz.  Js=Vmax/2
% */
%         j_sucrose = vcmax / 2. - rd;
% 
%         if(j_sucrose < aphoto)
%         aphoto = j_sucrose;
% 
%         cs = cca - aphoto / gb_mole;
% 
%         if(cs > 1000)
%         cs=input.co2air;
% 
% /*
%         Stomatal conductance for water vapor
% 
% 
%         forest are hypostomatous.
%         Hence we don't divide the total resistance
%         by 2 since transfer is going on only one side of a leaf.
% 
% 		alfalfa is amphistomatous...be careful on where the factor of two is applied
% 		just did on LE on energy balance
%     
% */
% 
%        gs_leaf_mole = (sfc_res.kballstr * rh_leaf * aphoto / cs) + bprime;
% 
% 
%         // convert Gs from vapor to CO2 diffusion coefficient
% 
% 
%         gs_co2 = gs_leaf_mole / 1.6;
% 
% /*
%         stomatal conductance is mol m-2 s-1
%         convert back to resistance (s/m) for energy balance routine
% */
% 
%         gs_m_s = gs_leaf_mole * tlk * met.pstat273;
% 
%         // need point to pass rstom out of subroutine 
%        
%         *rstompt = 1.0 / gs_m_s;
% 
% 
%         // to compute ci, Gs must be in terms for CO2 transfer
% 
% 
%         ci = cs - aphoto / gs_co2;
% 
% /*
%          if A < 0 then gs should go to cuticular value and recalculate A
%          using quadratic solution
% */
% 
% 
% 		// recompute wj and wc with ci
% 
%         
%         wj = j_photon * (ci - dd) / (4. * ci + b8_dd);
% 
%         wc = vcmax * (ci - dd) / (ci + bc);
%         
%  /* Collatz uses a quadratic model to compute a dummy variable wp to allow
%   for the transition between wj and wc, when there is colimitation.  this
%   is important because if one looks at the light response curves of the
%   current code one see jumps in A at certain Par values
%   
%    theta wp^2 - wp (wj + wc) + wj wc = 0
%    a x^2 + b x + c = 0
%    x = [-b +/- sqrt(b^2 - 4 a c)]/2a
% 
% 	*/
% 
% 		
% 
% 
% 
%        a=0.98;
%        b= -(wj +wc);
%        c=wj*wc;
%        
%        wp1=(-b + sqrt(b*b - 4*a*c))/(2*a);
%        wp2=(-b - sqrt(b*b - 4*a*c))/(2*a); 
% 	   
% 	   // wp = min (wp1,wp2);
% 
% 	   if(wp1 < wp2)
% 		    wp=wp1;
% 	   else
% 			wp=wp2;
% 
% 
% 
% // beta A^2 - A (Jp+Js) + JpJs = 0
% 
%        aa = 0.95;
%        bb= -(wp+ j_sucrose);
%        cc = wp* j_sucrose;
%        
%        
%        Aps1=(-bb + sqrt(bb*bb - 4*aa*cc))/(2*aa);
%        Aps2=(-bb - sqrt(bb*bb - 4*aa*cc))/(2*aa);
%        
%       // Aps=min(Aps1,Aps2);
% 
% 	   if(Aps1 < Aps2)
% 		   Aps=Aps1;
% 	   else
% 		   Aps = Aps2;
%        
%        if(Aps < aphoto && Aps > 0)
%            aphoto=Aps - rd;
%        
%        
% 
%         
%         if(aphoto <= 0.0)
%         goto quad;
% 
%         goto OUTDAT;
% 
%        
% 
%         // if aphoto < 0  set stomatal conductance to cuticle value   
% 
% quad:
% 
%        
%         gs_leaf_mole = bprime;
%         gs_co2 = gs_leaf_mole / 1.6;
% 
% /*
%         stomatal conductance is mol m-2 s-1
%         convert back to resistance (s/m) for energy balance routine
% */
% 
%         gs_m_s = gs_leaf_mole * tlk * (met.pstat273);
%         
%         // need point to pass rstom out of subroutine as a pointer
%         
%         *rstompt = 1.0 / gs_m_s;
% 
% 
% /*
%         a quadratic solution of A is derived if gs=ax, but a cubic form occurs
%         if gs =ax + b.  Use quadratic case when A is less than zero because gs will be
%         negative, which is nonsense
%    
% */
%             
%          ps_1 = cca * gb_mole * gs_co2;
%          delta_1 = gs_co2 + gb_mole;
%          denom = gb_mole * gs_co2;
% 
%          Aquad1 = delta_1 * E_ps;
%          Bquad1 = -ps_1 * E_ps - a_ps * delta_1 + E_ps * rd * delta_1 - B_ps * denom;
%          Cquad1 = a_ps * ps_1 - a_ps * dd * denom - E_ps * rd * ps_1 - rd * B_ps * denom;
% 
%          product=Bquad1 * Bquad1 - 4.0 * Aquad1 * Cquad1;
% 
%         if (product >= 0)
%         sqrprod= sqrt(product);
%          
%          aphoto = (-Bquad1 - sqrprod) / (2.0 * Aquad1);
% /*
%          Tests suggest that APHOTO2 is the correct photosynthetic root when
%          light is zero because root 2, not root 1 yields the dark respiration
%          value rd.
% */
% 
%         cs = cca - aphoto / gb_mole;
%         ci = cs - aphoto / gs_co2;
%        
% 
%         OUTDAT:
% 
% /*
%         compute photosynthesis with units of mg m-2 s-1 and pass out as pointers
%     
%         A_mg = APHOTO * 44 / 1000
% */
%       *A_mgpt = aphoto * .044;
% 
%       *resppt=rd;
% 
%       *cipnt=ci;  
% 
%       *wcpnt=wc;
% 
%       *wjpnt=wj;
%     
%    /*   
% 
%         printf(" cs       ci      gs_leaf_mole      CA     ci/CA  APS  root1  root2  root3\n");
%         printf(" %5.1f   %5.1f   %6.3f    %5.1f %6.3f  %6.3f %6.3f %6.3f  %6.3f\n", cs, ci, gs_leaf_mole, cca, ci / cca,aphoto,root1, root2, root3 );
% 
%    */
% 
%         return;
% }
% 
   