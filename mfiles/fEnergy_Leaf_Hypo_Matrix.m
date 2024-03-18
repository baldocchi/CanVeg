
function [energy]=fEnergy_Leaf_Hypo_Matrix(boundary_layer_res,Q_In,met,prof,gs,prm)


% Energy Balance subroutine
% leaf version, Hypostomatous
% Matlab CanVeg, Matrix form

% Energy and heat exchange on both sides of leaf. Transpiration and
% photosynthesis are on only one side of leaf

% D Baldocchi
% Nov 28,2020


% 	Initialize Parameters
	
% 	Incoming short and longwave radiation

  % Q_In is net incoming minus outgoin short wave and incoming Longwave
  % level above and below each layer
  
  energy.LE=zeros(prm.nn,prm.jtot);
  energy.H=zeros(prm.nn,prm.jtot);
  energy.Rnet=zeros(prm.nn,prm.jtot);
  energy.Lout=zeros(prm.nn,prm.jtot);
  Acoef=zeros(prm.nn,prm.jtot);
  Bcoef=zeros(prm.nn,prm.jtot);
  Ccoef=zeros(prm.nn,prm.jtot);
       
     
   gb=1./boundary_layer_res.vapor;
        
  % gw=1./rwater;  
  
 
  gw = (gb .* gs) ./(gb + gs);
  
    
   met.P_Pa=1000 .* met.P_kPa;   % air pressure, Pa
   
%    tk2=met.T_air_K .* met.T_air_K;
%    tk3=tk2 .* met.T_air_K;
%    tk4=tk3 .* met.T_air_K;

    tk2=prof.Tair_K(:,1:prm.jtot) .* prof.Tair_K(:,1:prm.jtot);
    tk3=tk2 .* prof.Tair_K(:,1:prm.jtot);
    tk4=tk3 .* prof.Tair_K(:,1:prm.jtot);

   
   % Using derivation for 2 sided leaves
    	       
	llout = prm.epsigma .* tk4;
   
	 
    % Hypostomatous

     d2est=fd2ESdT(prof.Tair_K(:,1:prm.jtot));  
     dest=fdESdT(prof.Tair_K(:,1:prm.jtot)); 
     est=fES(prof.Tair_K(:,1:prm.jtot));
     vpd_Pa=est-prof.eair_Pa(:,1:prm.jtot);
     llambda=fLambda(prof.Tair_K(:,1:prm.jtot));
     air_density = met.P_kPa .* prm.Mair ./ (prm.rugc .* prof.Tair_K(:,1:prm.jtot));  % kg m-3
     
     vpd_Pa(vpd_Pa < 0)=0;
     vpd_Pa(vpd_Pa > 5000)=5000;
    
    % gas, heat and thermodynamic coeficients
    
	 lecoef=air_density .*.622 .* llambda .*gw ./ met.P_Pa ;
     hcoef=air_density.*prm.Cp ./boundary_layer_res.heat;
     hcoef2 = 2 * hcoef;
     
     repeat = hcoef2 + prm.epsigma8 .* tk3;
     
     llout2=2*llout;  % longwave energy times 2
     
    % coefficients analytical solution
    
     Acoef1=lecoef.*d2est./(2.*repeat);
     Acoef=Acoef1/2;
    
     Bcoef=-(repeat) - lecoef .* dest  + Acoef .* (-2*Q_In + 4 .*llout);

     
      Ccoef=Acoef .* (Q_In .* Q_In -4*Q_In .*llout + 4 .*llout .*llout) +...
      lecoef.*(met.vpd_Pa .* repeat +dest .*(Q_In- 2.*llout));
     
           
    
    %  solve for LE
% 
%  a LE^2 + bLE + c = 0

% solve for both roots, but LE tends to be second root, le2

le1 = (-Bcoef + (Bcoef .*Bcoef - 4. * Acoef .* Ccoef).^.5) ./ (2. .* Acoef);
le2 = (-Bcoef - (Bcoef .* Bcoef - 4. * Acoef .* Ccoef).^.5) ./ (2. .* Acoef);


energy.LE=real(le2);


del_Tk = (Q_In-energy.LE-llout2) ./repeat;

Tsfc_K= prof.Tair_K(:,1:prm.nlayers) + (del_Tk);


% I was not getting good energy balance closure solving for dT separately
% as there were 2nd order non linear terms that are not well considered.
% Using analytical solurion for LE above and using that to solve for dT
% works better
% quadratic equation



%  H is sensible heat flux density from both sides of leaf

      energy.H = hcoef2 .* real(del_Tk);


% Lout is longwave emitted energy from both sides of leaves

        energy.Lout =2 .* prm.epsigma* Tsfc_K .^4.;
        energy.Qin=Q_In;
        energy.Rnet=Q_In-energy.Lout;  % net radiation as a function of longwave energy emitted
        energy.Tsfc=Tsfc_K;
       [energy.esTsfc]=fES(Tsfc_K);
        energy.vpd_Pa=vpd_Pa;
        energy.closure=energy.Rnet-energy.H-energy.LE;  % test energy balance closure
  
        %used to debug
        
%         figure(5)
%         clf
%         plot(energy.Rnet(:,1:prm.nlayers),energy.LE(:,1:prm.nlayers)+energy.H(:,1:prm.nlayers),'.')
%         xlabel('Rnet W m-2')
%         ylabel('H + LE W m-2')
        

end

	
