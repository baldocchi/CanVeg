
function [energy]=fEnergy_Leaf_Amphi_Matrix(boundary_layer_res,Q_In,met,prof,gs,prm)


% Energy Balance subroutine
% leaf version, Amphistomatous
% Matlab CanVeg, Matrix form

% D Baldocchi
% Aug 31, 2021


% 	Initialize Parameters
	
% 	Incoming short and longwave radiation

  % Q_In is net incoming minus outgoin short wave and incoming Longwave
  
  energy.LE=zeros(prm.nn,prm.jtot);
  energy.H=zeros(prm.nn,prm.jtot);
  energy.Rnet=zeros(prm.nn,prm.jtot);
  energy.Lout=zeros(prm.nn,prm.jtot);
  Acoef=zeros(prm.nn,prm.jtot);
  Bcoef=zeros(prm.nn,prm.jtot);
  Ccoef=zeros(prm.nn,prm.jtot);
       
     
  % trying to resolve issues of single and double sided stomata, from what
  % was measured and how applied to the model.
  
   gb=1./boundary_layer_res.vapor;
   
   gs_1side = gs/2;   % convert parallel conductance to one sided and then re-compute parallel total conductance with boundary layer conductances
   
   rs_1side=2 ./gs;
   
  % rtop = rstomtop + rb, stomata and boundary layer resistances in series
  % rbottom=rstombot + rb
  
  % top and bottom resistances are in parallel
  
  % rwaterleaf=rtop * rbottom/(rtop + rbottom)
  
  rtop=boundary_layer_res.vapor+rs_1side;
  rbottom=rtop;
  
  rwater=rtop .*rbottom ./(rtop+rbottom);
  
  gw=1./rwater;
   
  %  gw= 2* gs_1side .*gb./(gb+gs_1side);    % identical to the resistance
  %  version
  
 %  gw=gb.*gs./(gb+ 0.5 gs); water conductance,  two sided leaf Licor 6400
  
 % gw = (gb .* gs) ./(gb + 0.5 * gs);  % redid derivation of gw from
 % scratch based on Licor formula.. also the resistance and conductances
 % are identical
 
    
   met.P_Pa=1000 .* met.P_kPa;   % air pressure, Pa
   
    tk2=prof.Tair_K(:,1:prm.jtot) .* prof.Tair_K(:,1:prm.jtot);
    tk3=tk2 .* prof.Tair_K(:,1:prm.jtot);
    tk4=tk3 .* prof.Tair_K(:,1:prm.jtot);

   
   % Using derivation for 2 sided leaves
    	       
	llout = prm.epsigma .* tk4;
   
	 
    % Amphistomatous

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
    
     Acoef=lecoef .* d2est ./(2 *repeat);
     
     Bcoef=-lecoef.*dest-repeat-(Q_In ./repeat) .* (2*Acoef) +2 *Acoef .* (2.*llout./repeat); 
     
     Ccoef=repeat.*lecoef.* vpd_Pa+lecoef.* dest.*(Q_In-llout2)...
         +lecoef .* d2est/2.*(Q_In.*Q_In-4*Q_In.*llout+4.*llout.*llout)./repeat;
      
    
    %  solve for LE
% 
%  a LE^2 + bLE + c = 0

le1 = (-Bcoef + (Bcoef .*Bcoef - 4. * Acoef .* Ccoef).^.5) ./ (2. * Acoef);
le2 = (-Bcoef - (Bcoef .* Bcoef - 4. * Acoef .* Ccoef).^.5) ./ (2. * Acoef);


energy.LE=real(le2);

% solve for dT

% dT = (Q -LE - 2 ep sigma Ta^4)/(2 rho Cp gh + 8 ep sigma Ta^3)

del_Tk = (Q_In-energy.LE-llout2) ./repeat;

Tsfc_K= prof.Tair_K(:,1:prm.nlayers) + (del_Tk);


% Derivation for two sided
% 
%  aT^2 + bT + c =0


% solve for Ts using quadratic solution  
 
       
%  coefficients to the quadratic solution
%          
%          aT = prm.epsigma12 .* tk2 + d2est .* lecoef ./ 2.;
% 
%          bT = prm.epsigma8 .* tk3 + hcoef2 + lecoef .* dest;
%  
%          cT= -Q_In + 2 .* llout + lecoef .* vpd_Pa;
%  


% quadratic equation

% del_Tk=(-bT + (bT .* bT - 4. * aT .* cT).^.5) ./ (2. * aT);

% del_Tk=(-bT - (bT .* bT - 4. * aT .* cT).^.5) ./ (2. * aT);


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
        
     
        % used in program debug. Energy balance closes nicely
        
        figure(5)
        clf
        plot(energy.Rnet(:,1:prm.nlayers),energy.LE(:,1:prm.nlayers)+energy.H(:,1:prm.nlayers),'.')
        xlabel('Rnet W m-2')
        ylabel('H + LE W m-2')
        

end

	
