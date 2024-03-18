
function [Sun, Shade]=fEnergy_Carbon_Fluxes_Matrix(Sun,Shade,Qin, quantum,met,prof, prm)

% 8/31/2021 

%   The ENERGY_AND_CARBON_FLUXES routine to computes coupled fluxes
% 	of energy, water and CO2 exchange, as well as leaf temperature.  Computataions
% 	are performed for each layer in the canopy and on the sunlit and shaded fractions.
% 
% 	Analytical solution for leaf energy balance and leaf temperature is used.  The program
% 	is derived from work by Paw U (1986) and was corrected for errors with a re-derivation
% 	of the equations.  The quadratic form of the solution is used, rather than the quartic
% 	version that Paw U prefers.
% 
% 	Estimates of leaf temperature are used to adjust kinetic coefficients for enzyme kinetics,
% 	respiration and photosynthesis, as well as the saturation vapor pressure at the leaf surface.
% 
% 	The Analytical solution of the coupled set of equations for photosynthesis and stomatal
% 	conductance by Baldocchi (1994, Tree Physiology) is used.  This equation is a solution to
% 	a cubic form of the photosynthesis equation.  The photosynthesis algorithms are from the
% 	model of Farquhar.  Stomatal conductance is based on the algorithms of Ball-
% 	Berry and Collatz et al., which couple gs to A
% 
% 	Layer 1 is the soil, Layer 30 is top of the canopy

% 

% compute sun fluxes and ingest the right variables
   
    % compute for top and bottom of sun 
     
   % based on sunlit leaf temperature and air temperature of the layer
   % compute boundary layer resistances for heat, vapor and CO2
   
    Leaf.Tsfc=quantum.prob_beam(:,1:prm.jtot) .* Sun.Tsfc +quantum.prob_shade(:,1:prm.jtot) .* Shade.Tsfc;   
   
    [boundary_layer_res]=fBoundary_Resistance_Matrix(prof,met, Leaf.Tsfc, prm);
   
   [boundary_layer_res]=fBoundary_Resistance_Matrix(prof,met, Sun.Tsfc, prm);

   % Compute leaf photosynthesis
   
   switch(prm.stomata)
        case 'Hypo'
   
    [ps]=fLeafPsHypoMatrix(quantum.sun_abs,prof.co2(:,1:prm.jtot), Sun.Tsfc,...
       boundary_layer_res.co2,met.P_kPa,prof.eair_Pa(:,1:prm.jtot),prm);

    case 'Amphi'

      [ps]=fLeafPsAmphiMatrix(quantum.sun_abs(:,1:prm.jtot),prof.co2(:,1:prm.jtot), Sun.Tsfc,...
           boundary_layer_res.co2,met.P_kPa,prof.eair_Pa(:,1:prm.jtot),prm);
       
        otherwise
    end
   
   
   Sun.Ps=ps.aphoto;
   Sun.gs=ps.gs_m_s;
   Sun.Resp=ps.rd;
   
   
   % Compute energy balance on the top of sunlit leaves
  
   % pass and use prm if leaf is amphistomatous or hypostomatous
   
   disp('sun')
   
    switch(prm.stomata)
        case 'Hypo'
     [ans]=fEnergy_Leaf_Hypo_Matrix(boundary_layer_res,Qin.sun_abs,met,prof,Sun.gs,prm);
        case 'Amphi'
    [ans]=fEnergy_Leaf_Amphi_Matrix(boundary_layer_res,Qin.sun_abs,met,prof,Sun.gs,prm);
   
    
    % gives bad energy balancy
    %[ans]=fEnergy_Leaf_Amphi_MatrixV2(boundary_layer_res,Qin.sun_abs,met,prof,Sun.gs,prm);
        otherwise
    end
  
    % Compute energy balance for the bottom of sunlit leaves
    % if assuming hypostomatous assign gs a low value, eg 0.01 m/s
    
      
    Sun.LE=ans.LE ;
    Sun.H=ans.H ;
    Sun.Tsfc=ans.Tsfc;    % Radiation weighted not sure what the average is
    Sun.Lout=ans.Lout;
    Sun.Rnet=Qin.sun_abs-Sun.Lout;
    Sun.vpd_Pa=ans.vpd_Pa;
    Sun.closure=ans.closure;
  

    % redo for Shade fraction
    
      [boundary_layer_res]=fBoundary_Resistance_Matrix(prof,met, Shade.Tsfc, prm);
    
    switch(prm.stomata)
        case 'Hypo'
       
        [ps]=fLeafPsHypoMatrix(quantum.sh_abs(:,1:prm.jtot),prof.co2(:,1:prm.jtot),Shade.Tsfc,...
           boundary_layer_res.co2,met.P_kPa,prof.eair_Pa(:,1:prm.jtot),prm);
       
        case 'Amphi'
      [ps]=fLeafPsAmphiMatrix(quantum.sh_abs(:,1:prm.jtot),prof.co2(:,1:prm.jtot),Shade.Tsfc,...
           boundary_layer_res.co2,met.P_kPa,prof.eair_Pa(:,1:prm.jtot),prm);
   
     otherwise
    end
   
   Shade.Ps=ps.aphoto;
   Shade.gs=ps.gs_m_s;     % make sure I use conductance with m/s and it is for water
   Shade.Resp=ps.rd;

   
   disp('shade')
    
   
     switch(prm.stomata)
        case 'Hypo'
     [ans]=fEnergy_Leaf_Hypo_Matrix(boundary_layer_res,Qin.shade_abs,met,prof,Shade.gs,prm);
        case 'Amphi'
    [ans]=fEnergy_Leaf_Amphi_Matrix(boundary_layer_res,Qin.shade_abs,met,prof,Shade.gs,prm);
    
   % gives bad energy balance
   % [ans]=fEnergy_Leaf_Amphi_MatrixV2(boundary_layer_res,Qin.shade_abs,met,prof,Shade.gs,prm);
        otherwise
    end
  
    % radiation weighted leaf temperature
    Shade.Tsfc=ans.Tsfc;  
    Shade.LE=ans.LE;
    Shade.H=ans.H;
    Shade.Lout=ans.Lout;
    Shade.Rnet=Qin.shade_abs - Shade.Lout;
    Shade.vpd_Pa=ans.vpd_Pa;
    Shade.closure=ans.closure;
  
  
end
 

