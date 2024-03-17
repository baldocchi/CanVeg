function [solar] = fDiffuse_Direct_Radiation(rglobal,parin,press_kpa, sine_beta)
%
% /*
%  This subroutine uses the Weiss-Norman ( 1985, Agric. forest Meteorol. 34: 205-213)
%  routine tocompute direct and diffuse PAR from total par

% it was converted from my C code to Matlab
% 
%  revised 8/3/2023

% convert to matrix form

% There were two typos in Weiss and Norman (1985).
% Equation (3) should be
% Rdv=0.4(600 cos(theta) - RDV)
% Equation (5) should be
% Rdn=0.6(720 - RDN/cos(theta) -w) cos(theta).
% 
% Weiss and Normam assume a solar constant of 1320 W m-2, 
% which is much lower than 1360 W m-2

% revisiting the paper I see they excluded uv at the top of the atmosphere
% 
% 
%        fractions of NIR and PAR (visible)
% */


        solar.constant=1360.8;  % Kopp and Lean 2011 GRL
        
% Table 2.2: The distribution of extraterrestrial solar radiation.
% Ultra Violet	200-400 nm	8.7%
% Visible	400-700 nm	38.3%
% Near Infra Red	700-3500nm	51.7%
% http://www.itacanet.org/the-sun-as-a-source-of-energy/part-2-solar-energy-reaching-the-earths-surface/
        
        fuv=0.087;        
        fnir = .517;
        fvis = .383;

        ru = press_kpa ./ (101.3 * sine_beta);  % thickness of atmosphere for sites at sea level, like delta 
        
        ru(ru<=0)=NaN;  % for negative sun angles getting complex numbers
        
              
        Rgpotential= solar.constant .* exp(-.185 .* ru) .* sine_beta;
        
        TF=isnan(Rgpotential);
        
        Rgpotential(TF==1)=0;
       

% //         visible direct PAR

        viscoef=fvis*solar.constant;
        
        rdvis = viscoef * exp(-.185 .* ru) .* sine_beta;


% //      potential diffuse PAR
% 
% //        rsvis = .4 * (600.0 - rdvis) * solar.sine_beta;
% 
% // corrected version

        rsvis = 0.4 * (viscoef * sine_beta -rdvis); 



% /*
%         solar constant was assumed to be: 1320 W m-2. They excluded uv
% 
%         
% 
%         water absorption in NIR for 10 mm precip water

%        Reitan 1963 JApplMet has a regression fit on precipitable water
%        and Tdew

%       ln(W) (cm) = -0.981 + 0.034 Tdew(F)
% */

        wa = solar.constant .* .077 .* (2. * ru) .^ 0.3;

       
        %direct beam NIR
        
        nircoef=solar.constant*fnir;

        rdir = (nircoef .* exp(-.06 .* ru) - wa) .* sine_beta;
      
        TF=isnan(rdir);
        rdir(TF == 1)=0;
        

% /*
%         potential diffuse NIR
% */
% 
%  //       rsdir = .6 * (720 - wa - rdir) * solar.sine_beta;  // Eva asks if we should correct twice for angles?
%  
% // corrected version, Rdn=0.6(720 - RDN/cos(theta) -w) cos(theta).


		rsdir = 0.6 .* (nircoef -rdvis./sine_beta-wa) .* sine_beta;

        TF=isnan(rsdir);
        rsdir(TF == 1)=0;
        

        rvt = rdvis + rsvis;
        rit = rdir + rsdir;
        
        TF=isnan(rvt);
        rvt(TF == 1)=0;
        
         TF=isnan(rit);
        rit(TF == 1)=0;


        % check for negative values and values above 0.9
        % seem to be having large ratrad and need to double think about
        % night
        
        ratrad = rglobal ./ (rvt + rit);

        TF=isnan(ratrad);
        ratrad(TF == 1)=1;
        
        ratrad(ratrad==inf)=1;
        
% /*
%         ratio is the ratio between observed and potential radiation
% 
%         NIR flux density as a function of PAR
% 
%         since NIR is used in energy balance calculations
%         convert it to W m-2: divide PAR by 4.6
% */
% 

        nirx = rglobal - (parin ./ 4.6);  % W m-2


% //   ratio = (PARIN / 4.6 + NIRX) / (rvt + rit)

        
        ratrad(ratrad >= .9) = .89;
        ratrad(ratrad <= 0)=0.22;
       


% //     fraction PAR direct and diffuse

        xvalue=(0.9-ratrad)/.70;

        fvsb = rdvis ./ rvt .* (1. - power(xvalue,.67));

        fvsb(fvsb < 0)=0;
        

        fvsb(fvsb > 1)=1.0;
        
        fvd = 1. - fvsb;


% //      note PAR has been entered in units of uE m-2 s-1

        solar.par_beam = fvsb .* parin;
        solar.par_diffuse = fvd .* parin;

        if(solar.par_beam <= 0)
        
        solar.par_beam = 0;
        solar.par_diffuse = input.parin;
        
        end

        if(parin == 0)
        
        solar.par_beam=0.001;
        solar.par_diffuse=0.001;
        
        end

        xvalue=(0.9-ratrad)/.68;
        fansb = rdir ./ rit .* (1. - xvalue .^0.67);
        
         TF = isnan(fansb);
        fansb(TF)=0;
        fansb(fansb < 0) = 0;
        
        fansb(fansb > 1)=1.0;
       


        fand = 1. - fansb;


% //      NIR beam and diffuse flux densities

        solar.nir_beam = fansb .* nirx;
        solar.nir_diffuse = fand .* nirx;

        if(solar.nir_beam <= 0)
        
        solar.nir_beam = 0;
        solar.nir_diffuse = nirx;
        
        end

        if (nirx == 0)
        
        solar.nir_beam=0.1;
        solar.nir_diffuse=0.1;
        
        end

	    solar.nir_beam= nirx-solar.nir_diffuse;
        solar.par_beam = parin-solar.par_diffuse; 

        TF = isnan(solar.par_beam);
        solar.par_beam(TF)=0;
        
        TF = isnan(solar.par_diffuse);
        solar.par_diffuse(TF)=0;
        
        solar.ratrad=ratrad;
        
        solar.Rgpotential=Rgpotential;

end

