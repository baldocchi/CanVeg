
function [parcel]=DispCanveg_v2a(prm)

% DispCanveg


% /* =================================================================
% 
%   3-31-2021
% 
% 	Dennis Baldocchi
% 	Ecosystem Science Division
% 	Department of Environmental Science, Policy and Management
% 	345 Hilgard Hall
% 	University of California, Berkeley
% 	Berkeley, CA
% 	Baldocchi@berkeley.edu
% 
% -------------------------------------------
%          PROGRAM     

%        This version was converted to Matlab 
% 
%        This program computes the dispersion matrix, according to
%        Raupach (1989).  This model is coupled later to CANOAK
%        the oak photosynthesis model, to compute CO2 profiles.
% 
%        A number of particles is released at 40 levels
%        in a horizontally, homogeneous canopy.
% 
%        The vertical velocity w is computed with a Markov sequence.
%        The algorithm after the approach of Thomson (1987)
%        is used.
% 
%               dw(t) = a(w,z) dt + b(w,z) du
% 
%        where du is a random number with random increment with
%        mean of zero, a variance equal one and a Gaussian probability distribution.
% 
%        
%        This model computes the random flight field for a 1-dimensional
%        canopy and continuous release.
% 
%        Dispersion is only in the vertical and the canopy is
%        assumed to be horizontally homogeneous.
% 
%       
%        The system studied here is analogous to a
%        volume averaged point source, stacked one atop another.
%        Concentrations are computed on the principle of superposition.
% 
%        Since the canopy is horizontally homogeneous we can assume
%        an infinite extent in the y direction, so the concentration
%        does not vary in y and the source can be expressed in terms of
%        unit y.
% 
%       
%        In this program the particles released for a prescribe
%        time duration--TIMEMAX [s]. Tests in the past show that TIMEMAX is different
%        for tall forests and short crops.
% 
%        The simulation is for alfalfa
%        The met sensors are at 5 m and canopy height is about 3 m.
%        The dz is then different in the canopy vs above it
% 
% 	     For atmospheric stability effects, different Dispersion matrices are produced for different z/L values    
% 
%        the origional version took about 132 seconds to run for 100000
%        particles

%        pre allocated space for random numbers and set an array of random
%        numbers that I call with parcel.move. sped up computations to 
%        seconds.

%       Can I save Extra Time by creating lookup tables for TL, sigmaW and
%       dSigWdz

%       created on for TL and cut down time from 70 s for that call to
%       0.002, total run time to 70 s

%       created for sigmaw, total run time 53 s

%       created dvardz and total run time is 36s

% new review by Brunet 2020 BLM on parameters of turbulence



tic;
%profile on;

ht_atmos=prm.meas_ht-prm.veg_ht;
dht_canopy=prm.dht_canopy;
dht_atmos=prm.dht_atmos;
nlayers_atmos=prm.nlayers_atmos;

% consider look up tables for functional calls

zht=prm.zht;  % heights of each layer

sze=prm.nlayers ;                 % number of canopy layers plus one
domain.nlevel=sze;                 % number of canopy layers

npart = prm.npart; %; %1000000.;             % Number of parcels 1000000
timemax = 1000.0;                      % time of parcel run


 input.ustar = 1.00;               %  friction velocity (m s-1) 
 input.HH = prm.veg_ht;            %  canopy height (m) 
 input.DD=prm.dht;
  

 domain.delta_z=prm.delz;

% in this version domain delta_z is different in canopy or above. it
% changes in the code if z is within or above the canopy

domain.upper_boundary = prm.meas_ht;           % upper bound, ? times canopy height

%ihigh = domain.upper_boundary * domain.nlevel;    % number of vertical layers 

ihigh=prm.nlayers_atmos;

ihigh1 =ihigh + 1; 
	


% 		// Kaimal and Finnigan
% 
% 		// sigma_w/u* = 1.25(1+ 3 abs(z/L))^.333  for z/L < 0
% 
% 		// sigma_w/u* = 1.25 (1 + 0.2 z/L) for z/L 0 to 1
% 
% 		//   Thomson dispersion matrix  
% 
        turb.sigma_h = 1.25;      %  sigma w > HH  */

        %	turb.sigma_h= 1.506 ;           % // z/L = -0.25
        % 	turb.sigma_h = 1.696 ;          % z/L = -0.50
        %	turb.sigma_h = 1.984 ;         % z/L = -1.00
        %	turb.sigma_h = 2.692 ;       % z/L = -3.00
        %	turb.sigma_h = 1.35;     % z/L = 0.5
        	turb.sigma_h= 1.5;       %  z/L =1.0 


		turb.sigma_zo = 0.25;         % sigmaw over u* at z equal zero for exponential profile, Brunet 2020 BLM review

		turb.sigma_sur= turb.sigma_zo * turb.sigma_h * input.ustar;  % sigma w at zero for linear profile

        turb.del_sigma = (turb.sigma_h * input.ustar - turb.sigma_sur) / input.HH; % difference in sigma w


% seed random number with time 
rng('shuffle')   
 

% /*
% ***************************************************************
%        Time step length (* TL) and Lagrangian Length scale
% ****************************************************************
% */


fract=0.1;  % was 0.1 for short crops

parcel.delta_t = fract * TL(input.HH, input, turb);                         % time step 
turb.laglen = turb.sigma_h .* input.ustar .* TL(input.HH, input, turb);     % Lagrangian length scale */

% /*
% '       nn is the number of time increments that the particle travels
% '       TIMEMAX = t [s]
% '       Actual time of particle trave is TIMEMAX*TL.
% '       Number of travel steps is t/dt.
% */

        nn = fix((timemax / parcel.delta_t));

% pre-allocate random number and see if code is faster

random.rndij=zeros(nn+1,1);
random.rndij=normrnd(0,1,[nn+1,1]);

parcel.sumn = 0;

parcel.wwmean=zeros(ihigh1+1,1);
%parcel.DIJ=zeros(ihigh1,prm.jtot);
parcel.DIJ=zeros(ihigh,prm.jtot);

%parcel.partadd=zeros(ihigh1+1,1);


%         ******************************************************
%         number of particles per layer, redistribute npart #
%         with height according to the relative source strength
%         *****************************************************

 domain.nlev(1:domain.nlevel) = npart / domain.nlevel;
 parcel.sumn = sum(domain.nlev);

% '    ****************************************************************
% '       Start the release of particles:
% '    ****************************************************************
% */
       % parcel.sum_random = 0;
       parcel.sum_w = 0.;
       % parcel.sum_var_rand = 0;
        parcel.move = 1;

       % IT = 1;                          % /* particle counter */

         timescale = TL(input.HH, input, turb);

% /*
%         Release of particles is carried out separately for each layer,
%         starting with the lowest level.
% */
%parcel.sum_w = zeros(floor(nn),1);
 
  sigmaw=zeros(ihigh1,1);
 
 % calling sigma and dvarwdz once for each height grid
 % sigw=   SIGMA(zht,input,turb);
 
 sigmaw=Sigma_w_zht(zht,input,turb);
 
 dvarwdz=dw2dz_w_zht(zht,input,turb);
 
 tldz= tl_dz(zht,input);


    for ilevel = 1:domain.nlevel
     
%     ****************************************************************
%        1-D case. Particles are released from each level z.  We have
%        a continuous source and a horizontally homogeneous canopy.
%     ****************************************************************


%     ****************************************************************
%        1-D case. Particles are released from each level z.  We have
%        a continuous source and a horizontally homogeneous canopy.
%     ****************************************************************

    parcel.consum=zeros(ihigh1,1);


% ****************************************************************
%       at each level NLEV(LEVEL) particles are released
% ****************************************************************

        for part=1:domain.nlev(ilevel)
           
	parcel.z = double(ilevel) * prm.dht_canopy;   %// initial height at the level

% // the initial vertical velocity
  
       % random.random=normrnd(0,1);    % mean of zero, std dev 1

      
        IZF =  min(fix(parcel.z / prm.dht_canopy + 1), ihigh);

        I=1; 

% //       vertical velocity, WW


        %parcel.w = SIGMA(parcel.z,input,turb) * random.random;
        parcel.w =sigmaw(IZF)*random.rndij(I);


%         number of particle movements
        parcel.move = parcel.move + 1;
        

%         /* compute mean and variance of w and random number */
 
        parcel.sum_w = parcel.sum_w+parcel.w;
       

         while (I <=nn && IZF <= ihigh)     % <= ihigh

      

%        // Compute the vertical position and reflect z if it is zero
% 	   // Need to reflect also at top in future, but need deeper domain
	   
         % 

         % if (parcel.z <= 0)  %// reflect particle if at ground 
         parcel.w(parcel.z <= 0) = -parcel.w;
         parcel.z(parcel.z <= 0) = -parcel.z;
         % end

          IZF = floor(min((parcel.z / prm.dht_canopy + 1), ihigh1));  %ihigh1

% /*
% 
%     Compute the concentration of material in the controlled volume.
% 
%     Here we use the algorithm of Raupach (1989).  The ensemble average
%     concentration considers the fact that we have an extensive,
%     horizontally homogeneous and continuous source.  Information from
%     every step from t=0 to t=T is used. It is identical to releasing a plane
%     source with length x or ut, as does Wilson et al.
% */
       
        parcel.consum(IZF) =  parcel.consum(IZF) + parcel.delta_t;

% /*     Compute the new vertical velocity.
%        Introduce the bias velocity for the case of inhomogeneous
%        turbulence.  This is needed to prevent the false accumulation
%        of material near the ground that would otherwise occur as
%        fast air is brought downward, yet slower air at lower levels
%        is less apt to leave. (see Wilson et al. Thompson etc.)
% */

               % random.random= normrnd(0,1);


%                 wnew = -wold dt/Tl) + 1/2 dvarw/dz (1+w^2/varw)+
%                  (2 varw dt/Tl)du

 
% try and speed computations by using look up table for Tl
                %timescale = TL(parcel.z, input, turb);
                timescale=tldz(IZF);
                
                dtT=parcel.delta_t / timescale;
                
               % parcel.std_w = SIGMA(parcel.z, input, turb);
                parcel.std_w=sigmaw(IZF);
                
                parcel.var_w = parcel.std_w * parcel.std_w;
                
                %parcel.term1 = -parcel.w * parcel.delta_t / timescale;
                parcel.term1 = -parcel.w * dtT;
                
                %parcel.term2 = .5 * DW2DZ(parcel.z, input, turb) * (1. + (parcel.w * parcel.w) / parcel.var_w) * parcel.delta_t;
                parcel.term2 = .5 * dvarwdz(IZF)* (1. + (parcel.w * parcel.w) / parcel.var_w) * parcel.delta_t;
                
                %parcel.term3 = power((2 * parcel.var_w * parcel.delta_t / timescale),.5) * random.random;
		        %parcel.term3 = sqrt(2 * parcel.var_w * dtT) * random.random;	
                 parcel.term3 = sqrt(2 * parcel.var_w * dtT) * random.rndij(I);
                 
                
                  
                parcel.w = parcel.w + parcel.term1 + parcel.term2 + parcel.term3;


% /*    ****************************************************************
% '       STATISTICAL CHECK OF RANDOM NUMBER AND MEAN VERTICAL VELOCITY
% '    ****************************************************************
% */
%                 /*
%                 number of occurences at height IZF and its
%                 mean vertical velocity
%                 */

                parcel.move = parcel.move + 1;
                parcel.sum_w = parcel.sum_w+ parcel.w;
                
                I = I + 1;
                
                parcel.z = double(parcel.z) + parcel.w * parcel.delta_t;

               
         end   % end while particle less than max steps or in domain
             
               
              % parcel.partadd(IZF) = parcel.partadd(IZF)+ 1;
              random.rndij=normrnd(0,1,[nn+1,1]);

        end %//  end while next particle 



%     Introduce computation of concentration at each level.
%     Use super-position principle to compute canopy
%     concentration profile.

%      values of parcel.consum(ihigh+1) include parcels that have gone past
%      the top, so just want to normalize to ihigh

        parcel.conc(1:ihigh) = parcel.consum(1:ihigh);
       
        
       %// Compute the dispersion matrix then reset concentrations 
        
        parcel.DIJ(1:ihigh,ilevel) = (parcel.conc(1:ihigh) - parcel.conc(ihigh)) / (domain.delta_z(ilevel) * double(domain.nlev(ilevel)));
       
       % parcel.DIJ(1:ihigh,ilevel) = (parcel.conc(1:ihigh) ) / (domain.delta_z(ilevel) * domain.nlev(ilevel));
   
% '    ****************************************************************
% '       Statistical check of random number: VARR=variance (should = 1),
% '       SUMN = mean (should = 0), sumw# = mean vertical velocity
% '       (should = 0)
% '    ****************************************************************

	
		parcel.wwmean(ilevel) = mean(parcel.sum_w);
		
        
        disp(ilevel);
		
  end   % next domain Ilevel

% plot of Dispersion Matrix

figure(1)
clf
%plot(parcel.DIJ, 1:ihigh1)
plot(parcel.DIJ, 1:ihigh)
ylabel('level')
xlabel('Dij')
title(prm.title)
xlim([0 25])

toc


