% DispAlfalfa

clear all;
close all;

cd 'C:\Users\Baldocchi\Documents\MATLAB\DispAlfalfa'


% /* =================================================================
% 
% 3-21-2019
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
%          PROGRAM     DispAlfalfa.m

%        This version was Compiled on Microsoft C++ and has been converted to Matlab 
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
%        The simulation volume extends up to 3 * h, the crop height.
%        This altitude is divided in 40 layers (DZ ihigh).
% 
% 	 For atmospheric stability effects, different Dispersion matrices are produced for different z/L values    
% 
%    need to speed up code. found in DispCanveg if I pre allocate space for
%    random numbers based on parcel movements and assign them with
%    normrnd(0,1,[1,#move]) the code runs much faster

tic;

PI = 3.14159;        % Pi
sze = 31 ;            % number of canopy layers plus one
sze3 = 151 ;          % number of  5h atmospheric layers plus one

npart = 10000; %1000000.;             % Number of parcels 1000000
timemax = 1000.0;          % time of parcel run


 input.ustar = 1.00;               %  friction velocity (m s-1) 
 input.HH = 0.55;                  %  canopy height (m) 
 input.DD = 0.33;				   % zero plane displacement	
    
 domain.nlevel=30;                 % number of canopy layers

 % identify variables 

domain.delta_z = input.HH / domain.nlevel;     % distance between vertical layers 

domain.upper_boundary = 5;           % upper bound,five times canopy height,  
ihigh = domain.upper_boundary * domain.nlevel;    % number of vertical layers 
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
% 		   
% 
% 	     // if(( err=fopen_s(&fptr1,"d:\\canalfalfa\\DIJ5000_zLm025.csv","w")) =! 0);  
% 	     // turb.sigma_h= 1.506 ;            // z/L = -0.25

% //		if(( err=fopen_s(&fptr1,"d:\\canalfalfa\\DIJ5000_zLm050.csv","w")) =! 0);
% //	    turb.sigma_h = 1.696 ;          // z/L = -0.50
% 	
% 	 //   if(( err=fopen_s(&fptr1,"d:\\canalfalfa\\DIJ5000_zLm100.csv","w")) =! 0);
% 	//	turb.sigma_h = 1.984 ;         // z/L = -1.00
% 
% 	//	if(( err=fopen_s(&fptr1,"d:\\canalfalfa\\DIJ5000_zLm200.csv","w")) =! 0);
% 	//    turb.sigma_h = 2.391 ;        // z/L = -2.00
	    
% 	//	if(( err=fopen_s(&fptr1,"d:\\canalfalfa\\DIJ5000_zLm300.csv","w")) =! 0);
% 	//	 turb.sigma_h = 2.692 ;       // z/L = -3.00
% 
% 		if(( err=fopen_s(&fptr1,"d:\\canalfalfa\\DIJ5000_zLp050.csv","w")) =! 0); 
% 		turb.sigma_h = 1.35;      // z/L = 0.5
% 
% 	//	if(( err=fopen_s(&fptr1,"d:\\canalfalfa\\DIJ5000_zLp100.csv","w")) =! 0); 
% 	//	turb.sigma_h= 1.5;        // z/L =1.0 


		turb.sigma_zo = 0.1496;         % sigmaw over u* at z equal zero for exponential profile

		turb.sigma_sur= turb.sigma_zo * turb.sigma_h * input.ustar;  % sigma w at zero for linear profile

        turb.del_sigma = (turb.sigma_h * input.ustar - turb.sigma_sur) / input.HH; % difference in sigma w


% seed random number with time 
rng('shuffle')   
   

% /*
% ***************************************************************
%        Time step length (* TL) and Lagrangian Length scale
% ****************************************************************
% */


parcel.delta_t = .1 * TL(input.HH, input, turb);                         % time step */
turb.laglen = turb.sigma_h * input.ustar * TL(input.HH, input, turb);     %// Lagrangian length scale */

parcel.sumn = 0;

parcel.wwmean=zeros(ihigh1+1,1);
%parcel.partadd=zeros(ihigh1+1,1);


%         ******************************************************
%         number of particles per layer, redistribute npart #
%         with height according to the relative source strength
%         *****************************************************

 domain.nlev(1:domain.nlevel) = npart / domain.nlevel;
 parcel.sumn = sum(domain.nlev);

% /*
% '       nn is the number of time increments that the particle travels
% '       TIMEMAX = t [s]
% '       Actual time of particle trave is TIMEMAX*TL.
% '       Number of travel steps is t/dt.
% */

        nn = (timemax / parcel.delta_t);
        
        
        
       
        
        
       
% /*
% '    ****************************************************************
% '       Start the release of particles:
% '    ****************************************************************
% */
       % parcel.sum_random = 0;
       % parcel.sum_w = 0.;
       % parcel.sum_var_rand = 0;
        parcel.move = 1;

       % IT = 1;                          % /* particle counter */


%        /*
%         assume values of a Gaussian distribution
%         random numbers with a mean of zero and a variance of one.
%         */
          

         timescale = TL(input.HH, input, turb);

% /*
%         Release of particles is carried out separately for each layer,
%         starting with the lowest level.
% */
 parcel.sum_w = zeros(floor(nn),1);

  for ilevel = 1:domain.nlevel
     
       
		%parcel.sum_random = zeros(nrand,1);
  
% /*
% 
%     ****************************************************************
%        1-D case. Particles are released from each level z.  We have
%        a continuous source and a horizontally homogeneous canopy.
%     ****************************************************************
% */


    %for I = 1:ihigh    %   Having problem end is at next Ilevel not next I
%     /*
% 
%     ****************************************************************
%        1-D case. Particles are released from each level z.  We have
%        a continuous source and a horizontally homogeneous canopy.
%     ****************************************************************
% */


%     for(I = 1:ihigh1)
%     parcel.consum(I) = 0;
%     end

    parcel.consum=zeros(ihigh1,1);

% /*
% ****************************************************************
%       at each level NLEV(LEVEL) particles are released
% ****************************************************************
% */
        for part=1:domain.nlev(ilevel)
            
            %(part = 1; part <= domain.nlev[ilevel]; part++)
        

      %  IT = IT+1;
        
            
% 	xmod= IT;
% 	zmod=mod(xmod,100);

	parcel.z = ilevel * domain.delta_z;   %// initial height at the level

% // the initial vertical velocity
      

        random.random=normrnd(0,1);  % do I want a Gaussian random


% //       vertical velocity, WW
% 

        parcel.w = SIGMA(parcel.z,input,turb) * random.random;

%         /*
%         number of particle movements
%         */

        parcel.move = parcel.move + 1;
        

%         /* compute mean and variance of w and random number */

        %parcel.sum_random = parcel.sum_random + random.random;
        %parcel.sum_random(parcel.move) = random.random;
        parcel.sum_w(parcel.move) = parcel.w;
        %parcel.wwmean(ilevel) = parcel.wwmean(ilevel) + parcel.w;
        %parcel.sum_var_rand = parcel.sum_var_rand + random.random*random.random;

%    // The particle starts its run

%/*         for (I=1; I <= nn;I++) */

        IZF =  min((parcel.z / domain.delta_z + 1), ihigh1);

        I=1; 

         while (I <=nn && IZF <= ihigh) 

      

%        // Compute the vertical position and reflect z if it is zero
% 	   // Need to reflect also at top in future, but need deeper domain
	   
         % 

          if (parcel.z <= 0)  %// reflect particle if at ground 
          parcel.w = -parcel.w;
          parcel.z = -parcel.z;
          end

          IZF = floor(min((parcel.z / domain.delta_z + 1), ihigh1));

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

                random.random= normrnd(0,1);

%                 /*
%                 wnew = -wold dt/Tl) + 1/2 dvarw/dz (1+w^2/varw)+
%                  (2 varw dt/Tl)du
%                 */

                timescale = TL(parcel.z, input, turb);
                
                dtT=parcel.delta_t / timescale;
                
                parcel.std_w = SIGMA(parcel.z, input, turb);
                parcel.var_w = parcel.std_w * parcel.std_w;
                
                %parcel.term1 = -parcel.w * parcel.delta_t / timescale;
                parcel.term1 = -parcel.w * dtT;
                
                parcel.term2 = .5 * DW2DZ(parcel.z, input, turb) * (1. + (parcel.w * parcel.w) / parcel.var_w) * parcel.delta_t;
                
                %parcel.term3 = power((2 * parcel.var_w * parcel.delta_t / timescale),.5) * random.random;
		        %parcel.term3 = power((2 * parcel.var_w * dtT),.5) * random.random;	
                parcel.term3 = sqrt(2 * parcel.var_w * dtT) * random.random;	
                  
                parcel.w = parcel.w + parcel.term1 + parcel.term2 + parcel.term3;


% /*    ****************************************************************
% '       STATISTICAL CHECK OF RANDOM NUMBER AND MEAN VERTICAL VELOCITY
% '    ****************************************************************
% */
%                 /*
%                 number of occurences at height IZF and its
%                 mean vertical velocity
%                 */

                %parcel.wwmean(IZF) = parcel.wwmean(IZF)+ parcel.w;

                parcel.move = parcel.move + 1;
                %parcel.sum_random = parcel.sum_random + random.random;
                %parcel.sum_random(parcel.move) = random.random;
                parcel.sum_w(parcel.move) =  parcel.w;
                %parcel.sum_var_rand = parcel.sum_var_rand + random.random*random.random;

                I = I + 1;
                
                parcel.z = parcel.z + parcel.w * parcel.delta_t;

               
         end   % end while
             
               
              % parcel.partadd(IZF) = parcel.partadd(IZF)+ 1;

        end %//  end while next particle 


%     /*
%     Introduce computation of concentration at each level.
%     Use super-position principle to compute canopy
%     concentration profile.
%     */

%         for I = 1:ihigh
%         parcel.conc(I) = parcel.consum(I);
%         end % next I

        parcel.conc(1:ihigh) = parcel.consum(1:ihigh);
        
       %// Compute the dispersion matrix then reset concentrations 
		
% 		for I=1:ihigh
%         parcel.DIJ(I,ilevel) = (parcel.conc(I) - parcel.conc(ihigh)) / (domain.delta_z * domain.nlev(ilevel));
%         end % next I

        
        parcel.DIJ(1:ihigh,ilevel) = (parcel.conc(1:ihigh) - parcel.conc(ihigh)) / (domain.delta_z * domain.nlev(ilevel));
       
        
% /*
% '    ****************************************************************
% '       Statistical check of random number: VARR=variance (should = 1),
% '       SUMN = mean (should = 0), sumw# = mean vertical velocity
% '       (should = 0)
% '    ****************************************************************
% */
	

		%parcel.sum_var_rand = (parcel.sum_var_rand - (power(parcel.sum_random,2.0) / parcel.move)) / (parcel.move - 1.);

	   % parcel.sum_var_rand(ilevel)=var(parcel.sum_random);
		parcel.wwmean(ilevel) = mean(parcel.sum_w);
		%parcel.mean_random(ilevel)=mean(parcel.sum_random);
        
        disp(ilevel);
		
  end   % next Ilevel

% 		
% 		printf(" \n");
% 		printf(" %6i \n", parcel.move);
% 		printf("mean r:  %f \n ", parcel.sum_random);
% 		printf("var. r:  %f \n ", parcel.sum_var_rand);
% 		printf("mean w:  %f \n ", parcel.sum_w);


% plot of Dispersion Matrix

figure(1)
clf
plot(parcel.DIJ, 1:150)
ylabel('level')
xlabel('Dij')
title('alfalfa dispersion matrix')
xlim([0 12])

toc;


