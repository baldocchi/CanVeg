function [leafang]=LeafAngle(sunang,prm)

% 3/14/2024

% debug
% missing multiplication of cos(A) and sin(A) for Gdiff and expdiff

% suspicious G function  is not right

 % estimate leaf angle for 50 classes between 0 and pi/2
    % at midpoint between each angle class.  This is a big improvement over
    % the older code that divided the sky into 9 classes.
       
    % thetaLeaf=zeros(50,1);   
    % 
    % thetaLeaf=(pi/100:pi/100:pi/2) - pi/200; % solar zenith angle  
    % thetaLeaf=thetaLeaf'; % transpose matrix

    % talking with Martin and he too recommended doing this and I upped it
    % to 90, he mentioned Norman went from 9 to 90 with better computers

    thetaLeaf=zeros(90,1);   
      
    thetaLeaf=(pi/180:pi/180:pi/2) - pi/360; % solar zenith angle  
    thetaLeaf=thetaLeaf'; % transpose matrix
            
    % chose the leaf angle probability distribution, pdf, 
    % use distributions from deWit
    
    % Wang, W.M., Li, Z.L. and Su, H.B., 2007. 
%     Comparison of leaf angle distribution functions: 
%     Effects on extinction coefficient and fraction of sunlit foliage.
%     Agricultural and Forest Meteorology, 143(1): 106-122.



         
    switch prm.leafangle
        case 'planophile'
            leafang.pdf=2.*(1-cos(2*thetaLeaf))./pi;
            
        case 'spherical'
            leafang.pdf=sin(thetaLeaf);
            
        case 'erectophile'
            leafang.pdf=2*(1+cos(2*thetaLeaf))./pi;
            
        case 'plagiophile'
            leafang.pdf=2*(1-cos(4*thetaLeaf))./pi; %plagiophile
            
        case 'extremophile'
            leafang.pdf=2*(1+cos(4*thetaLeaf))./pi; %extremophile
            
        case 'uniform'
           leafang.pdf=2/pi; %uniform

        otherwise
            leafang.pdf=sin(thetaLeaf); % spherical
    end
        
    
   % plot(180*thetaLeaf/pi,pdf)

    % evaluate G function, the direction cosine as a function of solar
    % zenith angle, sun.theta_rad
    
    % using the algorithm from Warren Wilson and Wang et al
    
%     Wang, W. M., Z. L. Li, and H. B. Su. 2007.
%     Comparison of leaf angle distribution functions: 
%     Effects on extinction coefficient and fraction of sunlit foliage.
%     Agricultural and Forest Meteorology 143:106-122.

% Youngryel wrote an Gfunc code for ESPM 129, circa 2008 that uses this simpler
% approach. In my old code I used a convoluted algorithm from Lemeur. So I
% am borrowing from Youngryels code for this subroutine.
    
      

    % call function for G function, the direction cosine, for sun zenith angle
    
    %compute a matrix of Gfunc for all the inputs
    
    thetaSun=sunang.theta_rad;   
  
     leafang.Gfunc=fGfunc_dir(thetaSun, thetaLeaf,leafang.pdf);
     
     leafang.Gfunc(leafang.Gfunc < 0)=0.5;
    
    % call function for Gfunction for each sky sector of the hemisphere
    
    leafang.thetaSky=thetaLeaf;
    
     DA = pi/(2*length(leafang.thetaSky));  % azimuth increment -> (pi/2) / # sectors (50)
         
  
    
    % compute G function for all sky sectors
    
    leafang.Gfunc_Sky=fGfunc_diff(leafang.thetaSky, thetaLeaf, leafang.pdf);
       
    % Integrate the transmission of light through gaps in the hemisphere
      
    % azimuthal symmetry for integrating around 2 pi
      
    % zenith asymmetry integrating for pi/2
      
   
    % compute probability of beam transfer with a markov function for clumped leaves
    
     %dff_Markov % LAI of layer times Markov clumping factor
    
     dff_Markov = prm.dff .* prm.markov;    % the new LAI profile data of Belane
             
     exp_diffuse = exp(-dff_Markov .* (leafang.Gfunc_Sky ./ cos(leafang.thetaSky))');

       
     XX = sum(exp_diffuse .* (cos(leafang.thetaSky) .*sin(leafang.thetaSky))',2);
     leafang.integ_exp_diff = 2. .* XX .* DA;
     
     
     
% from Bonan
     
      % td - exponential transmittance of diffuse radiation through a single leaf layer
% with thickness dlai, estimated for nine sky angles in increments of 10 degrees

% for p = 1:params.npts
%    for iv = canopy.nbot(p):canopy.ntop(p)
%       td(p,iv) = 0;
% 
%       for j = 1:9
% 
%          % Sky angles (5, 15, 25, 35, 45, 55, 65, 75, 85)
% 
%          angle = (5 + (j - 1) * 10) * pi / 180;
% 
%          % Relative projected area of leaf in the direction of sky angle
% 
%          gdirj = phi1(p) + phi2(p) * cos(angle);
% 
%          % Sum transmittance
% 
%          td(p,iv) = td(p,iv) ...
%          + exp(-gdirj / cos(angle) * canopy.dlai(p,iv) * canopy.clumpfac(p)) * sin(angle) * cos(angle);
% 
%       end
%       td(p,iv) = td(p,iv) * 2 * (10 * pi / 180);
% 
%    end
% end
%      
      
      

end