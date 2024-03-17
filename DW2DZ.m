function [ y ] = DW2DZ( Z, input, turb )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here



% 	double y, dsigw2dz;
% /*
%     ****************************************************************
%        COMPUTES ds2/dz FOR GIVEN s(z)
%     ****************************************************************
% */

	if (Z < input.HH)

y = 2 *Z*turb.del_sigma*turb.del_sigma*input.ustar*input.ustar+2 *turb.sigma_zo*turb.del_sigma*input.ustar; 
	 
% 	  // linear model */
%
%      y = turb.sigma_zo+ z*turb.del_sigma; %//linear model */

%      y2 = a2 + 2 a b z + b2 z2

%      d y2/dz = 2 a b + 2 b2 z
% 
% 		// first compute derivative of sigw^2/u*^2, need to convert to d(s2) /dz 
% 
% 		// sigw=turb.sigma_zo*exp(2.132 *Z/HH);
% 
% 		// sigw2=s(0)^2 * exp(4.264*Z/HH);
% 
%        // dsigw2dz=(s(0)^2/HH) * 4.264 * exp(4.264*Z/HH);
%        
% 		// y=sigw2*ustar*ustar;

% 	 // first compute derivative of sigw^2/u*^2
% 
%      //  dsigw2dz=(turb.sigma_zo*turb.sigma_zo*4.264/HH)*exp(4.264*Z/HH);
% 
% 	// need to convert to ds2/dz so multiplication by u* times u* is needed
% 
%      //   y=dsigw2dz*ustar*ustar;


	
	else
	y = 0.0;
	end






end

