function [ y ] = dw2dz_w_zht( Z, input, turb )

% /*
%     ****************************************************************
%        COMPUTES ds2/dz FOR GIVEN s(z)
%     ****************************************************************
% */

	y=zeros(1,length(Z));

y(Z < input.HH) = 2 *Z(Z < input.HH)*turb.del_sigma*turb.del_sigma*input.ustar*input.ustar+2 *turb.sigma_zo*turb.del_sigma*input.ustar; 
	 
% 	  // linear model */
% 
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

		
	y(Z >= input.HH) = 0.0;
	

end

