function [ y ] = SIGMA(z,input, turb)

% '    ****************************************************************
% '       This function gives the s(w) value for height z
% '       Use linear decrease in sigma with z/h, as Wilson et al
% '       show for a corn canopy.
% '    ****************************************************************
% '
% '   DELSIG=(SIGMAH-SIGMASUR)/HH
% '
% */

if (z < input.HH)

 
	y = turb.sigma_zo+ z*turb.del_sigma; %//linear model */


% // exponential model, computed as function of sigw/u* and z/h
% // need to convert to sigma w so final multiplication by u* is needed
% 
%      
% 		// sigw=turb.sigma_zo*exp(2.132 *Z/HH);
% 
%        // y=sigw*ustar;  // multiply by ustar to compute sigma w


	
	else
		y = turb.sigma_h * input.ustar; 

end




end

