function [sun] = fSunAngle(time_zone,lat_deg,long_deg, day, hhour)
% Sun angles

% 3/23/2018

% this version in matrix mode for the array of met inputs

% add code for extraterrestrial radiation from Whiteman and Allwine

% Case for Berkeley, CA
% time_zone=-8;
% lat_deg=37.866;
% long_deg=-122.25;


    lat_rad=lat_deg*pi/180;  % latitude, radians
    long_rad=long_deg*pi/180; % longitude, radians

    std_meridian = 0;
    delta_long=(long_deg - std_meridian)*pi/180;

    delta_hours=delta_long*12/pi;


    
    
    declin = -23.45*3.1415/180*cos(2*3.1415*(day+10)/365); % declination angle
    
    cos_hour=-tan(lat_rad)*tan(declin);
    
    sun.sunrise=12- 12* acos(cos_hour/pi); % time of sunrise
    
    sun.sunset=12 + 12* acos(cos_hour)/pi; % time of sunset
    
    sun.daylength=sun.sunset-sun.sunrise;  % hours of day length
    
    f=pi*(279.5+0.9856*day)/180;
    
   
        
       % equation of time, hours
       
        Et=(-104.7*sin(f)+596.2*sin(2*f)+4.3*sin(3*f)-12.7*sin(4*f)-...
            429.3*cos(f)-2.0*cos(2*f)+19.3*cos(3*f))/3600;
        
        % longitudinal correction
        
        Lc_deg = long_deg - time_zone*15; % degrees from local meridian
        
        Lc_hr=Lc_deg*4/60;  % hours, 4 minutes/per degree
        
        T0 = 12-Lc_hr-Et;
                
        hour=pi*(hhour-T0)/12;  % hour angle, radians
        
        % sine of solar elevation, beta
        
        sin_beta=sin(lat_rad)*sin(declin)+cos(lat_rad)*cos(declin)...
            .* cos(hour);
        
        % solar elevation, radians
        
        sun.beta_rad=asin(sin_beta);
        
        % solar elevation, degrees
        
        sun.beta_deg=sun.beta_rad*180/pi;
        
        % solar zenith angles
        
        sun.theta_rad=pi/2 - sun.beta_rad;
    
        sun.theta_deg=sun.theta_rad*180/pi;
        
        % extraterrestrial radiation
        %qsi = So (dbar/d)^2 cos B = So( Y ) [(sin ¢ cos h)(-cos a sin i)
        % - sin h (sln a sin i) + (cos ¢ cos h) cos i] cos 6
        % + [cos ¢ (.cos a sin i) + sin ¢ cos i] sin ~
        
        % eccentricity, 0.0167
        % omega, angular velocity, 360/365, degrees per day
        % need to convert omega to radians
        
        omega = 0.98630;
        dbard = (1 - 0.0167 * cos(omega*day*pi/180)).^-2;
        sun.extraterr= 1365 .* dbard .* cos(sun.theta_rad);
        
end
   



