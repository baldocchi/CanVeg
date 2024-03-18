function [fT] =fTBOLTZ (rate, eakin, topt, tl)
        
        % Boltzmann temperature distribution for photosynthesis 
        
        % matrix version
        
        param.Rstar = 8.3144;
        param.hkin = 200000.0;    % 200000 enthalpy term, J mol-1
        param.skin = 710.0;       % entropy term, J K-1 mol-1

        dtlopt = tl - topt;
        prodt = param.Rstar * topt .* tl;
        numm = rate * param.hkin * exp(eakin .* (dtlopt) ./ (prodt));
        denom = param.hkin - eakin * (1.0 - exp(param.hkin .* (dtlopt) ./ (prodt)));
        fT = numm ./ denom;
        
        
        
end

