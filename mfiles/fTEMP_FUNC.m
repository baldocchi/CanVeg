%/////////////////////////////////////////////////////////////      
        function [fT] =fTEMP_FUNC(rate,eact,tprime,tref,t_lk)

        %Arhennius temperature function  
        % Matrix version
        
        param.Rstar = 8.3144;
            
        fT = rate .* exp(tprime * eact ./ (tref * param.Rstar .*t_lk));
        

end

