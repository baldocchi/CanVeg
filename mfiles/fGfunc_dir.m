function [Gfunc] = fGfunc_dir(theta_rad, thetaLeaf,pdf)

    % evaluate G function, the direction cosine as a function of solar
    % zenith angle, theta_rad
    % March 26, 2019
    
    Gfunc=zeros(length(theta_rad),1);
    A=zeros(length(pdf),1);
    F=zeros(length(pdf),1);
       
   for i=1:length(theta_rad) 
    
    
    %%%% old code
    % using the algorithm from Warren Wilson
    
        use=find((abs(cot(theta_rad(i)).*cot(thetaLeaf))>1)==1);
        for j=1:length(use)
            A(use(j))=cos(theta_rad(i)).*cos(thetaLeaf(use(j)));
        end
        use=find((abs(cot(theta_rad(i)).*cot(thetaLeaf))>1)==0);
        for j=1:length(use)
            psi(use(j))=acos(cot(theta_rad(i)).*cot(thetaLeaf(use(j))));
            %  psi(use(j))=acos(tan(pi/2-the(i)).*cot(theL(use(j))));         
            A(use(j))=cos(theta_rad(i)).*cos(thetaLeaf(use(j))).*(1+2/pi*(tan(psi(use(j)))-psi(use(j))));
        end
        F=A .* pdf;
        
        Gfunc(i)=trapz(thetaLeaf,F);
        
   end
end

