function [ Gfunc_Sky ] = fGfunc_diff(thetaSky, thetaLeaf, pdf )
%
%   Detailed explanation goes here
 % evaluate Gfunc Sky for all the sky zones of the hemisphere and
        % assume azimuthal symmetry.
        
         %  Use Gfunc_sky to compute the exp Beers law probability of
        %  penetration from the sky for an incremental layer of leaves
        
         % Pdiff = 2 int 0 to pi/2 P0(theta) cos(theta) sin(theta) dtheta
        
        % think about this and is Gfunc sky 0 to pi/2 or pi
        
        % is it integrated 2 pi (360) by pi/2 (90)
        
        % or pi (180) by pi (180)??
        
        Adiff=zeros(length(thetaSky),length(thetaLeaf));
        F=zeros(length(thetaSky),length(thetaLeaf));
        Gfunc_Sky=zeros(length(thetaSky),1);
                        
        for i=1:length(thetaSky) % 0 to pi/2
            
        use=find((abs(cot(thetaSky(i)).*cot(thetaLeaf))>1)==1);
        for j=1:length(use)
            Adiff(i, use(j))=cos(thetaSky(i)).*cos(thetaLeaf(use(j)));
        end
        use=find((abs(cot(thetaSky(i)).*cot(thetaLeaf))>1)==0);
        for j=1:length(use)
            psi(use(j))=acos(cot(thetaSky(i)).*cot(thetaLeaf(use(j))));
            %     psi(use(j))=acos(tan(pi/2-the(i)).*cot(theL(use(j))));         
            Adiff(i, use(j))=cos(thetaSky(i)).*cos(thetaLeaf(use(j))).*(1+2/pi*(tan(psi(use(j)))-psi(use(j))));
        end
        F(i, :)=Adiff(i, :).* pdf';  % needed to transpose pdf matrix for multiplication
                
        Gfunc_Sky(i)=trapz(thetaLeaf,F(i, :));  % trapz is trapezoidal integration
        end
        
       

end

