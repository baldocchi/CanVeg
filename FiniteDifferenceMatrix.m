function [soil] = FiniteDifferenceMatrix(soil,prm)
    
    % Convert Soil Physics with Basic from Campbell for soil heat flux
    
   % March 2, 2021

   % Using finite difference approach to solve soil heat transfer
   % equations, but using geometric spacing due to the exponential decay in
   % heat and temperature with depth into the soil
   
   % Soil has N layers and we increment from 1 to N+1
    
   % old code from python and basic based computations on N layers
   % from levels 0 to N
    
   % Lower Boundary is level 1, upper Boundary is level N+1
    
   
   
   
dt=10;                           
timeStepMax=floor(1800/dt); 
area = 1 ;                   	
maxNrIterations = timeStepMax;

tolerance = 1.e-2  ;  

% pre allocate space for coefficients

a_soil=zeros(prm.nn,soil.n_soil_2);
b_soil=zeros(prm.nn,soil.n_soil_2);
c_soil=zeros(prm.nn,soil.n_soil_2);
d_soil=zeros(prm.nn,soil.n_soil_2);
mm=zeros(prm.nn,soil.n_soil_2);
    
    Fst = 0.7;   % (0: explicit, 1: implicit Euler),,C code used 0.6

    Gst = 1.0 - Fst  ;   
    
    energyBalance = 1.;
    
       

     %  in C code boundary condition was [0]. Matlab does not deal with
     %  zero index, so assume boundary at top of soil is (1)...so increment
     
     %%%%%%%% level 1, Upper Boundary
     
     % layer 1
     
     %%%%%%%% level 2
     
     %%%  .....
     
     % layer n
     
     %%%%%%%%% level n+1, lowest boundary
     
     %  revisit old Basic code and adapt
     
     % conditions A1 = 0
     % Cn=0
        
    nrIterations = 0;
    
    % looping to solve the Fourier heat transfer equation, dT/dt ~ dH/dz
    while ((energyBalance > tolerance) && (nrIterations < maxNrIterations))
       
       
              
       for I = 1: soil.n_soil   %     soil.n_soil   // define coef for each soil layer
        
        IM1=I-1;
        IP1=I+1;
      
        c_soil(:,I) = -soil.k_conductivity_soil(:,I) * Fst;
               
        if I == 1
            
%             b_soil(:,1) = Fst .* (soil.k_conductivity_soil(:,2) + soil.k_conductivity_soil_bound) + soil.cp_soil(:,1);  
%             
%             b_soil(:,I) =1; 
%             
%             d_soil(:,1) = soil.T_soil_up_boundary .*soil.k_conductivity_soil_bound  .* Gst  ...
%             + soil.T_soil(:,1) .* soil.cp_soil(:,1) ...                               
%             - soil.T_soil(:,1) .* soil.k_conductivity_soil(:,1) .* Gst ...
%             - soil.T_soil(:,1) .* soil.k_conductivity_soil_bound  .*Gst  ...           
%             + soil.T_soil(:,2) .* soil.k_conductivity_soil(:,1) .* Gst; 
                   
             % upper boundary condition
              
%             d_soil(:,1)=d_soil(:,1) + Fst .* soil.T_soil_up_boundary .* soil.k_conductivity_soil_bound ;
             
                      
            % when fluxes at the surface are computed use 
            % d_soil(:,1)=d_soil(:,1) + Fst .* soil.T_soil_up_boundary .* soil.k_conductivity_soil_bound - Rnet + LE ;
            
            % It works Finally
             
            a_soil(:,I)=0; 
            b_soil(:,I) =1; 
            c_soil(:,I)=0;  
            d_soil(:,I)=soil.T_soil_up_boundary ;
               
        else
            
            a_soil(:,I) = -soil.k_conductivity_soil(:,IM1) * Fst;
            
            b_soil(:,I) = Fst .* (soil.k_conductivity_soil(:,I) + soil.k_conductivity_soil(:,IM1)) + soil.cp_soil(:,I);
            
            d_soil(:,I) =  soil.T_soil(:,IM1).* soil.k_conductivity_soil(:,IM1) .*Gst ...
                  + soil.T_soil(:,I) .* soil.cp_soil(:,I )...
                   - soil.T_soil(:,I) .* soil.k_conductivity_soil(:,I) .* Gst ...
                  -  soil.T_soil(:,I) .* soil.k_conductivity_soil(:,IM1) .* Gst  ...
                  +  soil.T_soil(:,IP1) .* soil.k_conductivity_soil(:,I) .* Gst ;
                  
               
        end
       end
          
       
       % boundary conditions
            
         d_soil(:,soil.n_soil)= d_soil(:,soil.n_soil) + Fst .* soil.k_conductivity_soil(:,soil.n_soil).* soil.T_soil_low_bound;
       
       
    
       % Thomas Algorithm to solve for simulataneous equ
       
%     Sub TriDiagonal_Matrix_Algorithm(N%, A#(), B#(), C#(), D#(), X#())
%     Dim i%, W#
%     For i = 2 To N
%         W = A(i) / B(i - 1)
%         B(i) = B(i) - W * C(i - 1)
%         D(i) = D(i) - W * D(i - 1)
%     Next i
%     X(N) = D(N) / B(N)
%     For i = N - 1 To 1 Step -1
%         X(i) = (D(i) - C(i) * X(i + 1)) / B(i)
%     Next i
% End Sub

            
       
       
            c_soil(:,soil.n_soil)=0;     
            
            
            for i=2: soil.n_soil
            mm(:,i)=a_soil(:,i)./b_soil(:,i-1);
            b_soil(:,i)=b_soil(:,i) - mm(:,i).* c_soil(:,i-1);
            d_soil(:,i)=d_soil(:,i) - mm(:,i).* d_soil(:,i-1);
            end

                  
            soil.T_soil(:,soil.n_soil)=d_soil(:,soil.n_soil)./b_soil(:,soil.n_soil);
            
     
            % back substitution
            
            for i=soil.n_soil-1:-1:1  % 
            soil.T_soil(:,i)=(d_soil(:,i) - c_soil(:,i) .* soil.T_soil(:,i+1)) ./ b_soil(:,i) ; 
            end
       
          
       % test iterations
      
       % use sum function...
           
        dSum=sum(soil.cp_soil(:,2:soil.n_soil).*(soil.T_soil(:,2:soil.n_soil)-soil.T_soil_old(:,2:soil.n_soil)),2);
       
        dTest = ( soil.k_conductivity_soil(:,1) .*(soil.T_soil(:,1)- soil.T_soil(:,2)) .* dt ...
             + (soil.T_soil_up_boundary - soil.T_soil(:,soil.n_soil-1)) .* soil.k_conductivity_soil(:,soil.n_soil-1) .* dt);
         
         energyBalance=nanmean(abs(dSum-dTest));
            
%                 for i in range(2, n):
%             dSum += C_T[i]*(T[i]-oldT[i])
%         energyBalance = (abs(dSum - f[1]*(T[1]-T[2])*dt 
%                            + f[n-1]*(T[n-1]-boundaryT)*dt))

                
                
                nrIterations = 1 +nrIterations; 
                
        
        
        soil.T_soil_old=soil.T_soil;   % double check but if iterating need to reset Told

    end   % end While
    
    
    
    
    soil.gsoil = soil.k_conductivity_soil(:,1) .*(soil.T_soil(:,1)-soil.T_soil(:,2));       % soil heat flux
    
    % add storage term, too..
    
    
    

    
    
end

