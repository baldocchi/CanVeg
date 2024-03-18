function [fT]= fd2ESdT (T)

%Second derivative of the saturation vapor pressure curve

param.Rstar=8.314;

a3en=1.675;
a4en=0.01408;
a5en=0.0005818;


%fT = -2. * ES(T) * LAMBDA(T) * 18. / (param.Rstar * T * T * T * 1000.);
 tcel=T-273.16;

        fT=2*a3en+6*a4en*tcel+12*a5en*tcel.*tcel;


end


