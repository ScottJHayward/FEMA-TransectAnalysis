function [ Runup,logstring ] = Aces_Beach_Runup(Hs_dw,Tp,slope)
%ACES_RUNUP Summary
%   Hs_dw = deepwater Wave height
%   Tp    = Peak Wave Period
%   slope = average beach slope
%
%   Will calculate 5 runup values for a given case; 
%   R=[Rmax, R2%, R1/10, R1/3, Rmean]

a=[2.32 1.86 1.70 1.38 0.88];
b=[0.77 0.71 0.71 0.70 0.69];

%constants
g=32.17405; %ft/s^2

%calculate wavelength from period
L_0=g*Tp^2/(2*pi);

%irribarren number
Irb=slope/sqrt(Hs_dw/L_0);

%now calculate the runup heights
Runup=Hs_dw.*a.*Irb.^b;

logstring=sprintf('ACES IRREGULAR WAVE RUNUP ON BEACHES\n\n# Reference:\n# Leenknecht, David A., Andre Szuwaiski, and Ann Sherlock. 1992.\n# "Automated Coastal Engineering System Technical Reference",\n# Coastal Engineering Research Center, Department of the Army\n# Waterways Experiments Station, Corps of Eniggneers, 3909 Halls\n# Ferry Road, Vicksburg, Mississippi 39180-6199.\n\nINPUTS:\n\nAcceleration Due to Gravity,        g  = %8.3f \nDeepwater Significant Wave height,  Hs = %8.2f \nWave Period,                        T  = %8.2f \nBeach Slope,                        S  = %8.3f \n\nEQUATIONS:\n\nRunup,           R   = Hs * a * Irb^b\nIribarren,       Irb = S/sqrt(Hs/L0)\nWavelength,      L0  = g * T^2 / 2 / pi\n\nCOEFFICIENTS:\n(Mase, H. 1989, "Random Wave Runup Height on Gentle Slopes,"\nj. Waterway, Port, Coastal and Ocean Engineering Division,\nASCE, Vol 115, No. 5, pp 649-661.)\n\n               [Rmax,  R2%%, R-1/3, R-1/10, R-mean]\n           a = [2.32, 1.86,  1.70,   1.38,   0.88]\n           b = [0.77, 0.71,  0.71,   0.70,   0.69]\n\nRESULTS:\n\n       RUNUP = [%4.1f, %4.1f, %5.1f, %6.1f, %6.1f]\n',g,Hs_dw,Tp,slope,Runup(1),Runup(2),Runup(3),Runup(4),Runup(5))
 



end

