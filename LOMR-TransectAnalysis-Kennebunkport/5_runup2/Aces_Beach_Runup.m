function [ Runup ] = Aces_Beach_Runup(Hs_dw,Tp,slope)
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

end

