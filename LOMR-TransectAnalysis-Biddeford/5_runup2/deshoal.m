%% Function to de-shoal wave height
    function [H0_E,H0_H]=deshoal(H,T,D,g)
        % function [H0_E,H0_H]=deepwaterWaveheight(H,T,D,g);
        %
        %  H=wave height
        %  T=wave period
        %  D=depth where wave height was recorded
        %  g=gravitational acceleration
        %
        %  H0_E = deepwater wave height Eckart eqn
        %  H0_H = deepwater wave height Hunt eqn
        %
        %  computes deep water wave height from linear theory
        %  and conservation of energy assumptions
        %
        %  uses both the Eckart and Hunt equations
        %  for approximate solution of dispersion relation.
        %
        %  Hunt is probably more accurate.
        %
        %  Reference:
        %
        %  R.G. Dean and R.A. Dalrymple. 2000.  Water
        %  Wave Mechanics for Engineers and Scientists. World
        %  Scientific Publishing Company, River Edge New Jersy
        %
        %
        %  USACE (1985), Direct Methods for Calculating Wavelength, CETN-1-17
        %  US Army Engineer Waterways Experiment Station Coastel Engineering
        %  Research Center, Vicksburg, MS
        %
        %  also see CEM Part II-3 for discussion of shoaling coefficient
        %
        
        twopi=2*3.141592653589793;
        
        % get deep water celerity
        L0=g.*T.*T./twopi;
        C0=L0./T;
        
        % angular frequency
        sigma=twopi ./ T;
        sigmasqOg=sigma.*sigma./g;
        
        % local wave length  (Eckart, 1951)
        k=sigmasqOg./sqrt(tanh(sigmasqOg.*D));
        L=twopi./k;
        
        
        % get local celerity at depth D
        C1=L./T;
        
        % or use Hunt's (1979) approximation for Celerity
        y=sigma.*sigma.*D./g;
        C1H=sqrt( g.*D./(y+1./(1 + 0.6522.*y + 0.4622.*y.^2 + 0.0864.*y.^4 +0.0675.*y.^5)) );
        
        % shoaling coefficient
        Ks=sqrt(C0./C1);
        KsH=sqrt(C0./C1H);
        
        % get deep water wave heigh
        H0_E=H./Ks;
        H0_H=H./KsH;
    end
