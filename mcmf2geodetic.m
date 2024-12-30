% mcmf2geodetic.m

% Created 20 Dec 2024, JReagoso
% This script iteratively computes a vehicle's geodetic latitude on Mars
% based on its mars-centered-mars-fixed (planet-fixed) position (similar to
% earth-center-earth-fixed 'ECEF')

% Inputs: Rvector (3x1 vector)/ units: km
% Output: geodetic latitude,longitude and height
%         latgd(rad), longd(rad), h(km)

% Reference: www.oc.nps.edu/oc2902w/coord/coordcvt.pdf

function [latgd_return, longd_return, h] = mcmf2geodetic(Rvector)
    
    % Mars geocentric latitude:
    longc = atan2(Rvector(2),Rvector(1));
    latgc = atan(Rvector(3)/sqrt(Rvector(1)^2 + Rvector(2)^2));

    longd_return = longc;

    a = Mars_GenPhysCons.RE_EQ;
    mars_e = Mars_GenPhysCons.eccen;

    % initialize array for latgd for speed:
    latgd =zeros(25,1);

    % initialize first 'guess' for vehicle Mars geodetic latitude with Mars geocentric latitude:
    latgd(1) = latgc; 
    
    ii = 1; 
   
    z = Rvector(3);
    delta_latgd = 1;

    while delta_latgd > 1.0E-10  % will cycle through until computed latgd iterative delta < 1E-10 (rad)
    
        Rn = a/sqrt(1-mars_e^2 * sin(latgd(ii)));
    
        p = sqrt(Rvector(1)^2 + Rvector(2)^2);
    
        h = p/cos(latgd(ii)) - Rn;
    
        latgd(ii+1) = atan((z/p)*(1 - (mars_e^2)*Rn/(Rn+h))^-1);
    
        delta_latgd = abs(latgd(ii+1) - latgd(ii));

        ii = ii + 1;

    end

    latgd_return = latgd(ii);
end