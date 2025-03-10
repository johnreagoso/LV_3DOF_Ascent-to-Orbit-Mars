classdef Mars_GenPhysCons

    % Mars_GenPhysCons
    % Created 25 Oct 2024, JReagoso

    % Purpose: for use with flight dynamics tools for Mars ascent flight trajectory propagation, design and analysis

    properties (Constant)
   
        %% Mars Geophysical constants
        
        GRAV_ACCEL = 3.72076e0; 							%  Mars's gravitational acceleration at 45 deg lat. [m/s^2]
        GM_KM = 4.282837e4;								    %  Mars's gravitional parameter ("mu" or "k") [km^3/s^2]
        GM_M  = 4.282837e13;								%  Mars's gravitional parameter ("mu" or "k") [m^3/s^2]
        RE_EQ = 3396.2;									    %  Mars' equatorial radius [km]
        RE_PL = 3376.2;                                     %  Mars' polar radius [km]
        RE_VL = 3389.5;                                     %  Mars' mean (volumetric) radius [km]
        ALT_RE = 400000*12*2.54/100000.;					%  Navy standard reentry altitude, 400 kft [km]
        ROTRATE = 0.00417807462229;                         %  Mars' rotational rate (sidereal); [deg/s]
        OMEGA = 1.029*7.292115e-005;                        %  Mars' nominal mean rotational rate (sidereal) [rad/s]
        inv_f = 169.779286926995;                           % flattening value
        f = 1/169.779286926995;                             % inverse flattening
        J2 = 0.00196064291017307; % 0.48416685e-03 *sqrt(5)*1.811;  % Mars J2 parameter
        eccen = 0.10837577173889;
        
        %%	Astrophysical constants
        AU = 1.49597870700e08;                              %  Astronomical unit [km] (IAU, Jul 2012; now an exact value)

    end
    
    %=====================================================================================================
    methods (Static)
        % Created: 31 Oct 2024  J Reagoso

        function rho = mars_atm_density(altitude_km)  %alt in km(s)
        % Computes Mars generic atmospheric density
        % ref: https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mars-atmosphere-equation-english/

        % Input: altitude (km)
        % Output: rho (slug/ft^3)

            altitude_m = altitude_km*1e3;
        
            if altitude_m > 7000     
                temp  = -23.4 - 0.00222*altitude_m;
                press = 0.6990*exp(-0.00009*altitude_m);
            else 
                temp  = -31.0 - 0.000998*altitude_m;
                press = 0.699*exp(-0.00009*altitude_m);
            end

            denom_ft3 = 0.1921*temp + 52.46251; %(0.1921*(temp + 273.1))
            rho = press/denom_ft3;  %slugs/ft^3
        
        end

    end % End of methods

end  

%=======================================================================================================================
    % REFERENCES:
    %
    % http://physics.nist.gov/cuu/Constants/index.html              %  Fundamental constants (CODATA, 2014)
    % http://http://maia.usno.navy.mil/NSFA/NSFA_cbe.html#GME2009   
    %                                                               %  IAU 2009, Numerical Standards for Astronomy. See also 
    %                                                               %  references therein.
    % http://ssd.jpl.nasa.gov/?constants                            %  Defined, derived, and primary constants for solar system
	% https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html


