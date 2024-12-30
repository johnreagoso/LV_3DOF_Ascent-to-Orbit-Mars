function [aState] = Function_3DOF_MarsJ2_wThrust_STAR15G(~, aStateInput, ~)
%  State vector         --> [X   Y   Z   Vx  Vy  Vz  mass]
% (d/dt) State vector   --> [Vx  Vy  Vz  AccelX  AccelY  AccelZ  mass_dot] 

aX_UnitVec      = [1; 0; 0];    
aY_UnitVec      = [0; 1; 0];    
aZ_UnitVec      = [0; 0; 1];

vArea = MAV_LaunchVehicle_Config.vArea;
Cd    = MAV_LaunchVehicle_Config.Cd;

mass = aStateInput(7);                                  
aMarsRotVel = [0; 0; Mars_GenPhysCons.OMEGA];

Rvector = [aStateInput(1); aStateInput(2); aStateInput(3)];     %km
Vvector = [aStateInput(4); aStateInput(5); aStateInput(6)];     %km/sec

dx  = Vvector(1);   dy  = Vvector(2);   dz  = Vvector(3);

VelVectorUnit = [dx/norm(Vvector); dy/norm(Vvector); dz/norm(Vvector)];

%% Gravitational Acceleration Compute: 

% km/sec^2
J2 = 0.48416685e-03 *sqrt(5)*1.811;

%[latgd, long, alt_meters] = ecef2geodetic(Rvector(1)*1e3, Rvector(2)*1e3, Rvector(3)*1e3, referenceEllipsoid('Mars'));  
[latgd, long, ~] = mcmf2geodetic(Rvector);

% J2- geocentric x-comp:
ddxGeoJ2 = (-Mars_GenPhysCons.GM_KM/(norm(Rvector))^2)*(Mars_GenPhysCons.RE_EQ/(norm(Rvector)))^2 *(3*J2*sin(latgd)*cos(latgd));

% J2- geocentric y-comp:
ddyGeoJ2 = 0.0;

% J2- geocentric z-comp:
ddzGeoJ2 = (Mars_GenPhysCons.GM_KM/(norm(Rvector))^2)*(1 - (Mars_GenPhysCons.RE_EQ/norm(Rvector))^2 * 3/2 * J2*(3*sin(latgd)^2 - 1));

[ddxEciJ2, ddyEciJ2, ddzEciJ2] = ...
    ned2efg_vector(latgd,long,ddxGeoJ2, ddyGeoJ2, ddzGeoJ2);

ddx_mars = ddxEciJ2 - ...
    dot(cross(2*aMarsRotVel, Vvector), aX_UnitVec)...   
            - dot(cross(aMarsRotVel,(cross(aMarsRotVel, Rvector))), aX_UnitVec);       
        
ddy_mars = ddyEciJ2 - ...
    dot(cross(2*aMarsRotVel, Vvector), aY_UnitVec)...
            - dot(cross(aMarsRotVel,(cross(aMarsRotVel, Rvector))), aY_UnitVec);
        
ddz_mars = ddzEciJ2 -...
    dot(cross(2*aMarsRotVel, Vvector), aZ_UnitVec)...
            - dot(cross(aMarsRotVel,(cross(aMarsRotVel, Rvector))), aZ_UnitVec);

%% Atmospheric Drag Acceleration Compute:

alt_km = (norm(Rvector) - Mars_GenPhysCons.RE_EQ) + 1.9626; 

Vvecfor_msec = Vvector*1000;                        % km/sec --> meters/sec

atm_rho = mars_atm_density(alt_km)*515.3790;         %slug/ft^3 --> kg/m^3 

drag_accel_msec2  = 0.5*atm_rho*((norm(Vvecfor_msec))^2)*Cd*vArea/mass;    % m/sec^2
drag_accel_magntd = drag_accel_msec2*0.0010;  % m/sec^2 to km/sec^2                                          % km/sec^2

drag_accel_vector = -drag_accel_magntd*VelVectorUnit;

%% Thrusting Profile Acceleration Compute:

ISP  = MAV_LaunchVehicle_Config.ISP2;    
mdot = MAV_LaunchVehicle_Config.mdot2;

% burntime = time-basetime;
% adjust_coeff = 1.25;
% 
% if burntime <= 14.0
%     thrust_magntd1 =  adjust_coeff*(-0.0014*burntime^3 + 0.0077*burntime^2 + 0.6276*burntime + 5.7752); % kNewtons
% else
%     thrust_magntd1 =  adjust_coeff*(0.2392*burntime^2 - 9.61*burntime+ 100.95);                     % kNewtons
% end

%mdot1 = thrust_magntd1/(0.001*9.80665*ISP);

thrust_magntd = (0.001*9.80665*mdot*ISP);
accel_magntd = thrust_magntd/mass;

thrust_accel_vector = accel_magntd*VelVectorUnit;

%% Final Acceleration Compute:

ddx = ddx_mars + drag_accel_vector(1) + thrust_accel_vector(1); 
ddy = ddy_mars + drag_accel_vector(2) + thrust_accel_vector(2); 
ddz = ddz_mars + drag_accel_vector(3) + thrust_accel_vector(3);

aState = [dx; dy; dz; ...
    ddx; ddy; ddz; -mdot];

%disp(alt_km);

%fprintf('%.2f  |  %.4f  |  %.4f  |   %.4f |   %.4f  |  %.4f   |  %.4f \n', time, mass, mdot, mdot1, thrust_magntd, thrust_magntd1 ,alt_km);
end

%% Atmospheric Density Calculation (rho per alt) 

function rho = mars_atm_density(altitude_km)  %alt in km(s)

    altitude_ft = altitude_km*3280.84; % km to feet 
    
    if altitude_ft > 22960     
        temp  = -10.34 - 0.001217*altitude_ft;
        press = 14.62*exp(-0.00003*altitude_ft);
            
    elseif altitude_ft < 0
        temp  = -10.34 - 0.001217*altitude_ft;
        press = 0.0;
    else 
        temp  = -25.68 - 0.000548*altitude_ft;
        press = 14.62*exp(-0.00003*altitude_ft);
    end

rho = press/(1149*(temp + 459.7));  %slugs/ft^3

end

