
function [aState] = Function_3DOF_MarsJ2(~, aStateInput, ~)
%  State vector         --> [X   Y   Z   Vx  Vy  Vz  mass] (aStateInput)
% (d/dt) State vector   --> [Vx  Vy  Vz  AccelX  AccelY  AccelZ  mass_dot] 

Area = MAV_LaunchVehicle_Config.vArea;
Cd    = MAV_LaunchVehicle_Config.Cd;

aX_UnitVec      = [1; 0; 0];    
aY_UnitVec      = [0; 1; 0];    
aZ_UnitVec      = [0; 0; 1];

mass = aStateInput(7);                                  

aMarsRotVel = [0; 0; Mars_GenPhysCons.OMEGA];
%aMarsRotVel = [0; 0; 0];

Rvector = [aStateInput(1); aStateInput(2); aStateInput(3)];     %km
Vvector = [aStateInput(4); aStateInput(5); aStateInput(6)];     %km/sec

dx  = Vvector(1);   dy  = Vvector(2);   dz  = Vvector(3);

VelVectorUnit = [dx/norm(Vvector); dy/norm(Vvector); dz/norm(Vvector)];

%% Gravitational Acceleration Compute: 

% km/sec^2
J2 = 0.00196064291017307; % 0.48416685e-03 *sqrt(5)*1.811;

%[latgd1, long1, ~] = ecef2geodetic(Rvector(1)*1e3, Rvector(2)*1e3, Rvector(3)*1e3, referenceEllipsoid('Mars'));      
[latgd, long, ~] = mcmf2geodetic(Rvector);

% J2- geocentric x-comp:
ddxGeoJ2 = (-Mars_GenPhysCons.GM_KM/(norm(Rvector))^2)*(Mars_GenPhysCons.RE_EQ/(norm(Rvector)))^2 *(3*J2*sin(latgd)*cos(latgd));

% J2- geocentric y-comp:
ddyGeoJ2 = 0.0;

% J2- geocentric z-comp:
ddzGeoJ2 = (Mars_GenPhysCons.GM_KM/(norm(Rvector))^2)*(1 - (Mars_GenPhysCons.RE_EQ/norm(Rvector))^2 * 3/2 * J2*(3*sin(latgd)^2 - 1));

[ddxEciJ2, ddyEciJ2, ddzEciJ2] = ...
    ned2efg_vector(latgd,long,ddxGeoJ2, ddyGeoJ2, ddzGeoJ2);

ddx_earth = ddxEciJ2 - ...
    dot(cross(2*aMarsRotVel, Vvector), aX_UnitVec)...   
            - dot(cross(aMarsRotVel,(cross(aMarsRotVel, Rvector))), aX_UnitVec);       
        
ddy_earth = ddyEciJ2 - ...
    dot(cross(2*aMarsRotVel, Vvector), aY_UnitVec)...
            - dot(cross(aMarsRotVel,(cross(aMarsRotVel, Rvector))), aY_UnitVec);
        
ddz_earth = ddzEciJ2 -...
    dot(cross(2*aMarsRotVel, Vvector), aZ_UnitVec)...
            - dot(cross(aMarsRotVel,(cross(aMarsRotVel, Rvector))), aZ_UnitVec);

%% Atmospheric Drag Acceleration Compute:

f = Mars_GenPhysCons.f; 

RadiusEq_km = Mars_GenPhysCons.RE_EQ;

lamda = atan((1-f)^2 *tan(latgd));

Rpos_denom = (1 + (1/(1-f)^2 - 1)*sin(lamda)^2);

Rpos_quant = RadiusEq_km^2/Rpos_denom;
Rpos = sqrt(Rpos_quant);  %meters

%Rpos = sqrt(RadiusEq_km^2/(Rpos_denom));  %meters

alt_km = norm([Rvector(1); Rvector(2); Rvector(3)]) - Rpos;  % km

Vvecfor_msec = Vvector*1000;                        % km/sec --> meters/sec

atm_rho = mars_atm_density(alt_km)*515.3790;         %slug/ft^3 --> kg/m^3 

drag_accel_msec2  = 0.5*atm_rho*((norm(Vvecfor_msec))^2)*Cd*Area/mass;    % m/sec^2
drag_accel_magntd = drag_accel_msec2*0.0010;                                            % km/sec^2

drag_accel_vector = -drag_accel_magntd*VelVectorUnit;

%% Final Acceleration Compute:

ddx = ddx_earth + drag_accel_vector(1); % + thrust_accel_vector(1); 
ddy = ddy_earth + drag_accel_vector(2); % + thrust_accel_vector(2); 
ddz = ddz_earth + drag_accel_vector(3); % + thrust_accel_vector(3);

aState = [dx; dy; dz; ...
    ddx; ddy; ddz; 0];

% fprintf('%.2f  |  %.4f  |  %.4f  |   %.4f  |  %.4f  \n', time, mass, drag_accel_magntd, thrust_magntd, alt_km);
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

quant = 1149*temp + 528195.3;

%rho = press/(1149*(temp + 459.7));  %slugs/ft^3
rho = press/quant;  %slugs/ft^3
end
