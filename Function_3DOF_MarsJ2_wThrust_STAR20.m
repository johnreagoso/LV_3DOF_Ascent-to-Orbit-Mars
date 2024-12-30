function [aState] = Function_3DOF_MarsJ2_wThrust_STAR20(time, aStateInput, pitchdown_rate)
%  State vector         --> [X   Y   Z   Vx  Vy  Vz  mass]
% (d/dt) State vector   --> [Vx  Vy  Vz  AccelX  AccelY  AccelZ  mass_dot] 

Area = MAV_LaunchVehicle_Config.vArea;
Cd    = MAV_LaunchVehicle_Config.Cd;

aX_UnitVec      = [1; 0; 0];    
aY_UnitVec      = [0; 1; 0];    
aZ_UnitVec      = [0; 0; 1];

mass = aStateInput(7);                                  
aMarsRotVel = [0; 0; Mars_GenPhysCons.OMEGA];

Rvector = [aStateInput(1); aStateInput(2); aStateInput(3)];     %km
Vvector = [aStateInput(4); aStateInput(5); aStateInput(6)];     %km/sec

dx  = Vvector(1);   dy  = Vvector(2);   dz  = Vvector(3);

VelVectorUnit = [dx/norm(Vvector); dy/norm(Vvector); dz/norm(Vvector)];

%% Gravitational Acceleration Compute: 

% km/sec^2
J2 = 0.48416685e-03 *sqrt(5)*1.811;

%[latgd, long, alt_m] = ecef2geodetic(Rvector(1)*1e3, Rvector(2)*1e3, Rvector(3)*1e3, referenceEllipsoid('Mars'));  
[latgd, long, ~] = mcmf2geodetic(Rvector);

%grav_accel = (Mars_GenPhysCons.GM_KM/(norm(Rvector))^2);

% J2- geocentric x-comp:
ddxGeoJ2 = (-Mars_GenPhysCons.GM_KM/(norm(Rvector))^2)*(Mars_GenPhysCons.RE_EQ/(norm(Rvector)))^2 *(3*J2*sin(latgd)*cos(latgd));
%ddxGeoJ2 = (-grav_accel)*(Mars_GenPhysCons.RE_EQ/(norm(Rvector)))^2 *(3*J2*sin(latgd)*cos(latgd));

% J2- geocentric y-comp:
ddyGeoJ2 = 0.0;

% J2- geocentric z-comp:
ddzGeoJ2 = (Mars_GenPhysCons.GM_KM/(norm(Rvector))^2)*(1 - (Mars_GenPhysCons.RE_EQ/norm(Rvector))^2 * 3/2 * J2*(3*sin(latgd)^2 - 1));

% sinlatgd2 = 3*sin(latgd)^2 - 1;
% ddzGeoJ2 = grav_accel - grav_accel*(Mars_GenPhysCons.RE_EQ/norm(Rvector))^2 * 3/2 * J2*sinlatgd2;


%ddzGeoJ2 = (grav_accel - grav_accel*(Mars_GenPhysCons.RE_EQ/norm(Rvector))^2 * 3/2 * J2*(3*sin(latgd)^2 - 1));
%ddzGeoJ2 = (grav_accel - grav_accel*(Mars_GenPhysCons.RE_EQ/norm(Rvector))^2 * 3/2 * J2*3*sin(latgd)^2 - J2);


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

%alt_km = alt_m*0.001;

%alt_meters = alt_km*1000;                           % km --> meters
Vvecfor_msec = Vvector*1000;                        % km/sec --> meters/sec

atm_rho = mars_atm_density(alt_km)*515.3790;         %slug/ft^3 --> kg/m^3 

drag_accel_msec2  = 0.5*atm_rho*((norm(Vvecfor_msec))^2)*Cd*Area/mass;    % m/sec^2
drag_accel_magntd = drag_accel_msec2*0.0010;  % m/sec^2 to km/sec^2                                          % km/sec^2

drag_accel_vector = -drag_accel_magntd*VelVectorUnit;

%% Thrusting Profile Acceleration Compute:

ISP  = MAV_LaunchVehicle_Config.ISP1;    
mdot = MAV_LaunchVehicle_Config.mdot1;

%thrust_magntd1 =(0.000133*time^3 - 0.022978*time^2 + 0.9731*time + 17.729)*0.3192;%* 0.266925*1.1; %Newtons
%mdot1 = thrust_magntd1/(0.001*9.80665*ISP);

thrust_magntd = (0.001*9.80665*mdot*ISP);

accel_magntd = thrust_magntd/mass;

%% Thrust Direction compute:
input.E = Rvector(1);    input.F = Rvector(2);    input.G = Rvector(3);
input.dE = Vvector(1);   input.dF = Vvector(2);   input.dG = Vvector(3);

dcm_ecef2ned = [-sin(latgd)*cos(long)   -sin(latgd)*sin(long)  cos(latgd);...
                -sin(long)   cos(long)  0;...
                        -cos(latgd)*cos(long)  -cos(latgd)*sin(long)  -sin(latgd)];  

% rotate ECEF vel to NED vel
ned_vector = dcm_ecef2ned*Vvector;  
% compute az, fpa proxies (analagous body pitch/yaw)
[az, fpa, speed] = cart2sph(ned_vector(1), ned_vector(2), -ned_vector(3));

% compute 'new' vehicle 'pitch' using fpa proxy..

pitchdown_aoa = pitchdown_rate*time;
fpa_new = fpa - deg2rad(pitchdown_aoa); 

% compute updated ned vel:
[ned_dN, ned_dE, ned_dU] = sph2cart(az, fpa_new, speed);

% rotate ned vel back into ecef:
ecef_vel = inv(dcm_ecef2ned)*[ned_dN; ned_dE; -ned_dU];  dx = ecef_vel(1); dy = ecef_vel(2); dz = ecef_vel(3); 

% now compute thrust unit direction (in ecef): 
thrust_unit_dir = [dx; dy; dz]*(1/norm([dx; dy; dz]));

%% compute thrust components in ecef:
thrust_accel_vector = accel_magntd*thrust_unit_dir;
%thrust_accel_vector = accel_magntd*VelVectorUnit;

%% Final Acceleration Compute:

ddx = ddx_mars + drag_accel_vector(1) + thrust_accel_vector(1); 
ddy = ddy_mars + drag_accel_vector(2) + thrust_accel_vector(2); 
ddz = ddz_mars + drag_accel_vector(3) + thrust_accel_vector(3);

aState = [dx; dy; dz; ...
    ddx; ddy; ddz; -mdot];

% disp(alt_km);
% fprintf('%.2f  |  %.4f  |  %.4f  |   %.4f |   %.4f  |  %.4f   |  %.4f \n', time, mass, mdot, mdot1, thrust_magntd, thrust_magntd1 ,alt_km);

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

function [outputStruct] = ECEFtoRL_Convert(inputStruct)

%          inputStruct.E = inputStruct.E*1.0e3;          inputStruct.F = inputStruct.F*1.0e3;          inputStruct.G = inputStruct.G*1.0e3;   
%          inputStruct.dE = inputStruct.dE*1.0e3;        inputStruct.dF = inputStruct.dF*1.0e3;        inputStruct.dG = inputStruct.dG*1.0e3;   

         Rx = inputStruct.E;  Ry = inputStruct.F;  Rz = inputStruct.G;
         Vx = inputStruct.dE; Vy = inputStruct.dF; Vz = inputStruct.dG;
    
        %% Geocentric Latitude:
              inv_f = 169.779286926995;
              a     = Mars_GenPhysCons.RE_EQ*1000;
              b     = a*(1-1/inv_f);
              %ecc   = sqrt(1-(b/a)^2); 

%         marsSpheroid.Name = 'MarsGeodeticApprx';
%         marsSpheroid.LengthUnit = 'meter';
%         marsSpheroid.SemimajorAxis     = wgs84Ellipsoid('meters').SemimajorAxis*0.532;
%         marsSpheroid.SemiminorAxis     = wgs84Ellipsoid('meters').SemiminorAxis*0.531;
%         marsSpheroid.InverseFlattening = wgs84Ellipsoid('meters').InverseFlattening*1.76;
%         marsSpheroid.Eccentricity      = 0.1105; %wgs84Ellipsoid('meters').Eccentricity*

        [lat_rad, long_rad, ~] = ecef2geodetic(inputStruct.E*1e3, inputStruct.F*1e3, inputStruct.G*1e3, referenceEllipsoid('Mars'));      
        lat_deg = lat_rad*180/pi;
        
        %% Longitude:
        long_deg = long_rad*180/pi;
        
        HighVals = find(long_deg>180);        long_deg(HighVals) = long_deg(HighVals)-360;
        
        for ii = 1:1:numel(inputStruct.E)
            DCMecef2ned{ii} = [cos(lat_rad(ii))*cos(long_rad(ii))     cos(lat_rad(ii))*sin(long_rad(ii))      sin(lat_rad(ii));...
                        -sin(long_rad(ii))      cos(long_rad(ii))       0;...
                            -sin(lat_rad(ii))*cos(long_rad(ii))     -sin(lat_rad(ii))*sin(long_rad(ii))     cos(lat_rad(ii));];
                                
            VelNED{ii} = DCMecef2ned{ii}*[Vx(ii); Vy(ii); Vz(ii)];   %ecefVel = transpose(TewariDCM)*[Vn; Ve; Vd];
            PosNED{ii} = DCMecef2ned{ii}*[Rx(ii); Ry(ii); Rz(ii)];
        end
    
        f = 1/169.779286926995; 
        RadiusEq_km = Mars_GenPhysCons.RE_EQ;
        lamda = atan((1-f)^2. *tan(lat_rad));
        Rpos = sqrt(RadiusEq_km^2./(1 + (1/(1-f)^2 - 1)*sin(lamda).^2));  %meters

        for ii = 1:numel(inputStruct.E)
            alt(ii,1) = norm([Rx(ii); Ry(ii); Rz(ii)]) - Rpos(ii);  % km

            Velmag(ii,1) = norm([Vx(ii); Vy(ii); Vz(ii)]);
                
            VelNEDsingle = [VelNED{ii}(1);  VelNED{ii}(2);  VelNED{ii}(3)];
            PosNEDsingle = [PosNED{ii}(1);  PosNED{ii}(2);  PosNED{ii}(3)];

            vmag(ii,1) = norm([VelNED{ii}(1); VelNED{ii}(2); VelNED{ii}(3)]);
            
            Vaz(ii,1) = atan2d(sqrt(VelNED{ii}(1)^2 + VelNED{ii}(2)^2),VelNED{ii}(3)); 
            Vqe(ii,1) = 90 - atan2d(VelNED{ii}(2),VelNED{ii}(1));                   
    
        end

        %% Output

        outputStruct = inputStruct;
        outputStruct.speed = Velmag;
        outputStruct.latitude = lat_deg;
        outputStruct.longitude = long_deg;
        outputStruct.Vqe = Vqe;
        outputStruct.Vaz = Vaz;
        outputStruct.altitude = alt;

        %outputStruct = RL_prep(outputStruct);

end

