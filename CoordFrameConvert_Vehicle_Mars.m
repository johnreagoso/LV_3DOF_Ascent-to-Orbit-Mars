%% ConvertStructure using LV_3DOF Output

% Date created: 
% 24 Oct 2024   - Updated/replaced Earth specific constants and logic with
%                 Mars specific data

% Note: this conversion script origin is for Earth centered trajectories (ECEF, ECI etc.). Any use of 'ECI', 'ECEF' etc. is a holdover 
% as this script is used for MCI (Mars Centered Inertial) or MCMF (Mars Centered Mars Fixed) coordinate frames. Future versions of this script will
% be modified accordingly. 

% Reference for Mars physical/gravity parameters:  https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html 

function [vehicleObj]= CoordFrameConvert_Vehicle_Mars(vehicleObjInput, frame_input)

    if strcmp(frame_input, 'RL') == 1 ||  strcmp(frame_input, 'rl') == 1
    
        mars_const = Mars_GenPhysCons();

        vehicleObj = RL_prep(vehicleObjInput);

%     %% Mars Geodetic/Geocentric Constants: 
        
           inv_f = mars_const.inv_f;
           a     = mars_const.RE_EQ*1.0e3; 
%           %a    = Mars_GenPhysCons.RE_VL*1.0e3;              % aligned with LV_3DOF Mars gravity model (semi-major axis- equatorial radius- km)
            b    = a*(1-1/inv_f);      % Mars semi-minor axis (km)
           ecc   = sqrt(1-(b/a)^2);    % used for a non-spherical Mars model 
   
%% Call RL_to_ECEF Conversion subroutine:
        vehicleObj = RLtoECEF_Convert(vehicleObj);
    
    elseif strcmp(frame_input, 'ECEF') == 1 ||  strcmp(frame_input, 'ecef') == 1
        vehicleObj = ECEFtoRL_Convert(vehicleObjInput);
    end

%% Call ECEF_to_ECI Conversion subroutine:
    vehicleObj = ECEFtoECI_Convert(vehicleObj);             
    
%% Call NEDtoBody Conversion subroutine:
    vehicleObj = NEDtoBody_Convert(vehicleObj);  
    
%% Metric to US standard:
    vehicleObj.altitude_ft = vehicleObj.altitude*3280.84; % converts from km to feet
    vehicleObj.speed_ftsec = vehicleObj.speed*3280.84;    % converts from km/sec to ft/sec        

%% Dynamic Pressure/ Max-Q Compuation:

  for mmm = 1:numel(vehicleObj.time)

      atmRho_slugft3(mmm)= mars_atm_density_kg_m3(vehicleObj.altitude(mmm)); %/515.4 ;

      if vehicleObj.altitude(mmm) > 100
        atmRho_slugft3(mmm) = 1.35066e-8;
      end

      if vehicleObj.altitude(mmm)>300
        atmRho_slugft3(mmm) = 0.00;
      end
        
      vehicleObj.Q_psf(mmm,1) = 0.5*atmRho_slugft3(mmm)*(vehicleObj.speed_ftsec(mmm))^2; % converts km/sec to ft/sec
  end 
      
  vehicleObj.MaxQ = max(vehicleObj.Q_psf); % lbf/ft2
  
  end


%% Plotting:
  %outputStatus = PlotCall(vehicleObj, 'yes'); 

function [status] = AddPitchoverAndStagingTimes(~, ~, ~)
        
    line([75 75],   [0 100], 'Color', 'k', 'LineStyle', '--');   
    line([2 2],         [0 100], 'Color', 'b', 'LineStyle', '--'); 
    line([5.5 5.5],     [0 100], 'Color', 'b', 'LineStyle', '--'); 
     
    text(2.2, 72, 'Pitchover Start', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
    text(5.7, 71, 'Pitchover End', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
    status = 'complete';
    
end


function [outputStruct] = RL_prep(vehicleObjInput)

    %% Geodetic/Geocentric Constants: 
          inv_f = 169.779286926995;
          a     = Mars_GenPhysCons.RE_EQ*1.0e3; 
          %a    = Mars_GenPhysCons.RE_VL*1.0e3;             
           b    = a*(1-1/inv_f);      % Mars semi-minor axis (km)
          ecc   = sqrt(1-(b/a)^2);    % used for a non-spherical Mars model 
    
          vehicleObj = vehicleObjInput;

    %% Are we using rotating Mars or not?
%           vehicleObj.rotatingMars = vehicleObjInput.rotatingMars;
%           vehicleObj.time      = vehicleObjInput.time;  
%           vehicleObj.speed     = vehicleObjInput.speed;  
%           vehicleObj.latitude  = vehicleObjInput.latitude;
%           vehicleObj.longitude = vehicleObjInput.longitude;
%           vehicleObj.altitude  = vehicleObjInput.altitude;
%           vehicleObj.Vqe       = vehicleObjInput.Vqe;              % vertical FPA
%           vehicleObj.Vaz       = vehicleObjInput.Vaz;
          vehicleObj.aoa(1:numel(vehicleObj.time),:) = 0.0;     % this AoA (pitch-down) data is a holdover/placeholder from the Ascent-to-Orbit 
%                                                                 % analysis using a modified version of RotatingEart that provides inputs to 3rd stage burn thrust-angle(s). 
%                                                                 
    %% Compute Range using Sodano's Method:
         vehicleObj.range_km = zeros(numel(vehicleObj.time),1);
         
         for kkk = 1:numel(vehicleObj.time)
            %vehicleObj.range_Sod_km(kkk) = sodanosInverseMethod_Mars(vehicleObj.latitude(1), vehicleObj.longitude(1), vehicleObj.latitude(kkk), vehicleObj.longitude(kkk));
            vehicleObj.range_km(kkk) = great_circ_dist(vehicleObj.latitude(1), vehicleObj.longitude(1), vehicleObj.latitude(kkk), vehicleObj.longitude(kkk));
            
            if isnan(vehicleObj.range_km(kkk)) 
                vehicleObj.range_km(kkk) = 0.00;
            end   
         end
         vehicleObj.range_nm = vehicleObj.range_km*0.539957;
         
    %% RL to NED Cartesian for Velocity and Acceleration:
           [vehicleObj.ned_dN, vehicleObj.ned_dE, vehicleObj.ned_dU] = sph2cart(vehicleObj.Vaz.*pi/180, vehicleObj.Vqe.*pi/180, vehicleObj.speed);   %km/sec
           vehicleObj.ned_dD = -vehicleObj.ned_dU;
           
           %   Using Matlab 'diff' functionality, we take the time derivative of NED velocity:
           vehicleObj.ned_ddN = diff(vehicleObj.ned_dN)./diff(vehicleObj.time);        %km/sec    
           vehicleObj.ned_ddE = diff(vehicleObj.ned_dE)./diff(vehicleObj.time);        %km/sec  
           vehicleObj.ned_ddD = diff(vehicleObj.ned_dD)./diff(vehicleObj.time);        %km/sec  
                  
           vehicleObj.ned_accel_norm = (vehicleObj.ned_ddN.^2 + vehicleObj.ned_ddE.^2 + vehicleObj.ned_ddD.^2).^0.5 ;
    outputStruct = vehicleObj;

end

%% RL to ECEF Conversion: 
function [outputStruct] = RLtoECEF_Convert(inputStruct)
 global a  
 global ecc
 % global b
 % global inv_f

  % ECEF to NED Conversion matrix derived from TAOS manual and matches with a ballistic TAOS comparison. 
  % Also matches with Matlab technical page: 
  % https://www.mathworks.com/help/aeroblks/directioncosinematrixeceftoned.html
  % code results below matches with Matlab geodetic2ecef functionality:

  inv_f = Mars_GenPhysCons.inv_f; 
  a     = Mars_GenPhysCons.RE_EQ*1000;
  b     = a*(1-1/inv_f);
  ecc   = sqrt(1-(b/a)^2); 

  r = 1.0E-3 *(a./sqrt(1-ecc.^2*sind(inputStruct.latitude).^2)); 
  inputStruct.E = (r + inputStruct.altitude).*cosd(inputStruct.latitude).*cosd(inputStruct.longitude);
  inputStruct.F = (r + inputStruct.altitude).*cosd(inputStruct.latitude).*sind(inputStruct.longitude);
  inputStruct.G = ((1 - ecc^2).*r + inputStruct.altitude).*sind(inputStruct.latitude);

%   spheroid = referenceEllipsoid('WGS 84');
%   [inputStruct.E, inputStruct.F, inputStruct.G] = geodetic2ecef(inputStruct.latitude*pi/180, inputStruct.longitude*pi/180,... 
%     inputStruct.altitude, spheroid);
    
  for jjj = 1:numel(inputStruct.time)  

      lat_rad  = inputStruct.latitude(jjj).*pi/180;   %convert deg inputs to radians
      long_rad = inputStruct.longitude(jjj).*pi/180;

      % Below DCM is a 3-2-1 sequence conversion matrix to rotate from
      % ECEF to NED (or vice versa via computing the inverse of this DCM):
      dcm_ecef2ned = [-sin(lat_rad)*cos(long_rad)   -sin(lat_rad)*sin(long_rad)  cos(lat_rad);...
                        -sin(long_rad)   cos(long_rad)  0;...
                                -cos(lat_rad)*cos(long_rad)  -cos(lat_rad)*sin(long_rad)  -sin(lat_rad)];  

      ecefVel = inv(dcm_ecef2ned)*[inputStruct.ned_dN(jjj); inputStruct.ned_dE(jjj); inputStruct.ned_dD(jjj)]; %#ok<MINV>        

      inputStruct.dE(jjj, 1) = ecefVel(1); % km/sec
      inputStruct.dF(jjj, 1) = ecefVel(2); % km/sec
      inputStruct.dG(jjj, 1) = ecefVel(3); % km/sec      

           if jjj == numel(inputStruct.time)
              break;
           end
  end

  % Time derivative method- this produces the best results that match with other tools (TAOS etc.) 
  inputStruct.ddE = diff(inputStruct.dE)./diff(inputStruct.time);                  %km/sec2  
  inputStruct.ddF = diff(inputStruct.dF)./diff(inputStruct.time);          
  inputStruct.ddG = diff(inputStruct.dG)./diff(inputStruct.time);

  outputStruct = inputStruct;
end



function [outputStruct] = ECEFtoRL_Convert(inputStruct)

%          inputStruct.E = inputStruct.E*1.0e3;          inputStruct.F = inputStruct.F*1.0e3;          inputStruct.G = inputStruct.G*1.0e3;   
%          inputStruct.dE = inputStruct.dE*1.0e3;        inputStruct.dF = inputStruct.dF*1.0e3;        inputStruct.dG = inputStruct.dG*1.0e3;   

         Rx = inputStruct.E;  Ry = inputStruct.F;  Rz = inputStruct.G;
         Vx = inputStruct.dE; Vy = inputStruct.dF; Vz = inputStruct.dG;
    
        %% Geodetic Latitude:
              inv_f = Mars_GenPhysCons.inv_f; 
              a     = Mars_GenPhysCons.RE_EQ*1000;
              b     = a*(1-1/inv_f);
              %ecc   = sqrt(1-(b/a)^2); 

        [lat_rad, long_rad, ~] = ecef2geodetic(inputStruct.E*1e3, inputStruct.F*1e3, inputStruct.G*1e3, referenceEllipsoid('Mars'));      
        lat_deg = lat_rad*180/pi;
        
        %% Longitude:
        long_deg = long_rad*180/pi;
        
        HighVals = find(long_deg>180);        long_deg(HighVals) = long_deg(HighVals)-360;
        
        for ii = 1:1:numel(inputStruct.E)
            DCMecef2ned{ii} = [cos(lat_rad(ii))*cos(long_rad(ii))     cos(lat_rad(ii))*sin(long_rad(ii))      sin(lat_rad(ii));...
                        -sin(long_rad(ii))      cos(long_rad(ii))       0;...
                            -sin(lat_rad(ii))*cos(long_rad(ii))     -sin(lat_rad(ii))*sin(long_rad(ii))     cos(lat_rad(ii));];
                                
            VelNED{ii} = DCMecef2ned{ii}*[Vx(ii); Vy(ii); Vz(ii)];  
            PosNED{ii} = DCMecef2ned{ii}*[Rx(ii); Ry(ii); Rz(ii)];
        end
    
        f = Mars_GenPhysCons.f;  
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

        outputStruct = RL_prep(outputStruct);

end

%% ECEF to ECI Conversion:
function [outputStruct] = ECEFtoECI_Convert(inputStruct)

    if inputStruct.rotatingMars == 1 
        marsRotRate = Mars_GenPhysCons.OMEGA;  %rad/sec
    else
        marsRotRate = 0.00;  %rad/sec  No Mars rotation modeled!
    end
      
      elapsedTOF = inputStruct.time;
      
      omega_init = 0;   % !!! for our SLBM trajectory purposes- at beginning of a trajectory at t = 0, ECI and ECEF systems are aligned. !!!
                        % !!! for other applications, this is not suitable. !!!
      
      omega_Mars = omega_init + abs(marsRotRate)*elapsedTOF;  %quantity expressing angular distance (in rad) that the Mars rotates/moves around the z-axis.
      
      % Position Conversion - below represents ECI = ECEF_pos*[Rz(t)]- rotating about the z-axis..  
      % Rz = [cos(omega) -sin(omega);   sin(omega)  cos(omega)];
      inputStruct.X = cos(omega_Mars).*inputStruct.E - sin(omega_Mars).*inputStruct.F;  % km
      inputStruct.Y = sin(omega_Mars).*inputStruct.E + cos(omega_Mars).*inputStruct.F;  % km
      inputStruct.Z = inputStruct.G;                                                     % km
      
      % Velocity Conversion- must use calculus chain-rule due to changing
      % quantity of Mars's rotation (omega) within the DCM:
      % ECI_vel = d/dt[DCM*ECEF_position] -->  d/dt(DCM)*ECEF_position] + DCM* d/dt(ECEF_position)
      
      inputStruct.dX = cos(omega_Mars).*inputStruct.dE - sin(omega_Mars).*inputStruct.dF - ...
        abs(marsRotRate).*(sin(omega_Mars).*inputStruct.E + cos(omega_Mars).*inputStruct.F);    % km/sec
      
      inputStruct.dY = sin(omega_Mars).*inputStruct.dE + cos(omega_Mars).*inputStruct.dF + ...
        abs(marsRotRate).*(cos(omega_Mars).*inputStruct.E - sin(omega_Mars).*inputStruct.F);    % km/sec
      
      inputStruct.dZ = inputStruct.dG;                                                             % km/sec
      
      inputStruct.eci_vel_norm = sqrt(inputStruct.dX.^2 + inputStruct.dY.^2 + inputStruct.dZ.^2);  

      % Below equations directly converts ECEF acceleration to ECI acceleration. Matches with TAOS math specification. 
      
      % Vehicle ECI acceleration compute:
      % ddEi = inputStruct.ddE - 2.*abs(marsRotRate).*inputStruct.dF(1:end-1) - ((abs(marsRotRate))^2).*inputStruct.E(1:end-1);
      % ddFi = inputStruct.ddF + 2.*abs(marsRotRate).*inputStruct.dE(1:end-1) - ((abs(marsRotRate))^2).*inputStruct.F(1:end-1);
       
      % inputStruct.ddX = ddEi.*cos(omega_Mars(1:end-1)) - ddFi.*sin(omega_Mars(1:end-1));
      % inputStruct.ddY = ddEi.*sin(omega_Mars(1:end-1)) + ddFi.*cos(omega_Mars(1:end-1));
      % inputStruct.ddZ = inputStruct.ddG;      

      outputStruct = inputStruct;
end

%% NEDtoBody Conversion:
function outputStruct = NEDtoBody_Convert(inputStruct)

    % NEDtoBody- computes Specific Load Correction- we remove gravity's contribution 
     % removes gravitation acceleration quantity from overall vehicle acceleration quantity using NED frame
     % for simplicity- we just remove the downward, positive gravitation
     % acceleration of 3.72076 m/sec2. Please note- this only works for a spherical Mars model.. if using a non-spherical Mars, need to modify !!
      
      % grav_accel = earth_grav*pow(earthRadius_m/(earthRadius_m + altitude[ii]),2);
      % grav_accel_alt_corrected = Mars_GenPhysCons.GRAV_ACCEL*0.001* (Mars_GenPhysCons.RE_EQ/(Mars_GenPhysCons.RE_EQ  + altitude[ii]))^2 ;

      ddD_no_grav = inputStruct.ned_ddD - Mars_GenPhysCons.GRAV_ACCEL*0.001; %km/sec2

      inputStruct.spec_load_accel = sqrt(inputStruct.ned_ddN.^2 + inputStruct.ned_ddE.^2 + ddD_no_grav.^2);
      
      for kkk = 1:numel(inputStruct.time)-1
          accelNED_no_grav = [inputStruct.ned_ddN(kkk); inputStruct.ned_ddE(kkk); ddD_no_grav(kkk)];
               
          % Convert Euler angles to the Direction Cosine Matrix to rotate: 
          dcm_ned2body = euler2DCM(inputStruct.Vaz(kkk), inputStruct.Vqe(kkk),  0.0);  %yaw, pitch, roll... 
          %dcm_ned2body = euler2DCM(inputStruct.Vaz(kkk), inputStruct.Vqe(kkk) - inputStruct.aoa(kkk),  0.0);  %yaw, pitch, roll... 
          
          % Here we rotate the NED coordinate frame into the missile-body
          % frame. This is using a 3-2-1 sequence rotation:
          accelBody_no_grav = (dcm_ned2body)*accelNED_no_grav;
          
          inputStruct.Body_X_Accel(kkk,1)      = accelBody_no_grav(1); % km/sec2
          inputStruct.Body_LateralAccel(kkk,1) = accelBody_no_grav(2); % km/sec2
          inputStruct.Body_NormalAccel(kkk,1)  = accelBody_no_grav(3); % km/sec2
          
      end
      
      outputStruct = inputStruct;
end

function rng = great_circ_dist(lat1, lon1, lat2, lon2)

%     if nargin == 4
%     end

    r = Mars_GenPhysCons.RE_EQ; % km
    % Convert degrees to radians
    d2r = pi/180;
    lat1 = lat1*d2r;
    lon1 = lon1*d2r;
    lat2 = lat2*d2r;
    lon2 = lon2*d2r;
    
    % Calculate great circle distance between points on a sphere using the
    % Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
    % length and has the same units as the radius of the sphere, R.  (If R is
    % 1, then RNG is effectively arc length in radians.)
    
    a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;
    rng = r * 2 * atan2(sqrt(a),sqrt(1 - a));

end



function rho = mars_atm_density(altitude_ft)  %alt in km(s)

    if altitude_ft > 22960     
        temp  = -10.34 - 0.001217*altitude_ft;
        press = 14.62*exp(-0.00003*altitude_ft);
    else 
        temp  = -25.68 - 0.000548*altitude_ft;
        press = 14.62*exp(-0.00003*altitude_ft);
    end

%rho = press/(1149*(temp + 459.7));  %slugs/ft^3
quant = 1149*temp + 1149*459.7;

%rho = press/(1149*(temp + 459.7));  %slugs/ft^3
rho = press/quant;  %slugs/ft^3
end


function rho = mars_atm_density_kg_m3(altitude_km)  %alt in km(s)

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
rho = press/quant;  %slugs/ft^3

end
