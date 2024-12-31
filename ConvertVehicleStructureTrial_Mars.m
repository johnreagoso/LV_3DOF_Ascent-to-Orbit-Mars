%% ConvertStructure
% Reference for Mars physical/gravity parameters:  https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html 

% Note: this conversion script origin is for Earth centered trajectories (ECEF, ECI etc.). Any use of 'ECI', 'ECEF' etc. is a holdover 
% and refers to MCI (Mars Centered Inertial) or MCMF (Mars Centered Mars Fixed) coordinate frames. Future versions of this script will
% be modified accordingly. 

function [vehicleObj]= ConvertVehicleStructureTrial_Mars(vehicleObj, frame, units)
    
    if strcmp(frame, 'ecef') == 1 || strcmp(frame, 'ECEF') == 1 || strcmp(frame, 'edm') == 1
        
%         if strncmp(units, 'km', 2) == 1 || strncmp(units, 'KM', 2) == 1
%             vehicleObj.E = vehicleObj.E*1.0e3;          vehicleObj.F = vehicleObj.F*1.0e3;          vehicleObj.G = vehicleObj.G*1.0e3; 
%             vehicleObj.dE = vehicleObj.dE*1.0e3;        vehicleObj.dF = vehicleObj.dF*1.0e3;        vehicleObj.dG = vehicleObj.dG*1.0e3;   
%         end

         Rx = vehicleObj.E;  Ry = vehicleObj.F;  Rz = vehicleObj.G;
         Vx = vehicleObj.dE; Vy = vehicleObj.dF; Vz = vehicleObj.dG;
    
        %% Geocentric Latitude:
              inv_f = 169.779286926995;
              % a     = Mars_GenPhysCons.RE_EQ*1000;
              % b     = a*(1-1/inv_f);
              % ecc   = sqrt(1-(b/a)^2); 

%         marsSpheroid.Name = 'MarsGeodeticApprx';
%         marsSpheroid.LengthUnit = 'meter';
%         marsSpheroid.SemimajorAxis     = wgs84Ellipsoid('meters').SemimajorAxis*0.532;
%         marsSpheroid.SemiminorAxis     = wgs84Ellipsoid('meters').SemiminorAxis*0.531;
%         marsSpheroid.InverseFlattening = wgs84Ellipsoid('meters').InverseFlattening*1.76;
%         marsSpheroid.Eccentricity      = 0.1105; %wgs84Ellipsoid('meters').Eccentricity*

        [lat_rad, long_rad, ~] = ecef2geodetic(vehicleObj.E, vehicleObj.F, vehicleObj.G, referenceEllipsoid('Mars'));

%         r = sqrt(vehicleObj.E.^2 + vehicleObj.F.^2);
%         ab2 = a^2-b^2;
%         ar1 = 1./(a*r);
% 
%         borkE = b*abs(vehicleObj.G) - ab2*ar1;
%         borkF = b*abs(vehicleObj.G) + ab2*ar1;
% 
%         borkE = b*abs(vehicleObj.G) - a^2-b^2/a*r;
%         borkF = b*abs(vehicleObj.G) + (a^2-b^2))/(a*r);
%         
%         P = (4/3)*(borkE.*borkF + 1);
%         Q = 2*(borkE.^2 - borkF.^2);
%         D = P.^3 + Q.^2;
%         
%         v = -(sqrt(D)+Q).^(1/3) + (sqrt(D)-Q).^(1/3);
%         borkG = (sqrt(borkE.^2+v)+borkE)/2;
%         
%         t = (sqrt(borkG.^2 + (borkF-v.*borkG))./(2*borkG-borkE)) - borkG;
%         lat_rad = sign(vehicleObj.G).*atan2(a*(1-t.^2),(2*b*t));       
        lat_deg = lat_rad*180/pi;
        
        %% Longitude:
%        long_rad = 2*atan2(sqrt(vehicleObj.E^2 + vehicleObj.F^2)- vehicleObj.E, vehicleObj.F);
        long_deg = long_rad*180/pi;
        
        HighVals = find(long_deg>180);
        long_deg(HighVals) = long_deg(HighVals)-360;
        
        for ii = 1:1:numel(vehicleObj.E)
            DCMecef2ned{ii} = [cos(lat_rad(ii))*cos(long_rad(ii))     cos(lat_rad(ii))*sin(long_rad(ii))      sin(lat_rad(ii));...
                        -sin(long_rad(ii))      cos(long_rad(ii))       0;...
                            -sin(lat_rad(ii))*cos(long_rad(ii))     -sin(lat_rad(ii))*sin(long_rad(ii))     cos(lat_rad(ii));];
                                
            VelNED{ii} = DCMecef2ned{ii}*[Vx(ii); Vy(ii); Vz(ii)];   %ecefVel = transpose(TewariDCM)*[Vn; Ve; Vd];
            PosNED{ii} = DCMecef2ned{ii}*[Rx(ii); Ry(ii); Rz(ii)];
        end
    
        f = 1/169.779286926995; 
        RadiusEq = Mars_GenPhysCons.RE_EQ; %*1000;
        lamda = atan((1-f)^2. *tan(lat_rad));
        Rpos = sqrt(RadiusEq^2./(1 + (1/(1-f)^2 - 1)*sin(lamda).^2));

        for ii = 1:numel(vehicleObj.E)
            alt(ii,1) = norm([Rx(ii); Ry(ii); Rz(ii)]) - Rpos(ii);

            Velmag(ii,1) = norm([Vx(ii); Vy(ii); Vz(ii)]);
                
            VelNEDsingle = [VelNED{ii}(1);  VelNED{ii}(2);  VelNED{ii}(3)];
            PosNEDsingle = [PosNED{ii}(1);  PosNED{ii}(2);  PosNED{ii}(3)];

            vmag(ii,1) = norm([VelNED{ii}(1); VelNED{ii}(2); VelNED{ii}(3)]);
            
            Vaz(ii,1) = atan2d(sqrt(VelNED{ii}(1)^2 + VelNED{ii}(2)^2),VelNED{ii}(3)); 
            Vqe(ii,1) = 90 - atan2d(VelNED{ii}(2),VelNED{ii}(1));                   
    
        end

            %% Output
        vehicleObj.speed = Velmag;
        vehicleObj.latitude = lat_deg;
        vehicleObj.longitude = long_deg;
        vehicleObj.Vqe = Vqe;
        vehicleObj.Vaz = Vaz;
        vehicleObj.altitude = alt;

      %  vehicleObj = ECEFtoECI_Convert_Mars_StandAlone(vehicleObj);

    elseif strcmp(frame, 'rl')
              
%       lat_rad  = vehicleObj.latitude.*pi/180;
%       long_rad = vehicleObj.long.*pi/180;
% 
%       inv_f = 169.779286926995;
%       a     = Mars_GenPhysCons.RE_EQ;
%       b     = a*(1-1/inv_f);     
%       ecc   = sqrt(1-(b/a)^2); 
%       
%       r = a/sqrt(1-ecc^2*sind(vehicleObj.latitude)^2);
%       vehicleObj.E = (r + vehicleObj.alt)*cosd(vehicleObj.latitude)*cosd(vehicleObj.long);
%       vehicleObj.F = (r + vehicleObj.alt)*cosd(vehicleObj.latitude)*sind(vehicleObj.long);
%       vehicleObj.G = ((1 - ecc^2)*r + vehicleObj.alt)*sind(vehicleObj.latitude);
%       
%       %Matlab Toolbox Geodetic2ECEF Conversion. Works too!! .. WHY???!!!
%        mars_spheroid.Name = 'MarsGeodeticApprx';
%        mars_spheroid.LengthUnit = 'meter';
%        mars_spheroid.SemimajorAxis = ;
%        mars_spheroid.SemiminorAxis = ;
%        mars_spheroid.InverseFlattening = 
%        
%        [vehicleObj.E, vehicleObj.F, vehicleObj.G] = geodetic2ecef(lat_rad, long_rad, vehicleObj.alt, spheroid);
% 
%       % ECEF to NED Conversion matrix derived from TAOS manual. Also matches with Matlab technical page: 
%       % https://www.mathworks.com/help/aeroblks/directioncosinematrixeceftoned.html
%       ECEF2NED_DCM = [-sin(lat_rad)*cos(long_rad)   -sin(lat_rad)*sin(long_rad)  cos(lat_rad);...
%                     -sin(long_rad)   cos(long_rad)  0;...
%                         -cos(lat_rad)*cos(long_rad)  -cos(lat_rad)*sin(long_rad)  -sin(lat_rad)];  
%                     
%       [dN, dE, dD] = sph2cart(vehicleObj.Vaz, vehicleObj.Vqe, vehicleObj.vmag); VelNED = [dN; dE; dD];             
%                     
%       ecefVel = ECEF2NED_DCM^-1*VelNED;

      %% Output:
      vehicleObj.dE = ecefVel(1);
      vehicleObj.dF = ecefVel(2);
      vehicleObj.dG = ecefVel(3);
    
    end    
end

