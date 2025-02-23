%% ConvertStructure
% Reference for Mars physical/gravity parameters:  https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html 

% Note: this conversion script origin is for Earth centered trajectories (ECEF, ECI etc.). Any use of 'ECI', 'ECEF' etc. is a holdover 
% and refers to MCI (Mars Centered Inertial) or MCMF (Mars Centered Mars Fixed) coordinate frames. Future versions of this script will
% be modified accordingly. 

function [vehicleObj]= ConvertVehicleStructureTrial_Mars(vehicleObj, frame, units)
    
    if strcmp(frame, 'ecef') == 1 || strcmp(frame, 'ECEF') == 1 || strcmp(frame, 'edm') == 1
        

         Rx = vehicleObj.E;  Ry = vehicleObj.F;  Rz = vehicleObj.G;
         Vx = vehicleObj.dE; Vy = vehicleObj.dF; Vz = vehicleObj.dG;
    
        %% Geocentric Latitude:
              inv_f = 169.779286926995;
              % a     = Mars_GenPhysCons.RE_EQ*1000;
              % b     = a*(1-1/inv_f);
              % ecc   = sqrt(1-(b/a)^2); 


        [lat_rad, long_rad, ~] = ecef2geodetic(vehicleObj.E, vehicleObj.F, vehicleObj.G, referenceEllipsoid('Mars'));
    
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
             

      %% Output:
      vehicleObj.dE = ecefVel(1);
      vehicleObj.dF = ecefVel(2);
      vehicleObj.dG = ecefVel(3);
    
    end    
end

