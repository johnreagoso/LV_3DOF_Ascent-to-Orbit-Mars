classdef MAV_LaunchVehicle_Config

%     properties
%         time
%         launchAngle                     % Angle (wrt local horz) at kickover [deg]
%         launchLatitude  = 18.38;        % Launch site latitude [deg]
%         launchLongitude = 77.58;        % Launch site longitude [deg]
%         launchAltitude  = 0.0045;       % Launch site altitude [km]
%         launchAzimuth                   % Azimuth (wrt true north) at kickover [deg]
%         launchSpeed     = 0.0050;       % Initial speed at ignition [m/s]
%         coast2 
%         rotatingMars = 1;
%         Vqe
%         Vaz
%         aoa
%         range_km
%         range_nm 
%         ned_dN
%         ned_dE
%         ned_dD
%         ned_dU
%         ned_ddN
%         ned_ddE
%         ned_ddD
%         ned_ddU
%         ned_accel_norm
%         E
%         F
%         G
%         dE
%         dF
%         dG
%         ddE
%         ddF
%         ddG
%         altitude
%         latitude
%         longitude
%         speed
%         X
%         Y
%         Z
%         dX
%         dY
%         dZ
%         ddX
%         ddY
%         ddZ
%         eci_vel_norm
%         spec_load_accel
%         altitude_ft
%         speed_ftsec
%         Q_psf
%         MaxQ
%     end

    properties (Constant)

        Cd = 0.6;
        vArea = 0.196648; %m^2  0.25 meter radius/ 0.50 meter diameter

        %Dia = [0.5  0.5];      %  Outside diameter of each stage [m] Area = pi*Dia.^2/4;                   %  Cross-sectional area [m^2]
        %Area = pi*Dia.^2/4;
        
        %% Payload Metrics
        payload  = 16;
        m_shroud = 0.0
                
        %% Stage Mass Metrics
        % Stage-1 STAR-20G   Stage-2 STAR-15G 

        % 450kg GLOM:
        stagemass = [373.2   60.8];         % Mass total of each stage [kg]
        
        % 425kg GLOM:
        % stagemass = [339.1   69.9];       % Mass total of each stage [kg]

        % 400kg GLOM:
        % stagemass = [303.6   80.4];       % Mass total of each stage [kg]
        
        m_prop    = [213    47];            % Propellant mass of each stage [kg]
        % MFstage = m_prop./stagemass;      %  Mass fraction of each stage [0]  

        
        
        %% Propulsion Metrics:
        stage1_burntime = 75.0 % sec
        stage2_burntime = 20.0 % sec

        ISP1 = 288.5;   % sec
        ISP2 = 289.00;  % sec

    end

    methods (Static)

        function dum = mdot1()
             dum = MAV_LaunchVehicle_Config.m_prop(1)/MAV_LaunchVehicle_Config.stage1_burntime;  
        end

        function dum = mdot2()
             dum = MAV_LaunchVehicle_Config.m_prop(2)/MAV_LaunchVehicle_Config.stage2_burntime;  
        end

        function dum = m_init()
            dum = sum(MAV_LaunchVehicle_Config.stagemass) + MAV_LaunchVehicle_Config.m_shroud +  MAV_LaunchVehicle_Config.payload;  %kg
        end

    end


%     methods
% 
%         function LV_obj = propTraj(obj)
% 
% 
% 
%         end
% 
%     end


end