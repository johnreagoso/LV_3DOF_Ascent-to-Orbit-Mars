classdef MAV_LaunchVehicle_Traj

    properties

        launchAngle           % Angle (wrt local horz) at kickover [deg]
        launchLatitude        % Launch site latitude [deg]
        launchLongitude       % Launch site longitude [deg]
        launchAltitude        % Launch site altitude [km]
        launchAzimuth         % Azimuth (wrt true north) at kickover [deg]
        launchSpeed = 0;      % Initial speed at ignition [m/s]
        coast2 
    end

    properties (Constant)

        Cd = 0.65;
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

        %m_init = sum(stagemass) + m_shroud + payload;  %kg
        
        %% Propulsion Metrics:
        stage1_burntime = 75.0 % sec
        stage2_burntime = 20.0 % sec

        ISP1 = 288.5;   % sec
        ISP2 = 289.00;  % sec


    end


    methods (Static)

        function dum = mdot1()
             dum = MAV_LaunchVehicle_Config.m_prop(1)/MAV_LaunchVehicle_Config.stage1_burntime;  % lbm --> kg
        end

        function dum = mdot2()
             dum = MAV_LaunchVehicle_Config.m_prop(2)/MAV_LaunchVehicle_Config.stage2_burntime;  % lbm --> kg
        end

    end

end