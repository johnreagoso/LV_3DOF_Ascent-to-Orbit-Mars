function [output_array] = TRAJ_3DOF_TOOL_ATO_EXEC(x0) 

% x0(1) = 84.500829;
% x0(2) = 71.721288;
% x0(3) = 465.756155;
% x0(4) = 0.0;

%% Vehicle Config:
m_shroud =  MAV_LaunchVehicle_Config.m_shroud;

% Stage-1 STAR-20G  
% Stage-2 STAR-15G 

% 450kg GLOM:
stagemass = [MAV_LaunchVehicle_Config.stagemass(1)   MAV_LaunchVehicle_Config.stagemass(2)];  % Mass total of each stage [kg]
m_prop    = [MAV_LaunchVehicle_Config.m_prop(1) MAV_LaunchVehicle_Config.m_prop(2)];            % Propellant mass of each stage [kg]

%% PL and mass total:
payload = MAV_LaunchVehicle_Config.payload;
m_init  = sum(stagemass) + m_shroud + payload;  %kg

%% Initial Cond(s)
vehicleObjInput.time         = 0.0;  
vehicleObjInput.speed        = 0.0050;  % ~5-10 m/sec MAV gas ejection from ESA SRL  
vehicleObjInput.latitude     = 18.38;
vehicleObjInput.longitude    = 77.58;
vehicleObjInput.altitude     = 0.0045;  %4.5 meter (~15ft) ignition post-eject from SRL
vehicleObjInput.Vqe          = x0(1); 
vehicleObjInput.rotatingMars = 1;
vehicleObjInput.Vaz          = x0(2);

vehicleObj = CoordFrameConvert_Vehicle_Mars(vehicleObjInput, 'RL');

%% New Inputs: ECEF (Earth-Centered-Earth-Fixed)
Rvector = [vehicleObj.E;  vehicleObj.F; vehicleObj.G];
Vvecfor = [vehicleObj.dE; vehicleObj.dF; vehicleObj.dG];

%% Execute Runge-Kutta for trajectory propagation
coast2 = x0(3) + 1.0e-6;

x01  = [Rvector(1) Rvector(2) Rvector(3)...
       Vvecfor(1) Vvecfor(2) Vvecfor(3) m_init];

options_coast = odeset('events', 'myEvent3DOF_LMO', 'RelTol',1e-6,'InitialStep', 0.0001,'MaxStep', 1.000);
options       = odeset('events', 'myEvent3DOF_LMO', 'RelTol',1e-6, 'InitialStep',0.0001,'MaxStep', 1.000);

%% Vehicle Propagation:
prop = 1;
while prop ~= 0
    
    %tic
    [epoch1, state_vector1] = ode45(@Function_3DOF_MarsJ2_wThrust_STAR20, [0.0, MAV_LaunchVehicle_Config.stage1_burntime], x01, options, x0(4));
    %toc
        [row1,~] = find(isnan(state_vector1));
        if isempty(row1) ~= 1
            state_vector4 = state_vector1(1:row1(1)-1,:); epoch4 = epoch1(1:row1(1)-1);
            break;
        end

    state_vector1(end,7) = state_vector1(end,7) - (stagemass(1) - m_prop(1)) - m_shroud;
    
    [epoch2, state_vector2] = ode45(@Function_3DOF_MarsJ2, [epoch1(end), epoch1(end)+coast2], state_vector1(end,:), options_coast);
    
        [row2,~] = find(isnan(state_vector2));
        if isempty(row2) ~= 1
            state_vector4 = state_vector2(1:row2(1)-1,:); epoch4 = epoch2(1:row2(1)-1);
            break;
        end
    
    [epoch3, state_vector3] = ode45(@Function_3DOF_MarsJ2_wThrust_STAR15G, [epoch2(end), epoch2(end)+MAV_LaunchVehicle_Config.stage2_burntime],...
        state_vector2(end,:), options); %, [mdot2, ISP2, epoch2(end)]);
      
        [row3,~] = find(isnan(state_vector3));
        if isempty(row3) ~= 1
            state_vector4 = state_vector3(1:row3(1)-1,:); epoch4 = epoch3(1:row3(1)-1);
            break;
        end

    [epoch4, state_vector4] = ode45(@Function_3DOF_MarsJ2, [epoch3(end), epoch3(end) + 1.0], state_vector3(end,:), options_coast);
        
        [row4,~] = find(isnan(state_vector4));
        if isempty(row4) ~= 1
            state_vector4 = state_vector4(1:row4(1)-1,:); epoch4 = epoch4(1:row4(1)-1);
            break;
        end

    prop = 0;
end

vehicleObj.time = epoch4(end);
vehicleObj.E = state_vector4(end,1);  vehicleObj.F = state_vector4(end,2);  vehicleObj.G = state_vector4(end,3);
vehicleObj.dE = state_vector4(end,4); vehicleObj.dF = state_vector4(end,5); vehicleObj.dG = state_vector4(end,6);
 
epochCombine = cat(1, epoch1, epoch2(2:end),epoch3(2:end),epoch4(2:end));
massCombine  = cat(1, state_vector1(:,7), state_vector2(2:end,7), state_vector3(2:end,7), state_vector4(2:end,7));
 
ECombine = cat(1, state_vector1(:,1), state_vector2(2:end,1), state_vector3(2:end,1), state_vector4(2:end,1));
FCombine = cat(1, state_vector1(:,2), state_vector2(2:end,2), state_vector3(2:end,2), state_vector4(2:end,2));
GCombine = cat(1, state_vector1(:,3), state_vector2(2:end,3), state_vector3(2:end,3), state_vector4(2:end,3));

dECombine = cat(1, state_vector1(:,4), state_vector2(2:end,4), state_vector3(2:end,4), state_vector4(2:end,4));
dFCombine = cat(1, state_vector1(:,5), state_vector2(2:end,5), state_vector3(2:end,5), state_vector4(2:end,5));
dGCombine = cat(1, state_vector1(:,6), state_vector2(2:end,6), state_vector3(2:end,6), state_vector4(2:end,6));

%% backup:

vehicleObjOut_single.time = epochCombine(end);
vehicleObjOut_single.E  = ECombine(end);     vehicleObjOut_single.F  = FCombine(end);     vehicleObjOut_single.G  = GCombine(end);
vehicleObjOut_single.dE = dECombine(end);    vehicleObjOut_single.dF = dFCombine(end);    vehicleObjOut_single.dG = dGCombine(end);
vehicleObjOut_single.rotatingMars = 1;

vehicleObj_endStep = ECEFtoECI_Convert_Mars_StandAlone(vehicleObjOut_single);

vehicleObj_endStep   = state_elements_Mars_Struct(vehicleObj_endStep); 
 
%output_array = [vehicleObj_endStep.time, vehicleObj_endStep.e, vehicleObj_endStep.perigee_alt, vehicleObj_endStep.incl, vehicleObj_endStep.RAAN, vehicleObj_endStep.apogee_alt, vehicleObj_endStep.SMA];

%% primary: comment out below here for GA runs.. 
% vehicleObjOut.time = epochCombine;
% vehicleObjOut.E = ECombine;     vehicleObjOut.F = FCombine;     vehicleObjOut.G = GCombine;
% vehicleObjOut.dE = dECombine;   vehicleObjOut.dF = dFCombine;   vehicleObjOut.dG = dGCombine;
%  
% vehicleObjOut.rotatingMars = 1;
% 
% vehicleObj = CoordFrameConvert_Vehicle_Mars(vehicleObjOut, 'ecef');
%  
% % PlotCall_Mars(vehicleObj, 'yes');
% % PlotCall_Standard_Mars(vehicleObj, 'yes');
% % vehicleObj.time = epoch4(end);
% % 
% vehicleLast.time = epoch4(end);
% vehicleLast.X = vehicleObj.X(end);  vehicleLast.Y = vehicleObj.Y(end);  vehicleLast.Z = vehicleObj.Z(end);
% vehicleLast.dX = vehicleObj.dX(end); vehicleLast.dY = vehicleObj.dY(end); vehicleLast.dZ = vehicleObj.dZ(end);
% %output_array = [vehicleObj_endStep.time, vehicleObj_endStep.e, vehicleObj_endStep.perigee_alt, vehicleObj_endStep.incl, vehicleObj_endStep.RAAN, vehicleObj.apogee_alt, vehicleObj.SMA];    

output_array = [vehicleObj_endStep.time, vehicleObj_endStep.e, vehicleObj_endStep.perigee_alt, vehicleObj_endStep.incl, vehicleObj_endStep.RAAN, vehicleObj_endStep.apogee_alt, vehicleObj_endStep.SMA];

return


