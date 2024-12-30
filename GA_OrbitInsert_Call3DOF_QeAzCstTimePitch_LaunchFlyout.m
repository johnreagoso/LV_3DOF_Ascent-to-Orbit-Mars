function output_array = GA_OrbitInsert_Call3DOF_QeAzCstTimePitch_LaunchFlyout(x0)
global target_LMO_alt;
global target_LMO_incl;

%% Main Body- 3DOF Call and Execution 
try 
   %run_output = TrajectoryTool_3DOF_GA_Exec_MSFC_LTGSTEER(x0);
   run_output = TrajectoryTool_3DOF_GA_Exec_MSFC(x0);
   
   eccen_orbit    = run_output(2);
   perigee_alt_km = run_output(3);
   incl_orbit     = run_output(4);  
   apogee_alt_km  = run_output(6);


catch
    % time            = 1e8;
    eccen_orbit     = 10;
    perigee_alt_km  = -3000;
    apogee_alt_km  = 100;
    incl_orbit      = 89;    
    disp('===========================================catch');
end

[row,~] = find(isnan(run_output));

if isempty(row) ~=1
    eccen_orbit     = 10;
    perigee_alt_km  = -3000;
    incl_orbit      = 89;
    apogee_alt_km  = 100;
    
    run_output(2) = eccen_orbit;
    run_output(3) = perigee_alt_km;
    run_output(4) = incl_orbit;
    disp('==============================================Inf');
end

%% Cost Function:
% trying to minimize these delta(s) here:

eccen_delta   = abs(eccen_orbit - 0.00010)           ;
perigee_delta = abs(perigee_alt_km-target_LMO_alt);
apogee_delta  = abs(apogee_alt_km - target_LMO_alt);
incl_delta    = abs(incl_orbit - target_LMO_incl);

%% we're here to minimize(J) !!:
% J = 1/abs(perigee_alt_km) + eccen_delta + incl_delta/100;
J = 0.10*perigee_delta/100 + 0.8*eccen_delta + 0.1*incl_delta/10;
%J = 0.2*perigee_delta + 0.2*apogee_delta + 0.6*incl_delta;

try 
    output_array(1) = J;            
    output_array(2) = run_output(2);
    output_array(3) = run_output(3);          
    output_array(4) = run_output(4);
    output_array(5) = run_output(1);    
    output_array(6) = run_output(6);
    output_array(7) = run_output(7);
catch
    disp('here on line-60');
end


end