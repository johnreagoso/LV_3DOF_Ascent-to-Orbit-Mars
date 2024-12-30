function output_array = GA_ATO_Call3DOF(x0)

global target_LMO_alt;
global target_LMO_incl;
global target_LMO_e

%% Main Body- 3DOF Call and Execution 
try 
   run_output = TRAJ_3DOF_TOOL_ATO_EXEC(x0);
   
   eccen_orbit    = run_output(2);
   perigee_alt_km = run_output(3);
   incl_orbit     = run_output(4);  

catch
    eccen_orbit     = 10;
    perigee_alt_km  = -3000;
    incl_orbit      = 89;    
    disp('===========================================catch');
end

[row,~] = find(isnan(run_output));

if isempty(row) ~=1
    eccen_orbit     = 10;
    perigee_alt_km  = -3000;
    incl_orbit      = 89;

    run_output(2) = eccen_orbit;
    run_output(3) = perigee_alt_km;
    run_output(4) = incl_orbit;
    disp('==============================================Inf');
end

%% Cost Function:
% trying to minimize these delta(s) here:
eccen_delta   = abs(eccen_orbit - target_LMO_e)           ;
perigee_delta = abs(perigee_alt_km-target_LMO_alt);
incl_delta    = abs(incl_orbit - target_LMO_incl);

%% we're here to minimize(J) !!:
J = 0.10*perigee_delta/100 + 0.8*eccen_delta + 0.1*incl_delta/10;

output_array(1) = J;            
output_array(2) = run_output(2);
output_array(3) = run_output(3);          
output_array(4) = run_output(4);
output_array(5) = run_output(1);    
output_array(6) = run_output(6);
output_array(7) = run_output(7);

end