function [J] = OPTIMIZE_MAIN_3DOF_ATO(x_inputs)

global target_LMO_alt 
global target_LMO_incl 
global target_LMO_e;

%% Scaled input conversion to 3DOF values
x0 = x_inputs;

%% Main Body- 3DOF Call and Execution 
try 
   run_output =  TRAJ_3DOF_TOOL_ATO_EXEC(x0);

catch
   disp(x0);
   disp('Matlab choked on the propagation'); 
   disp('check traj !!');
   pause;
end

eccen_orbit    = run_output(2);
perigee_alt_km = run_output(3);
incl_orbit     = run_output(4);

eccen_delta   = abs(eccen_orbit - target_LMO_e)           ;
perigee_delta = abs(perigee_alt_km-target_LMO_alt);
incl_delta    = abs(incl_orbit - target_LMO_incl);

%% Cost Function:
% we're here to minimize(J) !!:
J = perigee_delta/100 + eccen_delta + incl_delta/10;
% J = 0.1*perigee_delta/100 + 0.8*eccen_delta + 0.1*incl_delta/10;  

%%  Command Output:
ID1=1;
disp('====================================================================');
fprintf(ID1, 'Qe(deg): %f   ', x0(1)); fprintf(ID1, '           Azimuth (deg): %f\n', x0(2));
disp('====================================================================');
fprintf(ID1, '2nd Stage Coast time (sec): %f\n', x0(3)); 
disp('====================================================================');
fprintf(ID1, 'PitchDown Rate (deg/sec): %f\n', x0(4));
disp('====================================================================');
fprintf(ID1, 'LMO perigee alt(km): %f\n',   run_output(3));
fprintf(ID1, 'LMO apogee alt(km): %f\n',    run_output(6));
fprintf(ID1, 'LMO eccen : %f\n',            eccen_orbit);
fprintf(ID1, 'LMO incl (deg): %f\n',        incl_orbit);
disp('                                                                    ');
disp('====================================================================');
fprintf(ID1, 'LMO perigee alt delta (km): %f\n',        perigee_delta);
fprintf(ID1, 'LMO eccen delta:: %f\n',                  eccen_delta);
fprintf(ID1, 'LMO incl delta (deg) : %f\n',             incl_delta);
disp('====================================================================');
disp('                                                                    ');
%% Output Input Variables
fprintf(ID1, 'Cost Function Computed: %f\n', J);
disp('====================================================================');
disp(' ');

end
