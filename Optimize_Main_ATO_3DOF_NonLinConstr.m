 function [c,ceq] = Optimize_Main_ATO_3DOF_NonLinConstr(x_inputs)

global target_LMO_alt target_LMO_incl target_LMO_e;

%% Scaled input conversion to 3DOF values
x0 = x_inputs;

%% Main Body- 3DOF Call and Execution 
try 
   run_output =  TrajectoryTool_3DOF_GA_Exec_MSFC_LTGSTEER(x0);
catch
   disp(x0);
   disp('Matlab choked on the propagation- Line 21 in NonLinConstraints'); 
   disp('check traj !!');
   pause;
end

eccen_orbit    = run_output(2);
perigee_alt_km = run_output(3);
incl_orbit     = run_output(4);

%% Constraint Function(s):

c = [];
ceq = [ perigee_alt_km-target_LMO_alt, eccen_orbit - target_LMO_e, incl_orbit - target_LMO_incl];

end
