%% here we call the Matlab Fmincon tool:
fclose all;
clc;
clear all; %#ok<CLALL>

global payload
global target_LMO_alt target_LMO_incl target_LMO_e;

%% Inputs:
payload = 16;

%% Mars Target:
target_LMO_alt  = 343.0; % km 
target_LMO_incl = 25.0;  % deg
target_LMO_e    = 0.00; 

%% Main Body:
lb = [70.00,            65.0,           400,            0.0]; 

x0 = [87.999999 ,       68.556182,      460.915269,     0.104070]; 

ub = [88.0,             70.0,           550,            0.50];

%% Options:
options = optimoptions('fmincon','Display','iter','Algorithm','SQP'); %, 'DiffMaxChange', 1.0e6, 'DiffMinChange', 1e-7);

%options = optimoptions('fmincon','Display','iter','Algorithm','active-set'); %,'DiffMaxChange', 1.0e-1, 'DiffMinChange', 1e-6); %
%options = optimoptions(options, 'TolFun', 1e-9,'Display', 'iter', 'TolX', 1e-13);  

options = optimoptions("fmincon", "Algorithm","interior-point", "EnableFeasibilityMode",true,...
                          "SubproblemAlgorithm","cg"); %,'DiffMaxChange', 1e2, 'DiffMinChange', 1e-12);

%% Crucial Settings for potentially large solution-space problems:
%boptions = optimset('DiffMaxChange', 1e6, 'DiffMinChange', 1e-6);   % Default: Inf.. 0.01 earlier
% options = optimset('DiffMinChange', 1e-6);  % Default: 0.0
% options = optimset('DiffMinChange', 0.0);  % Default: 0.0

%% Main Body Functional Call:
%[x, fval, exitflag] = fmincon(@OPTIMIZE_MAIN_3DOF_ATO, ...
%  x0, [], [], [],[], lb, ub, [], options);

[x, fval, exitflag] = fmincon(@OPTIMIZE_MAIN_3DOF_ATO, ...
    x0, [], [], [],[], lb, ub, @Optimize_Main_Reagoso3DOF_QeAz_NonLinConstr); %, options);

%% Output
fprintf('%f\n',x(1));
fprintf('%f\n',x(2))
fprintf('%f\n',x(3))
fprintf('%f\n',fval);

disp(x);
disp(fval);
