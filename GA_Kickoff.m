%% Trajectory Shaping Genetic Algorithm 
clear all;
close all;
clc

for ii = 1:1:3
    for jj = 1:1:3
    
        for kk = 1:20
            try
                % OrbitInsertGA_TrajShape_LaunchFlyoutVer_VehConfig_v2(jj, ii);
                % GA_ATO_EXEC_LaunchFlyout(jj,ii);
                GA_ATO_EXEC_LaunchFlyoutwPtchDwn(jj, ii);
            catch
                disp('failure for attempt no. :'); 
                disp(kk);
            end
        end

    end
end