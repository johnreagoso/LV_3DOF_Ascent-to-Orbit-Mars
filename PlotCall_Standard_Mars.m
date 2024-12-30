

function outputStatus = PlotCall_Standard_Mars(vehicleObj, closeInReq)

%% Plotting Utilities

srm1bbo_epoch_idx = find(vehicleObj.time <= 75.0, 1, 'last');

srm2_startepoch_idx = find(vehicleObj.time <= 75+465.756155, 1, 'last');  %728
srm2_bbo_epoch_idx   = find(vehicleObj.time >= 75+465.756155+20, 1, 'first'); %833

figure;
%subplot(1, 2, 1)
plot(vehicleObj.range_km,    vehicleObj.altitude, 'linewidth', 1.5); grid on; title('MAV Baseline Trajectory Altitude vs. Range');
xlabel('km'); ylabel('km'); 
set(0,'defaultAxesFontSize',11);    ax = gca;    ax.FontSize = 11;    ax.FontWeight = 'bold';
%axis equal; xlim([0,250]);

line([0 0],         [0 400], 'Color', 'b', 'LineStyle', '--'); 
line([vehicleObj.range_km(srm1bbo_epoch_idx) vehicleObj.range_km(srm1bbo_epoch_idx)],     [0 400], 'Color', 'b', 'LineStyle', '--'); 
line([vehicleObj.range_km(srm2_startepoch_idx+1) vehicleObj.range_km(srm2_startepoch_idx+1)],     [0 400], 'Color', 'b', 'LineStyle', '--'); 
line([vehicleObj.range_km(srm2_bbo_epoch_idx) vehicleObj.range_km(srm2_bbo_epoch_idx)],     [0 400], 'Color', 'b', 'LineStyle', '--'); 

text(20, 200, 'SRM-1 Burnout', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'b');
text(170, 250, 'SRM-2 Ignition', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'b');    
text(205, 280, 'SRM-2 Burnout', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'b');    

%subplot(1, 2, 3)
figure;
plot(vehicleObj.time,  vehicleObj.Q_psf, 'linewidth', 1.5); grid on; title('MAV Baseline Trajectory Dynamic Press (Q) vs. TOF');
xlabel('sec'); ylabel('psf'); ylim([-0.1,25]); %axis square; 
set(0,'defaultAxesFontSize',11);    ax = gca;    ax.FontSize = 11;    ax.FontWeight = 'bold';

line([0 0],         [0 400], 'Color', 'b', 'LineStyle', '--'); 
line([vehicleObj.time(srm1bbo_epoch_idx) vehicleObj.time(srm1bbo_epoch_idx)],     [0 400], 'Color', 'b', 'LineStyle', '--'); 
line([vehicleObj.time(srm2_startepoch_idx+1) vehicleObj.time(srm2_startepoch_idx+1)],     [0 400], 'Color', 'b', 'LineStyle', '--'); 
line([vehicleObj.time(srm2_bbo_epoch_idx) vehicleObj.time(srm2_bbo_epoch_idx)],     [0 400], 'Color', 'b', 'LineStyle', '--'); 
          
text(75, 15, 'SRM-1 Burnout', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'b');
text(425, 10, 'SRM-2 Ignition', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'b');    
text(560, 15, 'SRM-2 Burnout', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'b'); 










       PlottingUtilitySSAG(gcf, 'MAV_Baseline_QeVsTime');   

    subplot(2, 4, 1)
        plot(vehicleObj.time,    vehicleObj.altitude, 'linewidth', 1.5); grid on; title('Altitude vs. TOF');
        xlabel('sec'); ylabel('km');
        set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
        %xlim([0, 0.01]);
        %ylim([0, 0.50]);














    subplot(2, 4, 1)
        plot(vehicleObj.time,    vehicleObj.altitude, 'linewidth', 1.5); grid on; title('Altitude vs. TOF');
        xlabel('sec'); ylabel('km');
        set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
        %xlim([0, 0.01]);
        %ylim([0, 0.50]);
        
     subplot(2, 4, 2)
        plot(vehicleObj.range_nm,    vehicleObj.altitude, 'linewidth', 1.5); grid on; title('Altitude vs. Range');
        grid on; 
        xlabel('km'); ylabel('nm');
        set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';   
       
        
    subplot(2, 4, 3)
       plot(vehicleObj.time(1:end-1),  vehicleObj.spec_load_accel/Mars_GenPhysCons.GRAV_ACCEL, 'linewidth', 1.5); grid on; title('g(s) vs. TOF');
       xlabel('sec'); ylabel('g(s)');
       set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
       ylim([0, 20]);
       
    subplot(2, 4, 4)
        plot(vehicleObj.time(1:end),    vehicleObj.speed, 'linewidth', 1.5); grid on; title('Inertial Velocity vs. TOF');
        xlabel('sec'); ylabel('km/s');  
        ylim([0 8]);
        set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';

    %% Extra Plots:               
    subplot(2, 4, 5)
       plot(vehicleObj.time,  vehicleObj.Vqe, 'linewidth', 1.5); grid on; title('FPA vs. TOF');
       xlabel('sec'); ylabel('deg'); %axis square;
       set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
   
    %PlottingUtilitySSAG(gcf, 'LV_3DOFoutConvertedData');
       
      
%% Close-In Plots- 'yes or no'. If yes, proceed... 
  
  if strcmp(closeInReq, 'yes') == 1
  
        figure;
        subplot(4, 3, 1)
            plot(vehicleObj.time,    vehicleObj.altitude, 'linewidth', 1.5); grid on; title('Altitude vs. TOF');
            xlabel('sec'); ylabel('km');
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
            xlim([0, 65.2]);

        subplot(4, 3, 2)
           plot(vehicleObj.time(1:end-1),  vehicleObj.spec_load_accel/Mars_GenPhysCons.GRAV_ACCEL, 'linewidth', 1.5); grid on; title('g(s) vs. TOF');
           xlabel('sec'); ylabel('g(s)');
           set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
           xlim([0, 65.2]);

        subplot(4, 3, 3)
            plot(vehicleObj.time(1:end),    vehicleObj.speed, 'linewidth', 1.5); grid on; title('Inertial Velocity vs. TOF');
            xlabel('sec'); ylabel('km/s');  
            ylim([0 8]);
            xlim([0, 65.2]);
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';

        subplot(4, 3, 4)
            plot(vehicleObj.time(1:end-1),  round(vehicleObj.ned_ddN,5), 'linewidth', 1.5); grid on; title('NED North Acceleration Comp vs. TOF');
            xlabel('sec'); ylabel('km/s^2');
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
            xlim([0, 65.2]);

        subplot(4, 3, 5)
            plot(vehicleObj.time(1:end-1),  round(vehicleObj.ned_ddE,5),'linewidth', 1.5); grid on; title('NED East Acceleration Comp vs. TOF');
            xlabel('sec'); ylabel('km/s^2');
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
            xlim([0, 65.2]);

        subplot(4, 3, 6)
            plot(vehicleObj.time(1:end-1),  round(vehicleObj.ned_ddD,5), 'linewidth', 1.5); grid on; title('NED Down Acceleration Comp vs. TOF');
            xlabel('sec'); ylabel('km/s^2');
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
            xlim([0, 65.2]);

        subplot(4, 3, 7)
            plot(vehicleObj.time(1:end-1),  round(vehicleObj.Body_X_Accel/Mars_GenPhysCons.GRAV_ACCEL,5), 'linewidth', 1.5);
            grid on; title('Body Longitudinal Accel Comp vs. TOF');
            xlabel('sec'); ylabel('g(s)');
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
            xlim([0, 65.2]);

        subplot(4, 3, 8)
            plot(vehicleObj.time(1:end-1),  round(vehicleObj.Body_LateralAccel/Mars_GenPhysCons.GRAV_ACCEL,5), 'linewidth', 1.5); 
            grid on; title('Body Lateral Accel Comp vs. TOF');
            xlabel('sec'); ylabel('g(s)');
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
            xlim([0, 65.2]);

        subplot(4, 3, 9)
            plot(vehicleObj.time(1:end-1),  round(vehicleObj.Body_NormalAccel/Mars_GenPhysCons.GRAV_ACCEL,5),'linewidth', 1.5); 
            grid on; title('Body Normal Accel Comp vs. TOF');
            xlabel('sec'); ylabel('g(s)'); 
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
            xlim([0, 65.2]);
            ylim([-1, 1]);

        %% Extra Plots:               
        subplot(4, 3, 10)
            plot(vehicleObj.range_km,    vehicleObj.altitude, 'linewidth', 1.5); grid on; title('Altitude vs. Range');
            grid on; 
            xlabel('km'); ylabel('km');
            xlim([0, 65.2]);
            set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';

        subplot(4, 3, 11)
           plot(vehicleObj.time,  vehicleObj.Vqe, 'linewidth', 1.5); grid on; title('FPA vs. TOF');
           xlabel('sec'); ylabel('deg'); %axis square;
           set(0,'defaultAxesFontSize',10);    ax = gca;    ax.FontSize = 10;    ax.FontWeight = 'bold';
           xlim([0, 65.2]);

       PlottingUtilitySSAG(gcf, 'LV_3DOFoutConvertedData_CloseIn');

       %% Separate Close-In Plots: 
       figure;
       plot(vehicleObj.time,  vehicleObj.Vqe, 'linewidth', 2.0); grid on; title('Flight Path Angle vs. TOF (FS Boost)');
       xlabel('sec'); ylabel('deg'); %axis square;
       set(0,'defaultAxesFontSize',16);    ax = gca;    ax.FontSize = 16;    ax.FontWeight = 'bold';
       xlim([0, 65.2]);
       ylim([50, 90]);
       line([2 2],         [0 100], 'Color', 'b', 'LineStyle', '--'); 
       line([5.5 5.5],     [0 100], 'Color', 'b', 'LineStyle', '--'); 
       text(2.2, 72, 'Pitchover Start', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
       text(5.7, 65, 'Pitchover End', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
       PlottingUtilitySSAG(gcf, 'LV_3DOFoutConvertedData_QeVsTime');   

       figure;
       plot(vehicleObj.time(1:end-1),  round(vehicleObj.Body_NormalAccel/Mars_GenPhysCons.GRAV_ACCEL,5), 'linewidth', 2.0); 
       grid on; title('Body Normal Accel Comp vs. TOF (FS Boost)');
       xlabel('sec'); ylabel('g(s)'); 
       title('Body-Normal Acceleration vs TOF (FS)');
       set(0,'defaultAxesFontSize',16);    ax = gca;    ax.FontSize = 16;    ax.FontWeight = 'bold';
       xlim([0, 65.2]); ylim([-0.2, 1]);

       line([2 2],         [-1 1], 'Color', 'b', 'LineStyle', '--'); 
       line([5.5 5.5],     [-1 1], 'Color', 'b', 'LineStyle', '--'); 

       text(2.2, 0.80, 'Pitchover Start', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
       text(5.7, 0.70, 'Pitchover End', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');

       PlottingUtilitySSAG(gcf, 'LV_3DOFoutConvertedData_BodyNormalVsTime');  
      
  end
  
    outputStatus = 'complete';
   
end