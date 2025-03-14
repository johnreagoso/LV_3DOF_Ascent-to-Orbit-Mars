function [outputStruct] = ECEFtoECI_Convert_Mars_StandAlone(inputStruct)

% Note: this conversion script origin is for Earth centered trajectories (ECEF, ECI etc.). Any use of 'ECI', 'ECEF' etc. is a description 
% holdover as this script refers to MCI (Mars Centered Inertial) or MCMF (Mars Centered Mars Fixed) coordinate frames. 
% Future versions of this script will be modified accordingly. 

    if inputStruct.rotatingMars == 1 
        marsRotRate = Mars_GenPhysCons.OMEGA;  %rad/sec
    else
        marsRotRate = 0.00;  %rad/sec  No Mars rotation modeled!
    end
      
      elapsedTOF = inputStruct.time;
      
      omega_init = 0;   % !!! for our LV trajectory purposes- at beginning of a trajectory at t = 0, ECI and ECEF systems are aligned. !!!
                        % !!! for other applications, this is not suitable. !!!
      
      omega_Mars = omega_init + abs(marsRotRate)*elapsedTOF;  %quantity expressing angular distance (in rad) that the Mars rotates/moves around the z-axis.
      
      % Position Conversion - below represents ECI = ECEF_pos*[Rz(t)]- rotating about the z-axis..  
      % Rz = [cos(omega) -sin(omega);   sin(omega)  cos(omega)];
      inputStruct.X = cos(omega_Mars).*inputStruct.E - sin(omega_Mars).*inputStruct.F;  % km
      inputStruct.Y = sin(omega_Mars).*inputStruct.E + cos(omega_Mars).*inputStruct.F;  % km
      inputStruct.Z = inputStruct.G;                                                     % km
      
      % Velocity Conversion- must use calculus chain-rule due to changing
      % quantity of Mars's rotation (omega) within the DCM:
      % ECI_vel = d/dt[DCM*ECEF_position] -->  d/dt(DCM)*ECEF_position] + DCM* d/dt(ECEF_position)
      
      inputStruct.dX = cos(omega_Mars).*inputStruct.dE - sin(omega_Mars).*inputStruct.dF - ...
        abs(marsRotRate).*(sin(omega_Mars).*inputStruct.E + cos(omega_Mars).*inputStruct.F);    % km/sec
      
      inputStruct.dY = sin(omega_Mars).*inputStruct.dE + cos(omega_Mars).*inputStruct.dF + ...
        abs(marsRotRate).*(cos(omega_Mars).*inputStruct.E - sin(omega_Mars).*inputStruct.F);    % km/sec
      
      inputStruct.dZ = inputStruct.dG;                                                             % km/sec  

      outputStruct = inputStruct;
end
