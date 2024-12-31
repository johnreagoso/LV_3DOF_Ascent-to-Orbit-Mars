function [outputStruct] = ECEFtoECI_Convert_Mars_StandAlone(inputStruct)

% Note: this conversion script origin is for Earth centered trajectories (ECEF, ECI etc.). Any use of 'ECI', 'ECEF' etc. is a holdover 
% and refers to MCI (Mars Centered Inertial) or MCMF (Mars Centered Mars Fixed) coordinate frames. Future versions of this script will
% be modified accordingly. 

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
      
      % Below equations below to a scheme used while at APL that directly converts ECEF 
      % acceleration to ECI acceleration. Matches with TAOS output. 
      
      % Vehicle ECI acceleration compute:
%        ddEi = inputStruct.ddE - 2.*abs(marsRotRate).*inputStruct.dF(1:end-1) - ((abs(marsRotRate))^2).*inputStruct.E(1:end-1);
%        ddFi = inputStruct.ddF + 2.*abs(marsRotRate).*inputStruct.dE(1:end-1) - ((abs(marsRotRate))^2).*inputStruct.F(1:end-1);
%        
%        inputStruct.ddX = ddEi.*cos(omega_Mars(1:end-1)) - ddFi.*sin(omega_Mars(1:end-1));
%        inputStruct.ddY = ddEi.*sin(omega_Mars(1:end-1)) + ddFi.*cos(omega_Mars(1:end-1));
%        inputStruct.ddZ = inputStruct.ddG;      

      outputStruct = inputStruct;
end
