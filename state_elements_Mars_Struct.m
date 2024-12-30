%{
    ##########################__79_Character_Header__###########################
    PROCEDURE NAME:    state_elements_Mars_Struct.m

    PURPOSE:    This m-file computes the orbital elements of a spacecraft utilizing only
    its position/velocity vectors required to be synched at a common epoch. 
    (Common epochs are assumed within this script).

	External References:
 	Name                                                        Purpose
	----                                                        -------
    n/a                                                         n/a

    INPUTS:
		1.  (1) vehicle object          -

	OUTPUTS:
		1. 	eccentricity                         (e)              -
        2.  inclination (rad)                    (Incl)          -
        3.  rt ascention of Ascending Node(rad)  (RAAN)          -
        4.  argument of periapse (rad)           (ArgPer)        -
        5.  true anomaly (rad)                   (TrueAnom)      -
        6.  true longitude (rad)                 (TrueLong)      -

  ##########################__79_Character_Header__###########################
   External References:
   NAME								PURPOSE
   ----                             -------  (none currently)
 
   Development History:
   NAME				DATE			DESCRIPTION OF CHANGE

   Reagoso, J.D.  	OCT 25 2024     Script initialized 
 ---------------------------------------------------------------------------------------------------
%}

%% Main Body 

function[vehicleObj] = state_elements_Mars_Struct(vehicleObjInput)

Mu_mars = Mars_GenPhysCons.GM_KM;

kVector = [0;0;1];          % kHat unit vector

pos_vector = [vehicleObjInput.X;  vehicleObjInput.Y;    vehicleObjInput.Z];
vel_vector = [vehicleObjInput.dX; vehicleObjInput.dY;   vehicleObjInput.dZ];

pos_magntd = norm(pos_vector);      % range Magnitude
vel_magntd = norm(vel_vector);      % velocity magnitude

h_vector   = cross(pos_vector, vel_vector) ;        % Spec angular momentum vector
h_norm     = norm(h_vector);                        % Spec angular momentum magnitude
incl       = acos(h_vector(3)/norm(h_vector)) ;     % inclination (rad)
node_vector = cross(kVector, h_vector);             % node vector

p = (h_norm^2)/Mu_mars; % semi-latus rectum
 
e_vector     = (1/Mu_mars)*(( vel_magntd^2 - Mu_mars/pos_magntd)*pos_vector - dot(pos_vector, vel_vector)*vel_vector); % eccen vector                                                                      );
e           = norm(e_vector) ;                       % eccentricity (unitless)

perigee = p/(1+e);  % perigee
apogee  = p/(1-e);  % apogee

%% Rigth Ascension of the Ascending Node(RAAN) (rad):

if (node_vector(2) < 0.00)
    RAAN = 2*pi - acos(node_vector(1)/norm(node_vector));    
else
    RAAN = acos(node_vector(1)/norm(node_vector));
end

%% Argument of Periapse (rad):

if (e_vector(3) < 0.0)
    arg_per = 2*pi - acos(dot(node_vector, e_vector)/(norm(node_vector)*norm(e_vector))); 
else
    arg_per = acos(dot(node_vector, e_vector)/(norm(node_vector)*norm(e_vector))); 
end

%% True Anomaly (rad)

if (dot(pos_vector, vel_vector) < 0.0)
    true_anom = 2*pi - acos(dot(e_vector, pos_vector)/(e*pos_magntd)) ;
else
    true_anom = acos(dot(e_vector, pos_vector)/(e*pos_magntd)) ; 
end

%% Orbital Paramter Field Output (in deg): 
vehicleObj.time     = vehicleObjInput.time;
vehicleObj.e        = e;
vehicleObj.incl     = rad2deg(incl);
vehicleObj.TrueAnom = rad2deg(true_anom);
vehicleObj.RAAN     = rad2deg(RAAN);
vehicleObj.ArgPer   = rad2deg(arg_per);
vehicleObj.TrueLong = rad2deg(RAAN + arg_per + true_anom);
vehicleObj.perigee  = perigee;
vehicleObj.apogee   = apogee;
vehicleObj.apogee_alt  = apogee - Mars_GenPhysCons.RE_EQ;
vehicleObj.perigee_alt = perigee - Mars_GenPhysCons.RE_EQ;
vehicleObj.SMA         = 0.5*(perigee + apogee);

end