function [E,F,G] = ned2efg_vector(latgd, long, n, e, d)

% function [E,F,G] = ned2efg_vector( latgd, long, n, e, d)
% ned2efg_vector rotates a NED vector about latgd and long to ECEF
% 
% -Input:
% latgd,long    Latitude, Longitude origin on Earth Surface in radians 
% n,e,d      Vector in NED frame
%
% -Output:
% E,F,G      Vector in ECEF frame in same units as input N,E,D     
%
% Examples
% function [E,F,G] = ned2efg_vector( latgd, long, n, e, d)
% function [E_veh, F_veh, G_veh] = ned2efg_vector( latgd_veh, long_veh, n_veh, e_veh, d_veh)

E = -n .* sin(latgd) .* cos(long) - e .* sin(long) - d .* cos(latgd) .* cos(long);
F = -n .* sin(latgd) .* sin(long) + e .* cos(long) - d .* cos(latgd) .* sin(long);
G =  n .* cos(latgd) - d .* sin(latgd);