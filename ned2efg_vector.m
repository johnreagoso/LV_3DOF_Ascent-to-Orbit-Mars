% ned2efg_vector rotates a North-East-Down (NED) Relative-Local(RL) vector 
% about latgd and long to Earth-Centered-Earth-Fixed (ECEF) frame

% -Input:
% latgd,long    Latitude, Longitude origin on Earth Surface in radians 
% n,e,d      Vector in NED frame
%
% -Output:
% E,F,G      Vector in ECEF frame in same units as input N,E,D     

function [E,F,G] = ned2efg_vector(latgd, long, n, e, d)

E = -n.* sin(latgd).* cos(long) - e.* sin(long) - d.* cos(latgd).* cos(long);
F = -n.* sin(latgd).* sin(long) + e.* cos(long) - d.* cos(latgd).* sin(long);
G =  n.* cos(latgd) - d.* sin(latgd);