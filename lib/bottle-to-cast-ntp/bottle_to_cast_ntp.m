function [p2, s2, t2] = level_neutral_buoyancy(S2, T2, P2, K2, s1, t1, p1, tol) %#codegen
%LEVEL_NEUTRAL_BUOYANCY  A fluid parcel's level of neutral buoyancy in a water column.
% 
% p2 = level_neutral_buoyancy(S2, T2, P2, K2, s1, t1, p1, tol)
% finds the pressure p2 (with precision tol) that is the level of neutral
% buoyancy for a fluid parcel of known salinity s1, potential temperature
% t1, and pressure p1 in a water column with salinity S2, potential
% temperature T2, and pressure P2 (having non-NaN data from 1 : K2). 
% Specifically, p2 satisfies
%   eos(s2, t2, p0) = eos(s1, t1, p0),
% where s2 = interp1(P2, S2, p2),
% and   t2 = interp1(P2, T2, p2),
% and   p0 = (p1 + p2) / 2 is the average of the fluid parcel's original
%                          pressure and final pressure.
% The equation of state is given by eos.m in the path.
%
% [p2, s2, t2] = level_neutral_buoyancy(...)
% also returns the salinity s2 and potential temperature t2 at the level of
% neutral buoyancy.
%
% For a Boussinesq model, pass P2 as Z * Z2P and p1 as z1 * Z2P, where
% Z is the depth of each grid point in the water column [m, positive], z1
% is the depth of the fluid parcel [m, positive], and Z2P is the Boussinesq
% conversion from depth to pressure, typically Z2P = grav * rho_c * 1e-4,
% where grav is the gravitational acceleration and rho_c the Boussinesq
% reference density.
%
%
% --- Input:
% S2 [nz, 1]: practical / Absolute salinity at neighbouring cast
% T2 [nz, 1]: potential / Conservative temperature at neighbouring cast
% P2 [nz, 1]: pressure at neighbouring cast
% K2 [1 , 1]: number of non-nan data points in neighbouring cast
% s1 [1 , 1]: practical / Absolute salinity of current parcel
% t1 [1 , 1]: potential / Conservative temperature of current parcel
% p1 [1 , 1]: pressure of current parcel
% tol [1, 1]: tolerance for solving the level of neutral buoyancy [dbar] or [m]
%
% Above, nz is the maximum number of data points in any water column. 
%
% --- Output:
% p2 [1, 1]: pressure of neighbouring water column at level of neutral 
%            buoyancy
% s2 [1, 1]: practical / Absolute salinity of neighbouring water column at
%            level of neutral buoyancy
% t2 [1, 1]: potential / Conservative temperature of neighbouring water
%            column at level of neutral buoyancy
%
%
% --- Requirements:
% eos
% bisectguess - https://www.mathworks.com/matlabcentral/fileexchange/69710
% interp1q2 - https://www.mathworks.com/matlabcentral/fileexchange/69713
%
%
% --- Code generation:
% codegen can be run as follows (given a value for nz): 
% >> mexconfig = coder.config('mex');
% >> mexconfig.ExtrinsicCalls = false;
% >> mexconfig.ResponsivenessChecks = false;
% >> mexconfig.IntegrityChecks = false;
% >> codegen('level_neutral_buoyancy', '-config', mexconfig, '-o', 'level_neutral_buoyancy_mex', '-args', ...
%     {zeros(nz,1), zeros(nz,1), zeros(nz,1), 0, 0, 0, 0, 0});

% --- Copyright:
% Copyright 2019 Geoff Stanley
%
% This file is part of Topobaric Surface.
% 
% Topobaric Surface is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
% 
% Topobaric Surface is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with Topobaric Surface.  If not, see
% <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au 
% Email     : geoffstanley@gmail.com
% Version   : 1.0
%
% Modified by : --
% Date        : --
% Changes     : --

p2 = bisectguess(@diff_fun, P2(1), P2(K2), tol, p1, S2(1:K2), T2(1:K2), P2(1:K2), s1, t1, p1);
[s2,t2] = interp1qn2(p2, P2(1:K2), S2(1:K2), T2(1:K2,:));
end

function out = diff_fun(p2, S2, T2, P2, s1, t1, p1)
[s2,t2] = interp1qn2(p2, P2, S2, T2);
p0 = (p1 + p2) / 2;
out = eos(s1, t1, p0) - eos(s2, t2, p0);
end