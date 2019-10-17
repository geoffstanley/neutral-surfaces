function [x, s, t, success] = ntp_bottle_to_cast(S, T, X, sB, tB, xB, tolx) %#codegen
%NTP_BOTTLE_TO_CAST  Find a bottle's level of neutral buoyancy in a water
%                    column, using the Neutral Tangent Plane relationship.
%
% [x, s, t] = ntp_bottle_to_cast(S, T, X, sB, tB, xB, tolx)
% finds (s, t, x), with precision in x of tolx, that is at the level of
% neutral buoyancy for a fluid bottle of (sB, tB, xB) in a water column of
% (S, T, X).  Specifically, s and t are given by
%   [s,t] = interp_firstdim_twovars(x, X, S, T)
% and x satisfies
%      eos(s, t, x') = eos(sB, tB, x')
% where eos is the equation of state given by eos.m in MATLAB's path,
% and   x' is in the range [x_avg - tolx/2, x_avg + tolx/2],
% and   x_avg = (xB + x) / 2 is the average of the fluid bottle's original
%                          and final pressure or depth.
%
% [x, s, t, success] = ntp_bottle_to_cast(...)
% returns a flag value success that is true if a valid solution was found,
% false otherwise.
%
% For a non-Boussinesq ocean, X, xB, and x are pressure.
% For a Boussinesq ocean, X, xB, and x are depth.
%
%
% --- Input:
% S [nk, 1]: practical / Absolute salinity in water column
% T [nk, 1]: potential / Conservative temperature in water column
% X [nk, 1]: pressure or depth in water column
% sB [1 , 1]: practical / Absolute salinity of current bottle
% tB [1 , 1]: potential / Conservative temperature of current bottle
% xB [1 , 1]: pressure or depth of current bottle
% tolx [1, 1]: tolerance for solving the level of neutral buoyancy (same
%             units as X and xB)
%
% Note: physical units for S, T, X, sB, tB, xB, x, s, t  are
% determined by eos.m.
%
% Note: X must increase monotonically along its first dimension.
%
%
% --- Output:
% x [1, 1]: pressure or depth in water column at level of neutral buoyancy
% s [1, 1]: practical / Absolute salinity in water column at level of neutral buoyancy
% t [1, 1]: potential / Conservative temperature in water column at level of neutral buoyancy
% success [1,1]: true if a valid solution was found, false otherwise.
%
%
% --- Requirements:
% eos, interp_firstdim_twovars
% bisectguess - https://www.mathworks.com/matlabcentral/fileexchange/69710
%
%
% --- Code generation:
% See ntp_bottle_to_cast_codegen.m

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2019  Geoff Stanley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com
% Version   : 2.0.0
%
% Modified by : --
% Date        : --
% Changes     : --

x = bisectguess(@diff_fun, X(1), X(end), tolx, xB, S, T, X, sB, tB, xB);

if isnan(x)
    s = nan;
    t = nan;
    success = false;
else
    [s,t] = interp_firstdim_twovars(x, X, S, T);
    success = true;
end
end

function out = diff_fun(x, S, T, X, sB, tB, xB)
% Evaluate difference between (a) eos at location on the cast
% (S, T, X) where the pressure is p, and (b) eos of the bottle
% (sB, tB, xB); here, density volume is always evaluated at the average
% pressure, (x + x0)/2.
[s,t] = interp_firstdim_twovars(x, X, S, T);
x_avg = (xB + x) / 2;
out = eos(sB, tB, x_avg) - eos(s, t, x_avg);
end