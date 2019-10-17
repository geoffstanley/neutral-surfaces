function [x, s, t] = obs_vertsolve(S, T, X, K, x, dfnb, dfnc, s0, t0, tolx) %#codegen
%OBS_UPDATE  Root finding of pressure or depth that matches equation of
%            state with single-valued function.
%
%
% [x, s, t] = obs_vertsolve(S, T, X, K, x, branch, vafnp, s0, t0, tolx)
% finds the pressure or depth x (within tolerance tolx) and its associated
% practical / Absolute salinity s and potential / Conservative temperature
% t, of a surface on which delta equals that determined by the
% single-valued function determined by dfnb (spline break points) and dfnc
% (spline coefficients), in an ocean of practical / Absolute salinity S and
% potential / Conservative temperature T at datasites whose pressure or
% depth is X. The number of valid data points in each water column is given
% by K. The equation of state is given by eos.m in the path, taking S, T,
% and X as its 3 inputs. delta is the in-situ density anomaly or specific
% volume anomaly, defined as eos(s,t,x) - eos(s0,t0,x) where s,t are S,T
% interpolated from X to x, and s0, t0 are reference values.
%
%
% --- Input:
% S [M, N]: practical/Absolute salinity at each data point on each cast
% T [M, N]: potential/Conservative temperature at each data point on each cast
% X [M, N]: pressure [dbar] or depth [m, positive] at each data point on each cast
% K [1, N]: number of valid data points on each cast
% x [1, N]: initial pressure [dbar] or depth [m] of the surface at each cast
% dfnb [1,B]  : break points for the spline giving delta as a function of x
% dfnc [B-1,O]: coefficient matrix for the spline giving delta as a function of x
% s0 [1, 1]: reference S value for delta
% t0 [1, 1]: reference T value for delta
% tolx [1, 1]: tolerance on pressure [dbar] or depth [m] for bisection solver
%
% Note: above, N is the number of water columns (possibly including land).
%              M the maximum number of data points on any cast.
%              B is the number of break points in the spline.
%              O is the order of the spline.
%
% Note: variables can actually be higher dimensional, e.g. N = [ni, nj],
%       and x can be any dimensional matrix, so long as it has N elements
%       in total.
%
% Note: X must increase along its first dimension.
%
%
% --- Output:
% x [same as input x]: pressure or depth of the updated surface
% s [same as input x]: practical / Absolute salinity of the updated surface
% t [same as input x]: potential / Conservative temperature of the updated surface%
%
%
% --- Requirements:
% eos, interp_firstdim_twovars, bisectguess
%
%
% --- Acknowledgements:
% The sub-function ppval1 is adapted from MATLAB's function PPVAL.
%
%
% --- Code generation:
% See obs_vertsolve_codegen.m

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

[M,N] = size(S);
[~,xj] = size(X);
MX = M * double(xj > 1);

% Loop over each valid water column
idx = 0;
for I = 1:N
    k = K(I);
    if k > 1 && isfinite(x(I))
        x(I) = bisectguess(@diff_fun, X(idx+1), X(idx+k), tolx, x(I), ...
            S(1:k,I), T(1:k,I), X((idx+1:idx+k).'), dfnb, dfnc, s0, t0);
    end
    idx = idx + MX;
end

% Interpolate S and T onto the updated surface
[s, t] = interp_firstdim_twovars(reshape(x, [1 size(x)]), X, S, T);

end


function out = diff_fun(x, S, T, X, dfnb, dfnc, s0, t0)
% The difference in delta between the single-valued function and the
% equation of state.
[s,t] = interp_firstdim_twovars(x, X, S, T);

out = ppval1(dfnb, dfnc, x) - ( eos(s, t, x) - eos(s0, t0, x) );
end


function v = ppval1(b,c,xx)
% PPVAL1: evaluate 1-dimensional piecewise polynomials.
% b: pp.breaks
% c: pp.coefs
% xx: vector of evaluation locations
% Assumes that d == 1, where d == pp.dim.
% The following code is adapted from MATLAB's PPVAL.

xx = xx(:).';   % Ensure row vector
lx = numel(xx);

% for each evaluation site, compute its breakpoint interval
% (mindful of the possibility that xx might be empty)
if lx > 0
    [~,index] = histc(xx,[-inf,b(2:end-1),inf]);
else
    index = ones(1,lx);
end

% adjust for NaN. (Inf's are handled naturally by histc).
index(isnan(xx)) = 1;

% now go to local coordinates ...
xx = xx-b(index);

% ... and apply nested multiplication:
v = c(index,1);
for i=2:size(c,2)
    v = xx(:).*v + c(index,i);
end

end
