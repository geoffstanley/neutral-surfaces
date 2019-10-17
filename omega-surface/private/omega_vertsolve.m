function [x, s, t] = omega_vertsolve(S, T, X, K, s, t, x, tolx, phi) %#codegen
%OMEGA_VERTSOLVE  Update surface so  in each water column for
%
%
% [s, t, x] = omega_vertsolve(S, T, X, K, s0, t0, x0, tolx, phi)
% determines pressures x that satisfy
% |eos(S_i(x(i)), T_i(x(i)), (x(i) + x0(i))/2) + phi(i) - eos(s0(i), t0(i), (x(i) + x0(i))/2)| < tolx
% where S_i(x') and T_i(x') are interpolants of the data
% (S(:,i),X(:,i)) and (T(:,i),X(:,i)).
% This also determines salinities s and temperatures t such that
% s(i) = S_i(x(i)) and t(i) = T_i(x(i)).
% The function eos.m, found in the same directory as this function,
% determines either the in-situ density or the specific volume.
%
%
% --- Input
% S [M, N]: practical / Absolute Salinity of the casts
% T [M, N]: potential / Conservative Temperature of the casts
% X [M, N]: pressure or depth of the casts
% K [1, N]: number of valid data points on each cast
% s [1, N]: practical / Absolute salinity on the initial surface
% t [1, N]: potential / Conservative temperature on the initial surface
% x [1, N]: pressure or depth on the initial surface
% tolx [1,1]: precision of solution in pressure or depth
% phi [1,N]: the desired in-situ density or specific volume change of the
%            surface
%
% Note: N is the number of water columns (possibly including land).
%       M the maximum number of data points on any cast.
%       Variables can actually be higher dimensional, e.g. N = [nx, ny],
%       and s, t, and x can be any dimensional matrix, so long as they have
%       N elements in total.
%       X can have size [M, 1], in which case it is used for each cast.
%       K should be given by K = squeeze(sum(isfinite(S), 1));
%
% --- Output
% x [same as input x]: pressure or depth of the updated surface
% s [same as input x]: practical / Absolute salinity of the updated surface
% t [same as input x]: potential / Conservative temperature of the updated surface
%
%
% --- Units
% The units of s, t, x, S, T, X, s0, t0, p0, tolx, and phi are determined
% by the function eos.m.
%
% --- Requirements:
% bisectguess.m, interp_firstdim_twovars.m, eos.m

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


% Inputs s0, t0, and p0 are named s, t, x so operations are done in-place.

[M,N] = size(S);
[~,Xj] = size(X);
MX = M * double(Xj > 1); % Used to handle when X is actually a vector

% Loop over each water column
idx = 0;
for i = 1:N
    d = phi(i);
    k = K(i);
    if ~isnan(d) && k > 0
        % The following attempts to limit the search direction based on the
        % sign of phi, but in practice this does not seem to work well.
        %{
        if phiI > 0  %, only search in one direction by setting lower bound (upper bound?)
            % Target specific volume is larger than current specific volume: move up, to lower pressures
            x(i) = bisectguess(@diff_fun, X(idx+1), x(i), tolx, x(i), ...
                S(1:k,i), T(1:k,i), X((idx+1:idx+k).'), s(i), t(i), x(i), err);
        else
            % Target specific volume is smaller than current specific volume: move down, to higher pressures
            x(i) = bisectguess(@diff_fun, x(i), X(idx+k), tolx, x(i), ...
                S(1:k,i), T(1:k,i), X((idx+1:idx+k).'), s(i), t(i), x(i), err);
        end
        %}
        
        % Bisect through the whole water column to find a solution to
        % eos_diff(x(i), S(:,i), T(:,i), X(:,i), s(i), t(i), x(i), d) == 0
        x(i) = bisectguess(@eos_diff, X(idx+1), X(idx+k), tolx, x(i), ...
            S(1:k,i), T(1:k,i), X((idx+1:idx+k).'), s(i), t(i), x(i), d);
    end
    idx = idx + MX;
end

% Interpolate S and T onto the updated surface
[s, t] = interp_firstdim_twovars(reshape(x, [1 size(x)]), X, S, T);

end


function out = eos_diff(x, S, T, X, s0, t0, x0, phi)
% Evaluate difference between (a) eos at location on the cast (S, T, X)
% where the pressure or depth is x, and (b) phi + eos of the bottle (s0,
% t0, x0); here, eos is always evaluated at the average pressure or depth,
% (x + x0)/2.

% Interpolate S and T to the current pressure or depth
[s,t] = interp_firstdim_twovars(x, X, S, T);

% Average the current pressure or depth and the original pressure or depth
x_avg = (x + x0) / 2;

% Calculate the specific volume error
out =  eos(s, t, x_avg) + phi - eos(s0, t0, x_avg) ;

end