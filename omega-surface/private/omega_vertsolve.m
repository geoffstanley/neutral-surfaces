function [x, s, t] = omega_vertsolve(SppX, TppX, X, BotK, s, t, x, tolx, phi) %#codegen
%OMEGA_VERTSOLVE  Root finding of a new surface with a specified a density 
%                 difference from the current surface
%
%
% [s, t, x] = omega_vertsolve(S, T, X, BotK, s0, t0, x0, tolx, phi)
% determines pressures x that satisfy
%   |eos(S_n(x(n)), T_n(x(n)), (x(n) + x0(n))/2) + phi(n) - eos(s0(n), t0(n), (x(n) + x0(n))/2)| < tolx
% where S_n(x') and T_n(x') are interpolants of the data
%   (S(:,n),X(:,n)) and (T(:,n),X(:,n)).
% This also determines salinities s and temperatures t such that
%   s(n) = S_n(x(n)) and t(n) = T_n(x(n)).
% The function eos.m, found in the same directory as this function,
% determines either the in-situ density or the specific volume.
%
%
% --- Input
% SppX [O, K-1, N]: coefficients for piecewise polynomial for practical 
%                   / Absolute Salinity in terms of X
% TppX [O, K-1, N]: coefficients for piecewise polynomial for potential 
%                   / Conservative Temperature in terms of X
% X [K, N]: knots for the pressure or depth of the casts
% BotK [1, N]: number of valid data points on each cast
% s [1, N]: practical / Absolute salinity on the initial surface
% t [1, N]: potential / Conservative temperature on the initial surface
% x [1, N]: pressure or depth on the initial surface
% tolx [1,1]: precision of solution in pressure or depth
% phi [1,N]: the desired in-situ density or specific volume change of the
%            surface
%
% Note: O is the order of the piecewise polynomials
%       K is the maximum number of knots in these piecewise polynomials, 
%           i.e. the maximum number of bottles in any cast
%       N is the number of water columns (possibly including land).
%
% Note: X must increase along its first dimension.
%
% Note: X can have size [K, 1], in which case it is used for each cast.
%
% Note: variables can actually be higher dimensional, e.g. N = [ni, nj],
%       and x can be any dimensional matrix, so long as it has N elements
%       in total.
%
% Note: BotK should be given by 
%           BotK = squeeze(sum(uint16(isfinite(S)), 1, 'native'));
%
%
% --- Output
% x [same as input x]: pressure or depth of the updated surface
% s [same as input x]: practical / Absolute salinity of the updated surface
% t [same as input x]: potential / Conservative temperature of the updated surface
%
%
% --- Units
% The units of s, t, x, S, T, X, tolx, and phi are determined by the
% function eos.m.

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

N = numel(x);
[K,XN] = size(X);
KX = K * double(XN > 1);

% Loop over each water column
nX = 0;
for n = 1:N
    d = phi(n);
    k = BotK(n);
    if ~isnan(d) && k > 0
        
        % Select this water column
        SppXn = SppX(:,1:k-1,n);
        TppXn = TppX(:,1:k-1,n);
        Xn = X((nX+1:nX+k).');
        
        % DEV: The following attempts to limit the search direction based
        % on phi's sign, but in practice this doesn't seem to work well.
        %{
        if phiI > 0  %, only search in one direction by setting lower bound (upper bound?)
            % Target specific volume is larger than current specific volume: move up, to lower pressures
            x(n) = bisectguess(@diff_fun, X(nX+1), x(n), tolx, x(n), ...
                SppXn, TppXn, Xn, s(n), t(n), x(n), err);
        else
            % Target specific volume is smaller than current specific volume: move down, to higher pressures
            x(n) = bisectguess(@diff_fun, x(n), X(nX+k), tolx, x(n), ...
                SppXn, TppXn, Xn, s(n), t(n), x(n), err);
        end
        %}
        
        % Bisect through the whole water column to find a solution to the
        % nonlinear root finding problem
        x(n) = bisectguess(@eos_diff, Xn(1), Xn(k), tolx, x(n), ...
            SppXn, TppXn, Xn, s(n), t(n), x(n), d);
        
        % Interpolate S and T onto the updated surface
        [s(n),t(n)] = ppc_val2(Xn, SppXn, TppXn, x(n));
        
    end
    nX = nX + KX;
end

end


function out = eos_diff(x, SppX, TppX, X, s0, t0, x0, d)
% Evaluate difference between (a) eos at location on the cast (S, T, X)
% where the pressure or depth is x, and (b) d + eos of the bottle (s0, t0,
% x0); here, eos is always evaluated at the average pressure or depth, (x +
% x0)/2.

% Interpolate S and T to the current pressure or depth
[s,t] = ppc_val2(X, SppX, TppX, x);

% Average the current pressure or depth and the original pressure or depth
x_avg = (x + x0) / 2;

% Calculate the density or specific volume difference
out =  eos(s, t, x_avg) + d - eos(s0, t0, x_avg) ;

end