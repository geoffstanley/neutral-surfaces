function [x, s, t] = tbs_vertsolve(S, T, X, K, x, branchmap, dfn, s0, t0, tolx, DX) %#codegen
%TBS_VERTSOLVE  Root finding of pressure or depth that matches equation of
%               state with multivalued function.
%
%
% [x, s, t] = tbs_vertsolve(S, T, X, K, x, branch, dfn, s0, t0, tolx, r)
% finds the pressure or depth x (within tolerance tolx) and its associated
% practical / Absolute salinity s and potential / Conservative temperature
% t, of a surface on which delta equals that determined by the multivalued
% function dfn, in an ocean of practical / Absolute salinity S and
% potential / Conservative temperature T at datasites whose pressure or
% depth is X.  The number of valid data points in each water column is
% given by K.  The equation of state is given by eos.m in the path, taking
% S, T, X as its 3 inputs. delta is the in-situ density anomaly or specific
% volume anomaly, defined as eos(s,t,x) - eos(s0,t0,x) where s,t are S,T
% interpolated from X to x, and s0, t0 are reference values.  In each water
% column I, the branch of the multivalued function for delta as a function
% of x is given by dfn(:,b) where b = branchmap(I).  This branch is a
% polynomial of x between dfn(1,b) and dfn(2,b) with coefficients given by
% dfn(3:end,b), and a linear extension of this polynomial outside this
% domain. Solutions are sought in the interval [A,B] where A = max(X(1,I),
% dfn(1,b) - DX) and B = min(X(K(I),I), dfn(2,b) + DX), i.e. within the
% valid X range of the local water column AND within DX of the polynomial
% domain of the local branch. (If X is a vector, the I indexing of X is
% dropped.)  Limiting this search domain with DX helps prevent the surface
% from jumping between multiple solutions in weakly stratified waters such
% as around Antarctica.
%
%
% --- Input:
% S [M, N]: practical/Absolute salinity at each data point on each cast
% T [M, N]: potential/Conservative temperature at each data point on each cast
% X [M, N]: pressure [dbar] or depth [m, positive] at each data point on each cast
% K [1, N]: number of valid data points on each cast
% x [1, N]: initial pressure [dbar] or depth [m] of the surface at each cast
% branchmap [1, N]: branch that each cast belongs to.
%                   Must have 1 <= branch(i) <= A, for all 1 <= i <= N.
% dfn [*, A]: domain and coefficients for each branch's polynomial function
%             for delta as a function of x, the b'th branch being defined
%             by dfn(:,b), which can be evaluated by pvaln or pvallin.
% s0 [1, 1]: reference S value for delta
% t0 [1, 1]: reference T value for delta
% tolx [1, 1]: tolerance on pressure [dbar] or depth [m] for bisection solver
% DX [1, 1]: seek solutions in the domain of the local dfn branch expanded
%            by DX in both directions.
%
% Note: above, N is the number of water columns (possibly including land).
%              M the maximum number of data points on any cast.
%              A is the number of branches for the multivalued function dfn.
%              * is a wildcard.
%
% Note: variables can actually be higher dimensional, e.g. N = [ni, nj],
%       and p can be any dimensional matrix, so long as it has N elements
%       in total.
%
% Note: X must increase along its first dimension.
%
%
% --- Output:
% x [same as input x]: pressure or depth of the updated surface
% s [same as input x]: practical / Absolute salinity of the updated surface
% t [same as input x]: potential / Conservative temperature of the updated surface
%
%
% --- Requirements:
% eos, interp_firstdim_twovars, bisectguess, pvallin
%
%
% --- Code generation:
% See tbs_vertsolve_codegen.m

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
    b = branchmap(I);
    if b > 0
        k = K(I);
        if k > 1
            
            % Bisect in the water column +/- DX metres or dbar from the
            % current surface to find a solution to
            % eos_diff(x(i), S(:,i), T(:,i), X(:,i), dfn(:,b), s0, t0) == 0
            x(I) = bisectguess(@diff_fun, max(X(idx+1), dfn(1,b) - DX), min(X(idx+k), dfn(2,b) + DX), tolx, x(I), ...
                S(1:k,I), T(1:k,I), X((idx+1:idx+k).'), dfn(:,b), s0, t0);
        end
    end
    idx = idx + MX;
end

% Interpolate S and T onto the updated surface
[s, t] = interp_firstdim_twovars(reshape(x, [1 size(x)]), X, S, T);

end


function out = diff_fun(x, S, T, X, dfn, s0, t0)
% The difference between delta evaluated (a) using the local branch of the
% multivalued function, and (b) using the equation of state with the local
% water properties.

% Evaluate water properties on the surface
[s,t] = interp_firstdim_twovars(x, X, S, T);

% Evaluate the delta difference
out = pvallin(dfn, x) - ( eos(s, t, x) - eos(s0, t0, x) );

end
