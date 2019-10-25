function [x, s, t] = tbs_vertsolve(SppX, TppX, X, BotK, s, t, x, branchmap, d_fn, s0, t0, tolx, DX) %#codegen
%TBS_VERTSOLVE  Root finding of pressure or depth that matches equation of
%               state with multivalued function.
%
%
% [x, s, t] = tbs_vertsolve(SppX, TppX, X, BotK, s, t, x, branch, d_fn, s0, t0, tolx, DX)
% finds the pressure or depth x (within tolerance tolx) and its associated
% practical / Absolute salinity s and potential / Conservative temperature
% t, of a surface on which delta equals that determined by the multivalued
% function d_fn, in an ocean whose practical / Absolute salinity and
% potential / Conservative temperature as functions of pressure or depth X
% are given by piecewise polynomials whose coefficients are SppX and TppX,
% and whose knots are X.  The number of valid data points in each water
% column is given by BotK.  The equation of state is given by eos.m in the
% path, taking S, T, X as its 3 inputs. delta is the in-situ density
% anomaly or specific volume anomaly, defined as eos(s,t,x) - eos(s0,t0,x)
% where s,t are S,T interpolated from X to x, and s0, t0 are reference
% values.  In each water column I, the branch of the multivalued function
% for delta as a function of x is given by d_fn(:,b) where b =
% branchmap(I).  This branch is a polynomial of x between d_fn(1,b) and
% d_fn(2,b) with coefficients given by d_fn(3:end,b), and a linear
% extension of this polynomial outside this domain. Solutions are sought in
% the interval [A,B] where A = max(X(1,I), d_fn(1,b) - DX) and B =
% min(X(BotK(I),I), d_fn(2,b) + DX), i.e. within the valid X range of the
% local water column AND within DX of the polynomial domain of the local
% branch. (If X is a vector, the I indexing of X is dropped.)  Limiting
% this search domain with DX helps prevent the surface from jumping between
% multiple solutions in weakly stratified waters such as around Antarctica.
% The inputs s and t are not used, but provided so these variables may be
% manipulated in-place.
%
%
% --- Input:
% SppX [O, K-1, N]: coefficients for piecewise polynomial for practical 
%                   / Absolute Salinity in terms of X
% TppX [O, K-1, N]: coefficients for piecewise polynomial for potential 
%                   / Conservative Temperature in terms of X
% X [K, N]: knots for the pressure or depth of the casts
% BotK [1, N]: number of valid data points on each cast
% s [1, N]: initial practical / Absolute salinity on the initial surface
% t [1, N]: initial potential / Conservative temperature on the initial surface
% x [1, N]: initial pressure [dbar] or depth [m] on the initial surface
% branchmap [1, N]: branch that each cast belongs to.
%                   Must have 1 <= branch(i) <= A, for all 1 <= i <= N.
% d_fn [D+3, A]: domain and coefficients for each branch's polynomial 
%             for delta as a function of x, the b'th branch being defined
%             by d_fn(:,b), which can be evaluated by pvaln or pvallin.
% s0 [1, 1]: reference S value for delta
% t0 [1, 1]: reference T value for delta
% tolx [1, 1]: tolerance on pressure [dbar] or depth [m] for bisection solver
% DX [1, 1]: seek solutions in the domain of the local d_fn branch expanded
%            by DX in both directions.
%
% Note: O is the order of the piecewise polynomials down each cast
%       K is the maximum number of knots in these piecewise polynomials, 
%           i.e. the maximum number of bottles in any cast
%       N is the number of water columns (possibly including land).
%       A is the number of branches for the multivalued function d_fn.
%       D is the degree of each branch's polynomial
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
% --- Output:
% x [same as input x]: pressure or depth of the updated surface
% s [same as input x]: practical / Absolute salinity of the updated surface
% t [same as input x]: potential / Conservative temperature of the updated surface

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

N = numel(x);
[K,XN] = size(X);
KX = K * double(XN > 1);

% Loop over each cast
nX = 0;
for n = 1:N
    b = branchmap(n);
    k = BotK(n);
    if b > 0 && k > 1
        
        % Select this water column
        SppXn = SppX(:,1:k-1,n);
        TppXn = TppX(:,1:k-1,n);
        Xn = X((nX+1:nX+k).');
        
        % Bisect in the water column +/- DX metres or dbar from the current
        % surface to find a solution to the nonlinear root finding problem
        x(n) = bisectguess(@diff_fun, max(Xn(1), d_fn(1,b) - DX), min(Xn(k), d_fn(2,b) + DX), tolx, x(n), ...
            SppXn, TppXn, Xn, d_fn(:,b), s0, t0);
        
        % Interpolate S and T onto the updated surface
        [s(n),t(n)] = ppc_val2(Xn, SppXn, TppXn, x(n));
        
    end
    nX = nX + KX;
end

end


function out = diff_fun(x, SppX, TppX, X, d_fn, s0, t0)
% The difference between delta evaluated (a) using the local branch of the
% multivalued function, and (b) using the equation of state with the local
% water properties.

% Evaluate water properties on the surface
[s,t] = ppc_val2(X, SppX, TppX, x);

% Evaluate the delta difference
out = pvallin(d_fn, x) - ( eos(s, t, x) - eos(s0, t0, x) );

end
