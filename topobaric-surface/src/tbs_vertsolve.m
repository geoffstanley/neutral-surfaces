function p = topobaric_update(S, T, P, K, p, branchmap, dfn, s0, t0, tol) %#codegen
%TOPOBARIC_UPDATE  Root finding of pressure that matches equation of state with multivalued function.
%
%
% p = topobaric_update(S, T, P, K, p, branch, dfn, s0, t0, tol)
% finds (within a tolerance of tol pressure units) the pressure p of a
% surface on which delta equals that determined by the multivalued function
% dfn, in an ocean of practical / Absolute salinity S, potential /
% Conservative temperature T, and pressure P. The number of valid data
% points in each water column is given by K. The equation of state is given
% by eos.m in the path. delta is the in-situ density minus a reference
% in-situ density profile, defined at any pressure p' as eos(s0,t0,p')
% where the reference S and T values are s0 and t0. If eos.m instead
% calculates the specific volume, delta is the specific volume anomaly,
% defined similarly. Each single-valued branch of the multivalued function
% dfn is given by a column of fn. The branch that applies at each water
% column is given by branchmap.
%
% 
% --- Input:
% S [M, N]: practical/Absolute salinity at each data point on each cast
% T [M, N]: potential/Conservative temperature at each data point on each cast
% P [M, N]: pressure at each data point on each cast [dbar]
% K [1, N]: number of valid data points on each cast
% p [1, N]: initial pressure of the surface at each cast [dbar]
% branchmap [1, N]: branch that each cast belongs to.
%                   Must have 1 <= branch(i) <= A, for all 1 <= i <= N.
% dfn [5, A]: domain and coefficients for each branch's quadratic function
%             for delta as a function of pressure. dfn(:,a) defines the
%             a'th branch of the multivalued function for delta as a
%             function of pressure, which can be evaluated by pvaln or
%             pvallin.
% s0 [1, 1]: reference S value for delta 
% t0 [1, 1]: reference T value for delta 
% tol [1, 1]: tolerance on pressure for bisection solver [dbar]
%
% Note: above, N is the number of water columns (possibly including land).
%              M the maximum number of data points on any cast.
%              A is the number of branches for the multivalued function dfn
%
% Note: variables can actually be higher dimensional, e.g. N = [nx, ny],
%       and p can be any dimensional matrix, so long as it has N elements
%       in total.
%
% Note: P must be increasing along its first dimension.
%
%
% --- Output:
% p [same as input p]: pressure of new surface [dbar]
%
%
% --- Requirements:
% eos, bisectguess, interp1q2, pvallin
%
%
% --- Code generation:
% See ../run_codegen.m 

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

[M,N] = size(S);
[~,pj] = size(P);
MP = M * double(pj > 1);

% Loop over each water column
idx = 0;
for I = 1:N
    b = branchmap(I);
    if b > 0
        k = K(I);
        if k > 1
            p(I) = bisectguess(@diff_fun, P(idx+1), P(idx+k), tol, p(I), ...
                S(1:k,I), T(1:k,I), P((idx+1:idx+k).'), dfn(:,b), s0, t0);
        end
    end
    idx = idx + MP;
end
end


function out = diff_fun(p, S, T, P, dfn, s0, t0)
% The difference in delta between the multivalued function and the equation
% of state.
[s,t] = interp1qn2(p, P, S, T);

out = pvallin(dfn, p) - ( eos(s, t, p) - eos(s0, t0, p) );

end
