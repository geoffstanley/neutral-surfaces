function [p,s,t] = delta_surf_vertsolve(Sppc, Tppc, P, BotK, p, s_ref, t_ref, d0, tolp) %#codegen
%DELTA_SURF_VERTSOLVE  Helper function for delta_surf, solving non-linear
%                      root finding problem in each water column.
%
%
% Note: to ensure the generated MEX function gives the same output as
% running this in native MATLAB, the MEX function must be generated with K
% specified as an integer class, not double: we use uint16.

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2020  Geoff Stanley
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


N = numel(p);
Pmat = ~isvector(P);

s = nan(size(p));
t = nan(size(p));

% Loop over each cast
for n = 1:N
    k = BotK(n);
    if k > 1
        
        % Select this water column
        Sppcn = Sppc(:,1:k-1,n);
        Tppcn = Tppc(:,1:k-1,n);
        if Pmat
          Pn = P(1:k,n);
        else
          Pn = P((1:k).'); % .' is for codegen, so X and (1:k).' both column vectors
        end
        
        % Search for a sign-change, expanding outward from an initial guess 
        [lb, ub] = fzero_guess_to_bounds(@myfcn, p(n), Pn(1), Pn(k), ...
          Sppcn, Tppcn, Pn, s_ref, t_ref, d0);
        
        if ~isnan(lb)
          % A sign change was discovered, so a root exists in the interval.
          % Solve the nonlinear root-finding problem using Brent's method
          p(n) = fzero_brent(@myfcn, lb, ub, tolp, ...
            Sppcn, Tppcn, Pn, s_ref, t_ref, d0);
          
          % Interpolate S and T onto the updated surface
          [s(n),t(n)] = ppc_val2(Pn, Sppcn, Tppcn, p(n));
        else
          p(n) = nan;
          s(n) = nan;
          t(n) = nan;
        end
        
    end
end

end

function out = myfcn(p, Sppc, Tppc, P, s_ref, t_ref, d0)
[s,t] = ppc_val2(P, Sppc, Tppc, p);
out = eos(s, t, p) - eos(s_ref, t_ref, p) - d0 ;
end