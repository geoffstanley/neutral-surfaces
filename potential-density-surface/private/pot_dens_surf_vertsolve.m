function [x,s,t] = pot_dens_surf_vertsolve(SppX, TppX, X, BotK, x, xref, val, tolx) %#codegen
%POT_DENS_SURF_VERTSOLVE  Helper function for pot_dens_surf, solving 
%                         non-linear root finding problem in each water column
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
%
% Modified by : --
% Date        : --
% Changes     : --

N = numel(x);
Xmat = ~isvector(X);


s = nan(size(x));
t = nan(size(x));

% Loop over each cast
for n = 1:N
    k = BotK(n);
    if k > 1
        
        % Select this water column
        SppXn = SppX(:,1:k-1,n);
        TppXn = TppX(:,1:k-1,n);
        if Xmat
          Xn = X(1:k,n);
        else
          Xn = X((1:k).'); % .' is for codegen, so X and (1:k).' both column vectors
        end
        
        % Search for a sign-change, expanding outward from an initial guess 
        [lb, ub] = fzero_guess_to_bounds(@myfcn, x(n), Xn(1), Xn(k), ...
          SppXn, TppXn, Xn, xref, val);
        
        if ~isnan(lb)
          % A sign change was discovered, so a root exists in the interval.
          % Solve the nonlinear root-finding problem using Brent's method
          x(n) = fzero_brent(@myfcn, lb, ub, tolx, ...
            SppXn, TppXn, Xn, xref, val);
          
          % Interpolate S and T onto the updated surface
          [s(n),t(n)] = ppc_val2(Xn, SppXn, TppXn, x(n));
        else
          x(n) = nan;
          s(n) = nan;
          t(n) = nan;
        end
        
    end
end

end

function out = myfcn(x, SppX, TppX, X, xref, val)
[s,t] = ppc_val2(X, SppX, TppX, x);
out = eos(s, t, xref) - val;
end