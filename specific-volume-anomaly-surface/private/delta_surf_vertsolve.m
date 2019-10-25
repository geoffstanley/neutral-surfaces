function [x,s,t] = delta_surf_vertsolve(SppX, TppX, X, BotK, x, s0, t0, d0, tolx) %#codegen
%DELTA_SURF_VERTSOLVE  Helper function for delta_surf, solving non-linear
%                      root finding problem in each water column.
%
%
% Note: to ensure the generated MEX function gives the same output as
% running this in native MATLAB, the MEX function must be generated with K
% specified as an integer class, not double: we use uint16.

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

s = nan(size(x));
t = nan(size(x));

% Loop over each cast
nX = 0;
for n = 1:N
    k = BotK(n);
    if k > 1
        
        % Select this water column
        SppXn = SppX(:,1:k-1,n);
        TppXn = TppX(:,1:k-1,n);
        Xn = X((nX+1:nX+k).');
        
        % Solve non-linear problem at each water column
        x(n) = bisectguess(@diff_fun, Xn(1), Xn(k), tolx, x(n), ...
            SppXn, TppXn, Xn, s0, t0, d0);
        
        % Interpolate S and T onto the updated surface
        [s(n),t(n)] = ppc_val2(Xn, SppXn, TppXn, x(n));
        
    end
    nX = nX + KX;
end

end

function out = diff_fun(x, SppX, TppX, X, s0, t0, d0)
[s,t] = ppc_val2(X, SppX, TppX, x);
out = eos(s, t, x) - eos(s0, t0, x) - d0 ;
end