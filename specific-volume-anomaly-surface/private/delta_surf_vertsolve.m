function [x,d0,s0,t0] = delta_surf_vertsolve(S, T, X, s0, t0, d0, tolx) %#codegen
%DELTA_SURF_VERTSOLVE  Helper function for delta_surf, solving non-linear
%                      root finding problem in each water column.

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

[nz, nx, ny] = size(S);
[~,Xj] = size(X);
MX = nz * double(Xj > 1);

if eos(34.5, 3, 1000) > 1
    sortdir = 'ascend'; % eosis in-situ density, increasing in 3rd dim
else
    sortdir = 'descend'; % eos is specific volume, decreasing in 3rd dim
end

% Count number of valid bottles per cast
BotK = squeeze(sum(isfinite(S),1));


% Loop over each cast
idx = 0;
x = nan(nx,ny);
for c = 1:nx*ny
    K = BotK(c);
    if K >= 2
        
        % Select cast
        Xc = X((idx+1:idx+K).');
        Sc = S(1:K,c);
        Tc = T(1:K,c);
        
        % Get started with the discrete version (and linear interpolation)
        Dc = sort(eos(Sc, Tc, Xc) - eos(s0, t0, Xc), 1, sortdir);
        x(c) = interp1qn(d0, Dc, Xc);
        
        % Now refine by solving the non-linear problem at each water column
        x(c) = bisectguess(@diff_fun, X(idx+1), X(idx+K), tolx, x(c), ...
            X((idx+1:idx+K).'), S(1:K,c), T(1:K,c), s0, t0, d0);
    end
    idx = idx + MX;
end

end

function out = diff_fun(x, X, S, T, s0, t0, d0)
[s,t] = interp_firstdim_twovars(x, X, S, T);
out = eos(s, t, x) - eos(s0, t0, x) - d0 ;
end