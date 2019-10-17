function x = pot_dens_surf_vertsolv(S, T, X, xref, val, tol) %#codegen
%POT_DENS_SURF_VERTSOLV  Helper function for pot_dens_surf, solving 
%                        non-linear root finding problem in each water column

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

[nk, ni, nj] = size(S);
[~,Xj] = size(X);
NX = nk * double(Xj > 1);

if eos(35,0,0) > 1
    sortdir = 'ascend'; % eos is in-situ density, increasing with 3rd argument
else
    sortdir = 'descend'; % eos is specific volume, decreasing with 3rd argument
end

% Count number of valid bottles per cast
BotK = squeeze(sum(isfinite(S),1));

% Loop over each cast
x = nan(ni,nj);
idx = 0;
for c = 1:ni*nj
    K = BotK(c);
    if K >= 2
        
        % Select cast
        Xc = X((idx+1:idx+K).');
        Sc = S(1:K,c);
        Tc = T(1:K,c);
        
        % Get started with the discrete version (and linear interpolation)
        Dc = sort(eos(Sc, Tc, xref), 1, sortdir);
        x(c) = interp1qn(val, Dc, Xc);
        
        % Now refine by solving the non-linear problem at each water column
        x(c) = bisectguess(@diff_fun, X(idx+1), X(idx+K), tol, x(c), ...
            X((idx+1:idx+K).'), S(1:K,c), T(1:K,c), xref, val);
    end
    idx = idx + NX;
end


end


function out = diff_fun(x, X, S, T, xref, val)
[s,t] = interp_firstdim_twovars(x, X, S, T);
out = eos(s, t, xref) - val;
end