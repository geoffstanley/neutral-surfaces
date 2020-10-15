function Y = hsap3(X, ATMP, ETAN, M, grav, varargin)
%HSAP3  Integrate hydrostatic balance to obtain acceleration potential at every data point.
%
%
% Y = hsap3(P, ATMP, ETAN, A, grav)
% integrates hydrostatic balance to obtain the acceleration potential Y --
% the depth times the gravitational acceleration -- at the data sites of A
% and P, in an ocean with atmospheric pressure ATMP, sea-surface height
% ETAN, specific volume A at the data sites of pressure P, and
% gravitational acceleration grav.
%
% Y = hsap3(Z, ATMP, ETAN, R, grav, rho_c)
% integrates hydrostatic balance to obtain the acceleration potential Y --
% the pressure divided by the Boussinesq reference density -- at the data
% sites of R and Z, in a Boussinesq ocean with atmospheric pressure ATMP,
% sea-surface height ETAN, in-situ density R at the data sites of depth Z,
% with gravitational acceleration grav, and Boussinesq reference density
% rho_c.
%
%
% --- Input:
% ATMP [ni, nj]: atmospheric pressure loading [dbar]
% ETAN [ni, nj]: sea-surface height [m]
% A [nk, ni, nj]: specific volume [m^3 kg^-1]
% R [nk, ni, nj]: in-situ density [kg m^-3]
% Z [nk, ni, nj] or [nk, 1]: depth [m, positive]
% P [nk, ni, nj] or [nk, 1]: pressure [dbar]
% grav [1, 1]: gravitational acceleration [m s^-2]
% rho_c [1, 1]: Boussinesq reference density [kg m^-3]
%
%
% --- Output:
% Y [nk, ni, nj]: acceleration potential from hydrostatic balance [m^2 s^-2]

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


narginchk(5,6);

db2Pa = 1e4;
lead1 = @(x) reshape(x, [1 size(x)]);

if nargin == 6
    % Boussinesq form
    rho_c = varargin{1};
    if isvector(X)
        Y = (db2Pa/rho_c) * lead1(ATMP) + grav * lead1(ETAN) ...
            + cumsum((grav/(2*rho_c)) * diff([0; X]) .* (M([1, 1:end-1],:,:) + M), 1);
    else
        Y = (db2Pa/rho_c) * lead1(ATMP) + grav * lead1(ETAN) ...
            + cumsum((grav/(2*rho_c)) * cat(1, X(1,:,:), diff(X,[],1)) .* (M([1, 1:end-1],:,:) + M), 1);
    end
    
else
    % Non-Boussinesq form
    Y = + grav * lead1(ETAN) ...
        - (db2Pa / 2) * cumsum(cat(1, X(1,:,:) - lead1(ATMP), diff(X,[],1)) .* (M([1 1:end-1],:,:) + M), 1);
    
end
