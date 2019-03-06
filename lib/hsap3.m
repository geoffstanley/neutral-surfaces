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
% ATMP [nx, ny]: atmospheric pressure loading [dbar]
% ETAN [nx, ny]: sea-surface height [m]
% A [nz, nx, ny]: specific volume [m^3 kg^-1]
% R [nz, nx, ny]: in-situ density [kg m^-3]
% Z [nz, nx, ny] or [nz, 1]: depth [m, positive]
% P [nz, nx, ny] or [nz, 1]: pressure [dbar]
% grav [1, 1]: gravitational acceleration [m s^-2]
% rho_c [1, 1]: Boussinesq reference density [kg m^-3]
%
%
% --- Output:
% Y [nz, nx, ny]: acceleration potential from hydrostatic balance [m^2 s^-2]

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


narginchk(5,6);

db2Pa = 1e4;
lead1 = @(x) reshape(x, [1 size(x)]);

if nargin == 6
    % Boussinesq form
    rho_c = varargin{1};
    if isvector(X)
        Y = (db2Pa/rho_c) * ATMP + grav * lead1(ETAN) ...
            + cumsum((grav/(2*rho_c)) * diff([0; X]) .* (M([1, 1:end-1],:,:) + M), 1);
    else
        Y = (db2Pa/rho_c) * ATMP + grav * lead1(ETAN) ...
            + cumsum((grav/(2*rho_c)) * cat(1, X(1,:,:), diff(X,[],1)) .* (M([1, 1:end-1],:,:) + M), 1);
    end
    
else
    % Non-Boussinesq form
    Y = + grav * lead1(ETAN) ...
        - (db2Pa / 2) * cumsum(cat(1, X(1,:,:) - ATMP, diff(X,[],1)) .* (M([1 1:end-1],:,:) + M), 1);
    
end
