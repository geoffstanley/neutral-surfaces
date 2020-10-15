function gsf = zhanghogg92(s, t, x, X, M, Y, s0, t0, x0, varargin)
%ZHANGHOGG92  The Zhang and Hogg (1992) geostrophic stream function.
%
%
% gsf = zhanghogg92(s, t, p, P, A, Y, s0, t0, p0)
% computes the Zhang and Hogg (1992) geostrophic stream function gsf on a
% surface with practical / Absolute salinity, potential / Conservative
% temperature, and pressure given by s, t, and p respectively, and with s,
% t, and p reference values given by s0, t0, and p0 respectively, in an
% ocean with pressure P, specific volume A, and hydrostatic acceleration
% potential Y, and with equation of state for the specific volume given by
% eos.m in the path. The specific volume should be calculated by A = eos(S,
% T, P) where S is the practical / Absolute salinity and T is the potential
% / Conservative temperature at the data sites of P. The hydrostatic
% acceleration potential should be calculated as Y = hsap3(P, ATMP, ETAN,
% A, grav), where ATMP is the atmospheric pressure, ETAN is the sea-surface
% height, and grav the gravitational acceleration. If s0, t0, or p0 are
% empty, they are taken as the mean of s, t, or p, respectively. Inputs s
% and t may instead be piecewise polynomial interpolants for S and T in
% terms of P in each water column, in which case s and t are determined by
% evaluating these interpolants onto p.
%
% gsf = zhanghogg92(s, t, z, Z, R, Y, s0, t0, z0, grav, rho_c)
% computes the Zhang and Hogg (1992) geostrophic stream function gsf on a
% surface with practical / Absolute salinity, potential / Conservative
% temperature, and depth given by s, t, and z respectively, and with s, t,
% and z reference values given by s0, t0, and z0 respectively, in a
% Boussinesq ocean with in-situ density R and hydrostatic acceleration
% potential Y given at depths Z, with gravitational acceleration grav and
% Boussinesq reference density rho_c, and with equation of state for the
% in-situ density given by eos.m in the path. The in-situ density should be
% calculated by R = eos(S, T, Z) where S is the practical / Absolute
% salinity and T is the potential / Conservative temperature at the depths
% Z. The hydrostatic acceleration potential should be calculated as Y =
% hsap3(Z, ATMP, ETAN, R, grav, rho_c), where ATMP is the atmospheric
% pressure and ETAN is the sea-surface height. If s0, t0, or z0 are empty,
% they are taken as the mean of s, t, or z, respectively. Inputs s
% and t may instead be piecewise polynomial interpolants for S and T in
% terms of Z in each water column, in which case s and t are determined by
% evaluating these interpolants onto z.
%
%
% --- Input:
% s [ni, nj] or [O, nk-1, ni, nj]: the practical / Absolute salinity,
%  on the surface, or a piecewise polynomial interpolant in each water column
% t [ni, nj] or [O, nk-1, ni, nj]: the potential / Conservative temperature
%  on the surface, or a piecewise polynomial interpolant in each water column
% p [ni, nj]: pressure on surface [dbar]
% z [ni, nj]: depth on surface [m, positive]
% P [nk, ni, nj] or [nk, 1]: pressure at all data sites [dbar]
% Z [nk, ni, nj] or [nk, 1]: depth at all data sites [m, positive]
% A [nk, ni, nj]: specific volume [m^3 kg^-1]
% R [nk, ni, nj]: in-situ density [kg m^-3]
% Y [nk, ni, nj]: acceleration potential from hydrostatic balance [m^2 s^-2]
% s0 [1, 1] or []: reference s value
% t0 [1, 1] or []: reference t value
% p0 [1, 1] or []: reference p value
% z0 [1, 1] or []: reference z value
% grav [1, 1]: gravitational acceleration [m s^-2]
% rho_c [1, 1]: Boussinesq reference density [kg m^-3]
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: physical units for s, t, s0, and t0 are determined by eos.m.
%
%
% --- Output:
% gsf [ni, nj]: geostrophic stream function [m^2 s^-2]
%
%
% --- References:
% Zhang, H.-M. & Hogg, N. G. Circulation and water mass balance in the
% Brazil Basin. Journal of marine research 50, 385?420 (1992).

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


% Input checking
narginchk(9,11);

BOUSSINESQ = length(varargin) == 2; % grav and rho_c provided

db2Pa = 1e4;  % Conversion from [dbar] to [Pa]

[ni,nj] = size(x);
nk = size(X,1);
is4D = @(F) ndims(F) == 4 && size(F,2) == nk-1 && size(F,3) == ni && size(F,4) == nj;
lead1 = @(x) reshape(x, [1 size(x)]);

if BOUSSINESQ
    grav = varargin{1};
    rho_c = varargin{2};
    fac = -grav / rho_c; % (-) because x is z we have been z>0 people!
else
    fac = db2Pa;
end

if is4D(s) % Evaluate interpolants for S and T onto the surface
    [s,t] = ppc_val2(X, s, t, lead1(x));
end

% If s0, t0, x0 were not provided, use the mean
if isempty(s0)
    s0 = nanmean(s(:));
end
if isempty(t0)
    t0 = nanmean(t(:));
end
if isempty(x0)
    x0 = nanmean(x(:));
end

% Calculate the easy terms involving integrals of hydrostatic balance:
[y , m ] = hsap2(s , t , x, X, M, Y, varargin{:});
[y0, m0] = hsap2(s0, t0, x,          varargin{:});

% Calculate the specific volume anomaly or the in-situ density anomaly
d = m - m0;

% Calculate the geostrophic stream function
gsf = y - y0 + fac * d .* (x - x0);

end
