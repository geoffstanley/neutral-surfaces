function gsf = mcdougallklocker10(s, t, x, X, M, Y, s0, t0, x0, varargin)
%MCDOUGALLKLOCKER10  The McDougall and Klocker (2010) geostrophic stream function.
%
%
% gsf = mcdougallklocker10(s, t, p, P, A, Y, s0, t0, p0)
% computes the McDougall and Klocker (2010) geostrophic stream function gsf
% on a surface with practical / Absolute salinity, potential / Conservative
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
% and t may instead be S and T, in which case s and t is determined by
% linear interpolation of S and T through P onto p.
% 
% gsf = mcdougallklocker10(s, t, z, Z, R, Y, s0, t0, z0, grav, rho_c)
% computes the McDougall and Klocker (2010) geostrophic stream function gsf
% on a surface with practical / Absolute salinity, potential / Conservative
% temperature, and depth given by s, t, and z respectively, and with s, t,
% and z reference values given by s0, t0, and z0 respectively, in a
% Boussinesq ocean with in-situ density R and hydrostatic acceleration
% potential Y given at depths Z, with gravitational acceleration grav and
% Boussinesq reference density rho_c, and with equation of state for the
% in-situ density given by eos.m in the path. The in-situ density should be
% calculated by R = eos(S, T, Z * 1e-4 * grav * rho_c) where S is the
% practical / Absolute salinity and T is the potential / Conservative
% temperature at the depths Z. The hydrostatic acceleration potential
% should be calculated as Y = hsap3(Z, ATMP, ETAN, R, grav, rho_c), where
% ATMP is the atmospheric pressure and ETAN is the sea-surface height. If
% s0, t0, or z0 are empty, they are taken as the mean of s, t, or z,
% respectively. Inputs s and t may instead be S and T, in which case s and
% t is determined by linear interpolation of S and T through Z onto z. Note
% z and Z must be positive and increase downwards.
%
%
% --- Input:
% s [nx, ny] or [nz, nx, ny]: the practical / Absolute salinity, 
%  on the surface, or at all data sites
% t [nx, ny] or [nz, nx, ny]: the potential / Conservative temperature
%  on the surface, or at all data sites
% p [nx, ny]: pressure on surface [dbar]
% z [nx, ny]: depth on surface [m, positive]
% P [nz, nx, ny] or [nz, 1]: pressure at all data sites [dbar]
% Z [nz, nx, ny] or [nz, 1]: depth at all data sites [m, positive]
% A [nz, nx, ny]: specific volume [m^3 kg^-1]
% R [nz, nx, ny]: in-situ density [kg m^-3]
% Y [nz, nx, ny]: acceleration potential from hydrostatic balance [m^2 s^-2]
% s0 [1, 1] or []: reference s value
% t0 [1, 1] or []: reference t value
% p0 [1, 1] or []: reference p value
% z0 [1, 1] or []: reference z value
% grav [1, 1]: gravitational acceleration [m s^-2]
% rho_c [1, 1]: Boussinesq reference density [kg m^-3]
%
% Note: nz is the maximum number of data points per water column,
%       nx is the number of data points in longitude,
%       ny is the number of data points in latitude.
%
% Note: physical units for s, t, s0, and t0 are determined by eos.m.
%
%
% --- Output:
% gsf [nx, ny]: geostrophic stream function [m^2 s^-2]
%
%
% --- Requirements:
% hsap2
% interp1qn2 - https://www.mathworks.com/matlabcentral/fileexchange/69713
%
%
% --- References:
% McDougall, T. J. & Klocker, A. An approximate geostrophic streamfunction
% for use in density surfaces. Ocean Modelling 32, 105?117 (2010).

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

Thermobaricity_div_rho = 2.7e-15; % K^-1 Pa^-2 m^2 s^-2

% Input checking
narginchk(9,11);
BOUSSINESQ = (nargin == 11);

db2Pa = 1e4;  % Conversion from [dbar] to [Pa]

[nx,ny] = size(x);
nz = size(X,1);
is3D = @(F) ndims(F) == 3 && all(size(F) == [nz,nx,ny]);
lead1 = @(x) reshape(x, [1 size(x)]);

if BOUSSINESQ
    grav = varargin{1};
    rho_c = varargin{2};
    fac = -grav / rho_c; % (-) because x is z we have been z>0 people!
    fac2 = (grav * rho_c).^2;
else
    fac = db2Pa;
    fac2 = fac^2;
end

if is3D(s) % Interpolate 3D s and t onto the surface
    try
        [s,t] = interp1qn2_mex(lead1(x), X, s, t);
    catch
        [s,t] = interp1qn2(lead1(x), X, s, t);
    end
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
gsf = y - y0 + (0.5 * fac) * d .* (x - x0) ... 
    - fac2 * Thermobaricity_div_rho / 12 .* (t - t0) .* (x - x0).^2 ;


% The non-Boussinesq form for the MK geostrophic stream function is
% ZH - (1/2) (p - p0) delta + \int_P^p - (1/12) (Tb/r) (t - t0) (p - p0)^2
% where ZH is the Zhang and Hogg (1992) geostrophic stream function.
% To transform MK to Boussinesq, first note how Montgomery transforms:
% Non-Boussinesq:           p delta  +  g z        + integral... where delta = spec vol anom
%     Boussinesq: (g/rho_c) z delta  +  p / rho_c  + integral... where delta = in-situ density anom
% So, p -> rho_c g z  as expected,
% but since delta also changes,  p delta -> (g/rho_c) z delta
% Now look at (55) in McDougall and Klocker (2010), which reads
%           (p - p_0) nabla delta  =(approx)= (T_b/r) (t - t_0)     *       (p - p_0) nabla (p - p_0)
% Transforming that to Boussinesq using the above rules, we get
%  (g/rho_c) (z - z0) nabla delta  =(approx)= (T_b/r) (t - t_0) (g rho_c)^2 (z - z_0) nabla (z - z_0)
% That tells us how the (1/12) term in the MK stream function transforms, namely
%  (1/12) (T_b/r) (g rho_c)^2 (t - t_0) (z - z_0)^2.