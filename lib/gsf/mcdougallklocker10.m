function gsf = mcdougallklocker10(s, t, p, P, A, Y, s0, t0, p0, varargin)
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
% calculated by R = eos(S, T, Z) where S is the practical / Absolute
% salinity and T is the potential / Conservative temperature at the depths
% Z. The hydrostatic acceleration potential should be calculated as Y =
% hsap3(Z, ATMP, ETAN, R, grav, rho_c), where ATMP is the atmospheric
% pressure and ETAN is the sea-surface height. If s0, t0, or z0 are empty,
% they are taken as the mean of s, t, or z, respectively. Inputs s and t
% may instead be S and T, in which case s and t is determined by linear
% interpolation of S and T through Z onto z. Note z and Z must be positive
% and increase downwards.
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
% McDougall, T. J. & Klocker, A. An approximate geostrophic streamfunction
% for use in density surfaces. Ocean Modelling 32, 105?117 (2010).

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


Thermobaricity_div_rho = 2.7e-15; % K^-1 Pa^-2 m^2 s^-2

% Input checking
narginchk(9,11);

BOUSSINESQ = length(varargin) == 2; % grav and rho_c provided

db2Pa = 1e4;  % Conversion from [dbar] to [Pa]

[ni,nj] = size(p);
nk = size(P,1);
is4D = @(F) ndims(F) == 4 && size(F,2) == nk-1 && size(F,3) == ni && size(F,4) == nj;
lead1 = @(p) reshape(p, [1 size(p)]);
nanmean = @(x) mean(x, 'omitnan');

if BOUSSINESQ
    grav = varargin{1};
    rho_c = varargin{2};
    fac = -grav / rho_c; % (-) because p is z we have been z>0 people!
    fac2 = (grav * rho_c).^2;
else
    fac = db2Pa;
    fac2 = fac^2;
end

if is4D(s) % Evaluate interpolants for S and T onto the surface
    [s,t] = ppc_val2(P, s, t, lead1(p));
end

% If s0, t0, p0 were not provided, use the mean
if isempty(s0)
    s0 = nanmean(s(:));
end
if isempty(t0)
    t0 = nanmean(t(:));
end
if isempty(p0)
    p0 = nanmean(p(:));
end

% Calculate the easy terms involving integrals of hydrostatic balance:
[y , a ] = hsap2(s , t , p, P, A, Y, varargin{:});
[y0, a0] = hsap2(s0, t0, p,          varargin{:});

% Calculate the specific volume anomaly or the in-situ density anomaly
d = a - a0;

% Calculate the geostrophic stream function
gsf = y - y0 + (0.5 * fac) * d .* (p - p0) ...
    - fac2 * Thermobaricity_div_rho / 12 .* (t - t0) .* (p - p0).^2 ;


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