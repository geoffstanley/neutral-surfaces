function gsf = orthobaric_montgomery(s, t, p, P, A, Y, s0, t0, OPTS, varargin)
%ORTHOBARIC_MONTGOMERY  The Montgomery potential with a functional offset.
%
% Zhang and Hogg (1992) upgraded the Montgomery (1937) potential by
% subtracting a constant reference pressure, to minimize errors when the
% surface in question is not a specific volume anomaly surface, on which
% the Montgomery potential is exact. Stanley (2018) extends this further,
% subtracting a reference pressure that is a function of specific volume
% anomaly.
%
%
% gsf = orthobaric_montgomery(s, t, p, P, A, Y, s0, t0, OPTS)
% computes the orthobaric Montgomery potential gsf on a surface with
% practical / Absolute salinity, potential / Conservative temperature, and
% pressure given by s, t, and p respectively, and with s and t reference
% values given by s0 and t0, respectively, in an ocean with pressure P,
% specific volume A, and hydrostatic acceleration potential Y, and with
% equation of state for the specific volume given by eos.m in the path. The
% specific volume should be calculated by A = eos(S, T, P) where S is the
% practical / Absolute salinity and T is the potential / Conservative
% temperature at the data sites of P. The hydrostatic acceleration
% potential should be calculated as Y = hsap3(P, ATMP, ETAN, A, grav),
% where ATMP is the atmospheric pressure, ETAN is the sea-surface height,
% and grav the gravitational acceleration. If s0 or t0 are empty, they are
% both taken as the values of s and t, respectively, where p is a maximum;
% this is done separately for each connected region of the surface. Inputs
% s and t may instead be S and T, in which case s and t is determined by
% linear interpolation of S and T through P onto p. Algorithmic parameters
% are determined by OPTS. On each connected region of the surface, the
% pressure is fit as a spline function of the specific volume anomaly (the
% specific volume minus that using the reference values s0 and t0) of order
% OPTS.SPLINE_ORDER and having at most OPTS.SPLINE_PIECES number of pieces.
% To determine the connected regions, the surface is treated as periodic in
% its first dimension (longitude) if OPTS.WRAP(1) is true, and likewise
% treated periodic in its second dimension (latitude) if OPTS.WRAP(2) is
% true.
%
% gsf = orthobaric_montgomery(s, t, z, Z, R, Y, s0, t0, OPTS, grav, rho_c)
% computes the orthobaric Montgomery potential gsf on a surface with
% practical / Absolute salinity, potential / Conservative temperature, and
% depth given by s, t, and z respectively, and with s and t reference
% values given by s0 and t0, respectively, in a Boussinesq ocean with
% in-situ density R and hydrostatic acceleration potential Y given at
% depths Z, with gravitational acceleration grav and Boussinesq reference
% density rho_c, and with equation of state for the in-situ density given
% by eos.m in the path. The in-situ density should be calculated by R =
% eos(S, T, Z) where S is the practical / Absolute salinity and T is the
% potential / Conservative temperature at the depths Z. The hydrostatic
% acceleration potential should be calculated as Y = hsap3(Z, ATMP, ETAN,
% R, grav, rho_c), where ATMP is the atmospheric pressure and ETAN is the
% sea-surface height. If s0 or t0 are empty, they are both taken as the
% values of s and t, respectively, where z is a maximum; this is done
% separately for each connected region of the surface. Inputs s and t may
% instead be S and T, in which case s and t is determined by linear
% interpolation of S and T through Z onto z. Algorithmic parameters are
% determined by OPTS. On each connected region of the surface, the depth is
% fit as a spline function of the in-situ density anomaly (the in-situ
% density minus that using the reference values s0 and t0) of order
% OPTS.SPLINE_ORDER and having at most OPTS.SPLINE_PIECES number of pieces.
% To determine the connected regions, the surface is treated as periodic in
% its first dimension (longitude) if OPTS.WRAP(1) is true, and likewise
% treated periodic in its second dimension (latitude) if OPTS.WRAP(2) is
% true. Note z and Z must be positive and increase downwards.
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
% OPTS [struct]: algorithmic parameters.
%   OPTS.WRAP [2 element array]: periodicity of the surface in each horizontal dimension
%   OPTS.SPLINE_ORDER [1, 1]: order of the spline (default: 4 for cubic splines)
%   OPTS.SPLINE_PIECES [1, 1]: maximum number of pieces of the spline (default: 12)
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
% Stanley, G. J.  The exact geostrophic stream function for neutral
% surfaces. Ocean Modelling, submitted.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


db2Pa = 1e4;  % Conversion from [dbar] to [Pa]

% Input checking
narginchk(9,11);

BOUSSINESQ = length(varargin) == 2; % grav and rho_c provided

[ni,nj] = size(p);
nk = size(P,1);
is1D = @(F) ismatrix(F) && all(size(F) == [nk,1]);
is2D = @(F) ismatrix(F) && all(size(F) == [ni,nj]);
is3D = @(F) ndims(F) == 3 && all(size(F) == [nk,ni,nj]);
is4D = @(F) ndims(F) == 4 && size(F,2) == nk-1 && size(F,3) == ni && size(F,4) == nj;
lead1 = @(p) reshape(p, [1 size(p)]);

assert((is2D(s) || is4D(s)), 'S must be [ni,nj] or [O,nk-1,ni,nj]');
assert(all(size(s) == size(t)), 'S and T must be the same size');
assert(is1D(P) || is3D(P), 'P must be [nk,1] or [nk,ni,nj]');
assert(is2D(p), 'p must be [ni,nj]');

assert(isstruct(OPTS), 'OPTS must be a struct')
assert(isfield(OPTS, 'WRAP') && isvector(OPTS.WRAP) && length(OPTS.WRAP) == 2, 'OPTS.WRAP must be a vector of length 2');
if isfield(OPTS, 'SPLINE_ORDER')
    assert(isscalar(OPTS.SPLINE_ORDER), 'OPTS.SPLINE_ORDER must be a scalar');
else
    OPTS.SPLINE_ORDER = 4; % Cubic
end
if isfield(OPTS, 'SPLINE_PIECES')
    assert(isscalar(OPTS.SPLINE_PIECES), 'OPTS.SPLINE_PIECES must be a scalar');
else
    OPTS.SPLINE_PIECES = 12; % Maximum number of knots, linearly spaced between domain limits
end

if BOUSSINESQ
    grav = varargin{1};
    rho_c = varargin{2};
    fac = -grav / rho_c; % (-) because p is z we have been z>0 people!
else
    fac = db2Pa;
end

if is4D(s) % Evaluate interpolants for S and T onto the surface
    [s,t] = ppc_val2(P, s, t, lead1(p));
end

% Find the connected regions
wet = isfinite(p);
neigh = grid_adjacency([ni,nj], 4, OPTS.WRAP);
[qu,qts,ncc] = bfs_conncomp(wet, neigh);

% The number of data points in the largest connected region
maxpix = max(diff(qts));


% Use s0 and t0 if provided, otherwise use the s and t values found where p
% is a maximum (deepest part of the surface) in each connected region.
if isscalar(s0) && isscalar(t0)
    [y0, a0] = hsap2(s0, t0, p, varargin{:});
else
    y0 = nan(ni,nj);
    a0 = nan(ni,nj);
    for l = 1 : ncc
        reg = qu(qts(l) : qts(l+1)-1);
        p_ = p(reg);
        s_ = s(reg);
        t_ = t(reg);
        [~,I] = max(p_);
        s0 = s_(I);
        t0 = t_(I);
        [y0(reg), a0(reg)] = hsap2(s0, t0, p_, varargin{:});
    end
end

% Calculate easy term involving an integral of hydrostatic balance
[y , a ] = hsap2(s, t, p, P, A, Y, varargin{:});

% Calculate the specific volume anomaly or the in-situ density anomaly
d = a - a0;

% Initialize the geostrophic stream function
gsf = y - y0 + fac * d .* p;

% For each connected region, add an integral of p as a function of d to the
% geostrophic stream function
for l = 1 : ncc
    
    reg = qu(qts(l) : qts(l+1)-1);
    d_ = d(reg);
    p_ = p(reg);
    
    % Number of spline pieces, by linearly scaling between 1 piece for
    % an empty region and OPTS.SPLINE_PIECES pieces for the largest region.
    pieces = ceil(OPTS.SPLINE_PIECES * length(reg) / maxpix);
    
    % Select the breakpoints linearly through the p data
    breaks = linspace(min(d_), max(d_), pieces+1);
    
    % Fit pressure (or depth) as a function of delta on the surface
    if breaks(1) == breaks(end)
        % All data has same p value. Simply choose the mean of the y data.
        pfnd = mkpp(breaks([1 end]), mean(p_));
    else
        % Fit a spline
        pfnd = splinefit(d_, p_, breaks, OPTS.SPLINE_ORDER);
    end
    
    % Finalize the geostrophic stream function
    gsf(reg) = gsf(reg) - fac * ppint(pfnd, pfnd.breaks(1), d_);
    
end


end