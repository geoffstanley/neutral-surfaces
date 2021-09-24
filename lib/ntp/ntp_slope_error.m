function [sx, sy] = ntp_slope_error(Sppc, Tppc, Z, z, tolz, dx, dy, wrap) %#codegen
%NTP_SLOPE_ERROR  Find the slope error between a given surface and the Neutral Tangent Plane
%
%
% [sx,sy] = ntp_slope_error(Sppc, Tppc, Z, z, tolz, dx, dy, wrap)
% finds the differences, sx and sy, between the slope of the surface of
% depth z, and the slope of the local neutral tangent plane.  Essentially,
% this runs ntp_slope and then subtracts the slope of the given surface.
%
%
% --- Input:
% Sppc [O, K-1, ni, nj]: coefficients for piecewise polynomial for practical 
%     / Absolute Salinity in terms of Z on the cast.  Here, O is the
%     polynomial order and K is the number of data points (including NaN's)
%     on the casts.
% Tppc [O, K-1, ni, nj]: as above but for potential / Conservative Temperature.
% Z [K, ni, nj] or [K, 1]: pressure or depth data points on casts
% z [ni, nj]: pressure or depth of a surface
% tolz [1, 1]: tolerance for solving the level of neutral buoyancy (same
%              units as Z)
% dx [ni, nj]:  distance between grid points (i,j) and (i-1,j)
% dy [ni, nj]:  distance between  grid points (i,j) and (i,j-1)
% wrap [2, 1]:  2 element vector specifying which dimensions are periodic
%
% Note: physical units for Sppc, Tppc, Z, z are determined by eos.m.  See
%       "Non-Boussinesq case" below.
%
% Note: Z must increase monotonically along its first dimension.
%
% Note: dx and dy can have any dimension replaced by a singleton dimension.
%
%
% --- Output:
% sx [ni, nj]: difference in slope of the neutral tangent plane centred at
%             (i-1/2,j) and the surface between casts (i,j) and (i-1,j)
% sy [ni, nj]: difference in slope of the neutral tangent plane centred at
%             (i,j-1/2) and the surface between casts (i,j) and (i,j-1)
%
%
% --- Non-Boussinesq case:
% in the non-Boussinesq case where Z and z are actually pressure then the
% units of sx and sy are [dbar / m].  To convert to [dbar / dbar] = [1],
% divide sx and sy by dp/dz which can be obtained from hydrostatic balance.
% Specifically, do the following.
%               Pa2db = 1e-4;
%               grav = 9.81; % [m s-2]  adjust as needed
%               [s,t] = ppc_val2(Z, Sppc, Tppc, z);
%               if eos(34.5, 3, 1000) > 1  % eos gives density
%                 dpdz = Pa2db * grav * eos(s,t,z);
%               else  % eos gives specific volume
%                 dpdz = Pa2db * grav ./ eos(s,t,z);
%               end
%               sx = sx * 2 ./ (dpdz + circshift(dpdz, [+1, 0]));
%               sy = sy * 2 ./ (dpdz + circshift(dpdz, [0, +1]));
% A better approach still would be to know pressure as a function of depth.
% This is not provided for here.
%
%
% --- See Also:
% ntp_midpoint_to_casts
% ntp_slope
% ppc_linterp, ppc_pchip

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% Input checks, set parameters
narginchk(4, 9);

if nargin < 5|| isempty(tolz)
  tolz = 1e-5;
end
if nargin < 6 || isempty(dx)
  dx = 1;
end
if nargin < 7 || isempty(dy)
  dy = 1;
end
if nargin < 8 || isempty(wrap)
  wrap = [true, true]; % doubly periodic
end

[ni,nj] = size(z);
nk = size(Z,1);

is2Ds = @(F) ismatrix(F) && ismember(size(F,1), [1,ni]) && ismember(size(F,2), [1,nj]);
is3Ds = @(F) ismember(size(F,2), [1,ni]) && ismember(size(F,3), [1,nj]);
is4D = @(F) ndims(F) == 4 && size(F,2) == nk-1 && size(F,3) == ni && size(F,4) == nj;
assert(is4D(Sppc), 'the dimensions of Sppc must be [?, nk-1, ni, nj] where nk = size(Z,1) and [ni, nj] = size(z)')
assert(is4D(Tppc), 'the dimensions of Tppc must be [?, nk-1, ni, nj] where nk = size(Z,1) and [ni, nj] = size(z)')
assert(is3Ds(Z)  , 'The last two dimensions of Z must match those of s (or be singletons).')
assert(is2Ds(dx) , 'the dimensions of dx must match those of z (or be singletons)');
assert(is2Ds(dy) , 'the dimensions of dy must match those of z (or be singletons)');
assert(isscalar(tolz), 'tolz must be a scalar');
assert(numel(wrap) == 2, 'wrap must be a 2 element (logical) vector');

% Just-in-time code generation:
ntp_slope_codegen(nk, ni, nj, isvector(Z));

% Solve nonlinear problem for NTP depth difference between each pair of adjacent casts
[dzi, dzj] = ntp_slope_mex(Sppc, Tppc, Z, z, tolz, 1, 1);

% slope error = grad_n z - grad_a z,  where z < 0.
% Since we've been z > 0 people, dzi and dzj got a negative sign added from
% ntp_slope, so here we add rather than subtract the grad_a z term.
im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
sx = (dzi + (z - im1(z)) ) ./ dx;
sy = (dzj + (z - jm1(z)) ) ./ dy;

% Handle non-periodic dimensions (ntp_slope assumed doubly-periodic).
if ~wrap(1)
  sx(1,:) = nan;
end
if ~wrap(2)
  sy(:,1) = nan;
end
