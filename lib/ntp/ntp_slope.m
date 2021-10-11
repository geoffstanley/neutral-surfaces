function [sx, sy] = ntp_slope(Sppc, Tppc, Z, z, tolz, dx, dy) %#codegen
%NTP_SLOPE  Find the slope of the Neutral Tangent Plane at all points along a surface
%
%
% [sx,sy] = ntp_slope(Sppc, Tppc, Z, z, tolz, dx, dy)
% finds the slopes in the two horizontal directions, sx and sy, of the
% neutral tangent plane at all points on a surface.  Essentially, this runs
% ntp_midpoint_to_casts at each point on a surface in the two directions,
% so see ntp_midpoint_to_casts for further documentation.
%
%
% --- Input:
% Sppc [O, K-1, ni, nj]: coefficients for piecewise polynomial for
%                        practical / Absolute Salinity in terms of Z
% Tppc [O, K-1, ni, nj]: coefficients for piecewise polynomial for
%                        potential / Conservative Temperature in terms of Z
% Z [K, ni, nj] or [K, 1]: pressure or depth in water column
% z [ni, nj]: pressure or depth of a surface
% tolz [1, 1]: tolerance for solving the level of neutral buoyancy (same
%              units as Z)
% dx [ni, nj]:  distance between grid points (i,j) and (i-1,j)
% dy [ni, nj]:  distance between  grid points (i,j) and (i,j-1)
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
% sx [ni, nj]: slope of the neutral tangent plane between (i,j) and (i-1,j)
% sy [ni, nj]: slope of the neutral tangent plane between (i,j) and (i,j-1)
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
% ntp_slope_error
% ppc_linterp, ppc_pchip

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% --- Input checks, set parameters
narginchk(4, 7);

if nargin < 5 || isempty(tolz)
  tolz = 1e-5;
end
if nargin < 6 || isempty(dx)
  dx = 1;
end
if nargin < 7 || isempty(dy)
  dy = 1;
end

[ni,nj] = size(z);

% --- Solve nonlinear problem for NTP depth difference between each pair of adjacent casts
sx = nan(ni, nj);
sy = nan(ni, nj);

Zmat = ~isvector(Z);
Zm = Z(:,1);
Zn = Z(:,1);

% K gives the number of valid bottles in each water column.
% Note that K >= 1.  Even where all bottles are invalid, i.e. land, K = 1.
K = squeeze(sum(isfinite(Sppc(1,:,:,:)))) + 1; 

% Loop over each water column.  
% Note this calculation assumes a doubly-periodic domain
% The NTP slope is grad_n z, where z < 0
% Negative sign added to ntp_midpoint_to_casts results because we've been z > 0 people.
for j = 1 : nj
  jm1 = mod(j-2, nj) + 1;
  for i = 1 : ni
    im1 = mod(i-2, ni) + 1;
    
    zm = z(i,j);
    km = K(i,j);
    if Zmat
      Zm = Z(:,i,j);
    end
    
    % --- NTP with neighbour in i dimension
    zn = z(im1,j);
    kn = K(im1,j);
    if Zmat
      Zn = Z(:,im1,j);
    end
    sx(i,j) = -ntp_midpoint_to_casts(Sppc(:,:,i,j), Tppc(:,:,i,j), Zm, km, Sppc(:,:,im1,j), Tppc(:,:,im1,j), Zn, kn, zm, zn, tolz);
    
    % --- NTP with neighbour in j dimension
    zn = z(i, jm1);
    kn = K(i,jm1);
    if Zmat
      Zn = Z(:,i,jm1);
    end
    sy(i,j) = -ntp_midpoint_to_casts(Sppc(:,:,i,j), Tppc(:,:,i,j), Zm, km, Sppc(:,:,i,jm1), Tppc(:,:,i,jm1), Zn, kn, zm, zn, tolz);
    
  end % i
end   % j

% Divide by horizontal distances
if nargin == 7 && ~isempty(dx) && ~isempty(dy) ...
    && ~all(dx(:) == 1) ...
    && ~all(dy(:) == 1)
  sx = sx ./ dx;
  sy = sy ./ dy;
end

end
