function [ex,ey,sx,sy] = ntp_epsilon_r_x(s, t, x, dx, dy, centre, wrap, grav, S, T, X)
%NTP_EPSILON_R_X  The epsilon neutrality error and slope error from the
%                 neutral tangent plane using density and pressure (or depth)
%
%
% [ex,ey] = ntp_epsilon_r_x(s,t,x,dx,dy,centre,wrap,grav)
% computes the zonal and meridional neutrality errors, ex and ey, for an
% approximately neutral surface on which the practical / Absolute salinity
% is s, the potential / Conservative temperature is t, and the pressure or
% depth is x. The zonal and meridional distances between neighbouring grid
% points are dx and dy. Results are calculated on the grid points of s, t
% and x if centre is true; otherwise, zonal and meridional errors are
% located zonally and meridionally, respectively, midway between s, t, and
% x grid points. Data is treated periodic in the i'th (i=1 or 2) dimension
% of x if and only if wrap(i) is true.  The equation of state is given by
% eos.m in the path which must accept s, t, and x as inputs.
%
% [ex,ey,sx,sy] = ntp_epsilon_r_x(s,t,x,dx,dy,centre,wrap,grav,S,T,X)
% also computes the zonal and meridional slope errors, sx and sy, for an
% ocean of practical / Absolute salinity S and potential / Conservative
% temperature T at datasites where the pressure or depth is X.
%
% For a non-Boussinesq ocean, X and x are pressure [dbar].  When sx, sy are
% requested, the gravitational acceleration grav must be given.
%
% For a Boussinesq ocean, X and x are depth [m, positive].  When sx, sy are
% requested, grav must be the empty vector [].
%
%
% --- Input:
% s     [ni, nj]: practical / Absolute salinity on the surface
% t     [ni, nj]: potential / Conservative temperature on the surface
% x     [ni, nj]: pressure [dbar] or depth [m, positive] of the surface
% dx    [ni, nj]: Zonal      grid distance [m] between S(:,i,j) and S(:,i-1,j)
% dy    [ni, nj]: Meridional grid distance [m] between S(:,i,j) and S(:,i,j-1)
% centre  [1, 1]: true to compute outputs on the central grid, false on the
%                 half grids [logical]
% wrap    [2, 1]: wrap(i) is true when the data is periodic in i'th
%                 dimension of x
% grav    [1, 1]: gravitational acceleration [m s^-2]
% S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% X [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: data can be arranged as latitude by longitude.  Simply swap dx and
%       dy, adjust wrap as needed, and note the outputs ex and ey will also
%       be swapped, as will be sx and sy.
%
% Note: physical units for s,t,x,S,T,X are determined by eos.m and eos_x.m.
%
% Note: X must be increasing along its first dimension.
%
% Note: dx and dy can have any dimension replaced by a singleton dimension.
%
%
% --- Output:
% ex [ni, nj]: zonal epsilon neutrality error [kg m^-4]
% ey [ni, nj]: meridional epsilon neutrality error [kg m^-4]
% sx [ni, nj]: zonal slope error [dimensionless]
% sy [ni, nj]: meridional slope error [dimensionless]
%
%
% --- Mathematical form of ex, ey, sx, sy:
% The neutrality error (usually denoted by epsilon) is here taken as
%   [ex,ey] = rs grad s + rt grad t
% where s and t are the values of S and T on the surface projected onto a
% sphere, and rs and rt are the partial derivatives of in-situ density with
% respect to S and T, respectively, taken on the surface. Note that grad s
% is the same as grad_a S evaluated on the surface, where grad_a is the
% projected non-orthogonal gradient. See discussion of the under-tilde
% operator in Stanley (2019).
%
% The above definition is equivalent to Stanley (2019) Eq. (28), but
% differs from the standard definition by a factor of r. The standard
% definition (McDougall and Jackett 1988, p. 162) is [ex,ey] / r.
%
% The slope error is related to the neutrality error by
%   [sx,sy] = g / (r N2) [ex,ey] = -1 / (d sigma / d z) [ex,ey]
% where g is the gravitational acceleration and N2 is the square of the
% vertical stratification, and sigma is the locally referenced potential
% density. This differs from the standard relation (see McDougall and
% Jackett 1988 Eq. (14), or Klocker and McDougall 2009, Eq. (10)) by a
% factor of r because of our definition for [ex,ey] above.
%
%
% --- References:
% Klocker, A., McDougall, T. J. & Jackett, D. R. A new method for forming
% approximately neutral surfaces. Ocean Science 5, 155?172 (2009).
%
% McDougall, T. J. & Jackett, D. R. On the helical nature of neutral
% trajectories in the ocean. Progress in Oceanography 20, 153?183 (1988).
%
% Stanley, G.J., 2019. Neutral surface topology. Ocean Modelling 138,
% 88–106. https://doi.org/10.1016/j.ocemod.2019.01.008

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
% Version   : 2.1.1
%
% Modified by : --
% Date        : --
% Changes     : --

% --- Input checks, set parameters
narginchk(7,11)
assert(~isempty(which('eos')), 'Cannot find eos.m in the path or current directory.');
assert(~isempty(which('eos_x')), 'Cannot find eos_x.m in the path or current directory.');
assert(nargout <= 2 || nargin == 11, 'If sx, sy are requested, must provide S, T, X');

Pa2db = 1e-4;
N_min = 2e-4;
N2rhogr_min = N_min^2 * 1000 / 10;
nonBOUSSINESQ = (nargin >= 8) && ~isempty(grav); % grav is given.


eos_is_specvol = (eos(34.5, 3, 1000) < 1);
lead1 = @(x) reshape(x, [1 size(x)]);


if centre
    % Get errors by centred differences on the T grid
    dx = dx * 2;
    dy = dy * 2;
    A_X = @(F) (im1_2d(F) + ip1_2d(F)) / 2;
    D_X = @(F)  ip1_2d(F) - im1_2d(F);
    A_Y = @(F) (jm1_2d(F) + jp1_2d(F)) / 2;
    D_Y = @(F)  jp1_2d(F) - jm1_2d(F);
else
    % Get errors by backward differences, and leave on the U, V grid.
    A_X = @(F) (F + im1_2d(F)) / 2;
    D_X = @(F)  F - im1_2d(F);
    A_Y = @(F) (F + jm1_2d(F)) / 2;
    D_Y = @(F)  F - jm1_2d(F);
end



r = eos(s, t, x);
if eos_is_specvol
    r = 1 ./ r;     % now r is density [kg m^3]
end

% Begin calculations for errors in x direction
sa = A_X(s);
ta = A_X(t);
xa = A_X(x);
rxa = eos_x(sa, ta, xa);
if eos_is_specvol
    % convert rxa into to in-situ density derivative w.r.t. x:
    rxa = -rxa ./ eos(sa, ta, xa).^2;
end

ex = (D_X(r) - rxa .* D_X(x)) ./ dx;
% Note, ex above is virtually identical to ex below (this equivalence is
% accurate to third order):
% >> [rs,rt] = eos_s_t(sa,ta,xa);
% >> ex = (rs .* D_X(s) + rt .* D_X(t)) ./ dx;
% In particular, note that for centred differences, this is (written with
% the original, unmodified dx):
% >> ex = (rs .* (ip1_2d(s) - im1_2d(s)) + rt .* (ip1_2d(t) - im1_2d(t))) ./ (2*dx);
% The above ex are NOT identical to the below form:
% >> rx = eos_x(s, t, x);
% >> ex = ((ip1_2d(r) - im1_2d(r)) - rx .* (ip1_2d(x) - im1_2d(x))) ./ (2*dx);
% Lesson: tempting though it is to use s,t and x in the local water column
% to calculate rx, the partial derivative of density w.r.t x, it's wrong!
if nargout < 2; return; end

% Begin calculations for errors in y direction.
sa = A_Y(s);
ta = A_Y(t);
xa = A_Y(x);
rxa = eos_x(sa, ta, xa);
if eos_is_specvol
    % convert rxa into to in-situ density derivative w.r.t. x:
    rxa = -rxa ./ eos(sa, ta, xa).^2;
end

ey = (D_Y(r) - rxa .* D_Y(x)) ./ dy;
if nargout < 3; return; end

% Start working on slope errors

% Calculate locally referenced potential density
sigma = eos(S, T, lead1(x));
if eos_is_specvol
    sigma = 1 ./ sigma;  % now sigma is density [kg m^3]
end

% Take derivative of sigma w.r.t pressure or depth (as a SPATIAL coordinate)
[~, N2rhogr] = pchipdqn(lead1(x), X, sigma);
if nonBOUSSINESQ
    % P is actually pressure, not depth analogue.
    % Convert d(sigma)/dp into d(sigma)/d|z| using hydrostatic balance
    N2rhogr = N2rhogr .* (Pa2db * grav * r);
end
% Now N2rhogr is d(sigma)/d|z| whether P is pressure or inverted depth
% N.B. don't multiply by -1, because P increases downwards even if it is depth.

% Apply threshold:
N2rhogr(N2rhogr < N2rhogr_min) = N2rhogr_min;


% Slope error in x direction
sx = ex ./ A_X(N2rhogr);


% Slope error in y direction
sy = ey ./ A_Y(N2rhogr);


    function out = im1_2d(in)
        out = circshift(in, [+1 0]);
        if ~wrap(1)
            out(1,:) = nan;
        end
    end

    function out = jm1_2d(in)
        out = circshift(in, [0 +1]);
        if ~wrap(2)
            out(:,1) = nan;
        end
    end

    function out = ip1_2d(in)
        out = circshift(in, [-1 0]);
        if ~wrap(1)
            out(end,:) = nan;
        end
    end

    function out = jp1_2d(in)
        out = circshift(in, [0 -1]);
        if ~wrap(2)
            out(:,end) = nan;
        end
    end

    function out = im1_3d(in)
        out = circshift(in, [0 +1 0]);
        if ~wrap(1)
            out(:,1,:) = nan;
        end
    end

    function out = jm1_3d(in)
        out = circshift(in, [0 0 +1]);
        if ~wrap(2)
            out(:,:,1) = nan;
        end
    end

end
