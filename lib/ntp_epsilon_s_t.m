function [ex,ey] = ntp_epsilon_s_t(s, t, x, dx, dy, centre, wrap)
%NTP_EPSILON_S_T  The epsilon neutrality error using salinity and temperature
%
%
% [ex,ey] = ntp_epsilon_s_t(s,t,x,dx,dy,centre,wrap)
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
% For a non-Boussinesq ocean, x is pressure [dbar].
%
% For a Boussinesq ocean, x is depth [m, positive].
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
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: data can be arranged as latitude by longitude.  Simply swap dx and
%       dy, adjust wrap as needed, and note the outputs ex and ey will also
%       be swapped, as will be sx and sy.
%
% Note: physical units for s,t,x,S,T,X are determined by eos.m and eos_s_t.m.
%
% Note: X must be increasing along its first dimension.
%
% Note: dx and dy can have any dimension replaced by a singleton dimension.
%
%
% --- Output:
% ex [ni, nj]: zonal epsilon neutrality error [kg m^-4]
% ey [ni, nj]: meridional epsilon neutrality error [kg m^-4]
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
% Version   : 2.0.0
%
% Modified by : --
% Date        : --
% Changes     : --

% --- Input checks, set parameters
assert(~isempty(which('eos')), 'Cannot find eos.m in the path or current directory.');
assert(~isempty(which('eos_s_t')), 'Cannot find eos_s_t.m in the path or current directory.');

eos_is_specvol = (eos(34.5, 3, 1000) < 1);


if centre
    % Get errors by centred differences
    
    % Begin calculations for errors in x direction.
    sa = (im1_2d(s) + ip1_2d(s)) / 2;
    ta = (im1_2d(t) + ip1_2d(t)) / 2;
    xa = (im1_2d(x) + ip1_2d(x)) / 2;
    [rs,rt] = eos_s_t(sa,ta,xa);
    ex = (rs .* (ip1_2d(s) - im1_2d(s)) + rt .* (ip1_2d(t) - im1_2d(t))) ./ (2*dx);
    if eos_is_specvol
        va = eos(sa, ta, xa);   % specific volume
        ex = -ex ./ (va .* va); % convert to units of density gradient
    end
    if nargout < 2; return; end
    
    % Begin calculations for errors in y direction.
    sa = (jm1_2d(s) + jp1_2d(s)) / 2;
    ta = (jm1_2d(t) + jp1_2d(t)) / 2;
    xa = (jm1_2d(x) + jp1_2d(x)) / 2;
    [rs,rt] = eos_s_t(sa,ta,xa);
    ey = (rs .* (jp1_2d(s) - jm1_2d(s)) + rt .* (jp1_2d(t) - jm1_2d(t))) ./ (2*dy);
    if eos_is_specvol
        va = eos(sa, ta, xa);   % specific volume
        ey = -ey ./ (va .* va); % convert to units of density gradient
    end
    
else
    % Get errors by backward differences, and leave on the U, V grid.
    
    % Begin calculations for errors in x direction.
    sa = (im1_2d(s) + s) / 2;
    ta = (im1_2d(t) + t) / 2;
    xa = (im1_2d(x) + x) / 2;
    [rs,rt] = eos_s_t(sa,ta,xa);
    ex = (rs .* (s - im1_2d(s)) + rt .* (t - im1_2d(t))) ./ dx;
    if eos_is_specvol
        va = eos(sa, ta, xa);   % specific volume
        ex = -ex ./ (va .* va); % convert to units of density gradient
    end
    if nargout < 2; return; end
    
    % Begin calculations for errors in y direction.
    sa = (jm1_2d(s) + s) / 2;
    ta = (jm1_2d(t) + t) / 2;
    xa = (jm1_2d(x) + x) / 2;
    [rs,rt] = eos_s_t(sa,ta,xa);
    ey = (rs .* (s - jm1_2d(s)) + rt .* (t - jm1_2d(t))) ./ dy;
    if eos_is_specvol
        va = eos(sa, ta, xa);   % specific volume
        ey = -ey ./ (va .* va); % convert to units of density gradient
    end
    
end

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

end
