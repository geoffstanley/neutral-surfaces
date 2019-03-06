function [ex,ey,sx,sy] = error_neutral(S, T, P, p, dx, dy, centre, wrap, grav, rho_c, eos_is_specvol)
%ERROR_NEUTRAL  The neutrality error and slope error from the neutral tangent plane.
%
%
% [ex,ey,sx,sy] = error_neutral(S,T,P,p,dx,dy,centre,wrap)
% computes the zonal and meridional neutrality errors, ex and ey, and the
% zonal and meridional slope errors, sx and sy, for an approximately
% neutral surface whose pressure is p in an ocean whose equation of state
% for the in-situ density is given by eos.m in the path, and whose
% practical / Absolute salinity is S, potential / Conservative temperature
% is T, pressure is P. The zonal and meridional distances between
% neighbouring grid points are dx and dy. Results are calculated on the
% grid points of S, T and P if centre is true; otherwise, ex and sx are
% located zonally midway between S, T, and P grid points, and ey and sy are
% located meridionally midway between S, T, and P grid points. If wrap(1)
% is true, the data is periodic in longitude; if wrap(2) is true, the data
% is periodic in latitude.
%
% ... = error_neutral(S,T,Z,z,dx,dy,centre,wrap,grav,rho_c)
% as above, but for an approximately neutral surface of depth z in a
% Boussinesq ocean having S and T at depth Z, with grav the gravitational
% acceleration and rho_c the Boussinesq reference density.
%
% ... = error_neutral(..., eos_is_specvol)
% as above if eos_is_specvol is false. If eos_is_specvol is true, the
% specific volume, rather than the in-situ density, is given by eos.m.
%
%
% --- Input:
% S [nz, nx, ny]: practical / Absolute salinity
% T [nz, nx, ny]: potential / Conservative temperature
% P [nz, nx, ny]: pressure [dbar]
% Z [nz, nx, ny] or [nz, 1]: depth [m, positive]
% p     [nx, ny]: pressure of the surface [dbar]
% z     [nx, ny]: depth of the surface [m, positive]
% dx    [nx, ny]: Zonal      grid distance [m], between S(:,i,j) and S(:,i-1,j)
% dy    [nx, ny]: Meridional grid distance [m], between S(:,i,j) and S(:,i,j-1)
% centre  [1, 1]: true to compute outputs on the central grid, false on the half grids [logical]
% wrap    [2, 1]: wrap(i) is true when the data is periodic in i'th dimension of p
% grav    [1, 1]: gravitational acceleration [m s^-2]
% rho_c [] or [1, 1]: Boussinesq reference density
% eos_is_specvol [1, 1]: true when eos.m in the path gives specific volume [logical]
%
% Note: nz is the maximum number of data points per water column,
%       nx is the number of data points in longitude,
%       ny is the number of data points in latitude.
%
% Note: physical units for S and T are determined by eos.m and eosdp.m.
%
% Note: dx and dy can have any dimension replaced by 1.
%
%
% --- Output:
% ex [nx, ny]: zonal neutrality error [kg m^-4]
% ey [nx, ny]: meridional neutrality error [kg m^-4]
% sx [nx, ny]: zonal slope error [dimensionless]
% sy [nx, ny]: meridional slope error [dimensionless]
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
% operator in Stanley (2018).
%
% The above definition is equivalent to Stanley (2018) Eq. (28), but
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
% --- Equation of State:
% The MATLAB path must contain two files, eos.m and eosdp.m. If
% eos_is_specvol is false or not provided, R = eos(S,T,P) and 
% RP = eosdp(S,T,P). If eos_is_specvol is true, A = eos(S,T,P) and
% AP = eosdp(S,T,P). Here,
%   S is the practical / Absolute salinity,
%   T is the potential / Conservative temperature,
%   P is the pressure [dbar],
%   A is the specific volume [m^3 kg^-1], and
%   AP is the partial derivative of specific volume w.r.t. pressure [m^3 kg^-1 dbar^-1].
%   R is the in-situ density [kg m^-3], and
%   RP is the partial derivative of in-situ density w.r.t. pressure [kg m^-3 dbar^-1].
%
%
% --- Requirements:
% eos, eosdp
% interp1qn2, pchi1d1 - https://www.mathworks.com/matlabcentral/fileexchange/69713
%
%
% --- References:
% Klocker, A., McDougall, T. J. & Jackett, D. R. A new method for forming
% approximately neutral surfaces. Ocean Science 5, 155?172 (2009).
%
% McDougall, T. J. & Jackett, D. R. On the helical nature of neutral
% trajectories in the ocean. Progress in Oceanography 20, 153?183 (1988).
%
% Stanley, G. J.  Neutral surface topology. Ocean Modelling, submitted.

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

% --- Input checks, set parameters
narginchk(8,11)

assert(~isempty(which('eos')), 'Cannot find eos.m in the path or current directory.');

Pa2db = 1e-4;
N_min = 2e-4;
N2rhogr_min = N_min^2 * 1000 / 10;
BOUSSINESQ = (nargin >= 10) && isscalar(rho_c);
if BOUSSINESQ
    Z2P = Pa2db * rho_c * grav;
    p = p * Z2P; % Internally convert depths to Boussinesq pressure.
    P = P * Z2P;
end
if nargin < 11
    eos_is_specvol = false; % Default, assume eos is for in-situ density
end


if exist('pchi1d1_mex', 'file')
    pchi1d1fn = @pchi1d1_mex;
else
    pchi1d1fn = @pchi1d1;
end
lead1 = @(x) reshape(x, [1 size(x)]);

% --- Interpolate S and T onto the surface
[s,t] = interp1qn2(lead1(p), P, S, T);
s = squeeze(s);
t = squeeze(t);


if centre
    % --- Get errors by centred differences
    r = eos(s, t, p);
    if eos_is_specvol
        r = 1 ./ r;
    end
    
    
    sX = (im1_2d(s) + ip1_2d(s)) / 2;
    tX = (im1_2d(t) + ip1_2d(t)) / 2;
    pX = (im1_2d(p) + ip1_2d(p)) / 2;
    rpX = eosdp(sX, tX, pX);
    if eos_is_specvol
        rX = eos(sX, tX, pX);
        rpX = -rpX ./ rX.^2;
    end
    ex = ((ip1_2d(r) - im1_2d(r)) - rpX .* (ip1_2d(p) - im1_2d(p))) ./ (2*dx);
    % Note, ex above is virtually identical (third order accuracy) to exst below:
    % >> [rs,rt] = densjmd95_ds1dt1(sX,tX,pX);
    % >> ex = (rs .* (ip1_2d(s) - im1_2d(s)) + rt .* (ip1_2d(t) - im1_2d(t))) ./ (2*dx);
    % But the above two forms of ex are NOT identical to the below form:
    % >> rp = eosdp(s, t, p);
    % >> ex = ((ip1_2d(r) - im1_2d(r)) - rp .* (ip1_2d(p) - im1_2d(p))) ./ (2*dx);
    % Lesson: tempting though it is to use s,t and p in the local water
    % column to calculate the partial derivative of density w.r.t pressure
    % (rp), it's wrong!
    clear sX tX pX rX rpX
    if nargout < 2; return; end
    
    sY = (jm1_2d(s) + jp1_2d(s)) / 2;
    tY = (jm1_2d(t) + jp1_2d(t)) / 2;
    pY = (jm1_2d(p) + jp1_2d(p)) / 2;
    rpY = eosdp(sY, tY, pY);
    if eos_is_specvol
        rY = eos(sY, tY, pY);
        rpY = -rpY ./ rY.^2;
    end
    ey = ((jp1_2d(r) - jm1_2d(r)) - rpY .* (jp1_2d(p) - jm1_2d(p))) ./ (2*dy);
    clear sY tY pY rY rpY
    if nargout < 3; return; end
    
    % Calculate stratification using locally referenced potential density (sigma)
    % N2rhogr is rho N^2 / g
    sigma = eos(S, T, permute(p, [3 1 2]));
    if eos_is_specvol
        sigma = 1 ./ sigma;
    end
    [~, N2rhogr] = pchi1d1fn(lead1(p), P, sigma);
    clear sigma
    if BOUSSINESQ
        % Convert to d(sigma)/d|z|, as if we had used depth rather than P above
        N2rhogr = N2rhogr .* Z2P;
    else
        % P is actually pressure, not depth analogue.
        % Convert d(sigma)/dp into d(sigma)/d|z| using hydrostatic balance
        N2rhogr = N2rhogr .* (Pa2db * grav * r);
    end
    % Now N2rhogr is d(sigma)/d|z| whether P is pressure or inverted depth
    % N.B. don't multiply by -1, because P increases downwards even if it is depth.
    
    % Apply threshold:
    N2rhogr(N2rhogr < N2rhogr_min) = N2rhogr_min;
    
    % Slope error
    sx = ex ./ N2rhogr;
    if nargout < 4; return; end
    
    sy = ey ./ N2rhogr;
    
else
    % Get errors by backward differences, and leave on the U, V grid.
    r = eos(s, t, p);
    if eos_is_specvol
        r = 1 ./ r;
    end
    
    % --- Begin calculations for errors in x direction
    % rp is the derivative of density w.r.t. pressure as a thermodynamic coordinate, in the eos.
    sX = (im1_2d(s) + s) / 2;
    tX = (im1_2d(t) + t) / 2;
    pX = (im1_2d(p) + p) / 2;
    rpX = eosdp(sX, tX, pX);
    if eos_is_specvol || ~BOUSSINESQ
        rX = eos(sX, tX, pX);
    end
    if eos_is_specvol
        rX = 1 ./ rX;
        rpX = -rpX .* rX.^2;
    end
    clear sX tX
    
    
    ex = ((r - im1_2d(r)) - rpX .* (p - im1_2d(p))) ./ dx;
    % Note, ex above is virtually identical to ex below:
    % >> [rs,rt] = densjmd95_ds1dt1(sX,tX,pX);
    % >> ex = (rs .* (s - im1_2d(s)) + rt .* (t - im1_2d(t))) ./ dx;
    if nargout < 2; return; end
    
    % --- Begin calculations for errors in y direction.
    sY = (jm1_2d(s) + s) / 2;
    tY = (jm1_2d(t) + t) / 2;
    pY = (jm1_2d(p) + p) / 2;
    rpY = eosdp(sY, tY, pY);
    if eos_is_specvol || ~BOUSSINESQ
        rY = eos(sY, tY, pY);
    end
    if eos_is_specvol
        rY = 1 ./ rY;
        rpY = -rpY .* rY.^2;
    end
    clear sY tY s t
    
    ey = ((r - jm1_2d(r)) - rpY .* (p - jm1_2d(p))) ./ dy;
    if nargout < 3; return; end
    clear r
    
    % Start working on zonal slope error
    SX = (S + im1_3d(S)) / 2;
    TX = (T + im1_3d(T)) / 2;
    if isvector(P)
        PX = P;
    else
        PX = (P + im1_3d(P)) / 2;
    end
    
    % Compute locally referenced potential density
    sigma = eos(SX, TX, permute(pX, [3 1 2]));
    if eos_is_specvol
        sigma = 1 ./ sigma;
    end
    
    % Take derivative of sigma w.r.t pressure or depth (as a SPATIAL coordinate)
    [~, N2rhogr] = pchi1d1fn(lead1(pX), PX, sigma);
    clear SX TX PX sigma
    if BOUSSINESQ
        % Convert to d(sigma)/d|z|, as if we had used depth rather than P above
        N2rhogr = N2rhogr .* Z2P;
    else
        % P is actually pressure, not depth analogue.
        % Convert d(sigma)/dp into d(sigma)/d|z| using hydrostatic balance
        N2rhogr = N2rhogr .* (Pa2db * grav * rX);
    end
    % Now N2rhogr is d(sigma)/d|z| whether P is pressure or inverted depth
    % N.B. don't multiply by -1, because P increases downwards even if it is depth.
    
    % Apply threshold:
    N2rhogr(N2rhogr < N2rhogr_min) = N2rhogr_min;
    
    % Slope error
    sx = ex ./ N2rhogr ;
    
    
    % Repeat for meridional slope error
    SY = (S + jm1_3d(S)) / 2;
    TY = (T + jm1_3d(T)) / 2;
    if isvector(P)
        PY = P;
    else
        PY = (P + jm1_3d(P)) / 2;
    end
    sigma = eos(SY, TY, permute(pY, [3 1 2]));
    if eos_is_specvol
        sigma = 1 ./ sigma;
    end
    [~, N2rhogr] = pchi1d1fn(lead1(pY), PY, sigma);
    clear SY TY PY sigma
    if BOUSSINESQ
        N2rhogr = N2rhogr .* Z2P;
    else
        N2rhogr = N2rhogr .* (Pa2db * grav * rY);
    end
    N2rhogr(N2rhogr < N2rhogr_min) = N2rhogr_min;
    sy = ey ./ N2rhogr ;
    
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

    function out = im1_3d(in)
        out = circshift(in, [0 +1 0]);
        if ~wrap(1)
            out(1,:) = nan;
        end
    end

    function out = jm1_3d(in)
        out = circshift(in, [0 0 +1]);
        if ~wrap(2)
            out(:,1) = nan;
        end
    end

end
