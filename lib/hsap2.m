function [y,m] = hsap2(s, t, x, varargin)
%HSAP2  Integrate hydrostatic balance to obtain acceleration potential on a surface.
%
%
% y = hsap2D(s, t, p)
% with s and t scalars, computes 
% 1e4 * \int_0^p eos(s,t,p') dp'
% using trapezoidal integration with uniform interval of 0.1 dbar, where
% eos is a function in the path giving the specific volume.
%
% y = hsap2D(s, t, p, p0, P)
% computes
% 1e4 * \int_p0^p eos(s,t,p') dp'
% using trapezoidal integration on the intervals given by columns of P,
% where eos is a function in the path giving the specific volume. If s and
% t are the same size as P, A, and Y, they are linearly interpolated onto
% the surface.
%
% y = hsap2(s, t, p, P, A, Y)
% computes 
% grav * ETAN  +  1e4 * \int_ATMP^p A dp'
% using trapezoidal integration on the intervals given by columns of P.
% Here, grav is the gravitational acceleration, ETAN is the sea-surface
% height, and ATMP is the atmospheric pressure at the sea-surface, and A is
% the specific volume at the data sites of the pressure P. These are
% encapsulated in the acceleration potential from hydrostatic balance at
% the data sites of A and P, given by Y = hsap3(P, ATMP, ETAN, A, grav).
% The practical / Absolute salinity and the potential / Conservative
% temperature on the surface whose pressure is p are given by s and t,
% respectively; these are used to compute the specific volume on the
% surface itself a = eos(s, t, p), which is used in the final term in the
% summation that trapezoidally integrates between the surface itself and
% the data site just shallower than the surface. If s and t are the same
% size as P, A, and Y, they are linearly interpolated onto the surface.
%
% [y, a] = hsap2(...)
% also returns a = eos(s, t, p), the specific volume on the surface.
%
%
% For a Boussinesq ocean, instead use the following forms. Note that z > 0
% and Z > 0 and these increase downwards. 
%
% y = hsap2D(s, t, z, grav, rho_c)
% with s and t scalars, computes 
% (grav / rho_c) * \int_z0^z eos(s, t, z' * 1e-4 * grav * rho_c) dz'
% using trapezoidal integration with uniform interval of 0.1 m, where eos
% is a function in the path that gives the in-situ density.
%
% y = hsap2D(s, t, z, z0, Z, grav, rho_c)
% computes 
% (grav / rho_c) * \int_z0^z eos(s, t, z' * 1e-4 * grav * rho_c) dz'
% using trapezoidal integration on the intervals given by column(s) of Z,
% where eos is a function in the path giving the in-situ density. If s and
% t are the same size as R and Y, they are linearly interpolated onto the
% surface.
%
% y = hsap2D(s, t, z, Z, R, Y, grav, rho_c)
% computes 
% (1e-4 / rho_c) * ATMP + grav * ETAN + (grav / rho_c) * \int_0^z R dz'
% using trapezoidal integration on the intervals given by columns of Z.
% Here, grav is the gravitational acceleration, ETAN is the sea-surface
% height, and ATMP is the atmospheric pressure at the sea-surface, and R is
% the in-situ density at the depths Z. These are encapsulated in the
% acceleration potential from hydrostatic balance at the data sites of R,
% given by Y = hsap3(Z, ATMP, ETAN, R, grav, rho_c). The practical /
% Absolute salinity and the potential / Conservative temperature on the
% surface whose depth is z are given by s and t, respectively; these are
% used to compute the in-situ density on the surface itself r = eos(s, t, z
% * 1e-4 * grav * rho_c), which is used in the final term in the summation
% that trapezoidally integrates between the surface itself and the data
% site just shallower than the surface. If s and t are the same size as R
% and Y, they are linearly interpolated onto the surface.
%
% [y, r] = hsap2(..., grav, rho_c)
% also returns r = eos(s, t, z * 1e-4 * grav * rho_c), the in-situ density
% on the surface.
%
%
% --- Input:
% s [1, 1] or [nx, ny] or [nz, nx, ny]: the practical / Absolute salinity, 
%  as a reference value, the values on the surface, or at all data sites
% t [1, 1] or [nx, ny] or [nz, nx, ny]: the potential / Conservative temperature
%  as a reference value, the values on the surface, or at all data sites
% p [nx, ny]: pressure on surface [dbar]
% z [nx, ny]: depth on surface [m, positive]
% p0 [1, 1]: reference pressure [dbar]
% z0 [1, 1]: reference depth [m, positive]
% P [nz, nx, ny] or [nz, 1]: pressure [dbar]
% Z [nz, nx, ny] or [nz, 1]: depth [m, positive]
% A [nz, nx, ny]: specific volume [m^3 kg^-1]
% R [nz, nx, ny]: in-situ density [kg m^-3]
% Y [nz, nx, ny]: the acceleration potential from hydrostatic balance [m^2 s^-2]
% grav [1, 1]: the gravitational acceleration [m s^-2]
% rho_c [1, 1]: the Boussinesq reference density [kg m^-3]
%
% Note: nz is the maximum number of data points per water column,
%       nx is the number of data points in longitude,
%       ny is the number of data points in latitude.
%
% Note: physical units for s and t are determined by eos.m.
%
% Note: the Montgomery (1937) potential, and its upgraded forms owing to
% Zhang and Hogg (1992) and to McDougall and Klocker (2010), as well as the
% orthobaric Montgomery potential (Stanley 2018) all use s0 and t0 as
% constant scalar values. The Cunningham (2000) geostrophic stream function
% uses s0 and t0 as the S and T values on the surface.
%
%
% --- Output:
% y [nx, ny]: acceleration potential from hydrostatic balance [m^2 s^-2]
% a [nx, ny]: specific volume on the surface [m^3 kg^-1]
% r [nx, ny]: in-situ density on the surface [kg m^-3]
%
%
% --- Requirements:
% eos
% binsrchrightn
% interp1qn2 - https://www.mathworks.com/matlabcentral/fileexchange/69713
%
%
% --- References:
% Cunningham, S. A. Circulation and volume flux of the North Atlantic using
% synoptic hydrographic data in a Bernoulli inverse. Journal of marine
% research 58, 1?35 (2000).
%
% McDougall, T. J. & Klocker, A. An approximate geostrophic streamfunction
% for use in density surfaces. Ocean Modelling 32, 105?117 (2010).
%
% Montgomery, R. A suggested method for representing gradient flow in
% isentropic surfaces. Bull. Amer. Meteor. Soc 18, 210?212 (1937).
%
% Stanley, G. J.  The exact geostrophic stream function for neutral
% surfaces. Ocean Modelling, submitted.
%
% Zhang, H.-M. & Hogg, N. G. Circulation and water mass balance in the
% Brazil Basin. Journal of marine research 50, 385?420 (1992).

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

narginchk(3,8);

db2Pa = 1e4;  % Conversion from [dbar] to [Pa]
Pa2db = 1e-4; % Conversion from [Pa] to [dbar]
lead1 = @(x) reshape(x, [1 size(x)]);

[nx,ny] = size(x);
is2D = @(F) ismatrix(F) && all(size(F) == [nx,ny]);
assert(all(size(s) == size(t)), 'First and second inputs must be the same size');
assert(is2D(x), 'Third input must be [nx,ny]');

BOUSSINESQ = length(varargin) >= 2 && isscalar(varargin{end-1}) && isscalar(varargin{end});
if BOUSSINESQ
    grav = varargin{end-1};
    rho_c = varargin{end};
    varargin = varargin(1 : end-2);
    fac = grav / rho_c;
    Z2P = Pa2db * grav * rho_c;
else
    fac = -db2Pa;
end

if isscalar(s) && isscalar(t)
    % Either
    %   [y,a] = hsap2D(s0, t0, p)
    % or
    %   [y,r] = hsap2D(s0, t0, z, grav, rho_c)
    
    dx = .1;
    lo = min(0, floor((min(x(:)) / dx)) * dx - dx);
    hi = max(x(:)) + dx;
    X = (lo : dx : hi).';
    if BOUSSINESQ
        m = eos(s, t, x * Z2P);
        M = eos(s, t, X * Z2P);
    else
        m = eos(s, t, x);
        M = eos(s, t, X);
    end
    Y = cumsum(diff(X) .* (M(1:end-1) + M(2:end)), 1);
    
    k = ceil((x-lo) / dx); % X(k) < x <= X(k+1)
    k(isnan(k)) = 2; % Prevent out-of-bounds indexing. Where k should be nan, x is nan, so next line makes y nan at these spots.
    y = (fac/2) * (Y(k-1) + (M(k) + m) .* (x - X(k)));
    
else
    assert(length(varargin) >= 2, 'With non-constant s and t, must provide at least five inputs.')
    
    X = varargin{1 + isscalar(varargin{1})}; % either varargin{1} or varargin{2}
    
    nz = size(X,1);
    is1D = @(F) ismatrix(F) && all(size(F) == [nz,1]);
    is3D = @(F) ndims(F) == 3 && all(size(F) == [nz,nx,ny]);
    assert(is1D(X) || is3D(X), 'Fourth input must be [nz,1] or [nz,nx,ny]');
    
    if is3D(s) % Interpolate 3D s and t onto the surface
        try
            [s,t] = interp1qn2_mex(lead1(x), X, s, t);
        catch
            [s,t] = interp1qn2(lead1(x), X, s, t);
        end
    end
    
    if isscalar(varargin{1})
        % Either
        %   [y,a] = hsap2D(s, t, p, p0, P)
        % or
        %   [y,r] = hsap2D(s, t, z, z0, Z, grav, rho_c)
        
        x0 = varargin{1};
        
        if BOUSSINESQ
            m = eos(s, t, x * Z2P); % in-situ density on the surface
            M = eos(lead1(s), lead1(t), X * Z2P); % in-situ density profiles using S and T from surface
            Y = hsap3(X,0,0,M,grav,rho_c);
        else
            m = eos(s, t, x); % specific volume on the surface
            M = eos(lead1(s), lead1(t), X); % specific volume profiles using S and T from surface
            Y = hsap3(X,0,0,M,0);
        end
        
        % Find vertical grid index for cell centre just above the surface
        k = binsrchrightn(x, X); % X(k) <= x < X(k+1)
        
        [ii,jj] = ndgrid(1:nx, 1:ny);
        inds = sub2ind([nz,nx,ny], k, ii, jj);
        yk = Y(inds);     % Y at cell centre just shallower than the surface
        mk = M(inds);     % M at cell centre just shallower than the surface
        if isvector(X)
            xk = X(k);    % X at cell centre just shallower than the surface
        else
            xk = X(inds); % X at cell centre just shallower than the surface
        end
        
        % Compute y on the surface
        y = yk + (fac/2) * (m + mk) .* (x - xk);
        
        
        % Now subtract the integral from x0
        if BOUSSINESQ
            m0 = eos(s, t, x0 * Z2P);
        else
            m0 = eos(s, t, x0);
        end
        k0 = binsrchrightn(x0, X);
        if isvector(X) % then k0 is a scalar, since x0 is definitely a scalar
            yk = squeeze(Y(k0,:,:)); % Y at cell centre just shallower than x0
            mk = squeeze(M(k0,:,:)); % M at cell centre just shallower than x0
            xk = X(k0);              % X at cell centre just shallower than x0
        else
            inds = sub2ind([nz,nx,ny], k0, ii, jj);
            yk = Y(inds); % Y at cell centre just shallower than x0
            mk = M(inds); % M at cell centre just shallower than x0
            xk = X(inds); % X at cell centre just shallower than x0
        end
        y = y - yk - (fac/2) * (mk + m0) .* (x0 - xk);
        
    else
        % Either
        %   [y,a] = hsap2D(s, t, p, P, A, Y)
        % or:
        %   [y,r] = hsap2D(s, t, z, Z, R, Y, grav, rho_c)
        
        assert(length(varargin) >= 3, 'With first three inputs non-constant, must provide at least six inputs');
        
        M = varargin{2}; % either density or specific volume (M for mass)
        Y = varargin{3};
        
        if BOUSSINESQ
            m = eos(s, t, x * Z2P); % in-situ density on the surface
        else
            m = eos(s, t, x); % specific volume on the surface
        end
        
        % Find vertical grid index for cell centre just above the surface
        k = binsrchrightn(x, X); % X(k) <= x < X(k+1)
        
        [ii,jj] = ndgrid(1:nx, 1:ny);
        inds = sub2ind([nz,nx,ny], k, ii, jj);
        yk = Y(inds);     % Y at cell centre just shallower than the surface
        mk = M(inds);     % M at cell centre just shallower than the surface
        if isvector(X)
            xk = X(k);    % X at cell centre just shallower than the surface
        else
            xk = X(inds); % X at cell centre just shallower than the surface
        end
        
        % Compute y on the surface
        y = yk + (fac/2) * (m + mk) .* (x - xk);
        
    end
    
end
end