function [y,m] = hsap2(s, t, p, varargin)
%HSAP2  Integrate hydrostatic balance to obtain acceleration potential on a surface.
%
%
% y = hsap2(s, t, p)
% with s and t scalars, computes
% 1e4 * \int_0^p eos(s,t,p') dp'
% using trapezoidal integration with uniform interval of 0.1 dbar, where
% eos is a function in the path giving the specific volume.
%
% y = hsap2(s, t, p, p0, P)
% computes
% 1e4 * \int_p0^p eos(s,t,p') dp'
% using trapezoidal integration on the intervals given by columns of P,
% where eos is a function in the path giving the specific volume. If s and
% t have an extra (vertical) dimension relative to p, they are interpolated
% onto the surface.
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
% the data site just shallower than the surface. If s and t have 2 more
% dimensions than p, they are interpreted as cofficients of piecewise
% polynomials for the 3D S and T as functions of P in each water column,
% and they are evaluated onto the surface.
%
% [y, a] = hsap2(...)
% also returns a = eos(s, t, p), the specific volume on the surface.
%
%
% For a Boussinesq ocean, instead use the following forms. Note that z > 0
% and Z > 0 and these increase downwards.
%
% y = hsap2(s, t, z, grav, rho_c)
% with s and t scalars, computes
% (grav / rho_c) * \int_z0^z eos(s, t, z') dz'
% using trapezoidal integration with uniform interval of 0.1 m, where eos
% is a function in the path that gives the in-situ density.
%
% y = hsap2(s, t, z, z0, Z, grav, rho_c)
% computes
% (grav / rho_c) * \int_z0^z eos(s, t, z') dz'
% using trapezoidal integration on the intervals given by column(s) of Z,
% where eos is a function in the path giving the in-situ density. If s and
% t are the same size as R and Y, they are interpolated onto the surface.
%
% y = hsap2(s, t, z, Z, R, Y, grav, rho_c)
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
% used to compute the in-situ density on the surface itself r = eos(s, t,
% z), which is used in the final term in the summation that trapezoidally
% integrates between the surface itself and the data site just shallower
% than the surface. If s and t have 2 more dimensions than p, they are
% interpreted as cofficients of piecewise polynomials for the 3D S and T as
% functions of P in each water column, and they are evaluated onto the
% surface.
%
% [y, r] = hsap2(..., grav, rho_c)
% also returns r = eos(s, t, z), the in-situ density
% on the surface.
%
%
% --- Input:
% s [1, 1] or [ni, nj] or [O, nk-1, ni, nj]: the practical / Absolute
%   salinity, as a reference value, the values on the surface, or a
%   piecewise polynomial interpolant in each water column
% t [1, 1] or [ni, nj] or [O, nk-1, ni, nj]: the potential / Conservative
%   temperature as a reference value, the values on the surface, or a
%   piecewise polynomial interpolant in each water column
% p [ni, nj]: pressure on surface [dbar]
% z [ni, nj]: depth on surface [m, positive]
% p0 [1, 1]: reference pressure [dbar]
% z0 [1, 1]: reference depth [m, positive]
% P [nk, ni, nj] or [nk, 1]: pressure [dbar]
% Z [nk, ni, nj] or [nk, 1]: depth [m, positive]
% A [nk, ni, nj]: specific volume [m^3 kg^-1]
% R [nk, ni, nj]: in-situ density [kg m^-3]
% Y [nk, ni, nj]: the acceleration potential from hydrostatic balance [m^2 s^-2]
% grav [1, 1]: the gravitational acceleration [m s^-2]
% rho_c [1, 1]: the Boussinesq reference density [kg m^-3]
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: physical units for s and t are determined by eos.m.
%
% Note: the Montgomery (1937) potential, and its upgraded forms owing to
% Zhang and Hogg (1992) and to McDougall and Klocker (2010), as well as the
% orthobaric Montgomery potential (Stanley 2019) all use s0 and t0 as
% constant scalar values. The Cunningham (2000) geostrophic stream function
% uses s0 and t0 as the S and T values on the surface.
%
%
% --- Output:
% y [ni, nj]: acceleration potential from hydrostatic balance [m^2 s^-2]
% a [ni, nj]: specific volume on the surface [m^3 kg^-1]
% r [ni, nj]: in-situ density on the surface [kg m^-3]
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
% Stanley, G.J., 2019. The exact geostrophic streamfunction for neutral
% surfaces. Ocean Modelling 138, 107–121.
%
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


narginchk(3,8);

lead1 = @(p) reshape(p, [1 size(p)]);

[ni,nj] = size(p);
assert(all(size(s) == size(t)), 'First and second inputs must be the same size');

% Check for gravity and Boussinesq reference density as last 2 arguments
BOUSSINESQ = length(varargin) >= 2 && isscalar(varargin{end-1}) && isscalar(varargin{end});
if BOUSSINESQ
    grav = varargin{end-1};
    rho_c = varargin{end};
    varargin = varargin(1 : end-2);
    fac = grav / rho_c;
else
    db2Pa = 1e4;  % Conversion from [dbar] to [Pa]
    fac = -db2Pa;
end

if isscalar(s) && isscalar(t)
    % Either
    %   [y,a] = hsap2(s0, t0, p)
    % or
    %   [y,r] = hsap2(s0, t0, z, grav, rho_c)
    
    dp = .1;
    lo = min(0, floor((min(p(:)) / dp)) * dp - dp);
    hi = max(p(:)) + dp;
    P = (lo : dp : hi).';
    m = eos(s, t, p);
    M = eos(s, t, P);
    Y = cumsum(diff(P) .* (M(1:end-1) + M(2:end)), 1);
    
    k = ceil((p-lo) / dp); % P(k) < p <= P(k+1)
    k(isnan(k)) = 2; % Prevent out-of-bounds indexing. Where k should be nan, p is nan, so next line makes y nan at these spots.
    y = (fac/2) * (Y(k-1) + (M(k) + m) .* (p - P(k)));
    
else
    assert(length(varargin) >= 2, 'With non-constant s and t, must provide at least five inputs.')
    
    P = varargin{1 + isscalar(varargin{1})}; % either varargin{1} or varargin{2}
    
    nk = size(P,1);
    is1D = @(F) ismatrix(F) && all(size(F) == [nk,1]);
    is3D = @(F) ndims(F) == 3 && all(size(F) == [nk,ni,nj]);
    is4D = @(F) ndims(F) == 4 && size(F,2) == nk-1 && size(F,3) == ni && size(F,4) == nj;
    assert(is1D(P) || is3D(P), 'Fourth input must be [nk,1] or [nk,ni,nj]');
    
    if is4D(s) % Evaluate interpolants for S and T onto the surface
        [s,t] = ppc_val2(P, s, t, lead1(p));
    end
    
    if isscalar(varargin{1})
        % Either
        %   [y,a] = hsap2(s, t, p, p0, P)
        % or
        %   [y,r] = hsap2(s, t, z, z0, Z, grav, rho_c)
        
        p0 = varargin{1};
        
        if BOUSSINESQ
            m = eos(s, t, p); % in-situ density on the surface
            M = eos(lead1(s), lead1(t), P); % in-situ density profiles using S and T from the surface
            Y = hsap3(P,0,0,M,grav,rho_c);
        else
            m = eos(s, t, p); % specific volume on the surface
            M = eos(lead1(s), lead1(t), P); % specific volume profiles using S and T from the surface
            Y = hsap3(P,0,0,M,0);
        end
        
        % Find vertical grid index for cell centre just above the surface
        k = binsrchrightn(p, P); % P(k) <= p < P(k+1)
        
        [ii,jj] = ndgrid(1:ni, 1:nj);
        inds = sub2ind([nk,ni,nj], k, ii, jj);
        yk = Y(inds);     % Y at cell centre just shallower than the surface
        mk = M(inds);     % M at cell centre just shallower than the surface
        if isvector(P)
            pk = P(k);    % P at cell centre just shallower than the surface
        else
            pk = P(inds); % P at cell centre just shallower than the surface
        end
        
        % Compute y on the surface
        y = yk + (fac/2) * (m + mk) .* (p - pk);
        
        
        % Now subtract the integral from p0
        m0 = eos(s, t, p0);
        k0 = binsrchrightn(p0, P);
        if isvector(P) % then k0 is a scalar, since p0 is definitely a scalar
            yk = squeeze(Y(k0,:,:)); % Y at cell centre just shallower than p0
            mk = squeeze(M(k0,:,:)); % M at cell centre just shallower than p0
            pk = P(k0);              % P at cell centre just shallower than p0
        else
            inds = sub2ind([nk,ni,nj], k0, ii, jj);
            yk = Y(inds); % Y at cell centre just shallower than p0
            mk = M(inds); % M at cell centre just shallower than p0
            pk = P(inds); % P at cell centre just shallower than p0
        end
        y = y - yk - (fac/2) * (mk + m0) .* (p0 - pk);
        
    else
        % Either
        %   [y,a] = hsap2(s, t, p, P, A, Y)
        % or:
        %   [y,r] = hsap2(s, t, z, Z, R, Y, grav, rho_c)
        
        assert(length(varargin) >= 3, 'With first three inputs non-constant, must provide at least six inputs');
        
        M = varargin{2}; % either density or specific volume (M for mass)
        Y = varargin{3};
        
        m = eos(s, t, p); % in-situ density or specific volume on the surface
        
        % Find vertical grid index for cell centre just above the surface
        k = binsrchrightn(p, P); % P(k) <= p < P(k+1)
        
        [ii,jj] = ndgrid(1:ni, 1:nj);
        inds = sub2ind([nk,ni,nj], k, ii, jj);
        yk = Y(inds);     % Y at cell centre just shallower than the surface
        mk = M(inds);     % M at cell centre just shallower than the surface
        if isvector(P)
            pk = P(k);    % P at cell centre just shallower than the surface
        else
            pk = P(inds); % P at cell centre just shallower than the surface
        end
        
        % Compute y on the surface
        y = yk + (fac/2) * (m + mk) .* (p - pk);
        
    end
    
end
end