function [p,val] = isopycnal(S, T, P, pref, var, tol, Z2P) %#codegen
%ISOPYCNAL  Potential density surface by nonlinear solution in each water column.
%
%
% p = isopycnal(S, T, P, pref, val, tol)
% finds the pressure p (with precision tol) of the isosurface val of
% potential density referenced to pref, in the ocean with practical /
% Absolute salinity S, potential / Conservative temperature T, and pressure
% P, and equation of state given by eos.m in the path.
%
% [p,val] = isopycnal(S, T, P, pref, [i0, j0, p0], tol)
% as above but finds the surface intersecting a reference cast given by
% grid indices (i0,j0) at pressure p0.
%
% z = isopycnal(S, T, Z, zref, val, tol, Z2P) 
% finds the depth z (with precision tol) of the isosurface val of potential
% density referenced to zref in a Boussinesq ocean, with Z is the depth of
% each grid point, in a Boussinesq ocean with practical / Absolute salinity
% S, potential / Conservative temperature T, equation of state given by
% eos.m in the path, and depth to pressure conversion Z2P (typically given
% by 1e-4 * rho_c * grav, where rho_c is the Boussinesq reference density
% and grav is the gravitational acceleration).
%
% [z, val] = isopycnal(S, T, Z, zref, [i0, j0, z0], tol, Z2P)
% as above but finds the surface intersecting a reference cast given by
% grid indices (i0,j0) at depth z0.
%
%
% --- Input:
% S [nz, nx, ny]: % S [nz, nx, ny]: practical / Absolute salinity
% T [nz, nx, ny]: potential / Conservative temperature
% P [nz, nx, ny]: pressure [dbar]
% Z [nz, nx, ny] or [nz, 1]: depth [m, positive]
% pref [1, 1]: reference pressure [dbar]
% zref [1, 1]: reference depth [m, positive]
% var [1, 1] or [1, 3]: isovalue of potential density, or location
% tol [1, 1]: precision of the pressure or depth of the surface [dbar or m]
% Z2P [1, 1]: Boussinesq conversion from depth to pressure [m dbar^-1]
%
% Note: nz is the maximum number of data points per water column,
%       nx is the number of data points in longitude,
%       ny is the number of data points in latitude.
%
% Note: physical units for S and T are determined by eos.m. 
%
% Note: P and Z increase monotonically along the first dimension. 
%
% --- Output: 
% p [nx, ny]: pressure of the surface [dbar]
% z [nx,ny]: depth of the surface [m, positive]
% val [1, 1]: potential density of the surface [kg m^-3]
%
%
% --- Requirements: 
% eos
% bisectguess - https://www.mathworks.com/matlabcentral/fileexchange/69710
% interp1qn, interp1qn2 - https://www.mathworks.com/matlabcentral/fileexchange/69713
%
%
% --- Code generation:
% codegen can be run as follows, for the Boussinesq version intersecting a
% given water column at a given depth (modify as needed for other input
% forms):
% >> mexconfig = coder.config('mex');
% >> mexconfig.ExtrinsicCalls = false;
% >> mexconfig.ResponsivenessChecks = false;
% >> mexconfig.IntegrityChecks = false;
% >> type_S = coder.typeof(0, [nz, nx, ny], [false, false, false]);
% >> type_P  = coder.typeof(0, [nz, nx, ny], [false, true, true]);
% >> codegen('isopycnal', '-args', {type_S, type_S, type_P, 0, [0, 0, 0], 0, 0}, '-config', mexconfig, '-o', 'isopycnal_ijp_mex');

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

[nz, nx, ny] = size(S);
[~,pj] = size(P);
NP = nz * double(pj > 1);


if nargin < 6 || isempty(tol)
    tol = 1e-4;
end
BOUSSINESQ = nargin == 7 && isscalar(Z2P);
if BOUSSINESQ
    P = P * Z2P;
    pref = pref * Z2P;
    tol = tol * Z2P;
end

if eos(35,0,0) > 1
    sortdir = 'ascend'; % eosis in-situ density, increasing in 3rd dim
else
    sortdir = 'descend'; % eos is specific volume, decreasing in 3rd dim
end

% Count number of valid bottles per cast
BotK = squeeze(sum(isfinite(S),1));


if isscalar(var)
    
    val = var;
    
else % var should be a 3 element vector
    
    i0 = var(1);
    j0 = var(2);
    p0 = var(3);
    if BOUSSINESQ
        p0 = p0 * Z2P;
    end
    
    % Get linear index to reference cast
    ij0 = sub2ind([nx,ny], i0, j0);
    idx = (ij0-1) * NP;
    
    K = BotK(ij0);
    
    % Select the reference cast
    Pc = P((idx+1:idx+K).');
    Sc = S(1:K,ij0);
    Tc = T(1:K,ij0);
    
    % Choose iso-value that will intersect (i0,j0,p0).
    [s0,t0] = interp1qn2(p0, Pc, Sc, Tc);
    val = eos(s0, t0, pref);
    
end

% Loop over each cast
p = nan(nx,ny);
idx = 0;
for c = 1:nx*ny
    K = BotK(c);
    if K >= 2
        
        % Select cast
        Pc = P((idx+1:idx+K).');
        Sc = S(1:K,c);
        Tc = T(1:K,c);
        
        % Get started with the discrete version.
        Dc = sort(eos(Sc, Tc, pref), 1, sortdir);
        p(c) = interp1qn(val, Dc, Pc);
        
        % Now refine by solving the non-linear problem at each water column
        p(c) = bisectguess(@diff_fun, P(idx+1), P(idx+K), tol, p(c), ...
            P((idx+1:idx+K).'), S(1:K,c), T(1:K,c), pref, val);
    end
    idx = idx + NP;
end

if BOUSSINESQ
    p = p / Z2P;
end

end


function out = diff_fun(p, P, S, T, pref, val)
[s,t] = interp1qn2(p, P, S, T);
out = eos(s, t, pref) - val;
end