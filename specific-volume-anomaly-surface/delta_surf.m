function [p,d0,s0,t0] = deltasurf(S, T, P, s0, t0, var, tol, Z2P) %#codegen
%DELTASURF Specific volume anomaly surface by nonlinear solution in each water column.
%
%
% p = deltasurf(S, T, P, s0, t0, d0, tol) 
% finds the pressure p (with precision tol) of the isosurface d0 of delta =
% eos(S,T,P) - eos(s0,t0,P) with reference practical / Absolute salinity s0
% and reference potential / Conservative temperature t0, in an ocean with
% practical / Absolute salinity S, potential / Conservative temperature T,
% and pressure P, and equation of state for the specific volume given by
% eos.m in the path.
%
% p = deltasurf(S, T, P, s0, t0, [i0, j0, p0], tol)
% as above but finds the surface intersecting a reference cast given by
% grid indices (i0,j0) at pressure p0.
%
% p = deltasurf(S, T, P, [], [], [i0, j0, p0], tol)
% selects the reference practical / Absolute salinity s0 and reference
% potential / Conservative temperature t0 by linearly interpolating S and T
% at the reference cast (i0,j0) to pressure p0.
%
% z = deltasurf(S, T, Z, s0, t0, d0, tol, Z2P) 
% finds the depth z (with precision tol) of the isosurface d0 of delta =
% eos(S,T,Z,Z2P) - eos(s0,t0,Z,Z2P) with reference practical / Absolute
% salinity s0 and reference potential / Conservative temperature t0, in a
% Boussinesq ocean with practical / Absolute salinity S and potential /
% Conservative temperature T at depth Z, equation of state for the in-situ
% density given by eos.m in the path, and depth to pressure conversion Z2P
% (typically given by 1e-4 * rho_c * grav, where rho_c is the Boussinesq
% reference density and grav is the gravitational acceleration).
%
% z = deltasurf(S, T, Z, s0, t0, [i0, j0, z0], tol, Z2P)
% as above but finds the surface intersecting a reference cast given by
% grid indices (i0,j0) at depth z0.
%
% z = deltasurf(S, T, Z, [], [], [i0, j0, z0], tol)
% selects the reference practical / Absolute salinity s0 and reference
% potential / Conservative temperature t0 by linearly interpolating S and T
% at the reference cast (i0,j0) to depth z0.
%
% [...,d0,s0,t0] = deltasurf(...)
% also returns the isovalue of the delta surface d0, the reference
% practical / Absolute salinity s0, and the reference potential /
% Conservative temperature t0.
%
%
% --- Input:
% S [nz, nx, ny]: practical / Absolute salinity
% T [nz, nx, ny]: potential / Conservative temperature
% P [nz, nx, ny]: pressure [dbar]
% Z [nz, nx, ny] or [nz, 1]: depth [m, positive]
% s0 [1, 1] or []: reference S
% t0 [1, 1] or []: reference T
% var [1, 1] or [1, 3]: isovalue of delta, or location to intersect
% tol [1, 1]: precision of the pressure or depth of the surface [dbar or m]
% Z2P [1, 1]: Boussinesq conversion from depth to pressure [m dbar^-1]
%
% Note: nz is the maximum number of data points per water column,
%       nx is the number of data points in longitude,
%       ny is the number of data points in latitude.
%
% Note: physical units for S, T, s0, and t0 are determined by eos.m. 
%
% Note: P and Z increase monotonically along the first dimension. 
%
%
% --- Output:
% p [nx, ny]: pressure of the delta surface [dbar]
% z [nx, ny]: depth of the delta surface [m, positive]
% d0 [1, 1]: isovalue of the delta surface [kg m^-3]
% s0 [1, 1]: reference S
% t0 [1, 1]: reference T
%
%
% --- Requirements:
% eos
% bisectguess - https://www.mathworks.com/matlabcentral/fileexchange/69710
% interp1qn, interp1qn2 - https://www.mathworks.com/matlabcentral/fileexchange/69713
%
%
% --- Code generation:
% codegen can be run as follows, for the Boussinesq version (modify as
% needed for other input forms):
% >> mexconfig = coder.config('mex');
% >> mexconfig.ExtrinsicCalls = false;
% >> mexconfig.ResponsivenessChecks = false;
% >> mexconfig.IntegrityChecks = false;
% >> type_S = coder.typeof(0, [nz, nx, ny], [false, false, false]);
% >> type_Z  = coder.typeof(0, [nz, 1], [false, false]);
% >> codegen('deltasurf', '-args', {type_S, type_S, type_P, 0, 0, [0, 0, 0], 0, 0}, '-config', mexconfig, '-o', 'deltasurf_mex');

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


if nargin < 7 || isempty(tol)
    tol = 1e-4;
end
BOUSSINESQ = nargin == 8 && isscalar(Z2P);
if BOUSSINESQ
    P = P * Z2P;
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
    % s0 and t0 should be provided
    d0 = var(1);
    
else
    % var should be a 3 element vector
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
    
    if nargin < 6 || isempty(s0) && isempty(t0)
        % Get reference salinity and temperature at the chosen location
        [s0, t0] = interp1qn2(p0, Pc, Sc, Tc);
    end
    
    % Choose iso-value that will intersect (i0,j0,p0).
    Dc = sort(eos(Sc, Tc, Pc) - eos(s0, t0, Pc), 1, sortdir);
    d0 = interp1qn(p0, Pc, Dc);
    
end


% Loop over each cast
idx = 0;
p = nan(nx,ny);
for c = 1:nx*ny
    K = BotK(c);
    if K >= 2
        
        % Select cast
        Pc = P((idx+1:idx+K).');
        Sc = S(1:K,c);
        Tc = T(1:K,c);
        
        % Get started with the discrete version.
        Dc = sort(eos(Sc, Tc, Pc) - eos(s0, t0, Pc), 1, sortdir);
        p(c) = interp1qn(d0, Dc, Pc);

        % Now refine by solving the non-linear problem at each water column
        p(c) = bisectguess(@diff_fun, P(idx+1), P(idx+K), tol, p(c), ...
            P((idx+1:idx+K).'), S(1:K,c), T(1:K,c), s0, t0, d0);
    end
    idx = idx + NP;
end

if BOUSSINESQ
    p = p / Z2P;
end

end

function out = diff_fun(p, P, S, T, s0, t0, d0)
[s,t] = interp1qn2(p, P, S, T);
out = eos(s, t, p) - eos(s0, t0, p) - d0 ;
end