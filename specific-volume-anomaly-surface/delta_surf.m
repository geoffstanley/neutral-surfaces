function [x,d0,s0,t0] = delta_surf(S, T, X, s0, t0, var, tolx, OPTS)
%DELTASURF Specific volume anomaly surface by nonlinear solution in each water column.
%
%
% x = delta_surf(S, T, X, s0, t0, d0, tolx)
% finds pressure or depth x (with precision tolx) of the isosurface d0 of
% delta = eos(S,T,X) - eos(s0,t0,X) with reference practical / Absolute
% salinity s0 and reference potential / Conservative temperature t0, in an
% ocean with practical / Absolute salinity S and potential / Conservative
% temperature T at data sites where the pressure or depth is X.  The
% equation of state is given by eos.m in MATLAB's path, which accepts S, T,
% and X as its 3 inputs.  For a non-Boussinesq ocean, x and X are pressure
% [dbar] and eos gives the specific volume.  For a Boussinesq ocean, x and
% X are depth [m] positive and increasing down, and eos gives the in-situ
% density.
%
% [x, d0] = delta_surf(S, T, X, s0, t0, [i0, j0, x0], tolx)
% as above but finds the delta isosurface, delta = d0, that intersects the
% reference cast at grid indices (i0,j0) at pressure or depth x0.
%
% [x, d0, s0, t0] = delta_surf(S, T, X, [], [], [i0, j0, x0], tolx)
% as above but also finds the reference practical / Absolute salinity s0
% and reference potential / Conservative temperature t0 by interpolating S
% and T at the reference cast (i0,j0) to pressure or depth x0.
%
% ... = delta_surf(..., OPTS)
% passes the OPTS struct to the code generation.  See
% delta_surf_vertsolve_codegen.m for details.
%
%
% --- Input:
% S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% X [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
% s0 [1, 1] or []: reference S
% t0 [1, 1] or []: reference T
% var [1, 1] or [1, 3]: isovalue of delta, or location to intersect
% tolx [1, 1]: precision of the pressure [dbar] or depth [m] of the surface
% OPTS [struct]: code generation options
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: physical units for S, T, X, x, d0, s0, t0, x0 are determined by eos.m.
%
% Note: X must increase monotonically along its first dimension.
%
%
% --- Output:
% x [ni, nj]: pressure [dbar] or depth [m] of the delta surface
% d0 [1, 1]: isovalue of the delta surface
% s0 [1, 1]: reference S
% t0 [1, 1]: reference T
%
%
% --- Requirements:
% eos, interp_firstdim_twovars
% bisectguess - https://www.mathworks.com/matlabcentral/fileexchange/69710
% interp1qn - https://www.mathworks.com/matlabcentral/fileexchange/69713

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

[nk, ni, nj] = size(S);
[~,xj] = size(X);
NX = nk * double(xj > 1);


if nargin < 7 || isempty(tolx)
    tolx = 1e-4;
end
if nargin < 8 || isempty(OPTS)
    OPTS = struct(); % Let the code generation function decide
end

% Run codegen to create MEX function handling the main computation
delta_surf_vertsolve_codegen(nk, ni, nj, isvector(X), OPTS);


if isscalar(var)
    % s0 and t0 should be provided
    d0 = var(1);
    
else
    % var should be a 3 element vector
    i0 = var(1);
    j0 = var(2);
    x0 = var(3);
    
    % Get linear index to reference cast
    ij0 = sub2ind([ni,nj], i0, j0);
    idx = (ij0-1) * NX;
    
    % Number of valid bottles in reference cast
    K = sum(isfinite(S(:,ij0)),1);
    
    % Select the reference cast
    Xc = X((idx+1:idx+K).');
    Sc = S(1:K,ij0);
    Tc = T(1:K,ij0);
    
    if nargin < 6 || isempty(s0) && isempty(t0)
        % Get reference salinity and temperature at the chosen location
        [s0, t0] = interp_firstdim_twovars(x0, Xc, Sc, Tc);
    end
    
    % Choose iso-value that will intersect (i0,j0,p0).
    Dc = eos(Sc, Tc, Xc) - eos(s0, t0, Xc);
    d0 = interp_firstdim_twovars(x0, Xc, Dc, Dc);
    
end

% Solve non-linear root finding problem in each cast
x = delta_surf_vertsolve_mex(S, T, X, s0, t0, d0, tolx);
