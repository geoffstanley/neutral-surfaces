function [x,val] = pot_dens_surf(S, T, X, xref, var, tolx, OPTS)
%POT_DENS_SURF  Potential density surface by nonlinear solution in each water column.
%
%
% x = pot_dens_surf(S, T, X, xref, val, tolx)
% finds the pressure or depth x (with precision tolx) of the isosurface val
% of potential density referenced to xref, in the ocean with practical /
% Absolute salinity S and potential / Conservative temperature T at data
% sites where the and pressure or depth is X.  The equation of state is
% given by eos.m in MATLAB's path, which accepts S, T, xref as its 3
% inputs. For a non-Boussinesq ocean, x, X, and xref should be pressure
% [dbar].  For a Boussinesq ocean, x, X, and xref should be depth [m],
% positive and increasing down.
%
% [p,val] = pot_dens_surf(S, T, X, xref, [i0, j0, x0], tolx)
% as above but finds the surface intersecting a reference cast given by
% grid indices (i0,j0) at pressure or depth x0.
%
% ... = pot_dens_surf(..., OPTS)
% passes the OPTS struct to the code generation.  See
% pot_dens_surf_vertsolve_codegen.m for details.
%
%
% --- Input:
% S [nk, ni, nj]: % S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% X [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
% xref [1, 1]: reference pressure [dbar] or depth [m]
% var [1, 1] or [1, 3]: isovalue of potential density, or target location
% tolx [1, 1]: precision of the pressure [dbar] or depth [m] of the surface
% OPTS [struct]: code generation options
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: physical units for S, T, X, xref, x are determined by eos.m.
%
% Note: X must increase monotonically along its first dimension.
%
% --- Output:
% x [ni, nj]: pressure [dbar] or depth [m] of the surface
% val [1, 1]: potential density of the surface [kg m^-3]
%
%
% --- Requirements:
% eos, interp_firstdim_twovars, bisectguess
%
%
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
[~,Xj] = size(X);
MX = nk * double(Xj > 1);


if nargin < 6 || isempty(tolx)
    tolx = 1e-4;
end
if nargin < 7 || isempty(OPTS)
    OPTS = struct(); % Let the code generation function decide
end

% Run codegen to create MEX function handling the main computation
pot_dens_surf_vertsolve_codegen(nk, ni, nj, isvector(X), OPTS);


% Decide on the isosurface value
if isscalar(var)
    
    val = var;
    
else % var should be a 3 element vector
    
    i0 = var(1);
    j0 = var(2);
    p0 = var(3);
    
    % Get linear index to reference cast
    ij0 = sub2ind([ni,nj], i0, j0);
    idx = (ij0-1) * MX;
    
    % Number of valid bottles in reference cast
    K = sum(isfinite(S(:,ij0)),1);
    
    % Select the reference cast
    Xc = X((idx+1:idx+K).');
    Sc = S(1:K,ij0);
    Tc = T(1:K,ij0);
    
    % Choose iso-value that will intersect (i0,j0,p0).
    [s0,t0] = interp_firstdim_twovars(p0, Xc, Sc, Tc);
    val = eos(s0, t0, xref);
    
end

% Solve non-linear root finding problem in each cast
x = pot_dens_surf_vertsolve_mex(S, T, X, xref, val, tolx);
