function [x,s,t,d0] = pot_dens_surf(S, T, X, xref, var, OPTS)
%POT_DENS_SURF  Potential density surface by nonlinear solution in each water column.
%
%
% x = pot_dens_surf(S, T, X, xref, d0)
% finds the pressure or depth x of the isosurface d0 of potential density
% referenced to xref, in the ocean with practical / Absolute salinity S and
% potential / Conservative temperature T at data sites where the and
% pressure or depth is X.  The equation of state is given by eos.m in
% MATLAB's path, which accepts S, T, xref as its 3 inputs. For a
% non-Boussinesq ocean, x, X, and xref should be pressure [dbar].  For a
% Boussinesq ocean, x, X, and xref should be depth [m], positive and
% increasing down.  Algorithmic parameters are provided in OPTS (see
% "Options" below for further details).
%
% [p,d0] = pot_dens_surf(S, T, X, xref, [i0, j0, x0])
% as above but finds the surface intersecting a reference cast given by
% grid indices (i0,j0) at pressure or depth x0.
%
% ... = pot_dens_surf(..., OPTS)
% sets algorithmic parameters (see "Options" below for further details).
%
%
% --- Input:
% S [nk, ni, nj]: % S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% X [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
% xref [1, 1]: reference pressure [dbar] or depth [m]
% var [1, 1] or [1, 3]: isovalue of potential density, or target location
% OPTS [struct]: options (see below)
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: physical units for S, T, X, xref, x, d0 are determined by eos.m.
%
% Note: X must increase monotonically along its first dimension.
%
% --- Output:
% x [ni, nj]: pressure [dbar] or depth [m] of the surface
% s [ni, nj]: practical / Absolute salinity on the surface
% t [ni, nj]: potential / Conservative temperature the surface
% d0 [1, 1]: potential density [kg m^-3] or specific volume [m^3 kg^-1] of
%            the surface
%
%
% --- Options:
% OPTS is a struct containing the following fields.
%   FILE_ID [1, 1]: 1 to write any output to MATLAB terminal, or a file
%       identifier as returned by fopen() to write to a file. Default: 1.
%   INTERPFN [function handle]: vertical interpolation function, used to
%       evaluate SppX and TppX if those are not provided.  E.g. INTERPFN =
%       @ppc_linterp.
%   SppX [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose 
%       knots are X that interpolate S as a function of X in each water 
%       column.  E.g. SppX = ppc_linterp(X, S);
%   TppX [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose 
%       knots are X that interpolate T as a function of X in each water 
%       column.  E.g. TppX = ppc_linterp(X, T);
%   VERBOSE [scalar]: 0 for no output; 1 for summary of each iteration;
%                     2 for detailed information on each iteration.
%   X_TOL [1, 1]: error tolerance, in the same units as X [dbar] or [m], 
%      when root-finding to update the surface.

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
[~,XN] = size(X);
MX = nk * double(XN > 1);

% Set up default options:
DEFS = struct();
DEFS.X_TOL = 1e-4;  % error tolerance in the vertical dimension
DEFS.INTERPFN = @ppc_linterp;  % linear interpolation in the vertical
DEFS.FILE_ID = 1;  % standard output to MATLAB terminal
DEFS.VERBOSE = 1;  % show a moderate level of information

% Override any options with user-specified OPTS
OPTS = catstruct(DEFS, OPTS); 

% Run codegen to create MEX function handling the main computation
pot_dens_surf_vertsolve_codegen(nk, ni, nj, isvector(X), OPTS);

% Interpolate S and T as piecewise polynomials of X, or use pre-computed interpolants in OPTS.
if isfield(OPTS, 'SppX') && isfield(OPTS, 'TppX')
    [~, kIS, iIS, jIS] = size(OPTS.SppX);
    [~, kIT, iIT, jIT] = size(OPTS.TppX);
end
if exist('kIS', 'var') && ...
        nk-1 == kIS && kIS == kIT && ...
        ni   == iIS && iIS == iIT && ...
        nj   == jIS && jIS == jIT
    SppX = OPTS.SppX;
    TppX = OPTS.TppX;
else
    SppX = OPTS.INTERPFN(X, S);
    TppX = OPTS.INTERPFN(X, T);
end

% Count number of valid bottles per cast
BotK = squeeze(sum(uint16(isfinite(S)), 1, 'native'));

% Decide on the isosurface value
if isscalar(var)
    
    d0 = var;
    
else
    % var should be a 3 element vector
    i0 = var(1);
    j0 = var(2);
    x0 = var(3);
    
    % Get linear index to reference cast
    ij0 = sub2ind([ni,nj], i0, j0);
    n = (ij0-1) * MX;
    
    % Select the reference cast
    X0 = X((n+1:n+nk).');
    SppX0 = SppX(:,:,ij0);
    TppX0 = TppX(:,:,ij0);
    
    % Choose iso-value that will intersect (i0,j0,x0).
    [s0,t0] = ppc_val2(X0, SppX0, TppX0, x0);
    d0 = eos(s0, t0, xref);
    
end

% Calculate 3D field for vertical interpolation
if eos(34.5,3,1000) > 1
    sortdir = 'ascend'; % eos is in-situ density, increasing with 3rd argument
else
    sortdir = 'descend'; % eos is specific volume, decreasing with 3rd argument
end
D = sort(eos(S, T, xref), 1, sortdir);

% Get started with the discrete version (and linear interpolation)
x = interp1qn(d0, D, X);

% Solve non-linear root finding problem in each cast
[x, s, t] = pot_dens_surf_vertsolve_mex(SppX, TppX, X, BotK, x, xref, d0, OPTS.X_TOL);
