function [x,s,t,d0,s0,t0] = delta_surf(S, T, X, s0, t0, var, OPTS)
%DELTASURF Specific volume anomaly surface by nonlinear solution in each water column.
%
%
% x = delta_surf(S, T, X, s0, t0, d0)
% finds pressure or depth x of the isosurface d0 of delta = eos(S,T,X) -
% eos(s0,t0,X) with reference practical / Absolute salinity s0 and
% reference potential / Conservative temperature t0, in an ocean with
% practical / Absolute salinity S and potential / Conservative temperature
% T at data sites where the pressure or depth is X.  The equation of state
% is given by eos.m in MATLAB's path, which accepts S, T, and X as its 3
% inputs.  For a non-Boussinesq ocean, x and X are pressure [dbar] and eos
% gives the specific volume.  For a Boussinesq ocean, x and X are depth [m]
% positive and increasing down, and eos gives the in-situ density.
%
% [x, d0] = delta_surf(S, T, X, s0, t0, [i0, j0, x0])
% as above but finds the delta isosurface, delta = d0, that intersects the
% reference cast at grid indices (i0,j0) at pressure or depth x0.
%
% [x, d0, s0, t0] = delta_surf(S, T, X, [], [], [i0, j0, x0])
% as above but also finds the reference practical / Absolute salinity s0
% and reference potential / Conservative temperature t0 by interpolating S
% and T at the reference cast (i0,j0) to pressure or depth x0.
%
% ... = delta_surf(..., OPTS)
% sets algorithmic parameters (see "Options" below for further details).
%
%
% --- Input:
% S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% X [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
% s0 [1, 1] or []: reference S
% t0 [1, 1] or []: reference T
% var [1, 1] or [1, 3]: isovalue of delta, or location to intersect
% OPTS [struct]: options (see below)
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
% x [ni, nj]: pressure [dbar] or depth [m] of the surface
% s [ni, nj]: practical / Absolute salinity on the surface
% t [ni, nj]: potential / Conservative temperature the surface
% d0 [1, 1]: isovalue of the delta surface
% s0 [1, 1]: reference S that defines delta
% t0 [1, 1]: reference T that defines delta
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
%
%
% --- Equation of State:
% The MATLAB path* must contain the function eos.m. This must accept 3
% inputs: S, T, and X. eos(S, T, X) returns the specific volume [m^3 kg^-1]
% or the in-situ density [kg m^-3].
% *Note: It is not sufficient to simply have these eos functions in the
% current working directory, because the compiled MEX functions will not be
% able to find them there.  They must be in the MATLAB path.  If they are
% in the current working directory, use `addpath(pwd)` to add the current
% working directory to the top of MATLAB's path.


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

[nk, ni, nj] = size(S);
[~,XN] = size(X);
KX = nk * double(XN > 1);

% Set up default options:
DEFS = struct();
DEFS.X_TOL = 1e-4;  % error tolerance in the vertical dimension
DEFS.INTERPFN = @ppc_linterp;  % linear interpolation in the vertical
DEFS.FILE_ID = 1;  % standard output to MATLAB terminal
DEFS.VERBOSE = 1;  % show a moderate level of information

% Override any options with user-specified OPTS
if nargin < 7 || isempty(OPTS)
  OPTS = DEFS;
else
  OPTS = catstruct(DEFS, OPTS); 
end


% Run codegen to create MEX function handling the main computation
ni_ = max(ni, 2048); % using variable size code generation and avoiding
nj_ = max(nj, 2048); % recompiling all the time
delta_surf_vertsolve_codegen(nk, ni_, nj_, isvector(X), OPTS);

% Interpolate S and T as functions of X, or use pre-computed interpolants in OPTS.
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
BotK = squeeze(sum(isfinite(S), 1));

% Decide on the isosurface value
if isscalar(var)
    % s0 and t0 should be provided
    assert(isscalar(s0) && isscalar(t0), 's0 and t0 must be provided as scalars if d0 is provided.');
    d0 = var(1);
    
else
    % var should be a 3 element vector
    i0 = var(1);
    j0 = var(2);
    x0 = var(3);
    
    % Get linear index to reference cast
    ij0 = sub2ind([ni,nj], i0, j0);
    n = (ij0-1) * KX;
    
    % Select the reference cast
    X0 = X((n+1:n+nk).');
    S0 = S(1:nk,ij0);
    T0 = T(1:nk,ij0);
    SppX0 = SppX(:,:,ij0);
    TppX0 = TppX(:,:,ij0);
    
    % Get reference salinity and temperature at the chosen location
    if isempty(s0) || isempty(t0)
        [s0,t0] = ppc_val2(X0, SppX0, TppX0, x0);
    end
    
    % Choose iso-value that will intersect (i0,j0,x0).
    D0 = eos(S0, T0, X0) - eos(s0, t0, X0);
    d0 = ppc_linterp(X0, D0, x0);

end

% Calculate 3D field for vertical interpolation
if eos(34.5,3,1000) > 1
    sortdir = 'ascend'; % eos is in-situ density, increasing with 3rd argument
else
    sortdir = 'descend'; % eos is specific volume, decreasing with 3rd argument
end
D = sort(eos(S, T, X) - eos(s0, t0, X), 1, sortdir);

% Get started with the discrete version (and linear interpolation)
x = interp1qn(d0, D, X);

% Solve non-linear root finding problem in each cast
[x, s, t] = delta_surf_vertsolve_mex(SppX, TppX, X, BotK, x, s0, t0, d0, OPTS.X_TOL);
