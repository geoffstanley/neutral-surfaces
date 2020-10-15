function [x,s,t,d0,s_ref,t_ref,diags] = delta_surf(S, T, X, s_ref, t_ref, var, OPTS)
%DELTASURF Specific volume anomaly surface by nonlinear solution in each water column.
%
%
% [x,s,t] = delta_surf(S, T, X, s_ref, t_ref, d0)
% finds pressure or depth x -- and its salinity s and temperature t -- of
% the isosurface d0 of delta = eos(S,T,X) - eos(s_ref,t_ref,X) with reference
% practical / Absolute salinity s_ref and reference potential / Conservative
% temperature t_ref, in an ocean with practical / Absolute salinity S and
% potential / Conservative temperature T at data sites where the pressure
% or depth is X.  The equation of state is given by eos.m in MATLAB's path,
% which accepts S, T, and X as its 3 inputs.  For a non-Boussinesq ocean, x
% and X are pressure [dbar] and eos gives the specific volume.  For a
% Boussinesq ocean, x and X are depth [m] positive and increasing down, and
% eos gives the in-situ density.
%
% [x, s, t, d0] = delta_surf(S, T, X, s_ref, t_ref, [i0, j0, x0])
% as above but finds the delta isosurface, delta = d0, that intersects the
% reference cast at grid indices (i0,j0) at pressure or depth x0.
%
% [x, s, t, d0, s_ref, t_ref] = delta_surf(S, T, X, [], [], [i0, j0, x0])
% as above but also chooses the reference practical / Absolute salinity s_ref
% and reference potential / Conservative temperature t_ref by interpolating S
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
% s_ref [1, 1] or []: reference S
% t_ref [1, 1] or []: reference T
% var [1, 1] or [1, 3]: isovalue of delta, or location to intersect
% OPTS [struct]: options (see below)
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: physical units for S, T, X, x, d0, s_ref, t_ref, x0 are determined by eos.m.
%
% Note: X must increase monotonically along its first dimension.
%
%
% --- Output:
% x [ni, nj]: pressure [dbar] or depth [m] of the surface
% s [ni, nj]: practical / Absolute salinity on the surface
% t [ni, nj]: potential / Conservative temperature the surface
% d0 [1, 1]: isovalue of the delta surface
% s_ref [1, 1]: reference S that defines delta
% t_ref [1, 1]: reference T that defines delta
% diags [struct]: diagnostics of the solution and computation time
%
%
% --- Options:
% OPTS is a struct containing a subset of the following fields.
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
%   DIST1_iJ [ni, nj]: Distance [m] in 1st dimension centred at (I-1/2, J)
%   DIST2_Ij [ni, nj]: Distance [m] in 2nd dimension centred at (I, J-1/2)
%   DIST2_iJ [ni, nj]: Distance [m] in 2nd dimension centred at (I-1/2, J)
%   DIST1_Ij [ni, nj]: Distance [m] in 1st dimension centred at (I, J-1/2)

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
DEFS.TOL_X_UPDATE = 1e-4;  % error tolerance in the vertical dimension
DEFS.INTERPFN = @ppc_linterp;  % linear interpolation in the vertical
DEFS.FILE_ID = 1;  % standard output to MATLAB terminal
DEFS.VERBOSE = 1;  % show a moderate level of information
DEFS.DIST1_iJ = 1;   % Distance [m] in 1st dimension centred at (I-1/2, J)
DEFS.DIST2_Ij = 1;   % Distance [m] in 2nd dimension centred at (I, J-1/2)
DEFS.DIST2_iJ = 1;   % Distance [m] in 2nd dimension centred at (I-1/2, J)
DEFS.DIST1_Ij = 1;   % Distance [m] in 1st dimension centred at (I, J-1/2)

% Override any options with user-specified OPTS
if nargin < 7 || isempty(OPTS)
  OPTS = DEFS;
else
  OPTS = catstruct(DEFS, OPTS); 
end

DIAGS = nargout >= 7;


% Run codegen to create MEX function handling the main computation
ni_ = max(ni, 4096); % using variable size code generation and avoiding
nj_ = max(nj, 4096); % recompiling all the time
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
    % s_ref and t_ref should be provided
    assert(isscalar(s_ref) && isscalar(t_ref), 's_ref and t_ref must be provided as scalars if d0 is provided.');
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
    SppX0 = SppX(:,:,ij0);
    TppX0 = TppX(:,:,ij0);
    
    % Evaluate salinity and temperature at the chosen location
    [s0,t0] = ppc_val2(X0, SppX0, TppX0, x0);
    if isempty(s_ref) || isempty(t_ref)
        s_ref = s0;
        t_ref = t0;
        d0 = 0; % iso-value that will intersect (i0,j0,x0).
    else
        d0 = eos(s0, t0, x0) - eos(s_ref, t_ref, x0); % iso-value that will intersect (i0,j0,x0).
    end
    
end

% Calculate 3D field for vertical interpolation
if eos(34.5,3,1000) > 1
    sortdir = 'ascend'; % eos is in-situ density, increasing with 3rd argument
else
    sortdir = 'descend'; % eos is specific volume, decreasing with 3rd argument
end
D = sort(eos(S, T, X) - eos(s_ref, t_ref, X), 1, sortdir);

% Get started with the discrete version (and linear interpolation)
x = interp1qn(d0, D, X);

% Start timer after all the setup has been done.
iter_tic = tic();

% Solve non-linear root finding problem in each cast
[x, s, t] = delta_surf_vertsolve_mex(SppX, TppX, X, BotK, x, s_ref, t_ref, d0, OPTS.TOL_X_UPDATE);

if DIAGS
    diags = struct();
    diags.clocktime = toc(iter_tic);
    [epsL2, epsL1] = eps_norms(s, t, x, true, OPTS.WRAP, OPTS.DIST1_iJ, OPTS.DIST2_Ij, OPTS.DIST2_iJ, OPTS.DIST1_Ij);
    diags.epsL1 = epsL1;
    diags.epsL2 = epsL2;
end

