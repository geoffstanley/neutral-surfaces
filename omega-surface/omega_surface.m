function [x, s, t, diags] = omega_surface(S, T, X, x, OPTS)
% OMEGA_SURFACE  Create an omega surface, minimizing error from the neutral tangent plane.
%
%
% [x, s, t] = omega_surface(S, T, X, x, OPTS)
% returns the pressure or depth x, practical / Absolute salinity s, and
% potential / Conservative temperature t on an omega surface, initialized
% from an approximately neutral surface of (input) pressure or depth x, in
% an ocean whose practical / Absolute salinity and potential / Conservative
% temperature are S and T located at datasites where the pressure or depth
% is X.  An omega surface attempts to minimize the L2 norm of the
% neutrality error. The density or specific volume (either may be used) and
% its partial derivatives with respect to S and T are given by the
% functions eos.m and eos_s_t.m in MATLAB's path.  Algorithmic parameters
% are provided in OPTS (see "Options" below for further details).  For
% units, see "Equation of State" below.
%
% --- Input:
%  S [nk, ni, nj]: practical / Absolute Salinity
%  T [nk, ni, nj]: potential / Conservative Temperature
%  X [nk, ni, nj] or [nk, 1]: pressure or depth
%  x     [ni, nj]: pressure or depth on initial surface
% OPTS [struct]: options (see "Options" below)
%
%
% --- Output:
%  x [ni, nj]: pressure or depth on omega surface
%  s [ni, nj]: practical / Absolute salinity on omega surface
%  t [ni, nj]: potential / Conservative temperature on omega surface
%  diags [struct]: diagnostics such as clock time and norms of neutrality
%                  errors.  See code for info. Programmable as needed.
%
% Note: physical units of S, T, X, and x are determined by eos.m.
%
%
% --- Equation of State:
% The MATLAB path* must contain two functions, eos.m and eos_s_t.m. Both
% accept 3 inputs: S, T, and X. eos(S, T, X) returns the specific volume
% [m^3 kg^-1] or the in-situ density [kg m^-3]. eos_s_t(S, T, X) returns,
% as its two outputs, the partial derivatives of eos with respect to S and
% T.
% *Note: It is not sufficient to simply have these eos functions in the
% current working directory, because the compiled MEX functions will not be
% able to find them there.  They must be in the MATLAB path.  If they are
% in the current working directory, use `addpath(pwd)` to add the current
% working directory to the top of MATLAB's path.
%
% For a non-Boussinesq ocean, x and X are pressure [dbar].
%
% For a Boussinesq ocean, x and X are depth [m].  It is essential that
% these, like pressure, are positive and increasing down.
%
% Various equation of state functions are found in ../lib/eos/.  Simply
% copy the desired functions to another location in the MATLAB path (such
% as this directory) and rename them eos.m and eos_s_t.m.  Note, the
% Boussinesq equation of state is often (but not always) just the regular
% equation of state but using a hydrostatic pressure (10^-4 * grav * rho_c
% * z) where grav [m s^-2] is the gravitational acceleration, rho_c [kg
% m^-3] is the Boussinesq reference density, and z [m, positive] is the
% depth. In such a case, simply make new eos.m and eos_x.m functions that
% accept depth as the third input by modifying the original functions that
% take pressure; this involves hard-coding the gravitational acceleration
% and Boussinesq reference density into the function.  An example of a
% Boussinesq eos.m and eos_s_t.m are given for the densjmd95 equation of
% state, in ../lib/eos/eoscg_densjmd95_bsq.m and
% ../lib/eos/eoscg_densjmd95_bsq_s_t.m.  Finally, note that eos.m and
% eos_s_t.m must be compatible with MATLAB's code generation (codegen),
% which may entail eliminating input checks and/or expansion of input
% variables (MATLAB's automatic expansion now handles this).
%
%
% --- Options:
% OPTS is a struct containing the following fields. Those marked * are
% required. See ./private/omega_defaults.m for default values.
%   FILE_ID [1, 1]: 1 to write any output to MATLAB terminal, or a file
%       identifier as returned by fopen() to write to a file. Default: 1.
%   FIGS_SHOW [scalar]: true to show figures of specific volume adjustment
%       during computation. Default: false.
%   FINAL_ROW_VALUES [scalar]: value with which to fill the final row of
%       the sparse matrix for the purpose of maintaining the mean density
%       on the surface. (The other values in the matrix have values of 1.)
%       The matrix problem is perfectly (not over, not under) determined so
%       this value doesn't matter in theory.  In practice, large values
%       (e.g. 100) degrade the numerical solution.  Values in the range of
%       1e-4 to 1 were tested on 1x1deg OCCA data, and all work well.
%       Default: 1e-2.
%   INTEGRATING_FACTOR: use this [ni,nj] matrix of a pre-computed
%       integrating factor (b) to modify the weights of the matrix problem.
%       Use the regular weights when this is []. Default: [].
%   ITER_MAX [1, 1]: maximum number of iterations. Default: 10
%   ITER_START_WETTING [scalar]: Do wetting at this and subsequent
%       iterations. Set to +inf to disable wetting. When wetting, outcrops
%       and incrops of the surface that are adjacent to valid parts of the
%       surface are added back into the surface if they are connected by a
%       neutral trajectory. Default: 1.
%   INTERPFN [function handle]: vertical interpolation function, used to
%       evaluate SppX and TppX if those are not provided.  Default:
%       INTERPFN = @ppc_linterp.
%   MLX []: do not remove the mixed layer (default)
%   MLX [struct]: calculate the mixed layer using these parameters in mixed_layer().
%   MLX [ni, nj]: use a pre-computed mixed layer pressure [dbar] or depth [m]
%   REF_IJ [1, 2]: pixel indices for reference water column. The surface
%       will be pinned to the depth of the initial surface at this
%       location. Regions that are disconnected from the region containing
%       this location are instead controlled by maintaining their mean
%       depth.  Pass [] to do no pinning, and hold all regions to their
%       mean depth.
%   SppX [O, nk-1, ni, nj]: Coefficients for piecewise polynomials, whose
%       knots are X, that interpolate S as a function of X in each water
%       column.  E.g. SppX = ppc_linterp(X, S);
%   TppX [O, nk-1, ni, nj]: Coefficients for piecewise polynomials, whose
%       knots are X, that interpolate T as a function of X in each water
%       column.  E.g. TppX = ppc_linterp(X, T);
%   TOL_LRPD_L1 [scalar]: Error tolerance in Locally Referenced Potential Density [kg m^-3].
%       Iterations stop when the L1 norm of the LRPD change of the surface
%       is below this value. Even if eos gives specific volume, specify
%       this with units of density; it will be converted. Set to 0 to
%       ignore this stopping criterion.
%       Default: 10^-7 kg m^-3 (chosen to give an uncertainty in pressure
%       of roughly +/- 0.01 dbar.)
%   TOL_X_CHANGE_L2 [scalar]: Error tolerance in change of pressure [dbar].
%       Iterations stop when the L2 norm of the change in pressure on the
%       surface is below this value.  Set to 0 to ignore this stopping
%       criterion.
%       Default: inf
%   TOL_LSQR_REL [scalar]: Relative tolerance for LSQR. Default: 10^-6.
%   VERBOSE [scalar]: 0 for no output; 1 for summary of each iteration;
%                     2 for detailed information on each iteration.
%                     Default: 1.
% * WRAP [2 element vector]: determines which dimensions are treated
%       periodic. Set WRAP(1) to true when periodic in the 1st dimension of
%       x (ni); set WRAP(2) to true when periodic in the 2nd dimension of x
%       (nj).
%
%
% --- References:
% Klocker et al 2009: A new method of forming approximately neutral
%  surfaces,., Ocean Science, 5, 155-172.

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
%
% Acknowledgements: Adapted from 'analyze_surface' by Andreas Klocker, and
%                subsequently modified by Paul Barker and Trevor McDougall.

% --- Notes on the code:
% Upper case letters, e.g. X, denote 3D scalar fields [nk,ni,nj]
% Lower case letters, e.g. x, denote 2D scalar fields    [ni,nj]
% Developmental things are marked with a comment "DEV"

%iter_tic = tic; % timer for the initialization and for each iterative loop

%% Simple checks and preparations:
S = double(S);
T = double(T);
X = double(X);
x = double(x);

% Process mandatory options
assert(isstruct(OPTS) && isfield(OPTS, 'WRAP') && isvector(OPTS.WRAP) && length(OPTS.WRAP) == 2, 'OPTS.WRAP must be provided as a vector of length 2');

% Get size of 3D hydrography
[nk,ni,nj] = size(S);
nij = nj * ni;

% Setup anonymous functions:
lead1 = @(x) reshape(x, [1 size(x)]);

% Pre-calculate things for Breadth First Search
qu = zeros(nij, 1); % queue storing linear indices to pixels
A5 = grid_adjacency([ni,nj], 5, OPTS.WRAP); % all grid points that are adjacent to all grid points, using 5-connectivity
A4 = A5([1,2,4,5],:); % all grid points that are adjacent to all grid points, using 4-connectivity
%A = grid_adjacency([ni,nj], 4, OPTS.WRAP); % all grid points that are adjacent to all grid points, using 4-connectivity

% Number of bottles per cast. BotK(n) > 0 if and only if pixel n is ocean.
BotK = squeeze(sum(isfinite(S), 1));


msg1 = 'Initial surface has ................................................   mean(x) %.2f; mean(eos) %.6e\n';
msg2 = 'Iter %2d (%6.2f sec) log_10(|eps|_2) = %9.6f --> %9.6f by |phi''|_1 = %.6e; %03d regions; %4d casts freshly wet; |Delta_x|_2 = %.6e; mean(x) = %.2f; mean(eos) = %.6f\n';

%% Process OPTS

% Load default options, then override with any user-specified OPTS.
OPTS = catstruct(omega_defaults(), OPTS);

% Dereference for speed
FIGS_SHOW = OPTS.FIGS_SHOW;
FILE_ID = OPTS.FILE_ID;
ITER_MIN = OPTS.ITER_MIN;
ITER_MAX = OPTS.ITER_MAX;
ITER_START_WETTING = OPTS.ITER_START_WETTING;
VERBOSE = OPTS.VERBOSE;
WRAP = OPTS.WRAP;
REF_IJ = OPTS.REF_IJ;
POISSON = OPTS.POISSON; 

PIN = ~isempty(REF_IJ);
if PIN
  I0 = sub2ind([ni, nj], REF_IJ(1), REF_IJ(2));
else
  I0 = [];
end

eos0 = eos(34.5, 3, 1000);
if eos0 > 1
  TOL_LRPD_L1 = OPTS.TOL_LRPD_L1;          % Density tolerance [kg m^-3]
else
  TOL_LRPD_L1 = OPTS.TOL_LRPD_L1 * 1000^2; % Specific volume tolerance [m^3 kg^-1]
end

TOL_X_CHANGE_L2 = OPTS.TOL_X_CHANGE_L2;

TOL_X_UPDATE = OPTS.TOL_X_UPDATE; % tolerance in pressure [dbar] during vertical solve.

DIAGS = (VERBOSE > 0) || (nargout > 3);

DX = OPTS.DX;
DY = OPTS.DY;

if ~POISSON
  % Set up indices used when building the matrix mat
  shift_im1 = @(F) circshift(F, [+1, 0]); % if x is [lat x lon], this shifts south
  shift_jm1 = @(F) circshift(F, [0, +1]); % if x is [lat x lon], this shifts west
  idx = reshape(1:nij, ni, nj); % Linear index to each pixel on the full grid (ocean or not)
  idx_im1 = shift_im1(idx);     % Linear index to the pixel above (^)   the local pixel
  idx_jm1 = shift_jm1(idx);     % Linear index to the pixel left (<) of the local pixel
end

%% Just in time code generation
Xvec = isvector(X);
ni_ = max(ni, 2048); % using variable size code generation and avoiding
nj_ = max(nj, 2048); % recompiling all the time
omega_vertsolve_codegen(nk, ni_, nj_, Xvec, OPTS);
bfs_conncomp_codegen(nk, ni_, nj_, Xvec, false, OPTS);
if ITER_START_WETTING <= ITER_MAX
  bfs_wet_codegen(nk, ni_, nj_, Xvec, OPTS);
end

%% Get MLX: the pressure or depth of the mixed layer
if ITER_MAX > 1
  if isempty(OPTS.MLX)
    % Do not remove the mixed layer
    MLX = [];
  elseif isstruct(OPTS.MLX)
    MLX = mixed_layer(S, T, X, OPTS.MLX);
  else
    % Use a pre-computed mixed layer
    MLX = OPTS.MLX;
  end
end

%% Interpolate S and T casts onto surface
[~, K, N] = size(OPTS.SppX);
if K > 0
  assert(K == nk-1 && N == nij, 'size(SppX) should be [O, nk-1, ni, nj] == [?, %d, %d, %d]', nk-1, ni, nj);
  SppX = OPTS.SppX;
else
  SppX = OPTS.INTERPFN(X, S);
end

[~, K, N] = size(OPTS.TppX);
if K > 0
  assert(K == nk-1 && N == nij, 'size(TppX) should be [O, nk-1, ni, nj] == [?, %d, %d, %d]', nk-1, ni, nj);
  TppX = OPTS.TppX;
else
  TppX = OPTS.INTERPFN(X, T);
end

[s, t] = ppc_val2(X, SppX, TppX, lead1(x));


%% Prepare diagnostics
if DIAGS
  
  diags = struct();
  diags.mean_x      = nan(ITER_MAX + 1, 1);
  diags.mean_eos    = nan(ITER_MAX + 1, 1);
  diags.epsL2       = nan(ITER_MAX + 1, 1);
  
  diags.phiL1         = nan(ITER_MAX, 1);
  diags.x_change_L1   = nan(ITER_MAX, 1);
  diags.x_change_L2   = nan(ITER_MAX, 1);
  diags.x_change_Linf = nan(ITER_MAX, 1);
  diags.num_regions   = nan(ITER_MAX, 1);
  diags.freshly_wet   = nan(ITER_MAX, 1);
  diags.relaxation    = nan(ITER_MAX, 1);
  diags.epsL2_before  = nan(ITER_MAX, 1);
  diags.clocktime     = nan(ITER_MAX, 1);
  
  diags.timer_wetting   = nan(ITER_MAX, 1);
  diags.timer_ntperrors = nan(ITER_MAX, 1);
  diags.timer_solver    = nan(ITER_MAX, 1);
  diags.timer_update    = nan(ITER_MAX, 1);
  
  % Diagnostics about state BEFORE this (first) iteration
  [eps_i, eps_j] = ntp_errors(s, t, x, DX, DY, true, false, WRAP);
  epsL2 = nanrms([eps_i(:); eps_j(:)]);
  mean_x = nanmean(x(:));
  mean_eos = nanmean(eos(s(:), t(:), x(:)));
  diags.epsL2(1) = epsL2;
  diags.mean_x(1) = mean_x;
  diags.mean_eos(1) = mean_eos;
  %diags.clocktime(1) = toc(iter_tic);
  
  if VERBOSE > 0
    fprintf(FILE_ID, msg1, mean_x, mean_eos);
  end
  
end

%% Begin iterations
for iter = 1 : ITER_MAX
  iter_tic = tic;
  
  % --- Remove the Mixed Layer
  % But keep it for the first iteration, which may be initialized from a
  % not very neutral surface
  if iter > 1 && ~isempty(MLX)
    x(x < MLX) = nan;
  end
  
  % --- Wetting via Breadth First Search
  mytic = tic;
  if iter >= ITER_START_WETTING
    [s, t, x, freshly_wet, qu] = bfs_wet_mex(SppX, TppX, X, s, t, x, TOL_X_UPDATE, A4, BotK, qu);
  else
    freshly_wet = 0;
  end
  
  % --- Accumulate regions via Breadth First Search
  [qu, qts, ncc, ~, L] = bfs_conncomp_all_mex(isfinite(x), A4, [], qu);
  
  if DIAGS
    timer_wetting = toc(mytic);
  end
  
  
  
  
  % --- Build the compact representation of the matrix optimization problem
  mytic = tic;
  if POISSON
    [G, H, epsL2] = linprob_dens(x, s, t, [], [], 5);
  else
    [eps_i, eps_j] = ntp_errors(s, t, x, 1, 1, true, false, WRAP);
    epsL2 = nanrms([eps_i(:); eps_j(:)]);
  end
  if DIAGS
    timer_ntperrors = toc(mytic);
  end
  
  % --- Build and solve sparse matrix problem
  mytic = tic;
  if POISSON
    [phi, G, H] = omega_matsolve_poisson(G, H, A5, L, I0, qu, qts);
  else
    phi = omega_matsolve_grad(eps_i, eps_j, idx_im1, idx_jm1, L, I0, qu, qts);
  end
  if DIAGS
    timer_solver = toc(mytic);
  end
  
  
  % --- Update the surface
  mytic = tic();
  x_old = x;      % Record old surface for pinning or diagnostic purposes.
  [x, s, t] = omega_vertsolve_mex(SppX, TppX, X, BotK, s, t, x, TOL_X_UPDATE, phi);
  
  if PIN && ~isnan(phi(I0))
    % Force x(i0,j0) to stay constant at the reference column,
    % identically. This avoids any intolerance from the vertical solver.
    x(I0) = x_old(I0);
  end
  
  if DIAGS
    timer_update = toc(mytic);
  end
  
  % --- Closing Remarks
  phiL1 = nanmean(abs(phi(:)));
  if DIAGS || isfinite(TOL_X_CHANGE_L2)
    x_change = x - x_old;
    x_change_L2   = nanrms(x_change(:));
  end
  
  if DIAGS
    
    diags.clocktime(iter)   = toc(iter_tic);
    
    x_change_L1   = nanmean(abs(x_change(:)));
    x_change_Linf = max(abs(x_change(:)));
    
    [eps_i, eps_j] = ntp_errors(s, t, x, DX, DY, true, false, WRAP);
    epsL2_after = nanrms([eps_i(:); eps_j(:)]);
    
    mean_x = nanmean(x(:));
    mean_eos = nanmean(eos(s(:),t(:),x(:)));
    
    % Diagnostics about state BEFORE this iteration, but AFTER wetting
    diags.epsL2_before(iter) = epsL2; % <-- Taking advantage that epsilon had to be computed anyways
    
    % Diagnostics about what THIS iteration did
    diags.phiL1(iter)         = phiL1;
    diags.x_change_L1(iter)   = x_change_L1;
    diags.x_change_L2(iter)   = x_change_L2;
    diags.x_change_Linf(iter) = x_change_Linf;
    diags.num_regions(iter)   = ncc;
    diags.freshly_wet(iter)   = freshly_wet;
    
    diags.timer_wetting(iter)   = timer_wetting;
    diags.timer_ntperrors(iter) = timer_ntperrors;
    diags.timer_solver(iter)    = timer_solver;
    diags.timer_update(iter)    = timer_update;
    
    % Diagnostics about the state AFTER this iteration
    diags.epsL2(iter+1)     = epsL2_after;
    diags.mean_x(iter+1)    = mean_x;
    diags.mean_eos(iter+1)  = mean_eos;
    
    
    if VERBOSE > 0
      fprintf(FILE_ID, msg2, iter, diags.clocktime(iter), log10(epsL2), log10(epsL2_after), phiL1, ncc, freshly_wet, x_change_L2, mean_x, mean_eos);
    end
  end
  
  % --- Show Figures
  if FIGS_SHOW
    if mod(iter,12) == 1
      hf = figure('Position',[20, 20, 1000, 800], 'Name', 'phi''');
    end
    ax = subplot(4, 3, mod(iter-1,12)+1, 'Parent', hf );
    if nj > ni % x, and phi, are probably lon x lat
      %imagesc(ax, phi)
      imagesc(ax, x_change)
    else       % x, and phi, are probably lat x lon
      %imagesc(ax, phi.')
      imagesc(ax, x_change')
    end
    ax.YDir = 'normal';
    shading('flat')
    %title(sprintf('%d: |\\phi''|_1 = %.2e', iter, phiL1), 'fontsize',10);
    %caxis(prctile(phi(:), [1 99]));
    title(sprintf('%d: |\\Delta x''|_1 = %.2e', iter, x_change_L1), 'fontsize',10);
    caxis(prctile(x_change(:), [1 99]));
    colorbar();
    drawnow()
  end
  
  % --- Check for convergence
  if (phiL1 < TOL_LRPD_L1 || x_change_L2 < TOL_X_CHANGE_L2) && iter >= ITER_MIN
    break
  end
  
end



if DIAGS
  % Trim output
  fields = fieldnames(diags);
  for i = 1 : length(fields)
    f = fields{i};
    diags.(f) = diags.(f)( isfinite( diags.(f) ));
  end
end