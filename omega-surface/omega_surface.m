function [x, s, t] = omega_surface(S, T, X, x, OPTS)
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
%
% Note: physical units of S, T, X, and x are determined by eos.m.
%
%
% --- Equation of State:
% The MATLAB path must contain two functions, eos.m and eos_s_t.m. Both
% accept 3 inputs: S, T, and X. eos(S, T, X) returns the specific volume
% [m^3 kg^-1] or the in-situ density [kg m^-3]. eos_s_t(S, T, X) returns,
% as its two outputs, the partial derivatives of eos with respect to S and
% T.
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
%       As this is lowered, the solution to the matrix equation cares more
%       about neutrality and less about maintaining the mean density.
%       Default: 1e-2 (chosen empirically from tests on 1x1deg OCCA data).
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
%   SppX [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose 
%       knots are X that interpolate S as a function of X in each water 
%       column.  E.g. SppX = ppc_linterp(X, S);
%   TppX [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose 
%       knots are X that interpolate T as a function of X in each water 
%       column.  E.g. TppX = ppc_linterp(X, T);
%   TOL_DENS [scalar]: Target error tolerance in density [kg m^-3].
%       Iterations stop when the L1 norm of the density change of the
%       surface is below this value. Even if eos gives specific volume,
%       specify this with units of density; it will be converted.
%       Default: 10^-7 kg m^-3 (chosen to give an uncertainty in pressure
%       of roughly +/- 0.01 dbar.)
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
% Version   : 2.0.0
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
shift_im1 = @(F) circshift(F, [+1, 0]); % if x is [lat x lon], this shifts south
shift_jm1 = @(F) circshift(F, [0, +1]); % if x is [lat x lon], this shifts west
flat = @(x) x(:);
lead1 = @(x) reshape(x, [1 size(x)]);

% Set up indices used when building the matrix mat
idx = reshape(1:nij, ni, nj); % Linear index to each pixel on the full grid (ocean or not)
idx_im1 = shift_im1(idx);     % Linear index to the pixel above (^)   the local pixel
idx_jm1 = shift_jm1(idx);     % Linear index to the pixel left (<) of the local pixel

% The density or specific volume change in each water column
phi = nan(ni,nj);

% Pre-calculate things for Breadth First Search
qu = zeros(nij, 1); % queue storing linear indices to pixels
A = grid_adjacency([ni,nj], OPTS.WRAP, 4); % all grid points that are adjacent to all grid points, using 4-connectivity

% Number of bottles per cast. BotK(n) > 0 if and only if pixel n is ocean.
BotK = squeeze(sum(uint16(isfinite(S)), 1, 'native'));



%% Process OPTS

% Load default options, then override with any user-specified OPTS.
OPTS = catstruct(omega_defaults(), OPTS);

USE_INTEGRATING_FACTOR = ~isempty(OPTS.INTEGRATING_FACTOR);

% Hard-coded parameters:
drdp = 4e-03; % Approximate derivative of in-situ density w.r.t. pressure [kg m^-3 dbar^-1]
X_TOL = OPTS.TOL_DENS / (drdp * 2); % tolerance in pressure [dbar] during vertical solve.  Factor of 2 for good measure
eos0 = eos(34.5, 3, 1000);
if eos0 > 1
    PHI_TOL = OPTS.TOL_DENS;          % Density tolerance [kg m^-3]
    msg1 = 'Initial surface has ................................................ mean(x) = %.4f, mean(dens) = %.2f\n';
    msg2 = 'Iter %02d had |eps|_2 = %.6e to start; done by %.2f sec with mean(x) = %.2f, mean(dens) = %.4f, |phi''|_1 = %.6e in %d regions, %d casts freshly wet\n';
else
    PHI_TOL = OPTS.TOL_DENS * eos0^2; % Specific volume tolerance [m^3 kg^-1]
    msg1 = 'Initial surface has ................................................   mean(x) %.4f; mean(specvol) %.6e\n';
    msg2 = 'Iter %02d had |eps|_2 = %.6e to start; done by %.2f sec with mean(x) = %.2f, mean(specvol) = %.6e, |phi''|_1 = %.6e in %d regions, %d casts freshly wet\n';
end

% Experimenting with decreasing tolerance:
% tol_rel = 1e-4; % Initial relative tolerance  DEV
% fac = 100; % Final relative tolerance will be smaller than initial by this factor  DEV
% tol_geofac = fac^(-1/OPTS.ITER_MAX);   % geometric factor by which tol_rel decreases on each iteration.  DEV

%% Just in time code generation
mytic = tic;
Xvec = isvector(X);
omega_vertsolve_codegen(nk, ni, nj, Xvec, OPTS);
bfs_conncomp_codegen(nk, ni, nj, Xvec, false, OPTS);
if OPTS.ITER_START_WETTING <= OPTS.ITER_MAX
    bfs_wet_codegen(nk, ni, nj, Xvec, false, OPTS);
end
if OPTS.VERBOSE > 1
    fprintf(OPTS.FILE_ID, '%.4f sec for Code Generation\n', toc(mytic));
end

%% Get MLX: the pressure or depth of the mixed layer
if OPTS.ITER_MAX > 1
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
mytic = tic;
SppX = OPTS.INTERPFN(X, S);
TppX = OPTS.INTERPFN(X, T);
[s,t] = ppc_val2(X,SppX,TppX,lead1(x));

if OPTS.VERBOSE > 0
    fprintf(OPTS.FILE_ID, msg1, nanmean(x(:)), nanmean(eos(s(:),t(:),x(:))));
    if OPTS.VERBOSE > 1
        fprintf(OPTS.FILE_ID, '  %.4f sec to interpolate S and T \n', toc(mytic));
    end
end

%% Begin iterations
total_tic = tic;
for iter = 1 : OPTS.ITER_MAX
    
    % --- Remove the Mixed Layer
    % But keep it for the first iteration, which may be initialized from a
    % not very neutral surface
    if iter > 1 && ~isempty(MLX)
        x(x < MLX) = nan;
    end
    
    % --- Wetting via Breadth First Search
    if iter >= OPTS.ITER_START_WETTING
        mytic = tic;
        [qu, qt, s, t, x, freshly_wet] = bfs_wet_all_mex(SppX, TppX, X, s, t, x, X_TOL, A, BotK, [], qu);
        wet = false(ni,nj);
        wet(qu(1:qt)) = true;
        if OPTS.VERBOSE > 1
            fprintf(OPTS.FILE_ID, '  %.4f sec to wet \n', toc(mytic));
        end
    else
        freshly_wet = 0; 
        wet = isfinite(x);
    end
    
    % --- Accumulate regions via Breadth First Search
    mytic = tic;
    [qu, qts, ncc] = bfs_conncomp_all_mex(wet, A, [], qu);
    if OPTS.VERBOSE > 1
        fprintf(OPTS.FILE_ID, '  %.4f sec to build regions \n', toc(mytic));
    end
    
    % --- Calculate specific volume neutrality errors on the current surface
    % Note that these are errors, not gradients, i.e. the density differences
    % have not been divided by distance.
    mytic = tic;
    if USE_INTEGRATING_FACTOR
        [no_b_eps_i, no_b_eps_j, eps_i, eps_j] = ntp_error_s_t_withb(s, t, x, OPTS.WRAP, OPTS.INTEGRATING_FACTOR);
    else
        [eps_i, eps_j] = ntp_error_s_t_withb(s, t, x, OPTS.WRAP);
    end
    if OPTS.VERBOSE > 1
        fprintf(OPTS.FILE_ID, '  %.4f sec to calculate neutrality errors \n', toc(mytic));
    end
    
    % --- Build and solve sparse matrix problem
    mytic = tic;
    phi(:) = nan;
    idx(:) = 0;    % Overwrite with linear indices to regions
    phiL1 = 0;     % For accumulating the L1 norm of phi
    nwc_total = 0; % For accumulating the L1 norm of phi
    epsL2 = 0;     % For accumulating the L2 norm of the epsilon error
    neq_total = 0; % For accumulating the L2 norm of the epsilon error
    for icc = 1 : ncc % icc = index of region
        
        
        % Collect and sort linear indices to all pixels in this region
        I = sort(qu(qts(icc) : qts(icc+1)-1));  % sorting here makes mat better structured; overall speedup.
        
        nwc = length(I);  % Number of water columns
        if nwc <= 1  % There are definitely no equations to solve
            phi(I) = 0; % Leave this isolated pixel at current pressure
            continue
        end
        
        % Label the water columns in this region alone by 1, 2, ...
        % No need to reset other elements of idx to 0.
        idx(I) = 1:nwc;
        
        % Build the RHS of the matrix problem, containing neutrality errors
        rhs = [eps_i(I); eps_j(I)];
        good_eq = ~isnan(rhs);
        rhs = [rhs(good_eq); 0]; % Ignore equations where the east or the north water column is invalid.
        neq = length(rhs) - 1; % Number of equations, excluding density conserving equation.  Note neq > 0 is guaranteed, because bfs_conncomp used 4-connectivity
        
        % Build indices for the columns of the sparse matrix
        jj = [[idx_im1(I); idx_jm1(I)], repmat(I, 2, 1)];
        jj = idx(jj);  % Remap global indices to local indices for this region, numbered 1, 2, ...
        jj = [flat(jj(good_eq,:)); (1:nwc).']; % Augment with the specific volume conserving equation
        
        % Build indices for the rows of the sparse matrix
        ii = [flat(repelem((1:neq).', 1, 2)); repelem(neq+1, nwc, 1)];
        
        % Build the values of the sparse matrix
        vv = [flat(repelem([-1, 1], neq, 1)); repelem(OPTS.FINAL_ROW_VALUES, nwc,1)];
        
        
        % Build the sparse matrix, with neq+1 rows and nwc columns
        mat = sparse( ii, jj, vv, neq+1, nwc );
        
        % Accumulate sum of squares for L2 norm
        if OPTS.VERBOSE > 0
            if USE_INTEGRATING_FACTOR
                B_no_b = [no_b_eps_i(I); no_b_eps_j(I)];
                B_no_b(isnan(B_no_b)) = 0;
                epsL2 = epsL2 + sum(B_no_b .* B_no_b);
            else
                epsL2 = epsL2 + sum(rhs .* rhs);
            end
            neq_total = neq_total + neq;
        end
        
        
        
        % Solve the LSQR problem
        [sol,flag] = omega_lsqr(mat, rhs, OPTS.TOL_LSQR_REL);
        if flag > 0
            warning('omega_surface:lsqr did not converge on iter %d, region %d', iter, icc);
        end
        
        % Save solution into phi and accumulate statistics for phi
        phi(I) = sol;
        phiL1 = phiL1 + sum(abs(sol));
        nwc_total = nwc_total + nwc;
    end
    phiL1 = phiL1 / nwc_total;
    if USE_INTEGRATING_FACTOR
        phi = phi ./ OPTS.INTEGRATING_FACTOR;
    end
    if OPTS.VERBOSE > 1
        fprintf(OPTS.FILE_ID, '  %.4f sec to build and solve matrices \n', toc(mytic));
    end
    
    % Force mean(phi) = 0 exactly. The LSQR solution only tries to keep
    % mean(phi) near zero in a least-squares sense, while also trying to
    % make the neutrality errors zero. When FINAL_ROW_VALUES is badly
    % chosen, |phi|_1 can wobble as iterations proceed, but forcing
    % mean(phi) = 0 seems to help |phi|_1 to decrease monotonically.
    phi = phi - nanmean(phi(:));
    
    % --- Update the surface
    mytic = tic();
    [x, s, t] = omega_vertsolve_mex(SppX, TppX, X, BotK, s, t, x, X_TOL, phi);
    if OPTS.VERBOSE > 1
        fprintf(OPTS.FILE_ID, '  %.4f sec for global vertical bisection \n', toc(mytic));
    end
    
    % --- Show Figures
    if OPTS.FIGS_SHOW
        if mod(iter,12) == 1
            hf = figure('Position',[20, 20, 1000, 800], 'Name', 'phi''');
        end
        ax = subplot(4, 3, mod(iter-1,12)+1, 'Parent', hf );
        if nj > ni % x, and phi, are probably lon x lat
            pcolor(ax, phi)
        else       % x, and phi, are probably lat x lon
            pcolor(ax, phi.')
        end
        shading('flat')
        title(sprintf('%d: |\\phi''|_1 = %.2e', iter, phiL1), 'fontsize',10);
        caxis(prctile(phi(:), [1 99]));
        colorbar();
        drawnow()
    end
    
    % --- Closing Remarks
    if OPTS.VERBOSE > 0
        epsL2 = sqrt(epsL2 / neq_total);
        fprintf(OPTS.FILE_ID, msg2, iter, epsL2, toc(total_tic), nanmean(x(:)), nanmean(eos(s(:),t(:),x(:))), phiL1, ncc, freshly_wet);
    end
    
    % --- Check for convergence
    if phiL1 < PHI_TOL
        break
    end
    
    % Reduce tolerance for next iteration
    %tol_rel = tol_rel * tol_geofac;  % DEV
    
end



