function [p, s, t, diags] = omega_surface(S, T, P, p, ref_cast, OPTS)
% OMEGA_SURFACE  Create an omega surface, minimizing error from the neutral tangent plane.
%
%
% [p, s, t] = omega_surface(S, T, P, p, ref_cast, OPTS)
% returns the pressure (or depth) p, practical / Absolute salinity s, and
% potential / Conservative temperature t on an omega surface, initialized
% from an approximately neutral surface of (input) pressure (or depth) p,
% in an ocean whose practical / Absolute salinity and potential /
% Conservative temperature are S and T located at datasites where the
% pressure (or depth) is P.  The pressure of the omega surface is pinned,
% unchanging through the iterations, at the reference cast indexed by
% ref_cast.  An omega surface attempts to minimize the L2 norm of the
% neutrality error. The density or specific volume (either may be used) and
% its partial derivatives with respect to S and T are given by the
% functions eos.m and eos_s_t.m in MATLAB's path. Algorithmic parameters
% are provided in OPTS (see "Options" below for further details).  For
% units, see "Equation of State" below.
%
%
% --- Input:
%  S [nk, ni, nj]: practical / Absolute Salinity
%  T [nk, ni, nj]: potential / Conservative Temperature
%  P [nk, ni, nj] or [nk, 1]: pressure (or depth)
%  p     [ni, nj]: pressure (or depth) on initial surface
%  ref_cast [1, 1] or [2, 1] : linear index or 2D index to the reference cast
%  OPTS [struct]: options (see "Options" below)
%
%
% --- Output:
%  p [ni, nj]: pressure (or depth) on omega surface
%  s [ni, nj]: practical / Absolute salinity on omega surface
%  t [ni, nj]: potential / Conservative temperature on omega surface
%  diags [struct]: diagnostics such as clock time and norms of neutrality
%                  errors.  See code for info. Programmable as needed.
%
% Note: physical units of S, T, P, and p are determined by eos.m.
%
%
% --- Equation of State:
% The MATLAB path* must contain two functions, eos.m and eos_s_t.m. Both
% accept 3 inputs: S, T, and P. eos(S, T, P) returns the specific volume
% [m^3 kg^-1] or the in-situ density [kg m^-3]. eos_s_t(S, T, P) returns,
% as its two outputs, the partial derivatives of eos with respect to S and
% T.
% *Note: It is not sufficient to simply have these eos functions in the
% current working directory, because the compiled MEX functions will not be
% able to find them there.  They must be in the MATLAB path.  If they are
% in the current working directory, use `addpath(pwd)` to add the current
% working directory to the top of MATLAB's path.
%
% For a non-Boussinesq ocean, p and P are pressure [dbar].
%
% For a Boussinesq ocean, p and P are actually depth [m].  It is essential
% that these, like pressure, are positive and increasing down.
%
% Various equation of state functions are found in ../lib/eos/.  Simply
% copy the desired functions to another location in the MATLAB path (such
% as this directory) and rename them eos.m and eos_s_t.m.  Note, the
% Boussinesq equation of state is often (but not always) just the regular
% equation of state but using a hydrostatic pressure (10^-4 * grav * rho_c
% * z) where grav [m s^-2] is the gravitational acceleration, rho_c [kg
% m^-3] is the Boussinesq reference density, and z [m, positive] is the
% depth. In such a case, simply make new eos.m and eos_p.m functions that
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
% OPTS is a struct containing the following fields. 
%   FILE_ID [1, 1]: 1 to write any output to MATLAB terminal, or a file
%       identifier as returned by fopen() to write to a file. Default: 1.
%   FIGS_SHOW [scalar]: true to show figures of specific volume adjustment
%       during computation. Default: false.
%   FINAL_ROW_VALUES [scalar]: value with which to fill the final row of
%       the sparse matrix for the purpose of selecting a unique solution
%       for th density perturbation. This value doesn't matter in theory,
%       but in practice, excessively large or small values may degrade the
%       numerical solution.  Values in the range of 1e-4 to 1 were tested
%       on 1x1deg OCCA data, and all work well. Default: 1e-2.
%   ITER_MAX [1, 1]: maximum number of iterations. Default: 10
%   ITER_START_WETTING [scalar]: Start wetting on iterations that are
%       >= ITER_START_WETTING. To disable wetting, set to +inf. Default: 1.
%   ITER_STOP_WETTING [scalar]: Do wetting for iterations that are
%       <= ITER_STOP_WETTING. To disable wetting, set to 0. Default: 5.
%   INTERPFN [function handle]: vertical interpolation function, used to
%       evaluate Sppc and Tppc if those are not provided.  Default:
%       INTERPFN = @ppc_linterp.
%   ML []: do not remove the mixed layer (default)
%   ML [struct]: calculate the mixed layer using these parameters in mixed_layer().
%   ML [ni, nj]: use a pre-computed mixed layer pressure [dbar] or depth [m]
%   Sppc [O, nk-1, ni, nj]: Coefficients for piecewise polynomials, whose
%       knots are at P, that interpolate S as a function of P in each water
%       column.  E.g. Sppc = ppc_linterp(P, S);
%   Tppc [O, nk-1, ni, nj]: Coefficients for piecewise polynomials, whose
%       knots are at P, that interpolate T as a function of P in each water
%       column.  E.g. Tppc = ppc_linterp(P, T);
%   TOL_LRPD_L1 [scalar]: Error tolerance in Locally Referenced Potential Density [kg m^-3].
%       Iterations stop when the Mean Absolute Value of the LRPD change of the surface
%       is below this value. Even if eos gives specific volume, specify
%       this with units of density; it will be converted. Set to 0 to
%       ignore this stopping criterion.
%       Default: 10^-7 kg m^-3 (chosen to give an uncertainty in pressure
%       of roughly +/- 0.01 dbar.)
%   TOL_P_CHANGE_L2 [scalar]: Error tolerance in change of pressure [dbar].
%       Iterations stop when the root-mean-square of the change in pressure on the
%       surface is below this value.  Set to 0 to ignore this stopping
%       criterion.
%       Default: inf
%   TOL_LSQR_REL [scalar]: Relative tolerance for LSQR. Default: 10^-6.
%   VERBOSE [scalar]: 0 for no output; 1 for summary of each iteration;
%                     2 for detailed information on each iteration.
%                     Default: 1.
%   WRAP [2 element array]: logical array.  WRAP(i) is true iff the domain 
%       is periodic in the i'th lateral dimension.  
%
%
% --- References:
% Klocker, McDougall, Jackett 2009: A new method of forming approximately
%  neutral surfaces, Ocean Science, 5, 155-172.
%
% Stanley, McDougall, Barker 2021: Algorithmic improvements to finding
%  approximately neutral surfaces, Journal of Advances in Earth System
%  Modelling, 13(5).


% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com

% Acknowledgements: Adapted from 'analyze_surface' by Andreas Klocker, and
%                subsequently modified by Paul Barker and Trevor McDougall.

% --- Notes on the code:
% Upper case letters, e.g. S, denote 3D scalar fields [nk,ni,nj]
% Lower case letters, e.g. s, denote 2D scalar fields    [ni,nj]
% Developmental things are marked with a comment "DEV"

%% Simple checks and preparations:
S = double(S);
T = double(T);
P = double(P);
p = double(p);

% Get size of 3D hydrography
[nk,ni,nj] = size(S);
nij = nj * ni;

% Setup anonymous functions:
lead1 = @(x) reshape(x, [1 size(x)]);
autoexp = @(x) repmat(x, ni / size(x,1), nj / size(x,2)); % automatic expansion to [ni,nj]

msg1 = 'Initial surface has  log_10(|eps|_2) == %9.6f .................. \n';
msg2 = 'Iter %2d (%6.2f sec) log_10(|eps|_2) == %9.6f by |phi''|_1 = %.6e; %4d casts freshly wet; |Delta_p|_2 = %.6e\n';

p_change_L2 = 0; % ensure this is defined; needed if OPTS.TOL_P_CHANGE_L2 == 0
%% Process OPTS

% Load default options, then override with any user-specified OPTS.
OPTS = catstruct(omega_defaults(), OPTS);

% Dereference for speed
FIGS_SHOW = OPTS.FIGS_SHOW;
FILE_ID = OPTS.FILE_ID;
ITER_MIN = OPTS.ITER_MIN;
ITER_MAX = OPTS.ITER_MAX;
ITER_START_WETTING = OPTS.ITER_START_WETTING;
ITER_STOP_WETTING  = OPTS.ITER_STOP_WETTING ;
VERBOSE = OPTS.VERBOSE;
POISSON = OPTS.POISSON; 
TOL_P_CHANGE_L2 = OPTS.TOL_P_CHANGE_L2;
TOL_P_UPDATE = OPTS.TOL_P_UPDATE; % tolerance in pressure [dbar] during vertical solve.
WRAP = OPTS.WRAP(:);

DIAGS = (VERBOSE > 0) || (nargout > 3);

if isscalar(ref_cast)
  assert(ref_cast >= 1 && ref_cast <= nij, 'Out of bounds Linear index for ref_cast.');
elseif numel(ref_cast) == 2
  assert(all(ref_cast >= 1) && all(ref_cast(:) <= [ni; nj]), 'ref_cast must index a cast within the domain.')
  ref_cast = sub2ind([ni, nj], ref_cast(1), ref_cast(2)); % Convert into linear index to the reference cast
else
  assert(false, 'ref_cast must be a 1 or 2 element vector');
end

% Pre-calculate things for Breadth First Search
qu = zeros(nij, 1); % queue storing linear indices to pixels
A4 = grid_adjacency([ni,nj], 4, WRAP); % all grid points that are adjacent to all grid points, using 4-connectivity

% Number of bottles per cast. BotK(n) > 0 if and only if pixel n is ocean.
BotK = squeeze(sum(isfinite(S), 1));

eos0 = eos(34.5, 3, 1000);
if eos0 > 1
  TOL_LRPD_L1 = OPTS.TOL_LRPD_L1;          % Density tolerance [kg m^-3]
else
  TOL_LRPD_L1 = OPTS.TOL_LRPD_L1 * 1000^2; % Specific volume tolerance [m^3 kg^-1]
end

% Soft notation, similar to that in MOM6: i = I - 1/2, j = J - 1/2
DIST1_iJ = OPTS.DIST1_iJ; % Distance [m] in 1st dimension centred at (I-1/2, J)
DIST2_Ij = OPTS.DIST2_Ij; % Distance [m] in 2nd dimension centred at (I, J-1/2)
DIST2_iJ = OPTS.DIST2_iJ; % Distance [m] in 2nd dimension centred at (I-1/2, J)
DIST1_Ij = OPTS.DIST1_Ij; % Distance [m] in 1st dimension centred at (I, J-1/2)

% Calculate the ratios of distances, and autoexpand
if POISSON
  DIST2on1_iJ = autoexp(DIST2_iJ ./ DIST1_iJ);
  DIST1on2_Ij = autoexp(DIST1_Ij ./ DIST2_Ij);
else
  sqrtDIST2on1_iJ = autoexp(sqrt(DIST2_iJ ./ DIST1_iJ));
  sqrtDIST1on2_Ij = autoexp(sqrt(DIST1_Ij ./ DIST2_Ij));
end

% auto expand to [ni,nj] sizes, for eps_norms()
DIST1_iJ = autoexp(DIST1_iJ); % Distance [m] in 1st dimension centred at (I-1/2, J)
DIST2_Ij = autoexp(DIST2_Ij); % Distance [m] in 2nd dimension centred at (I, J-1/2)
DIST2_iJ = autoexp(DIST2_iJ); % Distance [m] in 2nd dimension centred at (I-1/2, J)
DIST1_Ij = autoexp(DIST1_Ij); % Distance [m] in 1st dimension centred at (I, J-1/2)
AREA_iJ = autoexp(DIST1_iJ .* DIST2_iJ);   % Area [m^2] centred at (I-1/2, J)
AREA_Ij = autoexp(DIST1_Ij .* DIST2_Ij);   % Area [m^2] centred at (I, J-1/2)


%% Just in time code generation
Pvec = isvector(P);
ni_ = max(ni, 4096); % using variable size code generation and avoiding
nj_ = max(nj, 4096); % recompiling all the time
omega_vertsolve_codegen(nk, ni_, nj_, Pvec, OPTS);
% bfs_conncomp1_codegen(nk, ni_, nj_, Pvec, OPTS);
if ITER_START_WETTING <= ITER_MAX && ITER_STOP_WETTING > 0
  bfs_conncomp1_wet_codegen(nk, ni_, nj_, Pvec, OPTS)
end

%% Get ML: the pressure of the mixed layer
if ITER_MAX > 1
  if isempty(OPTS.ML)
    % Do not remove the mixed layer
    ML = [];
  elseif isstruct(OPTS.ML)
    ML = mixed_layer(S, T, P, OPTS.ML);
  else
    % Use a pre-computed mixed layer
    ML = OPTS.ML;
  end
end

%% Interpolate S and T casts onto surface
[~, K, N] = size(OPTS.Sppc);
if K > 0
  assert(K == nk-1 && N == nij, 'size(Sppc) should be [O, nk-1, ni, nj] == [?, %d, %d, %d]', nk-1, ni, nj);
  Sppc = OPTS.Sppc;
else
  Sppc = OPTS.INTERPFN(P, S);
end

[~, K, N] = size(OPTS.Tppc);
if K > 0
  assert(K == nk-1 && N == nij, 'size(Tppc) should be [O, nk-1, ni, nj] == [?, %d, %d, %d]', nk-1, ni, nj);
  Tppc = OPTS.Tppc;
else
  Tppc = OPTS.INTERPFN(P, T);
end

[s, t] = ppc_val2(P, Sppc, Tppc, lead1(p));
p(isnan(s)) = nan; % ensure same nan structure between s, t, and p. Just in case user gives, e.g., repmat(1000,ni,nj) for a 1000dbar isobaric surface

%% Prepare diagnostics
if DIAGS
  
  diags = struct();
  %diags.mean_p       = nan(ITER_MAX + 1, 1);
  %diags.mean_eos     = nan(ITER_MAX + 1, 1);
  diags.epsL1         = nan(ITER_MAX + 1, 1);
  diags.epsL2         = nan(ITER_MAX + 1, 1);
  
  diags.phi_L1        = nan(ITER_MAX, 1);
  diags.p_change_L1   = nan(ITER_MAX, 1);
  diags.p_change_L2   = nan(ITER_MAX, 1);
  diags.p_change_Linf = nan(ITER_MAX, 1);
  diags.freshly_wet   = nan(ITER_MAX, 1);
  diags.clocktime     = nan(ITER_MAX, 1);
  
  diags.timer_bfs     = nan(ITER_MAX, 1);
  diags.timer_solver  = nan(ITER_MAX, 1);
  diags.timer_update  = nan(ITER_MAX, 1);
  
  % Diagnostics about state BEFORE this (first) iteration
  [epsL2, epsL1] = eps_norms(s, t, p, true, WRAP, {}, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij);
  %mean_p = nanmean(p(:));
  %mean_eos = nanmean(eos(s(:), t(:), p(:)));
  diags.epsL1(1) = epsL1;
  diags.epsL2(1) = epsL2;
  %diags.mean_p(1) = mean_p;
  %diags.mean_eos(1) = mean_eos;
  
  if VERBOSE > 0
    %fprintf(FILE_ID, msg1, log10(epsL2), mean_p, mean_eos);
    fprintf(FILE_ID, msg1, log10(epsL2));
  end
  
end

%% Begin iterations
% Note: the surface exists wherever p is non-nan.  The nan structure of s
% and t is made to match that of p when the vertical solve step is done. 
for iter = 1 : ITER_MAX
  iter_tic = tic;
  
  % --- Remove the Mixed Layer
  % But keep it for the first iteration, which may be initialized from a
  % not very neutral surface
  if iter > 1 && ~isempty(ML)
    p(p < ML) = nan;
  end
  
  
  % --- Determine the connected component containing the reference cast, via Breadth First Search
  mytic = tic;
  if iter >= ITER_START_WETTING && iter <= ITER_STOP_WETTING
    [s, t, p, freshly_wet, qu, qt] = bfs_conncomp1_wet_mex(Sppc, Tppc, P, s, t, p, TOL_P_UPDATE, A4, BotK, ref_cast, qu);
  else
    [qu, qt] = bfs_conncomp1(isfinite(p), A4, ref_cast, qu);
    freshly_wet = 0;
  end
  timer_bfs = toc(mytic);
  assert(qt > 0, "Error: surface is NaN at the reference cast")
  
  
  % --- Solve global matrix problem for ...
  mytic = tic;
  if POISSON
    % ... the exactly determined Poisson equation
    phi = omega_matsolve_poisson(s, t, p, DIST2on1_iJ, DIST1on2_Ij, WRAP, A4, qu, qt, ref_cast);
  else
    % ... the overdetermined gradient equations
    phi = omega_matsolve_grad(s, t, p, sqrtDIST2on1_iJ, sqrtDIST1on2_Ij, WRAP, A4, qu, qt, ref_cast);
  end
  if DIAGS
    timer_solver = toc(mytic);
  end
  

  % --- Update the surface
  mytic = tic();
  p_old = p;      % Record old surface for pinning or diagnostic purposes.
  [p, s, t] = omega_vertsolve_mex(Sppc, Tppc, P, BotK, s, t, p, TOL_P_UPDATE, phi);
  
  % Force p to stay constant at the reference column, identically. This
  % avoids any intolerance from the vertical solver.
  p(ref_cast) = p_old(ref_cast);

  if DIAGS
    timer_update = toc(mytic);
  end

  
  % --- Closing Remarks
  phi_L1 = nanmean(abs(phi(:)));
  if DIAGS || TOL_P_CHANGE_L2 > 0
    p_change = p - p_old;
    p_change_L2 = nanrms(p_change(:));
  end
  
  if DIAGS
    
    diags.clocktime(iter)   = toc(iter_tic);
    
    p_change_L1   = nanmean(abs(p_change(:)));
    p_change_Linf = max(abs(p_change(:)));
    
    
    % Diagnostics about what THIS iteration did
    diags.phi_L1(iter)        = phi_L1;
    diags.p_change_L1(iter)   = p_change_L1;
    diags.p_change_L2(iter)   = p_change_L2;
    diags.p_change_Linf(iter) = p_change_Linf;
    diags.freshly_wet(iter)   = freshly_wet;
    
    diags.timer_bfs(iter)     = timer_bfs;
    diags.timer_solver(iter)  = timer_solver;
    diags.timer_update(iter)  = timer_update;
    
    % Diagnostics about the state AFTER this iteration
    [epsL2, epsL1] = eps_norms(s, t, p, true, WRAP, {}, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij);
    %mean_p = nanmean(p(:));
    %mean_eos = nanmean(eos(s(:),t(:),p(:)));
    diags.epsL1(iter+1) = epsL1;
    diags.epsL2(iter+1) = epsL2;
    %diags.mean_p(iter+1)    = mean_p;
    %diags.mean_eos(iter+1)  = mean_eos;
    
    
    if VERBOSE > 0
      %fprintf(FILE_ID, msg2, iter, diags.clocktime(iter), log10(epsL2), phi_L1, freshly_wet, p_change_L2, mean_p, mean_eos);
      fprintf(FILE_ID, msg2, iter, diags.clocktime(iter), log10(epsL2), phi_L1, freshly_wet, p_change_L2);
    end
  end
  
  % --- Show Figures
  if FIGS_SHOW
    if mod(iter,12) == 1
      hf = figure('Position',[20, 20, 1000, 800], 'Name', 'phi''');
    end
    ax = subplot(4, 3, mod(iter-1,12)+1, 'Parent', hf );
    if nj > ni % p, and phi, are probably lon p lat
      %imagesc(ax, phi)
      imagesc(ax, p_change)
    else       % p, and phi, are probably lat p lon
      %imagesc(ax, phi.')
      imagesc(ax, p_change')
    end
    ax.YDir = 'normal';
    shading('flat')
    %title(sprintf('%d: |\\phi''|_1 = %.2e', iter, phi_L1), 'fontsize',10);
    %caxis(prctile(phi(:), [1 99]));
    title(sprintf('%d: |\\Delta p''|_1 = %.2e', iter, p_change_L1), 'fontsize',10);
    caxis(prctile(p_change(:), [1 99]));
    colorbar();
    drawnow()
  end
  
  % --- Check for convergence
  if (phi_L1 < TOL_LRPD_L1 || p_change_L2 < TOL_P_CHANGE_L2) && iter >= ITER_MIN
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