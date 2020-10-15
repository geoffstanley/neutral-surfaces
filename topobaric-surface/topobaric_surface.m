function [x, s, t, RG, s0, t0, d_fn, diags] = topobaric_surface(S, T, X, x, OPTS)
%TOPOBARIC_SURFACE  Create a topobaric surface.
%
%
% x = topobaric_surface(S, T, X, x, OPTS)
% returns the pressure x (output) of a topobaric surface formed by an
% iterative procedure, initialized from an approximately neutral surface on
% which the pressure is x (input), in an ocean with practical / Absolute
% salinity is S and potential / Conservative temperature is T at datasites
% where the pressure is X. The equation of state for the specific volume
% and its partial derivative with respect to pressure are given by eos.m
% and eos_x.m in the path.  If topobaric_geostrf will not be called, eos.m
% and eos_x.m can instead determine the in-situ density and its derivative
% with respect to pressure. Only one connected component of the surface is
% processed.  Algorithmic options are given by OPTS (see below).  For
% physical units, see "Equation of State" below.
%
% [x, s, t, RG, s0, t0, d_fn] = topobaric_surface(...)
% also returns the Reeb Graph RG, the reference practical / Absolute
% salinity s0, the reference potential / Conservative temperature t0, and
% the multivalued function for delta in terms of pressure, d_fn. delta is
% the specific volume anomaly or the in-situ density anomaly, depending on
% which is given by eos.m.  The actual delta on the surface will closely
% match the delta values evaluated by d_fn on the surface, as follows.
% >> lead1 = @(x) reshape(x, [1, size(x)]);              % add leading singleton dimension
% >> SppX = ppc_linterp(X, S);                           % Interpolant for S in terms of X
% >> TppX = ppc_linterp(X, T);                           % Interpolant for T in terms of X
% >> [s,t] = ppc_val2(X, SppX, TppX, lead1(x));          % get S and T on the surface
% >> d = eos(s, t, x) - eos(s0, t0, x);                  % get delta on the surface
% >> x_ = x(RG.wet);                                     % get x just on the valid surface
% >> d_fn_at_x = nan(size(x));                           % prepare to build a 2D map of evaluations of the delta function
% >> for e = 1 : RG.nArcs                                % loop over all arcs of the Reeb graph
% >>     I = RG.arc_segment{e};                          % get linear indices to all pixels in the region associated with arc e
% >>     d_fn_at_x(I) = pvallin(d_fn(:,e), x_(I);        % evaluate local branch of the multivalued function for delta in terms of x
% >> end
% The smaller OPTS.TOL_X_UPDATE, the smaller OPTS.ITER_L2_CHANGE, and the larger
% OPTS.ITER_MAX, the closer (d_fn_at_x - d) will be to zero.
%
% [x, s, t, RG, s0, t0, d_fn1] = topobaric_surface(S, T, X, x, OPTS)
% with OPTS.REEB == false instead calculates an "orthobaric" surface, in
% which d is a single-valued function of the pressure or depth, d_fn1. This
% is done in just one connected ReGion of the surface, namely the true
% pixels in RG.wet. delta is the specific volume anomaly or in-situ density anomaly,
% and delta on the surface will (less-so, relative to with OPTS.REEB == true) closely
% match the delta values evaluated by d_fn1 on the surface, as follows:
% >> lead1 = @(x) reshape(x, [1, size(x)]);              % add leading singleton dimension
% >> SppX = ppc_linterp(X, S);                           % Interpolant for S in terms of X
% >> TppX = ppc_linterp(X, T);                           % Interpolant for T in terms of X
% >> [s,t] = ppc_val2(X, SppX, TppX, lead1(x));          % get S and T on the surface
% >> d = eos(s, t, x) - eos(s0, t0, x);                  % get delta on the surface
% >> x_ = x(RG.wet);                                     % get x just on the valid surface
% >> d_fn_at_x = ppval(d_fn1, x);                        % evaluate the function for delta in terms of x
%
%
% --- Sizes:  below,
% nk is the maximum number of data points per water column,
% ni is the number of data points in longitude,
% nj is the number of data points in latitude,
% nc is the number of water columns intersecting the surface,
% A is the number of arcs in the Reeb graph,
% N is the number of nodes in the Reeb graph,
% C is the number of cycles in the cycle basis of the Reeb graph.
% O is the order of piecewise polynomials for S and T in terms of X.
%
%
% --- Input:
% S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% X [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
% x [ni, nj]: pressure [dbar] or depth [m, positive] on approximately
%             neutral surface
% OPTS [struct]: options (see "Options" below)
%
% Note: X must increase monotonically along the first dimension.
%
%
% --- Output:
% x [ni, nj]: pressure [dbar] or depth [m, positive] on topobaric surface
% s [ni, nj]: practical / Absolute salinity on topobaric surface
% t [ni, nj]: potential / Conservative temperature on topobaric surface
% RG [struct]: the Reeb Graph and the single connected ReGion (see below)
% s0 [1, 1]: reference practical / Absolute salinity
% t0 [1, 1]: reference potential / Conservative temperature
% d_fn [5, A]: the multivalued function for delta in terms of x. Each column
%             is a single-valued branch, to be evaluated by pvaln or pvallin.
% d_fn1 [struct]: the single-valued function for delta in terms of x, to be
%                evaluated by ppval.
%
% Note: physical units of S, T, X, x, s, t, s0, t0 are determined by eos.m.
% If topobaric_geostrf is to be used, the units of x and X must be be
% [dbar] or [m, positive].
%
%
% --- Options:
% OPTS is a struct containing the following fields. Those marked * are
% required. See ./private/tbs_defaults.m for default values.
%   INTERPFN [function handle]: vertical interpolation function, used to
%       evaluate SppX and TppX if those are not provided.  Default:
%       INTERPFN = @ppc_linterp.
%   SppX [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose
%       knots are X that interpolate S as a function of X in each water
%       column.  E.g. SppX = ppc_linterp(X, S);
%   TppX [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose
%       knots are X that interpolate T as a function of X in each water
%       column.  E.g. TppX = ppc_linterp(X, T);
%   REEB [1, 1]: true to compute topobaric surfaces, false to compute
%     "orthobaric" surfaces.
%   GEOSTRF [1, 1]: true to create a modified topobaric surface (when
%       OPTS.REEB is true), which applies extra constraints to the
%       empirical delta function to ensure the exact geostrophic
%       streamfunction on the surface is well-defined.
%   SPLINE_BREAKS [vector]: knots of the single-valued function for delta
%     in terms of x [units the same as x], when OPTS.REEB is false.
%   SPLINE_ORDER [1, 1]: order of the spline of the single-valued function
%     for delta in terms of x when OPTS.REEB is false. [E: 4 for cubic splines.]
%   MLX []: do not remove the mixed layer
%   MLX [struct]: calculate the mixed layer using these parameters in mixed_layer().
%   MLX [ni, nj]: use a pre-computed mixed layer pressure [dbar] or depth [m]
%   REF_IJ [1, 2]: pixel indices for reference water column.
%   REF_X [1, 1]: override the pressure [dbar] or depth [m] at which the
%       topobaric surface intersects the reference water column
%   REF_S [1, 1]: override the reference practical / Absolute salinity.
%   REF_T [1, 1]: override the reference potential /Conservative temperature.
% * WRAP [1, 2]: determines which dimensions are treated periodic [logical].
%   Set WRAP(i) true when periodic in the i'th dimension of x (i=1,2).
%   FILL_PIX [1, 1]: fill all holes in the surface containing with fewer
%       pixels than this, before computing the Reeb graph.
%   FILL_IJ = [*, 2]: fill all holes in the surface except those containing
%       a pixel index in these rows, before computing the Reeb graph.
%   SIMPLIFY_ARC_REMAIN [1, 1]: number of arcs to remain after graph
%       simplification.
%   SIMPLIFY_WEIGHT_PERSIST [1, 1]: weighting (between 0 and 1) of
%       persistence (the difference between pressure values associated with
%       nodes) as opposed to area (number of pixels in each assocaited region)
%       during graph simplification.
%   TOL_X_UPDATE [1, 1]: error tolerance, in the same units as X [dbar] or
%      [m], when root-finding to update the surface.
%   X_EXPN [1, 1]: solutions to the root-finding problem in each water
%       column are sought in the domain of the local branch of the
%       multivalued function expanded outwards by this amount, in the same
%       units as X [dbar] or [m].
%   ITER_MAX [1, 1]: maximum number of iterations
%   TOL_X_CHANGE_L2 [1, 1]: quit when the change in the L2 norm of
%       the change in surface pressure [dbar] or depth [m] exceeds this
%       value. Set to 0 to deactivate.
%   ITER_START_WETTING [1, 1]: Iteration number at which to start wetting
%                              (use inf to disable wetting entirely).
%   VERBOSE [1, 1]: 0 for no output, 1 for basic information
%   FILE_ID [1, 1]: 1 to write any output to MATLAB terminal, or a file
%       identifier as returned by fopen() to write to that file.
%
%
% --- Reeb Graph:
% The Reeb graph contains N nodes, A arcs, and C cycles in the cycle basis.
% It is encapsulated in the struct RG, containing the following fields:
% wet [ni,nj]: logical map, true where the surface is in the ocean, false
%   where it has grounded or outcropped or has disconnected from the main
%   connected region of the surface. Where wet(i,j) is true, pixel (i,j)
%   has contributed to the emprifical fit for d_fn, and is a "vertex" in a
%   triangular mesh that underlies the Reeb graph. The linear index of the
%   vertex at pixel (i,j) is
%     ij2v(i,j);
%   where
%     ij2v = reshape(cumsum(RG.wet(:)), ni, nj);
%     ij2v(~RG.wet) = 0;
%   maps from pixel space to vertex space.
% n_casts [1, 1]: the number of pixels in the connected region (the number
%   of true values in wet).
% nNodes [1,1]: the number of nodes in the Reeb graph == N.
% nArcs  [1,1]: the number of arcs in the Reeb graph == A.
% node_v  [N,1]: node_v(n) is the index of the critical vertex associated
%   with node n. The linear index of the pixel associated with this vertex
%   is
%     v2I(RG.node_v(n));
%   where
%     v2I = find(ij2v);
%   maps from vertex space to pixel space with linear indices.
%   The 2D index of the pixel associated with this vertex is
%     [i,j] = v2ij(RG.node_v(n));
%   where
%     v2ij = @(v) sub2ind([ni, nj], v2I(v));
%   maps from vertex space to pixel space with 2D coordinates.
% node_fn [N,1]: node_fn(n) is the critical value (pressure or depth)
%   associated with node n. For any node n, 1 <= n <= N,
%     RG.node_fn(n) == p(v2I(RG.node_v(n)));
%   in the non-Boussinesq case, or
%     RG.node_fn(n) == z(v2I(RG.node_v(n)));
%   in the Boussinesq case.
% arc_from [A,1]: arc_from(a) is the lower  node incident to arc a.
% arc_to   [A,1]: arc_to(a)   is the higher node incident to arc a.
% 	Note, arc_from(a) == n_1 and arc_to(a) == n_2 means arc a is incident
%   to nodes n_1 and n_2, and that node_fn(n_1) < node_fn(n_2).
% node_next {N,1}: node_next{n} are the arcs a for which arc_from(a) == n.
% node_prev {N,1}: node_prev{n} are the arcs a for which arc_to(a) == n.
% node_type [N,1]: node_type(n) is the type of node n, being
%   1 for a minima, for which node_prev{n} == [] and node_next{n} == a with
%     arc_from(a) == n, or
%   3 for a maxima, for which node_next{n} == [] and node_prev{n} == a with
%     arc_to(a) == n, or
%   2 for a saddle, for which length(node_prev{n}) >= 1 and
%     length(node_next{n}) >= 1.
% arc_segment {A,1}: arc_segment{a} are the indices, among only those
%   pixels in the single connected region, to the vertices in the
%   associated region for arc a. If the Reeb graph was not simplified
%   (OPTS.SIMPLIFY_ARCS >= A), then the function value (pressure or depth)
%   of all vertices in arc_segment{a} lies between the function value of
%   the "lower" node incident to arc a and the function value of the
%   "upper" node incident to arc a. That is,
%     RG.node_fn(RG.arc_segment{a}) >= RG.node_fn(RG.arc_from(a));
%   and
%     RG.node_fn(RG.arc_segment{a}) <= RG.node_fn(RG.arc_to(a));
% cb_nodes {C,1}: cb_nodes{c} are the nodes in cycle c of the cycle basis.
% cb_arcs {C,1}: cb_arcs{c} are arcs in cycle c of the cycle basis.
%   Note, [cb_nodes{c}(1), cb_nodes{c}(1), ..., cb_nodes{c}(end),
%   cb_arcs{c}(end), cb_nodes{c}(1)] is the walk for cycle c in the cycle
%   basis.
% graph [N,N]: (sparse), graph(m,n) == a means arc a is incident to nodes m
%   and n.
% bfs_parent_node [N,1]: bfs_parent_node(n) == m means node n was
%   discovered from node m in the breadth-first search (bfs).
% bfs_topo_order [N,1]: is a "topological ordering" of all nodes, such that
%   bfs_topo_order(m) is discovered in the bfs before bfs_topo_order(n)
%   whenever m < n.
% bfs_missing_arc [C,1]: bfs_missing_arc(c) == a means arc a was not
%   traversed in the bfs, and thus arc a defines cycle c in the cycle
%   basis, namely the only cycle in the minimum spanning tree (the graph
%   containing only arcs traversed in the bfs) augmented by arc a. Note,
%   bfs_missing_arc(c) == cb_arcs{c}(1).
%
%
% --- Equation of State:
% The MATLAB path* must contain two functions, eos.m and eos_x.m. Both
% accept 3 inputs: S, T, and X. eos(S, T, X) is the equation of state and
% eos_x(S, T, X) is the partial derivative of the equation of state with
% respect to X (holding S and T constant).
% *Note: It is not sufficient to simply have these eos functions in the
% current working directory, because the compiled MEX functions will not be
% able to find them there.  They must be in the MATLAB path.  If they are
% in the current working directory, use `addpath(pwd)` to add the current
% working directory to the top of MATLAB's path.
%
% For a non-Boussinesq ocean, x and X are pressure [dbar], and if
% topobaric_geostrf will be called then the equation of state must return
% the specific volume [m^3 kg^-1].
%
% For a Boussinesq ocean, x and X are depth [m, positive], and if
% topobaric_geostrf will be called then the equation of state must return
% the in-situ density [kg m^-3].
%
% Various equation of state functions are found in ../lib/eos/.  Simply
% copy the desired functions to another location in the MATLAB path (such
% as this directory) and rename them eos.m and eos_x.m.  Note, the
% Boussinesq equation of state is often (but not always) just the regular
% equation of state but using a hydrostatic pressure (10^-4 * grav * rho_c
% * z) where grav [m s^-2] is the gravitational acceleration, rho_c [kg
% m^-3] is the Boussinesq reference density, and z [m, positive] is the
% depth. In such a case, simply make new eos.m and eos_x.m functions that
% accept depth as the third input by modifying the original functions that
% take pressure; this involves hard-coding the gravitational acceleration
% and Boussinesq reference density into the function.  An example of a
% Boussinesq eos.m and eos_x.m are given for the densjmd95 equation of
% state, in ../lib/eos/eoscg_densjmd95_bsq.m and
% ../lib/eos/eoscg_densjmd95_bsq_dz.m.  As the latter must return [kg m^-4]
% not [kg m^-3 dbar^-1], note the final multiplication by a conversion
% factor with units [dbar m^-1].  Finally, note that eos.m and eos_x.m must
% be compatible with MATLAB's code generation (codegen), which may entail
% eliminating input checks and/or expansion of input variables (MATLAB's
% automatic expansion now handles this).
%
%
% --- References:
% Stanley, G.J., 2019a. Neutral surface topology. Ocean Modelling 138,
% 88–106. https://doi.org/10.1016/j.ocemod.2019.01.008

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2020  Geoff Stanley
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
%
% Modified by : --
% Date        : --
% Changes     : --

% --- Notes on the code:
% 3D scalar fields [nk,ni,nj] are upper case letters, e.g. S
% 2D scalar fields [ni,nj] are lower case letters, e.g. x
% 1D scalar fields [1,nc] are lower case letters followed by an underscore,
% e.g. x_.
% The data-heavy 3D variables S and T (and possibly X) are not modified
% internally (except to ensure double precision).
% The specific volume anomaly or the in-situ density anomaly is denoted d
% (for delta).
% The partial derivative of d with respect to x is denoted dx.

%% Simple checks and preparations
warning('off', 'MATLAB:nargchk:deprecated') % splinefit uses nargchk
S = double(S);
T = double(T);
X = double(X);
x = double(x);

% Process mandatory options
assert(isstruct(OPTS) && isfield(OPTS, 'WRAP') && isvector(OPTS.WRAP) && length(OPTS.WRAP) == 2, 'OPTS.WRAP must be provided as a vector of length 2');

[ni,nj] = size(x);
nij = ni * nj;
nk = size(X,1);
Xvec = isvector(X);
is1D = @(F) ismatrix(F) && all(size(F) == [nk,1]);
is2D = @(F) ismatrix(F) && all(size(F) == [ni,nj]);
is3D = @(F) ndims(F) == 3 && all(size(F) == [nk,ni,nj]);

assert(all(size(S) == size(T)), 'S and T must be the same size.');
assert(is3D(S), 'S and T must have 3 dimensions [nk, ni, nj]');
assert(is2D(x), 'x must have 2 dimensions [ni, nj]');
assert(is1D(X) || is3D(X), 'X must be a vector of length nk, or a matrix of size [nk, ni, nj]');

% Simple anonymous functions
lead1 = @(x) reshape(x, [1 size(x)]); % augment with leading singleton dimension
nanrms = @(x) sqrt(nanmean(x(:) .* conj(x(:)))); % root mean square, ignoring nans
autoexp = @(x) repmat(x, ni / size(x,1), nj / size(x,2)); % automatic expansion to [ni,nj]

% Set up empty outputs to avoid error if OPTS.ITER_MAX < 1
RG = struct();
d_fn = [];

% Pre-calculate things for Breadth First Search
qu = zeros(nij, 1); % queue storing linear indices to pixels
A = grid_adjacency([ni,nj], 4, OPTS.WRAP); % all grid points that are adjacent to all grid points, using 4-connectivity

% Number of bottles per cast. BotK(n) > 0 if and only if pixel n is ocean.
BotK = squeeze(sum(isfinite(S), 1));

%% Process OPTS

% Load default options, then override with any user-specified OPTS.
OPTS = catstruct(tbs_defaults(ni,nj), OPTS);

% Whether the rectangular data contains all vertices in the simplical decomposition
ALL_VERTS = OPTS.DECOMP(1) == 'd' && OPTS.FILL_PIX == 0 && isempty(OPTS.FILL_IJ);
assert(ALL_VERTS || OPTS.ITER_MAX < 1, ...
  'To empirically fit the multivalued function, must use ''diagonal'' simplical decomposition and fill no holes during pre-processing');

% Soft notation, similar to that in MOM6: i = I - 1/2, j = J - 1/2
DIST1_iJ = autoexp(OPTS.DIST1_iJ); % Distance [m] in 1st dimension centred at (I-1/2, J)
DIST2_Ij = autoexp(OPTS.DIST2_Ij); % Distance [m] in 2nd dimension centred at (I, J-1/2)
DIST2_iJ = autoexp(OPTS.DIST2_iJ); % Distance [m] in 2nd dimension centred at (I-1/2, J)
DIST1_Ij = autoexp(OPTS.DIST1_Ij); % Distance [m] in 1st dimension centred at (I, J-1/2)
AREA_iJ = DIST1_iJ .* DIST2_iJ;   % Area [m^2] centred at (I-1/2, J)
AREA_Ij = DIST1_Ij .* DIST2_Ij;   % Area [m^2] centred at (I, J-1/2)



%% Just In Time code generation
ni_ = max(ni, 4096); % using variable size code generation and avoiding
nj_ = max(nj, 4096); % recompiling all the time
if OPTS.REEB
  tbs_vertsolve_codegen(nk, ni_, nj_, Xvec, OPTS);
else
  obs_vertsolve_codegen(nk, ni_, nj_, Xvec, OPTS);
end
bfs_conncomp_codegen(nk, ni_, nj_, Xvec, true, OPTS);
if OPTS.ITER_START_WETTING <= OPTS.ITER_MAX
  bfs_wet_codegen(nk, ni_, nj_, Xvec, OPTS);
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
if isfield(OPTS, 'SppX')
  SppX = OPTS.SppX;
else
  SppX = OPTS.INTERPFN(X, S);
end
if isfield(OPTS, 'TppX')
  TppX = OPTS.TppX;
else
  TppX = OPTS.INTERPFN(X, T);
end
[s,t] = ppc_val2(X,SppX,TppX,lead1(x));

%% Prepare Diagnostics
DIAGS = (OPTS.VERBOSE > 0) || (nargout >= 8);


diags = struct();
diags.clocktime     = nan(OPTS.ITER_MAX,1);
diags.x_change_L1   = nan(OPTS.ITER_MAX,1);
diags.x_change_L2   = nan(OPTS.ITER_MAX,1);
diags.x_change_Linf = nan(OPTS.ITER_MAX,1);
diags.freshly_wet   = nan(OPTS.ITER_MAX,1);

diags.epsL1 = nan(OPTS.ITER_MAX + 1,1);
diags.epsL2 = nan(OPTS.ITER_MAX + 1,1);

diags.timer_wetting   = nan(OPTS.ITER_MAX, 1);
diags.timer_recon     = nan(OPTS.ITER_MAX, 1);
diags.timer_reebgraph = nan(OPTS.ITER_MAX, 1);
diags.timer_fitting   = nan(OPTS.ITER_MAX, 1);
diags.timer_update    = nan(OPTS.ITER_MAX, 1);

% Compute the L1 and L2 norms of the neutrality (epsilon) errors
[diags.epsL2(1), diags.epsL1(1)] = eps_norms(s, t, x, false, OPTS.WRAP, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij);

%% Process the reference cast:
% Get linear index to reference cast
assert(numel(OPTS.REF_IJ) == 2, 'OPTS.REF_IJ must be a vector of length 2');
I0 = sub2ind([ni,nj], OPTS.REF_IJ(1), OPTS.REF_IJ(2));

% Set x0: the x value the surface will be pinned to at this ref cast
if isscalar(OPTS.REF_X)
  x0 = OPTS.REF_X;
else
  x0 = x(I0);
end

% Get (s0,t0): (S,T) values at the reference cast at the reference x.
s0 = s(I0);
t0 = t(I0);
d0 = eos(s0, t0, x0); % temporary, will be over-written

% Overwrite with reference values if provided
if isscalar(OPTS.REF_S)
  s0 = OPTS.REF_S;
end
if isscalar(OPTS.REF_T)
  t0 = OPTS.REF_T;
end

% Pre-compute delta (using the reference S and T) at the reference cast
% at the reference x. If reference S and T are not provided, this is 0.
d0 = d0 - eos(s0,t0,x0);


%% Begin iterations
for iter = 1 : OPTS.ITER_MAX
  iter_tic = tic();
  
  % --- Remove the Mixed Layer
  % But keep it for the first iteration, which may be initialized from a
  % not very neutral surface
  if iter > 1 && ~isempty(MLX)
    x(x < MLX) = nan;
  end
  
  
  % --- Wetting via Breadth First Search
  mytic = tic();
  if iter >= OPTS.ITER_START_WETTING
    [s, t, x, freshly_wet, qu] = bfs_wet_mex(SppX, TppX, X, s, t, x, OPTS.TOL_X_UPDATE, A, BotK, qu);
  else
    freshly_wet = 0;
  end
  
  % --- Breadth First Search to find connected region
  [qu, qts] = bfs_conncomp_one_mex(isfinite(x), A, I0, qu);
  qt = qts(2)-1;
  
  % Keep only the component of the surface connected to the reference cast
  wet = false(ni,nj);
  wet(qu(1:qt)) = true;
  x(~wet) = nan;
  
  if DIAGS
    timer_wetting = toc(mytic);
  end
  
  
  if OPTS.REEB
    mytic = tic();
    % --- Calculate the Reeb graph:
    % 1. Pre-process to select one region, possibly filling in certain holes;
    % 2. Calculate the Reeb graph;
    % 3. Post-process to undo any hole-filling and possibly simplify the graph.
    [x, arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, nNodes, ~, ~, timer_recon] = calc_reeb_graph(x, OPTS);
    
    % --- Prepare info about cycles
    [~, graph, ~, bfs_parent_node, bfs_topo_order, bfs_missing_arc, cb_arcs, cb_nodes] = ...
      cycle_analy_bfs(arc_from, arc_to, nNodes);
    
    if DIAGS
      timer_reebgraph = toc(mytic) - timer_recon;
    end
  else
    timer_recon = 0;
    timer_reebgraph = 0;
  end
  
  
  % --- Begin empirical fitting and graph integration
  mytic = tic();
  
  % --- Build x_, the vector containing x at each valid point on the surface
  x_ = x(wet).'; % Row vector
  n_casts = length(x_);
  
  
  if OPTS.REEB
    % --- Build vector associating each water column with an arc
    % Note some water columns (saddles) are associated with multiple arcs!
    % arc and arc_ simply take one of them.
    arc_ = zeros(1, n_casts);
    for e = 1 : nArcs
      arc_(arc_segment{e}) = e;
    end % Now all arc_ > 0
    arc = zeros(ni,nj);
    arc(wet) = arc_;
  end
  
  
  % --- Get derivative of delta on surface
  dx = eos_x(s, t, x) - eos_x(s0, t0, x);
  dx_ = dx(wet).';
  
  
  % --- Fit d(delta) / d(pressure), then integrate it
  if OPTS.REEB
    % --- Fit partial derivative of delta w.r.t. x as a multivalued function of x,
    dx_fn = branches_fit(x_, dx_, arc_from, arc_to, arc_segment, node_fn, cb_arcs, cb_nodes, OPTS.GEOSTRF);
    
    % --- Integrate to get delta as a multivalued function of x
    d_fn = int_graph_fun(dx_fn, arc_from, arc_to, node_fn, graph, bfs_parent_node, bfs_topo_order, bfs_missing_arc);
    
    % Adjust d_fn so that the new surface coincides with the old one on the starting water column
    % i.e. force d_fn(:,arc(I0)) evaluated at x0 to equal d0.
    d_fn(end,:) = d_fn(end,:) + (d0 - pvallin(d_fn(:,arc(I0)), x0));
    
    % Double check matching conditions at nodes of the Reeb Graph.
    %{
        for n = 1:nNodes
            neigharcs = [node_prev{n}; node_next{n}].';
            e = neigharcs(1);
            node_x = node_fn(n);
            d_fn_at_e = pvaln(d_fn(:,e), node_x);
            for ee = neigharcs(2:end)
                diff = d_fn_at_e - pvaln(d_fn(:,ee), node_x);
                if abs(diff) > 1e-9
                    fprintf(OPTS.FILE_ID, 'Bad Matching (diff=%e) at node %d, arcs %d and %d\n', diff, n, e, ee);
                end
            end
        end
        clear neigharcs node_x d_fn_at_e diff ee e
    %}
  else
    % --- Fit partial derivative of delta w.r.t. x as a single-valued function of x
    breaks = OPTS.SPLINE_BREAKS ;
    breaks = min(breaks, max(x_));
    breaks = max(breaks, min(x_));
    breaks = unique(breaks);
    dx_fn = splinefit(x_, dx_, breaks, OPTS.SPLINE_ORDER);
    
    % --- Integrate to get delta as a single-valued function of x
    d_fn = ppint(dx_fn);
    
    % Adjust d_fn so that the new surface coincides with the old one on the starting water column
    % i.e. force d_fn(:,arc(I0)) evaluated at x0 to equal d0.
    d_fn.coefs(:,end) = d_fn.coefs(:,end) + ( d0 - ppval(d_fn, x0) );
  end
  clear dx_fn dx_
  if DIAGS
    timer_fitting = toc(mytic);
  end
  
  % --- Solve for new pressures at which specific volume is that of the multivalued function
  mytic = tic();
  x_old = x; % Record old surface for diagnostic purposes. 
  if OPTS.REEB
    [x, s, t] = tbs_vertsolve_mex(SppX, TppX, X, BotK, s, t, x, arc, d_fn, s0, t0, OPTS.TOL_X_UPDATE, OPTS.X_EXPN);
  else
    [x, s, t] = obs_vertsolve_mex(SppX, TppX, X, BotK, s, t, x, d_fn.breaks, d_fn.coefs, s0, t0, OPTS.TOL_X_UPDATE);
  end
  if DIAGS
    timer_update = toc(mytic);
  end
  
  % --- Get ready for next iteration
  x_change = x - x_old;
  x_change_L2 = nanrms(x_change(:));
  
  % --- Closing remarks
  if DIAGS
    
    clocktime = toc(iter_tic);
    
    x_change_L1 = nanmean(abs(x_change(:)));
    x_change_Linf = max(abs(x_change(:)));
    if OPTS.VERBOSE >= 1
      fprintf(OPTS.FILE_ID, 'Iter %2d (%.2fsec):  %4d casts wet; %4d casts in/outcropped; x change has: L_1 %.6e, L_2 %.6e, L_inf %.6e\n',...
        iter, toc(iter_tic), freshly_wet, sum(isnan(x(wet))), x_change_L1, x_change_L2, x_change_Linf);
    end
    diags.clocktime(iter) = clocktime;
    
    
    % Diagnostics about what THIS iteration did
    diags.x_change_L1(iter) = x_change_L1;
    diags.x_change_L2(iter) = x_change_L2;
    diags.x_change_Linf(iter) = x_change_Linf;
    diags.freshly_wet(iter) = freshly_wet;
    
    diags.timer_wetting(iter)   = timer_wetting;
    diags.timer_recon(iter)     = timer_recon;
    diags.timer_reebgraph(iter) = timer_reebgraph;
    diags.timer_fitting(iter)   = timer_fitting;
    diags.timer_update(iter)    = timer_update;
    
    % Diagnostics about the state AFTER this iteration
    
    % Compute the L1 and L2 norms of the neutrality (epsilon) errors
    [diags.epsL2(iter+1), diags.epsL1(iter+1)] = eps_norms(s, t, x, false, OPTS.WRAP, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij);
    
  end
  
  % --- Check for convergence
  if (x_change_L2 < OPTS.TOL_X_CHANGE_L2) && iter >= OPTS.ITER_MIN
    break
  end
  
end

if nargout >= 4
  % Return Reeb Graph (if it was computed) and ReGion
  RG.wet = wet;
  RG.n_casts = n_casts;
  if OPTS.REEB
    RG.arc_from        = arc_from;
    RG.arc_to          = arc_to;
    RG.arc_segment     = arc_segment;
    RG.node_next       = node_next;
    RG.node_prev       = node_prev;
    RG.node_type       = node_type;
    RG.node_v          = node_v;
    RG.nArcs           = nArcs;
    RG.nNodes          = nNodes;
    RG.node_fn         = node_fn;
    
    % Add extra info to RG
    RG.cb_nodes        = cb_nodes;
    RG.cb_arcs         = cb_arcs;
    RG.graph           = graph;
    RG.bfs_parent_node = bfs_parent_node;
    RG.bfs_topo_order  = bfs_topo_order;
    RG.bfs_missing_arc = bfs_missing_arc;
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