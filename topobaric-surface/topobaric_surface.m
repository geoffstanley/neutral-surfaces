function [p, s, t, RG, s0, t0, d_fn, diags] = topobaric_surface(S, T, P, p, ref_cast, WRAP, OPTS)
%TOPOBARIC_SURFACE  Create a topobaric surface.
%
%
% p = topobaric_surface(S, T, P, p, ref_cast, WRAP, OPTS)
% returns the pressure p (output) of a topobaric surface formed by an
% iterative procedure, initialized from an approximately neutral surface on
% which the pressure is p (input), in an ocean with practical / Absolute
% salinity is S and potential / Conservative temperature is T at datasites
% where the pressure is P.  The depth or pressure of the topobaric surface
% is pinned, unchanging through the iterations, at the reference cast
% indexed by ref_cast.  The equation of state for the specific volume and
% its partial derivative with respect to pressure are given by eos.m and
% eos_p.m in the path.  If topobaric_geostrf will not be called, eos.m and
% eos_p.m can instead determine the in-situ density and its derivative with
% respect to pressure.  Only one connected component of the surface is
% processed.  The domain is periodic in the i'th horizontal dimension iff
% WRAP(i) is true.  Algorithmic options are given by OPTS (see below).  For
% physical units, see "Equation of State" below.
%
% [p, s, t, RG, s0, t0, d_fn] = topobaric_surface(...)
% also returns the Reeb Graph RG, the reference practical / Absolute
% salinity s0, the reference potential / Conservative temperature t0, and
% the multivalued function for delta in terms of pressure, d_fn. delta is
% the specific volume anomaly or the in-situ density anomaly, depending on
% which is given by eos.m.  The actual delta on the surface will closely
% match the delta values evaluated by d_fn on the surface, as follows.
% >> Sppc = ppc_linterp(P, S);                           % Interpolant for S in terms of P
% >> Tppc = ppc_linterp(P, T);                           % Interpolant for T in terms of P
% >> [s,t] = ppc_val2(P, Sppc, Tppc, p);                 % get S and T on the surface
% >> d = eos(s, t, p) - eos(s0, t0, p);                  % gegeoffstanley@gmail.com"t delta on the surface
% >> p_ = p(RG.wet);                                     % get p just on the valid surface
% >> d_fn_at_p = nan(size(p));                           % prepare to build a 2D map of evaluations of the delta function
% >> for e = 1 : RG.nArcs                                % loop over all arcs of the Reeb graph
% >>     I = RG.arc_segment{e};                          % get linear indices to all pixels in the region associated with arc e
% >>     d_fn_at_p(I) = pvallin(d_fn(:,e), p_(I);        % evaluate local branch of the multivalued function for delta in terms of p
% >> end
% The smaller OPTS.TOL_P_UPDATE, the smaller OPTS.ITER_L2_CHANGE, and the larger
% OPTS.ITER_MAX, the closer (d_fn_at_p - d) will be to zero.
%
% [p, s, t, RG, s0, t0, d_fn1] = topobaric_surface(S, T, P, p, ref_cast, WRAP, OPTS)
% with OPTS.REEB == false instead calculates an "orthobaric" surface, in
% which d is a single-valued function of the pressure or depth, d_fn1. This
% is done in just one connected ReGion of the surface, namely the true
% pixels in RG.wet. delta is the specific volume anomaly or in-situ density
% anomaly, and delta on the surface will (less-so, relative to with
% OPTS.REEB == true) closely match the delta values evaluated by d_fn1 on
% the surface, as follows:
% >> Sppc = ppc_linterp(P, S);                           % Interpolant for S in terms of P
% >> Tppc = ppc_linterp(P, T);                           % Interpolant for T in terms of P
% >> [s,t] = ppc_val2(P, Sppc, Tppc, p);                 % get S and T on the surface
% >> d = eos(s, t, p) - eos(s0, t0, p);                  % get delta on the surface
% >> p_ = p(RG.wet);                                     % get p just on the valid surface
% >> d_fn_at_p = ppval(d_fn1, p);                        % evaluate the function for delta in terms of p
%
% [..., diags] = topobaric_surface(...)
% also returns diagnostic information about the solution at each iteration.
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
% O is the order of piecewise polynomials for S and T in terms of P.
%
%
% --- Input:
% S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% P [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
% p [ni, nj]: pressure [dbar] or depth [m, positive] on approximately
%             neutral surface
% ref_cast [1, 1] or [2, 1] : linear index or 2D index to the reference cast
% WRAP [2 element array]: determines which dimensions are treated periodic
%                         [logical].  Set WRAP(i) to true when periodic in 
%                         the i'th lateral dimension(i=1,2).
% OPTS [struct]: options (see "Options" below)
%
% Note: P must increase monotonically along the first dimension.
%
%
% --- Output:
% p [ni, nj]: pressure [dbar] or depth [m, positive] on topobaric surface
% s [ni, nj]: practical / Absolute salinity on topobaric surface
% t [ni, nj]: potential / Conservative temperature on topobaric surface
% RG [struct]: the Reeb Graph and the single connected ReGion (see below)
% s0 [1, 1]: reference practical / Absolute salinity
% t0 [1, 1]: reference potential / Conservative temperature
% d_fn [5, A]: the multivalued function for delta in terms of p. Each column
%             is a single-valued branch, to be evaluated by pvaln or pvallin.
% d_fn1 [struct]: the single-valued function for delta in terms of p, to be
%                evaluated by ppval.
%
% Note: physical units of S, T, P, p, s, t, s0, t0 are determined by eos.m.
% If topobaric_geostrf is to be used, the units of p and P must be be
% [dbar] or [m, positive].
%
%
% --- Options:
% OPTS is a struct containing the following fields. See ./private/tbs_defaults.m for default values.
%   DIST1_iJ [ni, nj]: 
%     Distance [m] in 1st dimension centred at (I-1/2, J): the distance 
%     between cell centres at (I,J) and (I-1,J)
%   DIST2_iJ [ni, nj]:
%     Distance [m] in 2nd dimension centred at (I-1/2, J): the length of
%     the face between cells at (I,J) and (I-1,J)
%   DIST1_Ij [ni, nj]:
%     Distance [m] in 1st dimension centred at (I, J-1/2): the length of
%     the face between cells at (I,J) and (I,J-1)
%   DIST2_Ij [ni, nj]: 
%     Distance [m] in 2nd dimension centred at (I, J-1/2): the distance 
%     between cell centres at (I,J) and (I,J-1)
%   INTERPFN [function handle]: vertical interpolation function, used to
%       evaluate Sppc and Tppc if those are not provided.  Default:
%       INTERPFN = @ppc_linterp.
%   Sppc [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose
%       knots are P that interpolate S as a function of P in each water
%       column.  E.g. Sppc = ppc_linterp(P, S);
%   Tppc [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose
%       knots are P that interpolate T as a function of P in each water
%       column.  E.g. Tppc = ppc_linterp(P, T);
%   REEB [1, 1]: true to compute topobaric surfaces, false to compute
%     "orthobaric" surfaces.
%   GEOSTRF [1, 1]: true to create a modified topobaric surface (when
%       OPTS.REEB is true), which applies extra constraints to the
%       empirical delta function to ensure the exact geostrophic
%       streamfunction on the surface is well-defined.
%   SPLINE_BREAKS [vector]: knots of the single-valued function for delta
%     in terms of p [units the same as p], when OPTS.REEB is false.
%   SPLINE_ORDER [1, 1]: order of the spline of the single-valued function
%     for delta in terms of p when OPTS.REEB is false. [E: 4 for cubic splines.]
%   ML []: do not remove the mixed layer
%   ML [struct]: calculate the mixed layer using these parameters in mixed_layer().
%   ML [ni, nj]: use a pre-computed mixed layer pressure [dbar] or depth [m]
%   REF_P [1, 1]: override the pressure [dbar] or depth [m] at which the
%       topobaric surface intersects the reference water column
%   REF_S [1, 1]: override the reference practical / Absolute salinity.
%   REF_T [1, 1]: override the reference potential / Conservative temperature.
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
%   TOL_P_UPDATE [1, 1]: error tolerance, in the same units as P [dbar] or
%      [m], when root-finding to update the surface.
%   P_EXPN [1, 1]: solutions to the root-finding problem in each water
%       column are sought in the domain of the local branch of the
%       multivalued function expanded outwards by this amount, in the same
%       units as P [dbar] or [m].
%   TOL_P_CHANGE_L2 [1, 1]: quit when the change in the L2 norm of
%       the change in surface pressure [dbar] or depth [m] exceeds this
%       value. Set to 0 to deactivate.
%   ITER_MAX [1, 1]: maximum number of iterations
%   ITER_START_WETTING [scalar]: Start wetting on iterations that are
%       >= ITER_START_WETTING. To disable wetting, set to +inf. Default: 1.
%   ITER_STOP_WETTING [scalar]: Do wetting for iterations that are
%       <= ITER_STOP_WETTING. To disable wetting, set to 0. Default: 5.
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
% The MATLAB path* must contain two functions, eos.m and eos_p.m. Both
% accept 3 inputs: S, T, and P. eos(S, T, P) is the equation of state and
% eos_p(S, T, P) is the partial derivative of the equation of state with
% respect to P (holding S and T constant).
% *Note: It is not sufficient to simply have these eos functions in the
% current working directory, because the compiled MEX functions will not be
% able to find them there.  They must be in the MATLAB path.  If they are
% in the current working directory, use `addpath(pwd)` to add the current
% working directory to the top of MATLAB's path.
%
% For a non-Boussinesq ocean, p and P are pressure [dbar], and if
% topobaric_geostrf will be called then the equation of state must return
% the specific volume [m^3 kg^-1].
%
% For a Boussinesq ocean, p and P are depth [m, positive], and if
% topobaric_geostrf will be called then the equation of state must return
% the in-situ density [kg m^-3].
%
% Various equation of state functions are found in ../lib/eos/.  Simply
% copy the desired functions to another location in the MATLAB path (such
% as this directory) and rename them eos.m and eos_p.m.  Note, the
% Boussinesq equation of state is often (but not always) just the regular
% equation of state but using a hydrostatic pressure (10^-4 * grav * rho_c
% * z) where grav [m s^-2] is the gravitational acceleration, rho_c [kg
% m^-3] is the Boussinesq reference density, and z [m, positive] is the
% depth. In such a case, simply make new eos.m and eos_p.m functions that
% accept depth as the third input by modifying the original functions that
% take pressure; this involves hard-coding the gravitational acceleration
% and Boussinesq reference density into the function.  An example of a
% Boussinesq eos.m and eos_p.m are given for the densjmd95 equation of
% state, in ../lib/eos/eoscg_densjmd95_bsq.m and
% ../lib/eos/eoscg_densjmd95_bsq_dz.m.  As the latter must return [kg m^-4]
% not [kg m^-3 dbar^-1], note the final multiplication by a conversion
% factor with units [dbar m^-1].  Finally, note that eos.m and eos_p.m must
% be compatible with MATLAB's code generation (codegen), which may entail
% eliminating input checks and/or expansion of input variables (MATLAB's
% automatic expansion now handles this).
%
%
% --- References:
% Stanley, G.J., 2019a. Neutral surface topology. Ocean Modelling 138,
% 88–106. https://doi.org/10.1016/j.ocemod.2019.01.008

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% --- Notes on the code:
% 3D scalar fields [nk,ni,nj] are upper case letters, e.g. S
% 2D scalar fields [ni,nj] are lower case letters, e.g. p
% 1D scalar fields [1,nc] are lower case letters followed by an underscore,
% e.g. p_.
% The data-heavy 3D variables S and T (and possibly P) are not modified
% internally (except to ensure double precision).
% The specific volume anomaly or the in-situ density anomaly is denoted d
% (for delta).
% The partial derivative of d with respect to p is denoted dp.

%% Simple checks and preparations
warning('off', 'MATLAB:nargchk:deprecated') % splinefit uses nargchk
S = double(S);
T = double(T);
P = double(P);
p = double(p);

if nargin < 7 || isempty(OPTS)
  OPTS = struct();
end

[ni,nj] = size(p);
nij = ni * nj;
nk = size(P,1);
Pvec = isvector(P);
is1D = @(F) ismatrix(F) && all(size(F) == [nk,1]);
is2D = @(F) ismatrix(F) && all(size(F) == [ni,nj]);
is3D = @(F) ndims(F) == 3 && all(size(F) == [nk,ni,nj]);

assert(all(size(S) == size(T)), 'S and T must be the same size.');
assert(is3D(S), 'S and T must have 3 dimensions [nk, ni, nj]');
assert(is2D(p), 'p must have 2 dimensions [ni, nj]');
assert(is1D(P) || is3D(P), 'P must be a vector of length nk, or a matrix of size [nk, ni, nj]');

if isscalar(ref_cast)
  assert(ref_cast >= 1 && ref_cast <= nij, 'Out of bounds Linear index for ref_cast.');
elseif numel(ref_cast) == 2
  assert(all(ref_cast >= 1) && all(ref_cast(:) <= [ni; nj]), 'ref_cast must index a cast within the domain.')
  ref_cast = sub2ind([ni, nj], ref_cast(1), ref_cast(2)); % Convert into linear index to the reference cast
else
  assert(false, 'ref_cast must be a 1 or 2 element vector');
end

assert(length(WRAP) == 2, 'WRAP must be provided as a vector of length 2');

% Simple anonymous functions
lead1 = @(x) reshape(x, [1 size(x)]); % augment with leading singleton dimension
nanmean = @(x) mean(x, 'omitnan'); % mean, ignoring nans
nanrms = @(x) sqrt(nanmean(x(:) .* conj(x(:)))); % root mean square, ignoring nans
autoexp = @(x) repmat(x, ni / size(x,1), nj / size(x,2)); % automatic expansion to [ni,nj]

% Set up empty outputs to avoid error if OPTS.ITER_MAX < 1
RG = struct();
d_fn = [];

% Pre-calculate things for Breadth First Search
qu = zeros(nij, 1); % queue storing linear indices to pixels
A4 = grid_adjacency([ni,nj], 4, WRAP); % all grid points that are adjacent to all grid points, using 4-connectivity

% Number of bottles per cast. BotK(n) > 0 if and only if pixel n is ocean.
BotK = squeeze(sum(isfinite(S), 1));

%% Process OPTS

% Load default options, then override with any user-specified OPTS.
OPTS = catstruct(tbs_defaults(), OPTS);

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
  tbs_vertsolve_codegen(nk, ni_, nj_, OPTS);
else
  obs_vertsolve_codegen(nk, ni_, nj_, OPTS);
end
% bfs_conncomp1_codegen(nk, ni_, nj_, OPTS);
if OPTS.ITER_START_WETTING <= OPTS.ITER_MAX && OPTS.ITER_STOP_WETTING > 0
  bfs_conncomp1_wet_codegen(nk, ni_, nj_, OPTS)
end

%% Get ML: the pressure of the mixed layer
if OPTS.ITER_MAX > 1
  if isempty(OPTS.ML)
    % Do not remove the mixed layer
    REMOVE_MIXED_LAYER = false;
    ML = -inf(ni, nj);  % this deactivates mixed layer removal in bfs_conncomp1_wet
  elseif isstruct(OPTS.ML)
    REMOVE_MIXED_LAYER = true;
    ML = mixed_layer(S, T, P, OPTS.ML);
  else
    % Use a pre-computed mixed layer
    REMOVE_MIXED_LAYER = true;
    ML = OPTS.ML;
  end
end

%% Interpolate S and T casts onto surface
[~, K, N] = size(OPTS.Sppc);
if K > 0
  assert(K == nk-1 && N == nij, 'size(OPTS.Sppc) should be [O, nk-1, ni, nj] == [?, %d, %d, %d]', nk-1, ni, nj);
  Sppc = OPTS.Sppc;
else
  Sppc = OPTS.INTERPFN(P, S);
end

[~, K, N] = size(OPTS.Tppc);
if K > 0
  assert(K == nk-1 && N == nij, 'size(OPTS.Tppc) should be [O, nk-1, ni, nj] == [?, %d, %d, %d]', nk-1, ni, nj);
  Tppc = OPTS.Tppc;
else
  Tppc = OPTS.INTERPFN(P, T);
end
[s,t] = ppc_val2(P,Sppc,Tppc,lead1(p));
p(isnan(s)) = nan; % ensure same nan structure between s, t, and p. Just in case user gives, e.g., repmat(1000,ni,nj) for a 1000dbar isobaric surface

%% Prepare Diagnostics
DIAGS = (OPTS.VERBOSE > 0) || (nargout >= 8);


diags = struct();
diags.clocktime     = nan(OPTS.ITER_MAX,1);
diags.p_change_L1   = nan(OPTS.ITER_MAX,1);
diags.p_change_L2   = nan(OPTS.ITER_MAX,1);
diags.p_change_Linf = nan(OPTS.ITER_MAX,1);
diags.freshly_wet   = nan(OPTS.ITER_MAX,1);

diags.epsL1 = nan(OPTS.ITER_MAX + 1,1);
diags.epsL2 = nan(OPTS.ITER_MAX + 1,1);

diags.timer_wetting   = nan(OPTS.ITER_MAX, 1);
diags.timer_recon     = nan(OPTS.ITER_MAX, 1);
diags.timer_reebgraph = nan(OPTS.ITER_MAX, 1);
diags.timer_fitting   = nan(OPTS.ITER_MAX, 1);
diags.timer_update    = nan(OPTS.ITER_MAX, 1);

% Compute the L1 and L2 norms of the neutrality (epsilon) errors
[diags.epsL2(1), diags.epsL1(1)] = eps_norms(s, t, p, false, WRAP, {}, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij);

%% Process the reference cast:
% Get linear index to reference cast

% Set p0: the p value the surface will be pinned to at this ref cast
if isscalar(OPTS.REF_P)
  p0 = OPTS.REF_P;
else
  p0 = p(ref_cast);
end

% Get (s0,t0): (S,T) values at the reference cast at the reference p.
s0 = s(ref_cast);
t0 = t(ref_cast);
d0 = eos(s0, t0, p0); % temporary, will be over-written

% Overwrite with reference values if provided
if isscalar(OPTS.REF_S)
  s0 = OPTS.REF_S;
end
if isscalar(OPTS.REF_T)
  t0 = OPTS.REF_T;
end

% Pre-compute delta (using the reference S and T) at the reference cast
% at the reference p. If reference S and T are not provided, this is 0.
d0 = d0 - eos(s0,t0,p0);


%% Begin iterations
% Note: the surface exists wherever p is non-nan.  The nan structure of s
% and t is made to match that of p when the vertical solve step is done. 
for iter = 1 : OPTS.ITER_MAX
  iter_tic = tic();
  
  % --- Remove the Mixed Layer
  % But keep it for the first iteration, which may be initialized from a
  % not very neutral surface
  if REMOVE_MIXED_LAYER && iter > 1
    p(p < ML) = nan;
  end
  
  
  % --- Determine the connected component containing the reference cast,
  % via Breadth First Search, and do wetting while at it
  mytic = tic();
  if iter >= OPTS.ITER_START_WETTING && iter <= OPTS.ITER_STOP_WETTING
    [s, t, p, freshly_wet, qu, qt] = bfs_conncomp1_wet_mex(Sppc, Tppc, P, s, t, p, ML, OPTS.TOL_P_UPDATE, A4, BotK, ref_cast, qu);
  else
    [qu, qt] = bfs_conncomp1(isfinite(p), A4, ref_cast, qu);
    freshly_wet = 0;
  end
  

  % Keep only the component of the surface connected to the reference cast
  wet = false(ni,nj);
  wet(qu(1:qt)) = true;
  p(~wet) = nan;
  
  if DIAGS
    timer_wetting = toc(mytic);
  end
  
  
  if OPTS.REEB
    mytic = tic();
    % --- Calculate the Reeb graph:
    % 1. Pre-process to select one region, possibly filling in certain holes;
    % 2. Calculate the Reeb graph;
    % 3. Post-process to undo any hole-filling and possibly simplify the graph.
    [p, arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, nNodes, ~, ~, timer_recon] = calc_reeb_graph(p, WRAP, OPTS);
    
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
  
  % --- Build p_, the vector containing p at each valid point on the surface
  p_ = p(wet).'; % Row vector
  n_casts = length(p_);
  
  
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
  dp = eos_p(s, t, p) - eos_p(s0, t0, p);
  dp_ = dp(wet).';
  
  if iter == 6
    0;
  end
  
  % --- Fit d(delta) / d(pressure), then integrate it
  if OPTS.REEB
    % --- Fit partial derivative of delta w.r.t. p as a multivalued function of p,
    dp_fn = branches_fit(p_, dp_, arc_from, arc_to, arc_segment, node_fn, cb_arcs, cb_nodes, OPTS.GEOSTRF);
    
    % --- Integrate to get delta as a multivalued function of p
    d_fn = int_graph_fun(dp_fn, arc_from, arc_to, node_fn, graph, bfs_parent_node, bfs_topo_order, bfs_missing_arc);
    
    % Adjust d_fn so that the new surface coincides with the old one on the starting water column
    % i.e. force d_fn(:,arc(ref_cast)) evaluated at p0 to equal d0.
    d_fn(end,:) = d_fn(end,:) + (d0 - pvallin(d_fn(:,arc(ref_cast)), p0));
    
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
    % --- Fit partial derivative of delta w.r.t. p as a single-valued function of p
    breaks = OPTS.SPLINE_BREAKS ;
    breaks = min(breaks, max(p_));
    breaks = max(breaks, min(p_));
    breaks = unique(breaks);
    dp_fn = splinefit(p_, dp_, breaks, OPTS.SPLINE_ORDER);
    
    % --- Integrate to get delta as a single-valued function of p
    d_fn = ppint(dp_fn);
    
    % Adjust d_fn so that the new surface coincides with the old one on the starting water column
    % i.e. force d_fn(:,arc(ref_cast)) evaluated at p0 to equal d0.
    d_fn.coefs(:,end) = d_fn.coefs(:,end) + ( d0 - ppval(d_fn, p0) );
  end
  clear dp_fn dp_
  if DIAGS
    timer_fitting = toc(mytic);
  end
  
  % --- Solve for new pressures at which specific volume is that of the multivalued function
  mytic = tic();
  p_old = p; % Record old surface for diagnostic purposes. 
  if OPTS.REEB
    [p, s, t] = tbs_vertsolve_mex(Sppc, Tppc, P, BotK, s, t, p, arc, d_fn, s0, t0, OPTS.TOL_P_UPDATE, OPTS.P_EXPN);
  else
    [p, s, t] = obs_vertsolve_mex(Sppc, Tppc, P, BotK, s, t, p, d_fn.breaks, d_fn.coefs, s0, t0, OPTS.TOL_P_UPDATE);
  end
  if DIAGS
    timer_update = toc(mytic);
  end
  
  % --- Get ready for next iteration
  p_change = p - p_old;
  p_change_L2 = nanrms(p_change(:));
  
  % --- Closing remarks
  if DIAGS
    
    clocktime = toc(iter_tic);
    
    p_change_L1 = nanmean(abs(p_change(:)));
    p_change_Linf = max(abs(p_change(:)));
    if OPTS.VERBOSE >= 1
      fprintf(OPTS.FILE_ID, 'Iter %2d (%.2fsec):  %4d casts wet; %4d casts in/outcropped; p change has: L_1 %.6e, L_2 %.6e, L_inf %.6e\n',...
        iter, toc(iter_tic), freshly_wet, sum(isnan(p(wet))), p_change_L1, p_change_L2, p_change_Linf);
    end
    diags.clocktime(iter) = clocktime;
    
    
    % Diagnostics about what THIS iteration did
    diags.p_change_L1(iter) = p_change_L1;
    diags.p_change_L2(iter) = p_change_L2;
    diags.p_change_Linf(iter) = p_change_Linf;
    diags.freshly_wet(iter) = freshly_wet;
    
    diags.timer_wetting(iter)   = timer_wetting;
    diags.timer_recon(iter)     = timer_recon;
    diags.timer_reebgraph(iter) = timer_reebgraph;
    diags.timer_fitting(iter)   = timer_fitting;
    diags.timer_update(iter)    = timer_update;
    
    % Diagnostics about the state AFTER this iteration
    
    % Compute the L1 and L2 norms of the neutrality (epsilon) errors
    [diags.epsL2(iter+1), diags.epsL1(iter+1)] = eps_norms(s, t, p, false, WRAP, {}, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij);
    
  end
  
  % --- Check for convergence
  if (p_change_L2 < OPTS.TOL_P_CHANGE_L2) && iter >= OPTS.ITER_MIN
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