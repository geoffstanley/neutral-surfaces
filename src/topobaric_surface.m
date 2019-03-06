function [p, RG, s0, t0, dfn] = topobaric_surface(S, T, P, p, OPTS)
%TOPOBARIC_SURFACE  Create a topobaric surface.
%
%
% p = topobaric_surface(S, T, P, p)
% returns the pressure p (output) of a topobaric surface formed by an
% iterative procedure, initialized from an approximately neutral surface on
% which the pressure is p (input), in an ocean whose practical / Absolute
% salinity is S, potential / Conservative temperature is T, and pressure is
% P. The specific volume and its partial derivative with respect to
% pressure are determined by eos.m and eosdp.m in MATLAB's path. If
% topobaric_geostrf will not be called, eos.m and eosdp.m can instead
% determine the in-situ density and its derivative with respect to
% pressure. Only one connected component of the surface is processed.
%
% p = topobaric_surface(..., OPTS)
% overrides default options with those provided in OPTS (see below). Any
% option not given in OPTS will be set to its default value (found in
% private/set_defaults.m).
%
% [p, RG, s0, t0, dfn] = topobaric_surface(...)
% also returns the Reeb Graph RG, the reference practical / Absolute
% salinity s0, the reference potential / Conservative temperature t0,
% and the multivalued function for delta in terms of pressure, dfn.
% delta is the specific volume anomaly or the in-situ density anomaly,
% depending on which is given by eos.m and eosdp.m, and is determined by
%   d = eos(s, t, p) - eos(s0, t0, p)
% where s and t are S and T on the surface, given by
%   [s,t] = interp1qn2(reshape(p, [1, size(p)]), P, S, T).
% Now, d will closely match dfn evaluated at p, computed as follows:
%   dfn_at_p = nan(size(p));
%   p_ = p(RG.ocean);
%   for e = 1 : RG.nArcs
%     dfn_at_p(RG.arc_segment{e}) = pvallin(dfn(:,e), p_(RG.arc_segment{e});
%   end
%
% [p, RG, s0, t0, dfn1] = topobaric_surface(S, T, P, p, OPTS)
% with OPTS.REEB == false instead calculates an "orthobaric" surface, in
% which d is a single-valued function of the pressure, dfn1. This is done
% in just one connected ReGion of the surface, returned as the true pixels
% in RG.ocean. The full specific volume or in-situ density on the surface,
%   eos(s, t, p),
% will (less) closely match
%   ppval(dfn1, p) + eos(s0, t0, p).
%
% [z, RG, s0, t0, dfn] = topobaric_surface(S, T, Z, z, OPTS)
% as above but for a Boussinesq ocean with reference density OPTS.RHOB and
% with S and T data given at depths Z, and surfaces specified by their
% depth z. If topobaric_geostrf is to be called, eos.m and eosdp.m must
% define the in-situ density and its partial derivative w.r.t. pressure
% (not depth). Note dfn is a function for delta in terms of depth (not
% pressure). Also note Z and z are positive and increasing downwards.
%
%
% Below: nz is the maximum number of data points per water column,
%        nx is the number of data points in longitude,
%        ny is the number of data points in latitude,
%        nc is the number of water columns intersecting the surface,
%        A is the number of arcs in the Reeb graph,
%        N is the number of ndoes in the Reeb graph,
%        C is the number of cycles in the cycle basis of the Reeb graph,
%        * is a wildcard.
%
% --- Input:
% S [nz, nx, ny]: practical / Absolute salinity
% T [nz, nx, ny]: potential / Conservative temperature
% P [nz, nx, ny]: pressure [dbar]
% Z [nz, nx, ny] or [nz, 1]: depth [m, positive]
% p [nx, ny]: pressure on approximately neutral surface [dbar]
% z [nx, ny]: depth of approximately neutral surface [m, positive]
% OPTS [struct]: options (see below)
%
% Note: P and Z must increase monotonically along the first dimension.
%
%
% --- Output:
% p [nx, ny]: pressure on the topobaric surface [dbar]
% z [nx, ny]: depth of the topobaric surface [m, positive]
% RG [struct]: the Reeb Graph and the single connected ReGion (see below)
% s0 [1, 1]: reference practical / Absolute salinity
% t0 [1, 1]: reference potential / Conservative temperature
% dfn [5, A]: the multivalued function for delta in terms of pressure or
%             depth. Each column is a single-valued branch, to be evaluated
%             by pvaln or pvallin.
% dfn1 [struct]: the single-valued function for delta in terms of pressure
%                or depth, to be evaluated by ppval.
%
% Note: physical units of S, T, s0, and t0 are determined by eos.m.
% Physical units of P, p, Z, and z are also determined by eos.m, though if
% topobaric_geostrf is to be used, these must be [dbar] and [m, positive]).
%
%
% --- Options:
% OPTS is a struct containing some subset of the following fields:
% GRAV [1, 1]: gravitational acceleration [m s^-2].
% RHOB [] or [1, 1]: Empty for non-Boussinesq analysis, or a scalar
%   reference density [kg.m^-3] for Boussinesq analysis.
% REEB [1, 1]: true to compute topobaric surfaces, false to compute
%   "orthobaric" surfaces.
% SPLINE_BREAKS [vector]: knots of the single-valued
%   function, when OPTS.REEB is false. Provide SPLINE_BREAKS as [m,
%   positive] if OPTS.RHOB is provided, else [dbar].
% SPLINE_ORDER [1, 1]: order of the spline of the single-valued function,
%   when OPTS.REEB is false. E.g. 4 for cubic splines.
% MLP []: use default parameters from mixed_layer_pressure().
% MLP [struct]: override parameters from mixed_layer_pressure().
% MLP [nx, ny]: use a pre-computed mixed layer pressure [dbar].
% REF_IJ [1, 2]: pixel indices for reference water column.
% REF_P [1, 1]: override the pressure at which the topobaric surface
%   intersects the reference water column [dbar].
% REF_Z [1, 1]: override the depth    at which the topobaric surface
%   intersects the reference water column [m, positive].
% REF_S [1, 1]: override the reference practical / Absolute salinity .
% REF_T [1, 1]: override the reference potential /Conservative temperature.
% WRAP [1, 2]: determines which dimensions are treated periodic [logical].
%   Set WRAP(1) true when periodic in longitude; set WRAP(2) true when
%   periodic in latitude.
% FILL_PIX [1, 1]: fill all holes in the surface containing with fewer
%   pixels than this, before computing the Reeb graph.
% FILL_IJ = [*, 2]: fill all holes in the surface except those containing a
%   pixel index in these rows, before computing the Reeb graph.
% SIMPLIFY_ARC_REMAIN [1, 1]: number of arcs to remain after graph
%   simplification.
% SIMPLIFY_WEIGHT_PERSIST [1, 1]: weighting (between 0 and 1) of
%   persistence (the difference between pressure values associated with
%   nodes) as opposed to area (number of pixels in each assocaited region)
%   during graph simplification.
% TOL [1, 1]: error tolerance when root-finding to update surface. Provide
%   TOL as [m] if OPTS.RHOB is provided, else [dbar].
% DAMP [1, 1]: damping when updating to the new surface. No damping at
%   0. Set > 0 for damping. Set < 0 for overrelaxation.
% ITER_MAX [1, 1]: maximum number of iterations (see "Exit Points" below).
% ITER_L2_CHANGE [1, 1]: quit when the change in the root-mean-square of
%   the surface pressure or depth exceeds this value. Set to 0 to
%   deactivate. Provide ITER_L2_CHANGE as [m] if OPTS.RHOB is provided,
%   else [dbar].
% GEOSTRF [1, 1]: true to create a modified topobaric surface, which
%   applies extra constraints to the empirical function that ensure the
%   exact geostrophic stream function on the surface is well-defined.
% VERBOSE [1, 1]: 0 for no output, 1 for basic information
% FILE_ID [1, 1]: 1 to write any output to MATLAB terminal, or a file
%   identifier as returned by fopen() to write to that file.
%
%
% --- Reeb Graph
% The Reeb graph contains N nodes, A arcs, and C cycles in the cycle basis.
% It is encapsulated in the struct RG, containing the following fields:
% ocean [nx,ny]: a single connected region on the surface. ocean(i,j) is
%   true iff pixel (i,j) contributed to the empirical fit for dfn.
%   A true pixel of ocean is a "vertex" in a triangular mesh that underlies
%   the Reeb graph. The index of the vertex at pixel (i,j) is
%     ij2v(i,j);
%   where
%     ij2v = reshape(cumsum(RG.ocean(:)), nx, ny);
%     ij2v(~RG.ocean) = 0;
%   maps from pixel space to vertex space.
% n_casts [1, 1]: the number of pixels in the connected region (the number
%   of true values in ocean).
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
%     v2ij = @(v) sub2ind([nx, ny], v2I(v));
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
% --- Exit points:
% Let p_ = p(RG.ocean),
%     d_ = d(RG.ocean).
% When OPTS.ITER_MAX is not an integer, then
%   p_(RG.node_v(n)) == RG.node_fn(n)  for all 1 <= n <= N.
% When OPTS.ITER_MAX is an integer,
%   pvallin(dfn(:,e), p_(RG.arc_segment{e}))
%      very nearly equals
%   d_(RG.arc_segment{e})
%      for all 1 <= e <= A, with closer approximation for smaller OPTS.TOL.
%
%
% --- Equation of State:
% The MATLAB path must contain two files, eos.m and eosdp.m.
% In the non-Boussinesq case (OPTS.RHOB is not provided), A = eos(S,T,P)
% and AP = eosdp(S,T,P).
% In the Boussinesq case (OPTS.RHOB is provided), R = eos(S,T,P)
% and RP = eosdp(S,T,P).
% Here,
%   S is the practical/Absolute salinity,
%   T is the potential/Conservative temperature,
%   P is the pressure [dbar],
%   A is the specific volume [m^3 kg^-1], and
%   AP is the partial derivative of specific volume w.r.t. pressure [m^3 kg^-1 dbar^-1].
%   R is the in-situ density [kg m^-3], and
%   RP is the partial derivative of in-situ density w.r.t. pressure [kg m^-3 dbar^-1].
% Note, in the Boussinesq case, eos and eosdp will receive not the in-situ
% pressure but rather the Boussinesq pressure, P = Z * 1e-4 * OPTS.GRAV *
% OPTS.RHOB where Z is the depth [m, positive].
%
% The default equation of state is the JMD95 in-situ density. Others can be
% used by changing eos.m and eosdp.m. To use non-Boussinesq analysis and
% JMD95, copy << ../lib/eos/specvoljmd95_tb.m >> to << ../lib/eos.m >> and
% copy << ../lib/eos/specvoljmd95_tb_dp.m >> to << ../lib/eosdp.m >>. The 
% TEOS-10 equation of state can be used by copying the appropriate gsw_*.m
% functions from << ../lib/eos/ >> to << ../lib/eos.m >> and 
% << ../lib/eosdp.m >>.  Other equations of state may be used by constructing
% eos.m and eosdp.m similarly. Note that run_codegen.m must be re-run after
% changing eos.m and eosdp.m. You can use a different equation of state by 
% creating eos.m and eosdp.m yourself; note you may need to eliminate input 
% checks and expansion of input variables (MATLAB's automatic expansion, used
% throughout this code, handles this) for compatibility with MATLAB's codegen.
%
%
% --- Requirements (beyond the Topobaric Surface package):
% see ../README.md
%
%
% --- References:
% Stanley, G. J.  Neutral surface topology. Ocean Modelling, submitted.

% --- Copyright:
% Copyright 2019 Geoff Stanley
%
% This file is part of Topobaric Surface.
%
% Topobaric Surface is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% Topobaric Surface is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with Topobaric Surface.  If not, see
% <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com
% Version   : 1.0
%
% Modified by : --
% Date        : --
% Changes     : --

% --- Notes on the code:
% 3D scalar fields [nz,nx,ny] are upper case letters, e.g. P
% 2D scalar fields [nx,ny] are lower case letters, e.g. p
% 1D scalar fields [1,nc] are lower case letters followed by an underscore,
% e.g. p_.
% The data-heavy variables S, T, and P are not modified internally (except
% to ensure double precision, and if P is actually Z in the Boussinesq
% case, it is converted to P).
% The specific volume is denoted by a and A, for alpha (rather than v and V
% for nu, to avoid confusion with the meridional velocity).
% The specific volume anomaly or the in-situ density anomaly is denoted d
% (for delta).
% The partial derivative of d with respect to pressure is denoted dp.

%% Simple checks and preparations
narginchk(4,5);
P = double(P);
S = double(S);
T = double(T);
p = double(p);
if nargin < 5
    OPTS = struct();
end

[nx,ny] = size(p);
nz = size(P,1);
is1D = @(F) ismatrix(F) && all(size(F) == [nz,1]);
is2D = @(F) ismatrix(F) && all(size(F) == [nx,ny]);
is3D = @(F) ndims(F) == 3 && all(size(F) == [nz,nx,ny]);


% An order [depth, latitude, longitude] is also okay.
assert(is2D(p), 'p must have 2 dimensions, ordered: longitude, latitude');
assert(is1D(P) || is3D(P), 'P must be vector of length nz, or of size [nz, nx, ny]');

if exist('interp1qn2_mex', 'file')
    interp1qn2_ = @interp1qn2_mex;
else
    interp1qn2_ = @interp1qn2;
end
if exist('topobaric_update_mex', 'file')
    topobaric_update_ = @topobaric_update_mex;
else
    topobaric_update_ = @topobaric_update;
end
if exist('orthobaric_update_mex', 'file')
    orthobaric_update_ = @orthobaric_update_mex;
else
    orthobaric_update_ = @orthobaric_update;
end

lead1 = @(x) reshape(x, [1 size(x)]);
nanrms = @(x) sqrt(nanmean(x(:) .* conj(x(:))));

warning('off', 'MATLAB:nargchk:deprecated') % splinefit uses nargchk

%% Process OPTS

% Load default options, then override with any user-specified OPTS.
OPTS = catstruct(set_defaults(nx,ny), OPTS);

assert(isfield(OPTS, 'WRAP') && isvector(OPTS.WRAP) && length(OPTS.WRAP) == 2, 'OPTS.WRAP must be provided as a vector of length 2');

BOUSSINESQ = isfield(OPTS, 'RHOB');
if BOUSSINESQ
    assert(isscalar(OPTS.RHOB), 'OPTS.RHOB, if provided, must be a scalar');
    assert(isfield(OPTS, 'GRAV') && isscalar(OPTS.GRAV), 'With OPTS.RHOB provided, OPTS.GRAV must be provided as a scalar');
    % P and p are actually Z and z, as are some OPTS.  Internally convert
    % these (and related options) from [m] to [dbar] internally for use
    % with eos() and eosdp().
    Z2P = 1e-4 * OPTS.RHOB * OPTS.GRAV;
    p = p * Z2P;
    P = P * Z2P;
    OPTS.TOL = OPTS.TOL * Z2P;
    OPTS.ITER_L2_CHANGE = OPTS.ITER_L2_CHANGE * Z2P;
    OPTS.SPLINE_BREAKS = OPTS.SPLINE_BREAKS * Z2P;
    OPTS.REF_P = OPTS.REF_Z * Z2P;
end

% Whether the rectangular data contains all vertices in the simplical decomposition
OPTS.ALL_VERTS = OPTS.DECOMP(1) == 'd' && OPTS.FILL_PIX == 0 && isempty(OPTS.FILL_IJ);
assert(OPTS.ALL_VERTS || OPTS.ITER_MAX < 1, ...
    'To empirically fit the multivalued function, must use ''diagonal'' simplical decomposition and fill no holes during pre-processing');

if OPTS.ITER_MAX >= 1
    assert(all(size(S) == size(T)), 'S and T must be the same size.');
    assert(is3D(S), 'S and T must have 3 dimensions [depth, longitude, latitude]');
    
    % Get number of valid data points per cast.  Used in topobaric_update
    K = squeeze(sum(isfinite(S), 1));
    
    
    % --- Process the reference cast:
    % Get linear index to reference cast
    assert(numel(OPTS.REF_IJ) == 2, 'OPTS.REF_IJ must be a vector of length 2');
    I0 = sub2ind([nx,ny], OPTS.REF_IJ(1), OPTS.REF_IJ(2));
    
    % Select the reference pressure if provided, or use the pressure of the
    % surface at the reference cast
    if isscalar(OPTS.REF_P)
        p0 = OPTS.REF_P;
    else
        p0 = p(I0);
    end
    
    % Get reference values at the reference cast at the reference pressure.
    if isvector(P)
        [s0, t0] = interp1qn2(p0, P, S(:,I0), T(:,I0));
    else
        [s0, t0] = interp1qn2(p0, P(:,I0), S(:,I0), T(:,I0));
    end
    a0 = eos(s0, t0, p0);
    
    % Overwrite with reference values if provided
    if isscalar(OPTS.REF_S)
        s0 = OPTS.REF_S;
    end
    if isscalar(OPTS.REF_T)
        t0 = OPTS.REF_T;
    end
    
    % Pre-compute delta (using the reference S and T) at the reference cast
    % at the reference pressure. If reference S and T are not provided,
    % this is 0.
    d0 = a0 - eos(s0,t0,p0);
    
    
    % Get MLP: the pressure of the mixed layer
    if OPTS.ITER_MAX > 1
        if isempty(OPTS.MLP)
            % Do not remove the mixed layer
            MLP = [];
        elseif isstruct(OPTS.MLP)
            MLP = mixed_layer_pressure(S, T, P, OPTS.MLP);
        else
            % Use a pre-computed mixed layer
            MLP = OPTS.MLP;
        end
    end
    
else
    s0 = [];
    t0 = [];
    dfn = [];
end



%% Begin iterations
for iter = 1 : (OPTS.ITER_MAX+.5)
    itertic = tic();
    
    % --- Remove the Mixed Layer
    % But keep it for the first iteration, which may be initialized from a not very neutral surface,
    if iter > 1 && ~isempty(MLP)
        p(p < MLP) = nan;
    end
    
    
    if OPTS.REEB
        % --- Calculate the Reeb graph:
        % 1. Pre-process to select one region, possibly filling in certain holes;
        % 2. Calculate the Reeb graph;
        % 3. Post-process to undo any hole-filling and possibly simplify the graph.
        [p, ocean, arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, nNodes] ...
            = calc_Reeb_graph(p, OPTS);
        
        % --- Prepare info about cycles
        [~, graph, ~, bfs_parent_node, bfs_topo_order, bfs_missing_arc, cb_arcs, cb_nodes] = ...
            cycle_analy_bfs(arc_from, arc_to, nNodes);
    else
        % Select one connected component of the surface
        ocean = select_conn_cpt(isfinite(p), OPTS.WRAP, OPTS.REF_IJ);
    end
    
    % --- Build p_, the vector containing p at each valid point on the surface
    p_ = p(ocean).'; % Row vector
    n_casts = length(p_);
    
    % --- Possible exit point, if user just wants the RG, perhaps for
    % computing the topobaric geostrophic stream function.
    if iter > OPTS.ITER_MAX
        break
    end
    
    if OPTS.REEB
        % --- Build vector associating each water column with an arc
        % Note some water columns (saddles) are associated with multiple arcs!
        % arcvec simply takes one of them.
        arc_ = zeros(1, n_casts);
        for e = 1 : nArcs
            arc_(arc_segment{e}) = e;
        end % Now all arc_ > 0
        arc = zeros(nx,ny);
        arc(ocean) = arc_;
    end
    
    % --- Get derivative of delta on surface
    [s,t] = interp1qn2_(lead1(p), P, S, T);
    dp = eosdp(s, t, p) - eosdp(s0, t0, p);
    dp_ = dp(ocean).';
    clear s t dp
    
    % --- Fit d(delta) / d(pressure), then integrate it
    if OPTS.REEB
        % --- Fit partial derivative of delta w.r.t. pressure as a multivalued function of pressure,
        dpfn = branches_fit(p_, dp_, arc_from, arc_to, arc_segment, node_fn, cb_arcs, cb_nodes, OPTS.GEOSTRF);
        
        % --- Integrate to get delta as a multivalued function of pressure
        dfn = int_graph_fun(dpfn, arc_from, arc_to, node_fn, graph, bfs_parent_node, bfs_topo_order, bfs_missing_arc);
        clear dpfn
        
        % Adjust dfn so that the new surface coincides with the old one on the starting water column
        % i.e. force dfn(:,arc(I0)) evaluated at p0 to equal d0.
        dfn(end,:) = dfn(end,:) + (d0 - pvallin(dfn(:,arc(I0)), p0));
        
        % Double check matching conditions at nodes of the Reeb Graph.
        %{
        for n = 1:nNodes
            neigharcs = [node_prev{n}; node_next{n}].';
            e = neigharcs(1);
            node_p = node_fn(n);
            dfn_at_e = pvaln(dfn(:,e), node_p);
            for ee = neigharcs(2:end)
                diff = dfn_at_e - pvaln(dfn(:,ee), node_p);
                if abs(diff) > 1e-9
                    fprintf(OPTS.FILE_ID, 'Bad Matching (diff=%e) at node %d, arcs %d and %d\n', diff, n, e, ee);
                end
            end
        end
        clear neigharcs node_p dfn_at_e diff ee e
        %}
    else
        % --- Fit partial derivative of delta w.r.t. pressure as a single-valued function of pressure
        breaks = OPTS.SPLINE_BREAKS ;
        breaks = min(breaks, max(p_));
        breaks = max(breaks, min(p_));
        breaks = unique(breaks);
        dpfn = splinefit(p_, dp_, breaks, OPTS.SPLINE_ORDER);
        
        % --- Integrate to get delta as a single-valued function of pressure
        dfn = ppint(dpfn);
        clear dpfn
        
        % Adjust dfn so that the new surface coincides with the old one on the starting water column
        % i.e. force dfn(:,arc(I0)) evaluated at p0 to equal d0.
        dfn.coefs(:,end) = dfn.coefs(:,end) + ( d0 - ppval(dfn, p0) );
    end
    clear dp_
    
    % --- Solve for new pressures at which specific volume is that of the multivalued function
    if OPTS.REEB
        pnew = topobaric_update_(S, T, P, K, p, arc, dfn, s0, t0, OPTS.TOL);
    else
        pnew = orthobaric_update_(S, T, P, K, p, dfn.breaks, dfn.coefs, s0, t0, OPTS.TOL);
    end
    
    
    % --- Get ready for next iteration
    pdiff = pnew - p;
    L_2_p_diff = nanrms(pdiff(:));
    
    % --- Check for convergence or maximum iterations
    STOP = (iter >= OPTS.ITER_MAX) || (L_2_p_diff < OPTS.ITER_L2_CHANGE) ;
    
    if STOP || OPTS.DAMP == 0
        p = pnew; % Definitely want the exact solution, when STOPping
    else
        p = p + (1 - OPTS.DAMP) * pdiff;
    end
    
    % --- Closing remarks
    if OPTS.VERBOSE >= 1
        if BOUSSINESQ
            fprintf(OPTS.FILE_ID, 'Iter %02d (%.2fsec): %d casts unsolved; z change has: Mean %0.5g, L_2 %0.5g, L_inf %0.5g\n',...
                iter, toc(itertic), sum(isnan(pnew(ocean))), nanmean(pdiff(:)) / Z2P, L_2_p_diff / Z2P, max(abs(pdiff(:))) / Z2P);
        else
            fprintf(OPTS.FILE_ID, 'Iter %02d (%.2fsec): %d casts unsolved; p change has: Mean %0.5g, L_2 %0.5g, L_inf %0.5g\n',...
                iter, toc(itertic), sum(isnan(pnew(ocean))), nanmean(pdiff(:)), L_2_p_diff, max(abs(pdiff(:))));
        end
    end
    clear pdiff pnew
    
    
    % --- Primary exit point for generating topobaric surfaces. 
    % Presently, d approximately equals dfn evaluated at (with better
    % approximation for lower OPTS.TOL)
    if STOP
        break
    end
    
end

if BOUSSINESQ
    p = p / Z2P; % Convert back to depth
end

if nargout < 2
    return
end

% Return Reeb Graph (if it was computed) and ReGion
RG = struct();
RG.ocean = ocean;
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
    if BOUSSINESQ
        RG.node_fn     = node_fn / Z2P; % Convert back to depth
    else
        RG.node_fn     = node_fn;
    end
    
    % Add extra info to RG
    RG.cb_nodes        = cb_nodes;
    RG.cb_arcs         = cb_arcs;
    RG.graph           = graph;
    RG.bfs_parent_node = bfs_parent_node;
    RG.bfs_topo_order  = bfs_topo_order;
    RG.bfs_missing_arc = bfs_missing_arc;
end

if nargout < 5
    return
end

if BOUSSINESQ && OPTS.ITER_MAX >= 1
    % Convert the domain of dfn back to depth
    if OPTS.REEB
        dfn(1,:) = dfn(1,:) / Z2P;
        dfn(2,:) = dfn(2,:) / Z2P;
        dfn(3,:) = dfn(3,:) * Z2P^2;
        dfn(4,:) = dfn(4,:) * Z2P;
    else
        fac = 1;
        for i = size(dfn.coefs,2)-1 : -1 : 1
            fac = fac * Z2P;
            dfn.coefs(:,i) = dfn.coefs(:,i) * fac;
        end
    end
end

end