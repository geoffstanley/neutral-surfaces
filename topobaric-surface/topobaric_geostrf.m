function [gsf, gsfdiff] = topobaric_geostrf(s, t, p, P, A, Y, s0, t0, ref_cast, OPTS, varargin)
%TOPOBARIC_GEOSTRF  The topobaric geostrophic stream function.
%
%
% gsf = topobaric_geostrf(s, t, p, P, A, Y, s0, t0, ref_cast, OPTS)
% computes the topobaric geostrophic stream function gsf on a surface with
% practical / Absolute salinity, potential / Conservative temperature, and
% pressure given by s, t, and p respectively, and with s and t reference
% values given by s0 and t0, respectively, in an ocean with pressure P,
% specific volume A, and hydrostatic acceleration potential Y, and with
% equation of state for the specific volume given by eos.m in the path. The
% specific volume should be calculated by A = eos(S, T, P) where S is the
% practical / Absolute salinity and T is the potential / Conservative
% temperature at the data sites of P. The hydrostatic acceleration
% potential should be calculated as Y = hsap3(P, ATMP, ETAN, A, grav),
% where ATMP is the atmospheric pressure, ETAN is the sea-surface height,
% and grav the gravitational acceleration. Inputs s and t may instead be
% piecewise polynomial interpolants for S and T in terms of P in each water
% column, in which case s and t are determined by evaluating these
% interpolants onto p. Algorithmic parameters are determined by OPTS.
% Providing OPTS.RG and OPTS.d_fn defines a multivalued function for the
% specific volume anomaly (the specific volume minus that using s0 and t0)
% on the surface in terms of the pressure on the surface; in this case,
% OPTS.RG, s0, t0, and OPTS.d_fn should be the final four outputs of
% topobaric_surface.m called with OPTS.REEB true. If OPTS.RG and OPTS.d_fn
% are not provided, the Reeb graph is internally calculated and the
% specific volume anomaly on the surface is empirically fit as a
% multivalued function of pressure on the surface. The surface is treated
% as periodic in its first dimension (longitude) if OPTS.WRAP(1) is true,
% and likewise treated periodic in its second dimension (latitude) if
% OPTS.WRAP(2) is true.
%
% [gsf, gsfdiff] = topobaric_geostrf(s, t, p, P, A, Y, s0, t0, ref_cast, OPTS)
% also returns the difference, gsfdiff, between gsf and an alternative
% geostrophic stream function that is obtained like gsf but arcs in the
% Reeb graph that define cycles in the cycle basis (i.e. those arcs that
% are not traversed in a breadth-first search) have their corresponding
% branches of the multivalued function integrated in the opposite direction
% to how gsf integrated them. If OPTS.RG, s0, t0, and OPTS.d_fn are provided
% using output from topobaric_surface.m that was called with OPTS.GEOSTRF =
% true, or if OPTS.GEOSTRF here is true, then the topobaric geostrophic
% stream function is well-defined, meaning gsfdiff is zero (up to a
% numerical precision determined by LSQLIN).
%
% gsf = topobaric_geostrf(s, t, p, P, A, Y, s0, t0, ref_cast, OPTS) with OPTS.REEB false
% computes the orthobaric geostrophic stream function gsf, with s, t, p, P,
% A, Y, s0, t0, and eos.m as above. Providing OPTS.RG and OPTS.d_fn defines
% a single-valued function for the specific volume anomaly (as above) on
% the surface in terms of the pressure on the surface; in this case,
% OPTS.RG, s0, t0, and OPTS.d_fn should be the final four outputs of
% topobaric_surface.m called with OPTS.REEB false. If OPTS.RG and OPTS.d_fn
% are not provided, the specific volume anomaly is empirically fit as a
% single-valued function of pressure on the surface, namely as a spline of
% order OPTS.SPLINE_ORDER having breaks at pressures given by
% OPTS.SPLINE_BREAKS.
%
% For a Boussinesq ocean, do the following:
% - replace P with Z, the depth of the R and Y data;
% - replace p with z, the depth of the surface;
% - both z and Z must be positive and increase downwards;
% - provide two additional parameters after OPTS, grav and rho_c;
% - replace A with R, the in-situ density, pre-computed as
%   R = eos(S, T, Z);
% - pre-compute Y as Y = hsap3(Z, ATMP, ETAN, R, grav, rho_c);
% - ensure eos.m defines the in-situ density as a function of practical /
%   Absolute salinity, potential / Conservative temperature, and pressure
%   in dbar.
%
%
% --- Input:
% s [ni, nj] or [O, nk-1, ni, nj]: the practical / Absolute salinity,
%  on the surface, or a piecewise polynomial interpolant in each water column
% t [ni, nj] or [O, nk-1, ni, nj]: the potential / Conservative temperature
%  on the surface, or a piecewise polynomial interpolant in each water column
% p [ni, nj]: pressure [dbar] or depth [m] on surface
% P [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m] at all data sites [dbar]
% A [nk, ni, nj]: specific volume [m^3 kg^-1] or in-situ density [kg m^-3]
% Y [nk, ni, nj]: acceleration potential from hydrostatic balance [m^2 s^-2]
% s0 [1, 1] or []: reference s value
% t0 [1, 1] or []: reference t value
% ref_cast [1 or 2 element array]: index (linear or 2D) to the reference cast
% OPTS [struct]: algorithmic parameters.
%   OPTS.WRAP [2 element array]: periodicity of the surface in each horizontal
%     dimension
%   REEB [1, 1]: true to compute the topobaric geostrophic stream
%     function, or false to compute the orthobaric geostrophic stream
%     function
%   GEOSTRF [1, 1]: false to allow the geostrophic stream function
%     to be ill-defined; has no effect when called with OPTS.RG and
%     OPTS.d_fn
%   RG [struct]: the Reeb Graph or ReGion, as output from
%     topobaric_surface.m
%   d_fn [5, na]: each branch of the multivalued function for specific
%     volume anomaly or in-situ density anomaly, as output from
%     topobaric_surface.m
%   d_fn [struct]: a piecewise polynomial for the single-valued specific
%     volume anomaly or in-situ density anomaly, as may be evaluated by
%     ppval.m
%   SPLINE_BREAKS [vector]: break points of the single-valued function for
%     specific volume anomaly or in-situ density anomaly as a function of
%     pressure or depth; only matters when OPTS.REEB is false. [dbar or m,
%     positive]
%   SPLINE_ORDER [1, 1]: order of the single-valued function for specific
%     volume anomaly or in-situ density anomaly as a function of pressure
%     or depth; only matters when OPTS.REEB is false
% grav [1, 1]: gravitational acceleration [m s^-2]
% rho_c [1, 1]: Boussinesq reference density [kg m^-3]
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude,
%       na is the number of arcs in the Reeb graph.
%
% Note: physical units for s, t, s0, and t0 are determined by eos.m.
%
%
% --- Output:
% gsf [ni, nj]: geostrophic stream function [m^2 s^-2]
% gsfdiff [ni, nj]: difference between gsf and the alternate geostrophic
%                   stream function [m^2 s^-2]
%
%
% --- References:
% Stanley, G.J., 2019b. The exact geostrophic streamfunction for neutral
% surfaces. Ocean Modelling 138, 107–121.
% https://doi.org/10.1016/j.ocemod.2019.04.002

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
% MERCHANTABILITY or FITNESS FOR M PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


db2Pa = 1e4;

% Input checking
narginchk(10,12);
BOUSSINESQ = (nargin == 12);
assert(isstruct(OPTS), 'OPTS must be a struct');
assert(isfield(OPTS, 'WRAP') && numel(OPTS.WRAP) == 2, 'OPTS.WRAP must be a vector of length 2');

[ni,nj] = size(p);
nk = size(P,1);
is4D = @(F) ndims(F) == 4 && size(F,2) == nk-1 && size(F,3) == ni && size(F,4) == nj;
lead1 = @(x) reshape(x, [1 size(x)]);

% Set defaults:
DEFAULT = struct();
DEFAULT.REEB = true;
DEFAULT.GEOSTRF = true;
DEFAULT.SPLINE_BREAKS = [0 200 1500 1800 6000]; % knots at these pressures or (positive) depths
DEFAULT.SPLINE_ORDER = 4; % Cubic
DEFAULT.FILL_PIX = 0;
DEFAULT.FILL_IJ = [];
DEFAULT.SIMPLIFY_ARC_REMAIN = inf;
DEFAULT.DECOMP = 'diagonal';

% Override with any user-specified OPTS.
OPTS = catstruct(DEFAULT, OPTS);

% Get linear index to reference cast
if isscalar(ref_cast)
  assert(ref_cast >= 1 && ref_cast <= ni*nj, 'Out of bounds Linear index for ref_cast.');
elseif numel(ref_cast) == 2
  assert(all(ref_cast >= 1) && all(ref_cast(:) <= [ni; nj]), 'ref_cast must index a cast within the domain.')
  ref_cast = sub2ind([ni, nj], ref_cast(1), ref_cast(2)); % Convert into linear index to the reference cast
else
  assert(false, 'ref_cast must be a 1 or 2 element vector');
end

if OPTS.REEB
    assert(isvector(OPTS.SPLINE_BREAKS), 'OPTS.SPLINE_BREAKS must be a vector');
    assert(isscalar(OPTS.SPLINE_ORDER), 'OPTS.SPLINE_ORDER must be a scalar');
end

if is4D(s) % Evaluate interpolants for S and T onto the surface
    [s,t] = ppc_val2(P, s, t, lead1(p));
end


if isfield(OPTS, 'd_fn') && isfield(OPTS, 'RG')
    d_fn = OPTS.d_fn;
    RG = OPTS.RG;
    wet = RG.wet;
    if OPTS.REEB
        arc_from        = RG.arc_from;
        arc_to          = RG.arc_to;
        arc_segment     = RG.arc_segment;
        nArcs           = RG.nArcs;
        node_fn         = RG.node_fn;
        cb_nodes        = RG.cb_nodes;
        cb_arcs         = RG.cb_arcs;
        graph           = RG.graph;
        bfs_parent_node = RG.bfs_parent_node;
        bfs_topo_order  = RG.bfs_topo_order;
        bfs_missing_arc = RG.bfs_missing_arc;
    end
else
    
    % The function (d_fn) was not given.  Prepare for its later
    % computation, by computing the Reeb graph of p, and the specific
    % volume anomaly on the surface
    d_fn = [];
    
    % Find the connected component containing the reference cast, using 4-connectivity
    [~,~,~,wet] = bfs_conncomp(isfinite(p), grid_adjacency([ni,nj], 4, OPTS.WRAP), ref_cast, 4);
    p(~wet) = nan;
    
    % Calculate the Reeb graph
    [p, arc_from, arc_to, arc_segment, ~, ~, ~, node_fn, ~, nArcs, nNodes] = calc_reeb_graph(p, OPTS);
    
    % Prepare info about cycles
    [~, graph, ~, bfs_parent_node, bfs_topo_order, bfs_missing_arc, cb_arcs, cb_nodes] = ...
        cycle_analy_bfs(arc_from, arc_to, nNodes);
    
end

% Initialize the geostrophic stream function with the easy terms
[y , a ] = hsap2(s , t , p, P, A, Y, varargin{:});
[y0, a0] = hsap2(s0, t0, p, varargin{:});
gsf = y - y0;

if isempty(d_fn)
    % Calculate the specific volume anomaly or the in-situ density anomaly
    d = a - a0;
    d_ = d(wet).';
end
p_ = p(wet).'; % Row vector
n_casts = length(p_);

% Integrate the specific volume (in-situ density) as a function of pressure
% (depth) in a non-Boussinesq (Boussinesq) ocean.
if OPTS.REEB
    
    if isempty(d_fn)
        % --- Fit specific volume as a multivalued function of pressure,
        
        if ~OPTS.GEOSTRF
            % Requested to make the stream function ill-defined. Overwrite
            % cycle information so that all branches will be fit
            % independently, with no cycle constraints. Note, this does
            % nothing if d_fn was provided. But if OPTS.GEOSTRF was false
            % when calling topobaric_surface(), then d_fn produced by
            % topobaric will similarly yield an ill-defined stream
            % function.
            cb_nodes = {};
            cb_arcs = {};
        end
        
        d_fn = branches_fit(p_, d_, arc_from, arc_to, arc_segment, node_fn, cb_arcs, cb_nodes, false);
    end
    
    % Integrate the multivalued function
    if nargout < 2
        intfn          = int_graph_fun(d_fn, arc_from, arc_to, node_fn, graph, bfs_parent_node, bfs_topo_order, bfs_missing_arc);
    else
        [intfn,intfn2] = int_graph_fun(d_fn, arc_from, arc_to, node_fn, graph, bfs_parent_node, bfs_topo_order, bfs_missing_arc);
    end
    
    % Evaluate the integral everywhere, using the local branch
    intval_ = nan(n_casts, 1);
    for e = 1:nArcs
        seg = arc_segment{e};
        intval_(seg) = pvaln(intfn(:,e), p_(seg)); % With no graph simplification, we can use pvaln instead of pvallin
    end
    
else
    
    if isempty(d_fn)
        
        % --- Fit specific volume as a single-valued function of pressure
        breaks = OPTS.SPLINE_BREAKS;
        breaks = min(breaks, max(p_));
        breaks = max(breaks, min(p_));
        breaks = unique(breaks);
        d_fn = splinefit(p_, d_, breaks, OPTS.SPLINE_ORDER);
    end
    
    % Integrate the single-valued function
    intval_ = ppint(d_fn, p_(1), p_);
    intval_ = intval_(:);
    
end

if BOUSSINESQ
    grav = varargin{1};
    rho_c = varargin{2};
    if isfield(OPTS, 'GRAV') && OPTS.GRAV ~= grav
        warning('topobaric_geostrf:gravMismatch', 'grav and OPTS.GRAV are not equal.')
    end
    if isfield(OPTS, 'RHOB') && OPTS.RHOB ~= rho_c
        warning('topobaric_geostrf:rhobMismatch', 'rho_c and OPTS.RHOB are not equal.')
    end
    fac = -grav / rho_c; % (-) because p is depth > 0!
else
    fac = db2Pa;
end

% Finish computing the geostrophic stream function
gsf(wet) = gsf(wet) + intval_ * fac;

if nargout < 2
    return
end

if ~OPTS.REEB
    % No possibility for different gsf's if using single-valued functions
    gsfdiff = zeros(ni,nj);
    gsfdiff(isnan(gsf)) = nan;
    return
end

% Calculate the difference between the above gsf and the alternative
% topobaric gsf.  The difference is so slight that a difference is returned
% rather than the alternative topobaric gsf itself, to avoid machine
% rounding error.
gsfdiff_ = zeros(n_casts, 1);
for i = 1:length(bfs_missing_arc)
    e = bfs_missing_arc(i);
    seg = arc_segment{e};
    intval2_seg = pvaln(intfn2(:,i), p_(seg)); % With no graph simplification, we can use pvaln instead of pvallin
    gsfdiff_(seg) = (intval2_seg(:) - intval_(seg)) * fac;
end
gsfdiff = zeros(ni,nj);
gsfdiff(isnan(gsf)) = nan;
gsfdiff(wet) = gsfdiff_;

end