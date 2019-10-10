function [gsf, gsfdiff] = topobaric_geostrf(s, t, x, X, M, Y, s0, t0, OPTS, varargin)
%TOPOBARIC_GEOSTRF  The topobaric geostrophic stream function.
%
%
% gsf = topobaric_geostrf(s, t, p, P, A, Y, s0, t0, OPTS)
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
% and grav the gravitational acceleration. Inputs s and t may instead be S
% and T, in which case s and t is determined by linear interpolation of S
% and T through P onto p. Algorithmic parameters are determined by OPTS.
% Providing OPTS.RG and OPTS.dfn defines a multivalued function for the
% specific volume anomaly (the specific volume minus that using s0 and t0)
% on the surface in terms of the pressure on the surface; in this case,
% OPTS.RG, s0, t0, and OPTS.dfn should be the final four outputs of
% topobaric_surface.m called with OPTS.REEB true. If OPTS.RG and OPTS.dfn
% are not provided, the Reeb graph is internally calculated and the
% specific volume anomaly on the surface is empirically fit as a
% multivalued function of pressure on the surface. The surface is treated
% as periodic in its first dimension (longitude) if OPTS.WRAP(1) is true,
% and likewise treated periodic in its second dimension (latitude) if
% OPTS.WRAP(2) is true.
%
% [gsf, gsfdiff] = topobaric_geostrf(s, t, p, P, A, Y, s0, t0, OPTS)
% also returns the difference, gsfdiff, between gsf and an alternative
% geostrophic stream function that is obtained like gsf but arcs in the
% Reeb graph that define cycles in the cycle basis (i.e. those arcs that
% are not traversed in a breadth-first search) have their corresponding
% branches of the multivalued function integrated in the opposite direction
% to how gsf integrated them. If OPTS.RG, s0, t0, and OPTS.dfn are provided
% using output from topobaric_surface.m that was called with OPTS.GEOSTRF =
% true, or if OPTS.GEOSTRF here is true, then the topobaric geostrophic
% stream function is well-defined, meaning gsfdiff is zero (up to a
% numerical precision determined by LSQLIN).
%
% gsf = topobaric_geostrf(s, t, p, P, A, Y, s0, t0, OPTS) with OPTS.REEB false
% computes the orthobaric geostrophic stream function gsf, with s, t, p, P,
% A, Y, s0, t0, and eos.m as above. Providing OPTS.RG and OPTS.dfn defines
% a single-valued function for the specific volume anomaly (as above) on
% the surface in terms of the pressure on the surface; in this case,
% OPTS.RG, s0, t0, and OPTS.dfn should be the final four outputs of
% topobaric_surface.m called with OPTS.REEB false. If OPTS.RG and OPTS.dfn
% are not provided, the specific volume anomaly is empirically fit as a
% single-valued function of pressure on the surface, namely as a spline of
% order OPTS.SPLINE_ORDER having breaks at pressures given by
% OPTS.SPLINE_BREAKS.
%
% For a Boussinesq ocean, do the following:
% - replace P with Z, the depth of the R and Y data;
% - replace p with z, the depth of the surface;
% - both z and Z must be positive and increase downwards;
% - provide two additional parameters, grav and rho_c;
% - replace A with R, the in-situ density, pre-computed as 
%   R = eos(S, T, Z, 1e-4 * grav * rho_c); 
% - pre-compute Y as Y = hsap3(Z, ATMP, ETAN, R, grav, rho_c);
% - ensure eos.m defines the in-situ density as a function of practical /
%   Absolute salinity, potential / Conservative temperature, and pressure 
%   in dbar.
%
%
% --- Input:
% s [nx, ny] or [nz, nx, ny]: the practical / Absolute salinity, 
%  on the surface, or at all data sites
% t [nx, ny] or [nz, nx, ny]: the potential / Conservative temperature
%  on the surface, or at all data sites
% p [nx, ny]: pressure on surface [dbar]
% z [nx, ny]: depth on surface [m, positive]
% P [nz, nx, ny] or [nz, 1]: pressure at all data sites [dbar]
% Z [nz, nx, ny] or [nz, 1]: depth at all data sites [m, positive]
% A [nz, nx, ny]: specific volume [m^3 kg^-1]
% R [nz, nx, ny]: in-situ density [kg m^-3]
% Y [nz, nx, ny]: acceleration potential from hydrostatic balance [m^2 s^-2]
% s0 [1, 1] or []: reference s value
% t0 [1, 1] or []: reference t value
% OPTS [struct]: algorithmic parameters.
%   OPTS.WRAP [2, 1]: periodicity of the surface in each horizontal 
%     dimension
%   REEB [1, 1]: true to compute the topobaric geostrophic stream
%     function, or false to compute the orthobaric geostrophic stream
%     function
%   GEOSTRF [1, 1]: false to allow the geostrophic stream function
%     to be ill-defined; has no effect when called with OPTS.RG and 
%     OPTS.dfn
%   RG [struct]: the Reeb Graph or ReGion, as output from 
%     topobaric_surface.m
%   dfn [5, na]: each branch of the multivalued function for specific
%     volume anomaly or in-situ density anomaly, as output from
%     topobaric_surface.m
%   dfn [struct]: a piecewise polynomial for the single-valued specific
%     volume anbomaly or in-situ density anomaly, as may be evaluated by
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
% Note: nz is the maximum number of data points per water column,
%       nx is the number of data points in longitude,
%       ny is the number of data points in latitude,
%       na is the number of arcs in the Reeb graph.
%
% Note: physical units for s, t, s0, and t0 are determined by eos.m.
%
%
% --- Output:
% gsf [nx, ny]: geostrophic stream function [m^2 s^-2]
% gsfdiff [nx, ny]: difference between gsf and the alternate geostrophic
%                   stream function [m^2 s^-2]
%
%
% --- Requirements:
% topobaric_surface, pvaln, hsap2
% splinefit, ppint - https://www.mathworks.com/matlabcentral/fileexchange/13812
%
%
% --- References:
% Stanley, G. J.  The exact geostrophic stream function for neutral
% surfaces. Ocean Modelling, submitted.

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


db2Pa = 1e4;

% Input checking
narginchk(9,11);
BOUSSINESQ = (nargin == 11);

[nx,ny] = size(x);
nz = size(X,1);
is3D = @(F) ndims(F) == 3 && all(size(F) == [nz,nx,ny]);
lead1 = @(x) reshape(x, [1 size(x)]);

assert(isstruct(OPTS), 'OPTS must be a struct');
if ~isfield(OPTS, 'REEB')
    OPTS.REEB = true;
end
if ~isfield(OPTS, 'GEOSTRF')
    OPTS.GEOSTRF = true;
end
if ~isfield(OPTS, 'SPLINE_BREAKS') % Only used for orthobaric gsf
    OPTS.SPLINE_BREAKS = [0 200 1500 1800 6000]; % knots at these pressures or (positive) depths
end
if ~isfield(OPTS, 'SPLINE_ORDER') % Only used for orthobaric gsf
    OPTS.SPLINE_ORDER = 4; % Cubic
end
assert(isfield(OPTS, 'WRAP') && numel(OPTS.WRAP) == 2, 'OPTS.WRAP must be a vector of length 2');
assert(isvector(OPTS.SPLINE_BREAKS), 'OPTS.SPLINE_BREAKS must be a vector');
assert(isscalar(OPTS.SPLINE_ORDER), 'OPTS.SPLINE_ORDER must be a scalar');


if is3D(s) % Interpolate 3D s and t onto the surface
    try
        [s,t] = interp1qn2_mex(lead1(x), X, s, t);
    catch
        [s,t] = interp1qn2(lead1(x), X, s, t);
    end
end


if isfield(OPTS, 'dfn') && isfield(OPTS, 'RG')
    
    dfn = OPTS.dfn;
    RG = OPTS.RG;
    
else
    
    % The function (dfn) was not given.  Prepare for its later
    % computation, by computing the Reeb graph of p, and the specific
    % volume anomaly on the surface
    dfn = [];
    
    % Compute the Reeb graph.
    OPTS.ITER_MAX = 0.5; % Ensure this
    [x, RG] = topobaric_surface([], [], X, x, OPTS);
    
end

% Initialize the geostrophic stream function with the easy terms
[y , m ] = hsap2(s , t , x, X, M, Y, varargin{:});
[y0, m0] = hsap2(s0, t0, x, varargin{:});
gsf = y - y0;

if isempty(dfn)
    % Calculate the specific volume anomaly or the in-situ density anomaly
    d = m - m0;
    d_ = d(RG.ocean).';
end
x_ = x(RG.ocean).'; % Row vector

% Integrate the specific volume (in-situ density) as a function of pressure
% (depth) in a non-Boussinesq (Boussinesq) ocean.
if OPTS.REEB
    
    if isempty(dfn)
        % --- Fit specific volume as a multivalued function of pressure,
        
        if ~OPTS.GEOSTRF
            % Requested to make the stream function ill-defined. Overwrite
            % cycle information so that all branches will be fit
            % independently, with no cycle constraints. Note, this does
            % nothing if dfn was provided. But if OPTS.GEOSTRF was false
            % when calling topobaric_surface(), then dfn produced by
            % topobaric will similarly yield an ill-defined stream
            % function.
            RG.cb_nodes = {};
            RG.cb_arcs = {};
        end
        
        dfn = branches_fit(x_, d_, RG.arc_from, RG.arc_to, RG.arc_segment, RG.node_fn, RG.cb_arcs, RG.cb_nodes, false);
    end
    
    % Integrate the multivalued function
    if nargout < 2
        intfn          = int_graph_fun(dfn, RG.arc_from, RG.arc_to, RG.node_fn, RG.graph, RG.bfs_parent_node, RG.bfs_topo_order, RG.bfs_missing_arc);
    else
        [intfn,intfn2] = int_graph_fun(dfn, RG.arc_from, RG.arc_to, RG.node_fn, RG.graph, RG.bfs_parent_node, RG.bfs_topo_order, RG.bfs_missing_arc);
    end
    
    % Evaluate the integral everywhere, using the local branch
    intval_ = nan(RG.n_casts, 1);
    arc_segment = RG.arc_segment; % dereference for faster execution
    for e = 1:RG.nArcs
        seg = arc_segment{e};
        intval_(seg) = pvaln(intfn(:,e), x_(seg)); % With no graph simplification, we can use pvaln instead of pvallin
    end
    
else
    
    if isempty(dfn)
        
        % --- Fit specific volume as a single-valued function of pressure
        breaks = OPTS.SPLINE_BREAKS;
        breaks = min(breaks, max(x_));
        breaks = max(breaks, min(x_));
        breaks = unique(breaks);
        dfn = splinefit(x_, d_, breaks, OPTS.SPLINE_ORDER);
    end
    
    % Integrate the single-valued function
    intval_ = ppint(dfn, x_(1), x_);
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
    fac = -grav / rho_c; % (-) because x is z we have been z>0 people!
else
    fac = db2Pa;
end

% Finish computing the geostrophic stream function
gsf(RG.ocean) = gsf(RG.ocean) + intval_ * fac;

if nargout < 2
    return
end

if ~OPTS.REEB
    % No possibility for different gsf's if using single-valued functions
    gsfdiff = zeros(nx,ny);
    gsfdiff(isnan(gsf)) = nan;
    return
end

% Calculate the difference between the above gsf and the alternative
% topobaric gsf.  The difference is so slight that a difference is returned
% rather than the alternative topobaric gsf itself, to avoid machine
% rounding error.
gsfdiff_ = zeros(RG.n_casts, 1);
bfs_missing_arc = RG.bfs_missing_arc; % dereference for faster execution
for i = 1:length(bfs_missing_arc)
    e = bfs_missing_arc(i);
    seg = arc_segment{e};
    intval2_seg = pvaln(intfn2(:,i), x_(seg)); % With no graph simplification, we can use pvaln instead of pvallin
    gsfdiff_(seg) = (intval2_seg(:) - intval_(seg)) * fac;
end
gsfdiff = zeros(nx,ny);
gsfdiff(isnan(gsf)) = nan;
gsfdiff(RG.ocean) = gsfdiff_;

end