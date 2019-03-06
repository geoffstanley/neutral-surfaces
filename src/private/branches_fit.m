function yfnx = branches_fit(x, y, arc_from, arc_to, arc_segment, node_fn, cb_arcs, cb_nodes, GEOSTRF)
%BRANCHES_FIT  Fit y as a multi-valued linear function of x subject to graph constraints.
%
%
% yfnx = branches_fit(x, y, arc_from, arc_to, arc_segment, node_fn, cb_arcs, cb_nodes)
% fits a multivalued function yfnx, each branch i of which is an affine
% linear function that best (in a least-squares sense) fits the data y(seg)
% to x(seg), where seg = arc_segment{i} indexes all data points associated
% with the arc i in the graph. The graph has an arc i incident to two
% nodes, arc_from(i) and arc_to(i); node n is associated with a value
% node_fn(n). A cycle basis for the graph is detailed by cb_arcs and
% cb_nodes. The branches of yfnx are fit subject to a set of constraints,
% one for each cycle in the cycle basis, such that the integral of yfnx (as
% can be computed from int_graph_fun) going around each cycle is zero.
%
% yfnx = branches_fit(x, y, arc_from, arc_to, arc_segment, node_fn, cb_arcs, cb_nodes, true)
% adds a second set of constraints such that the second integral of yfnx
% (as can be computed by applying int_graph_fun to the output of
% int_graph_fun called on yfnx) going around a cycle is zero.
%
%
% --- Input:
% x [1, D]: independent variable (must be row vector)
% y [1, D]: dependent variable (must be row vector)
% arc_from [A, 1]: the "lower" node that each arc is incident to
% arc_to   [A, 1]: the "higher" node that each arc is incident to
% arc_segment {A,1}: element a provides the indices, in the simplical
%   mesh, to the vertices in the associated region for arc a
% node_fn  [N, 1]: the value associated with each node
% cb_arcs {C,1}, cb_arcs{c} is arcs in cycle c of the cycle basis
% cb_nodes {C,1}, cb_nodes{c} is the nodes in cycle c of the cycle basis
% GEOSTRF [1, 1]: true if the second integral should also be constrained to
%   integrate to zero around each cycle.
%
% 
% --- Output:
% yfnx [4, A]: the multivalued function. Each column defines a
%   single-valued branch as an affine linear function that may be evaluated
%   by pvaln
%
% Note, above A is the number of arcs, N the number of nodes, and C the
% number of cycles in the cycle basis of the graph. D the number of data
% points in the simplical mesh that arc_segment indexes.
%
%
% --- Requirements:
% lsqlin (Optimization Toolbox)
%
%
% --- References
% Stanley, G. J.  Neutral surface topology. Ocean Modelling, submitted.
%
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


nArcs = length(arc_from);
if nargin < 9
    GEOSTRF = false;
end

% Each column of yfnx describes the linear function for each arc.
% First two rows give the domain of the data
% Third rows gives the slope
% Fourth rows  gives the value at the minimum of its domain
yfnx = zeros(4, nArcs);
yfnx(1,:) = node_fn(arc_from);
yfnx(2,:) = node_fn(arc_to);

nCycles = length(cb_arcs);
if nCycles > 0
    
    % Collect all arcs in all cycles, and make them unique.
    % Remember how to map between unqiue cycle arcs, and the complete set.
    cba = vertcat(cb_arcs{:});
    [arcs, ~, map] = unique(cba); % So, cba == arcs(map)
    
    % Collect all data we'll need to fit for these arcs
    verts = arc_segment(arcs);
    segsz = cellfun('length', verts);
    
    
    verts = vertcat(verts{:});
    x_ = x(verts) - repelem(yfnx(1,arcs), 1, segsz);
    y_ = y(verts);
    
    % Build block diagonal matrix, size M by N
    % The j'th block has length(segsz(arcs(j)) rows, and each row is of the form [x, 1]
    M = length(verts);
    N = 2 * length(arcs); % Number of parameters in yfnx. Each arc has a linear function for y, hence 2 free params per arc.
    MAT = sparse(repmat((1:M), 1, 2), ...
        repelem((0:2:N-2).', segsz) + [1,2], ...
        horzcat(x_, ones(1,M)), ...
        M, N);
    
    % --- Build constraint matrix, for all cycles:
    % Cycles that share an arc must be solved together.
    % Might be better to split this problem into sub-problems, each having
    % a collection of cycles that are disjoint from all the other cycles.
    % However, we'll just let lsqlin() do that internally.
    if GEOSTRF
        % Make modified topobaric surfaces, which have a well-defined exact geostrophic stream function.
        % There are two constraints per cycle in the cycle basis:
        % One to make the modified topobaric surface unique-depth, 
        % and one to make its stream function well-defined.
        nCON = 2 * nCycles;       % number of constraint equations
        nnzCON = 4 * length(cba); % number of non-zero entries in constraint matrix
        
        ii = zeros(nnzCON, 1);
        jj = zeros(nnzCON, 1);
        vv = zeros(nnzCON, 1);
        i = 0;
        I = 0;
        for c = 1 : nCycles
            nodes = cb_nodes{c};   % Actual node indices
            n = length(nodes);     % Number of nodes in the cycle
            e = map(i+1 : i+n).';  % Index of the arcs in the cycle, in the list of all arcs in all cycles!
            e = e*2 + [-1; 0];     % Prepare for matrix building.
            
            con = zeros(2,n);
            
            pn = node_fn(nodes).';  % pressure at each node
            dp = diff([pn, pn(1)]); % pressure difference between nodes
            sgn = sign(dp);
            
            % Build |dp| powers.
            dpp = cumprod(repmat(abs(dp), 3, 1), 1, 'reverse');
            dpp(1,:) = dpp(1,:) / 6; % == |dp|^3 / (3!)
            dpp(2,:) = dpp(2,:) / 2; % == |dp|^2 / (2!)
            %                 dpp(1,:) == |dp|
            
            % First do constraint for geostrophic stream function
            flip = (sgn == -1);
            con(:,flip) = con(:,flip) - dp(flip) .* dpp(2:3,flip);
            dpp = dpp .* sgn;
            con = con + dpp(1:2,:);
            con(:,1:end-1) = con(:,1:end-1) + (pn(1) - pn(2:end)) .* dpp(2:3,1:end-1);
            
            inds = I+1 : I+2*n;
            ii(inds) = 2*c-1;
            jj(inds) = e(:);
            vv(inds) = con(:);
            I = inds(end);
            
            % Now do constraint for specific volume / density.
            con = dpp(2:3,:) ;
            
            inds = I+1 : I+2*n;
            ii(inds) = 2*c;
            jj(inds) = e(:);
            vv(inds) = con(:);
            I = inds(end);
            
            i = i + n;
        end
    else
        % There is one constraint per cycle in the cycle basis, 
        % to make the topobaric surface unique-depth, 
        % or to make the topobaric geostrophic stream function well-defined
        nCON = nCycles;           % number of constraint equations
        nnzCON = 2 * length(cba); % number of non-zero entries in constraint matrix
        ii = zeros(nnzCON, 1);
        jj = zeros(nnzCON, 1);
        vv = zeros(nnzCON, 1);
        i = 0;
        I = 0;
        for c = 1 : nCycles
            nodes = cb_nodes{c};   % Actual node indices
            n = length(nodes);     % Number of nodes in the cycle
            e = map(i+1 : i+n).';  % Index of the arcs in the cycle, in the list of all arcs in all cycles!
            e = e*2 + [-1; 0];     % Prepare for matrix building.
            
            pn = node_fn(nodes).';  % pressure at each node
            dp = diff([pn, pn(1)]); % pressure difference between nodes
            
            con = repmat(dp, 2, 1);
            con(1,:) = con(1,:) .* abs(con(1,:)) / 2; % == dp .* abs(dp) / 2
            % con(2,:)                                  == dp;
            
            inds = I+1 : I+2*n;
            ii(inds) = c;
            jj(inds) = e(:);
            vv(inds) = con(:);
            I = inds(end);
            
            i = i + n;
        end
    end
    % Essentially, we build this (for each c as above, when ~DO_GEOSTRF): CON(c,a(:)) = dp(:);
    CON = sparse(ii, jj, vv, nCON, N, nnzCON);
    
    % --- Solve for model coefficients COEF that satisfy
    % CON * COEF = 0, and minimize |MAT * COEF - y_|
    opts.StepTolerance = 1e-16;
    opts.OptimalityTolerance = 1e-12;
    opts.ConstraintTolerance = 1e-12;
    opts.Display = 'off';
    COEF = lsqlin(MAT, y_, [], [], CON, zeros(nCON,1), [], [], [], opts);
    
    yfnx(3,arcs) = COEF(1:2:end);
    yfnx(4,arcs) = COEF(2:2:end);
    
    % Build list of arcs that must still be processed
    notdone = true(1,nArcs);
    notdone(arcs) = false;
    todo = find(notdone);
else
    todo = 1:nArcs;
end

% --- Fit all regions not on any cycle. These are uncoupled problems.
for e = todo
    seg = arc_segment{e};
    N = length(seg);
    x_ = x(seg).';   % Make column vector
    y_ = y(seg).'; % Make column vector
    
    MAT = ones(N, 2);
    MAT(:,1) = x_ - yfnx(1,e);
    yfnx(3:4,e) = MAT \ y_; % [Slope, offset at min x]
end


end