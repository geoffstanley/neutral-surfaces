function [F, F2] = int_graph_fun(f, arc_from, arc_to, node_fn, graph, bfs_parent_node, bfs_topo_order, bfs_missing_arc)
%INT_GRAPH_FUN  Integrate a multivalued function with branches on a graph.
%
%
% F = int_graph_fun(f, arc_from, arc_to, node_fn)
% integrates the multivalued function f having a single-valued branch for
% each arc of the graph. The graph has an arc a incident to two nodes,
% arc_from(a) and arc_to(a); node n is associated with a value node_fn(n).
% The single-valued branch of f on arc a is defined by f(:,a), and may be
% evaluated at x in the domain of f(:,a) by pvaln(f(:,a), x). The same
% applies to F. Given a node n, F is such that pvaln(F(:,a),x) is the same
% for all arcs a incident to node n. That is, the branches of F, associated
% with arcs incident to a common node, meet continuously at the values
% associated with the common node.
%
% ... = int_graph_fun(f, arc_from, arc_to, node_fn, graph, bfs_parent_node, bfs_topo_order, bfs_missing_arc)
% provides further information about the cycles in the graph, as output
% from cycle_analy_bfs.
%
% [F, F2] = int_graph_fun(...)
% also gives F2, an alternative integral of f where arcs defining cycles in
% the cycle basis of the graph (namely the arcs bfs_missing_arc), are
% integrated in the opposite direction to how F integrated thme. Whereas F
% meets continuously at the "lower" node associated with these arcs, F2
% meets continuously at the "upper" node associated with these arcs. In
% fact, only the branches of F at these arcs bfs_missing_arc are provided
% in F2.
%
%
% --- Input:
% f [P+2, A]: the multivalued function. Each column defines a single-valued
% branch as a polynomial of order P that may be evaluated by pvaln
% arc_from [A, 1]: the "lower" node that each arc is incident to
% arc_to   [A, 1]: the "higher" node that each arc is incident to
% node_fn  [N, 1]: the value associated with each node
% graph [N, N]: sparse matrix representing the graph without duplicate arcs
% bfs_parent_node [N, 1]: the node from which this node was discovered in the breadth-first search
% bfs_topo_order [N, 1]: topological ordering of the nodes in the breadth-first search
% bfs_missing_arc [C, 1]: arcs that were not traversed during the breadth-first search
%
%
% --- Output
% F  [P+3,A]: the integral of f. 
% F2 [P+3,C]: an alternative integral of f for the branches associated with
%             the arcs bfs_missing_arc.
%
% Note, above A is the number of arcs, N the number of nodes, and C the
% number of cycles in the cycle basis of the graph.

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


narginchk(4,8);

nNodes = length(node_fn);

if nargin < 8  % Locally compute graph cycle information
    [~, graph, ~, bfs_parent_node, bfs_topo_order, bfs_missing_arc] = cycle_analy_bfs(arc_from, arc_to, nNodes);
end

% Integrate all functions from the lower node towards the higher node
% Evaluating F at the lower node gives 0; the constant of integration must
% be added individually (per arc).
F = pintn(f);

% Integrate f from node 1 to node 2 in the ordering
n2 = bfs_topo_order(2);
n1 = bfs_parent_node(n2); % == topoorder(1) == 1
e = full(graph(n1, n2));
if n1 == arc_to(e) % Reverse the fact that F had integrated from lower node
    F(end,e) = -pvaln(F(:,e), node_fn(n1));
end


if nNodes > 2
    % Integrate f from node 2 to node 3 in the ordering
    n2 = bfs_topo_order(3);
    n1 = bfs_parent_node(n2);
    if n1 == 1
        % Handle rare case where root node (1) has two neighbours:
        % bfs_topo_order(3) == n2 --> o   o <-- n0 == bfs_topo_order(2)
        %                              \ /
        %                      root --> o <-- n1 == 1
        n0 = bfs_topo_order(2);
    else
        % Handle the usual case, where the root has one neighbour:
        %         o   o <-- n2
        %          \ /
        %           o <-- n1
        %           |
        % root -->  o <-- n0
        n0 = bfs_parent_node(n1);
    end
    e = full(graph(n1, n2));
    
    if n1 == arc_to(e) % Reverse the fact that F had integrated from lower node
        F(end,e) = -pvaln(F(:,e), node_fn(n1));
    end

    % Connect this branch with its parent branch
    a01 = full(graph(n0, n1));
    F(end,e) = F(end,e)  +  pvaln(F(:,a01), node_fn(n1));

    
    for k = 4:nNodes
        % Integrate f to node k in the topological ordering
        n2 = bfs_topo_order(k);
        n1 = bfs_parent_node(n2);
        n0 = bfs_parent_node(n1);
        e = full(graph(n1, n2));

        if n1 == arc_to(e) % Reverse the fact that F had integrated from lower node
            F(end,e) = -pvaln(F(:,e), node_fn(n1));
        end

        % Connect this branch with its parent branch
        a01 = full(graph(n0, n1));
        F(end,e) = F(end,e)  +  pvaln(F(:,a01), node_fn(n1));

    end

end

% Now integrate the arcs that define elements of the cycle basis.
% If f is well-defined, the integration direction shouldn't matter,
% so just choose to integrate from lower node towards higher node
for e = bfs_missing_arc(:).'
    n1 = arc_from(e);
    n0 = bfs_parent_node(n1);
    a01 = full(graph(n0, n1));
    F(end,e) = pvaln(F(:,a01), node_fn(n1));
end


if nargout >= 2
    % Requested alternative integration where cycle arc is traversed in the
    % other direction (higher node towards lower node)
    F2 = F(:,bfs_missing_arc); % Pre-allocate space
    for i = 1 : length(bfs_missing_arc)
        e = bfs_missing_arc(i);
        n1 = arc_to(e);
        n0 = bfs_parent_node(n1);
        a01 = full(graph(n0, n1));
        F2(end,i) = F2(end,i) + pvaln(F(:,a01), node_fn(n1)) - pvaln(F2(:,i), node_fn(n1));
    end
end
