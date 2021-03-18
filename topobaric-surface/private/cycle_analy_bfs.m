function [dup, G, bfs_dist, bfs_parent_node, bfs_topo_order, bfs_missing_arc, cb_arcs, cb_nodes] = cycle_analy_bfs(arc_from, arc_to, nNodes)
%CYCLE_ANALY_BFS  Analyse cycles in a graph by means of a breadth-first search.
%
%
% dup = cycle_analy_bfs(arc_from, arc_to, nNodes)
% finds all duplicate arcs in the graph, i.e. those arcs that are incident
% to a pair of nodes that are already adjacent by another arc. The number
% of nodes in the graph is nNodes, the number of arcs in the graph is
% length(arc_from), and there is an arc incident to node arc_from(a) and
% node arc_to(a), for each arc index a.
% 
% [dup, G] = cycle_analy_bfs(...)
% also returns a sparse matrix G of size nNodes by nNodes, representing the
% graph exclusive of its duplicate arcs. Specifically, G(m,n) = a when arc
% a is incident to nodes m and n.

% [dup, G, bfs_dist, bfs_parent_node] = cycle_analy_bfs(...)
% performs a breadth-first search on the graph (with 1 the root node). If
% node n is reachable from the root node, then bfs_dist(n) is the length of
% the shortest walk between n and the root node, and bfs_parent_node(n) is
% the penultimate node in this walk; otherwise, bfs_dist(n) is -1 and
% bfs_parent_node(n) is 0. of nodes between the root note and each node is
% given by bfs_dist. The number
% 
% [dup, G, bfs_dist, bfs_parent_node, bfs_topo_order] = cycle_analy_bfs(...)
% also provides a topological ordering of the nodes by sorting bfs_dist.
% That is, for all nodes m < n, the breadth-first search discovers node
% bfs_topo_order(m) before node bfs_topo_order(n).
% 
% [dup, G, bfs_dist, bfs_parent_node, bfs_topo_order, bfs_missing_arc] = cycle_analy_bfs(...)
% also returns the vector of arcs that were not traversed in the
% breadth-first search. Each of these arcs defines a cycle in the graph,
% for example the cycle consisting of this missing arc and the shortest
% walk in the graph exclusive of this missing arc that connects the two
% nodes incident to this missing arc.
%
% [dup, G, bfs_dist, bfs_parent_node, bfs_topo_order, bfs_missing_arc, cb_arcs, cb_nodes] = cycle_analy_bfs(...)
% also computes a cycle basis, the i'th element of which is a cycle whose
% walk is cb_nodes{i}(1), cb_arcs{i}(1), ... cb_nodes{i}(end),
% cb_arcs{i}(end), cb_nodes{i}(1). Also note cb_arcs{i}(1) ==
% bfs_missing_arc(i).
%
%
% --- Input:
% arc_from [A, 1]: the "lower" node that each arc is incident to
% arc_to   [A, 1]: the "higher" node that each arc is incident to
% nNodes [1, 1]: the number of nodes in the graph
%
%
% --- Output:
% dup [vector]: index for duplicate arcs
% G [N, N]: sparse matrix representing the graph without duplicate arcs
% bfs_dist [N, 1]: distance from the root node (1) to each node in a breadth-first search
% bfs_parent_node [N, 1]: the node from which this node was discovered in the breadth-first search
% bfs_topo_order [N, 1]: topological ordering of the nodes in the breadth-first search
% bfs_missing_arc [C, 1]: arcs that were not traversed during the breadth-first search
% cb_arcs  {C, 1}: cell array, each element an array of arcs  defining a cycle in the cycle basis
% cb_ndoes {C, 1}: cell array, each element an array of nodes defining a cycle in the cycle basis
%
% Note, above A is the number of arcs, N the number of nodes, and C the
% number of cycles in the cycle basis of the graph.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% --- Find duplicate arcs (between the same pair of nodes):
[~, iA, iC] = unique([arc_from, arc_to], 'rows', 'stable');
nondup = (iA(iC) - (1:length(iC))') == 0;
dup = find(~nondup);
if nargout == 1; return; end

% --- Build adjacency matrix (for non-duplicate arcs) holding the arc id
% G(i,j) = a where arc_from(a) == i and arc_to(a) == j.
nArcs = length(arc_from);
arc_id = 1:nArcs;
G = sparse(arc_from(nondup), arc_to(nondup), arc_id(nondup), nNodes, nNodes);
G = G + G.'; % Make symmetric.
if nargout == 2; return; end

% --- Breadth first search from the root (1) node
[bfs_dist, ~, bfs_parent_node] = bfs(G,1);
bfs_parent_node = bfs_parent_node(:); 
if nargout < 5; return; end

% --- Build a Topological Ordering of the nodes.
[~, bfs_topo_order] = sort(bfs_dist);
if nargout == 5; return; end


% --- Build list of arcs that define each cycle in the cycle basis
% For each (non-root) node, get the arc that led to it in the bfs.
e = full(G(sub2ind(size(G), 2:nNodes, bfs_parent_node(2:nNodes).')));

% All the other arcs define the cycles (ie. removing that arc from G removes that cycle)
bfs_missing_arc = true(nArcs, 1);
bfs_missing_arc(e) = false;
bfs_missing_arc = find(bfs_missing_arc);
if nargout == 6; return; end

% --- Build cycle basis
nCyclesInBasis = length(bfs_missing_arc);
cb_arcs = cell(nCyclesInBasis,1);
cb_nodes = cell(nCyclesInBasis,1);
assert(nCyclesInBasis == nArcs - nNodes + 1, 'Incorrect number of cycles in basis. Does graph have multiple components?');
if nCyclesInBasis == 0
    return
end

% Build the Minimum Spanning Tree
GMST = G;
for e = bfs_missing_arc.'
    if nondup(e)
        % Only remove the arc if it is a non-duplicate arc.
        % If e is a duplicate arc, it does not actually exist in G, though
        % its duplicate does, so do not remove it!
        i = arc_from(e);
        j = arc_to(e);
        GMST(i,j) = 0;
        GMST(j,i) = 0;
    end
end

% Pre-process for speedier bfs()
[rp, ci] = sparse_to_csr(GMST);
GMST = struct('rp', rp, 'ci', ci);

for c = 1 : nCyclesInBasis
    e = bfs_missing_arc(c);
    n1 = arc_from(e);
    n2 = arc_to(e);
    
    if nondup(e)
        
        % Find shortest path from n1 to n2, using bfs
        [dist,~,pred] = bfs(GMST, n1, n2);
        n1_to_n2 = path_from_pred(pred, dist, n2);
        
        % Collect the nodes and arcs on this cycle.
        cb_nodes{c} = circshift(n1_to_n2,+1);
        cb_arcs{c} = full(G(sub2ind(size(G), cb_nodes{c}, n1_to_n2)));

    else
        
        cb_nodes{c} = [n2; n1];
        cb_arcs{c} = [e; full(G(n1,n2))];
        
    end
    %assert(cb_arcs{c}(1) == bfs_missing_arc(c), 'Something is wrong');
    %assert(all(cb_arcs{c} > 0), 'Something is wrong');
end
