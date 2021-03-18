function [qu, qt] = bfs_conncomp1(G, A, r, qu) %#codegen
%BFS_CONNCOMP1  Find one connected component using Breadth First Search
%
%
% [qu, qt] = bfs_conncomp(G, A, r)
% performs a breadth first search (BFS) through nodes of a graph, starting
% from the root node r.  The graph is best thought of as laid out on a
% grid, with nodes at the true locations of the multi-dimensional array G,
% and numbered by their linear index in G.   Node m is adjacent to node n =
% A(m,j) provided (n <= numel(G) && G(n)), for each j = 1 : D, where D =
% size(A,2) is the maximal degree in the graph. The outputs are: the search
% queue for the BFS, qu; the tail index of qu for the BFS, qt.
% Specifically, qu(1 : qt) are the linear indices to true locations of G
% searched by the BFS.
%
% [...] = bfs_conncomp(..., qu)
% works on qu using memory in-place, where the fourth argument qu is a
% vector of length N and class double.
%
%
% --- Input:
% G [array with N elements]: true where there are nodes, false elsewhere
% A [D,N]: adjacency, where D is the most neighbours possible
% r [1,1]: perform one BFS from this root node
% qu [N,1]: vector to work in-place (optional)
%
%
% --- Output:
% qu [N,1]: the nodes visited by the BFS's in order from 1 to qt
% qt [1,1]: the queue tail index
%
%
% See also GRID_ADJACENCY

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


N = numel(G);

if nargin < 4 || isempty(qu)
  qu = zeros(N, 1); % pre-allocate queue storing linear indices to nodes
end


qt = 0; % Queue Tail
qh = 0; % Queue Head

D = size(A,1); % maximal degree

% Initialize BFS from root node
qt = qt + 1; % Add r to queue
qu(qt) = r;
G(r) = false; % mark r as discovered

while qt > qh
  qh = qh + 1; % advance head of the queue
  m = qu(qh);  % me node; pop from head of queue
  for d = 1 : D
    n = A(d,m); % neighbour node
    if n <= N && G(n) % First condition checks n is not a neighbour across a non-periodic boundary
      % n is good, and undiscovered
      qt = qt + 1;  % Add n to queue
      qu(qt) = n;
      G(n) = false; % mark n as discovered
    end
  end
end


