function [qu, qts, ncc, G, L] = bfs_conncomp(G, A, r, qu) %#codegen
%BFS_CONNCOMP  Find connected components using Breadth First Search
%
%
% [qu, qts, ncc] = bfs_conncomp(G, A, r)
% performs a breadth first search (BFS) through nodes of a graph, starting
% from the root node r.  The graph is best thought of as laid out on a
% grid, with nodes at the true locations of the multi-dimensional array G,
% and numbered by their linear index in G.   Node m is adjacent to node n =
% A(m,j) provided (n > 0 && G(n)), for each j = 1 : D, where D = size(A,2)
% is the maximal degree in the graph. The outputs are: the search queue for
% the BFS, qu; the tail indices of qu for the BFS, qts; and the number of
% connected components found, ncc = 1. Specifically, qu(qts(1) : qts(2)-1)
% are the linear indices to true locations of G searched by the BFS.  Note,
% qu(1) = r, and qts(1) = 1.
%
% [qu, qts, ncc] = bfs_conncomp(G, A)
% as above, but scans the nodes from first to last, starting a separate BFS
% on any node previously undiscovered by a BFS.  Each BFS discovers an
% entire connected component. The outputs are: the combined search queue
% for all the BFS, qu: the tail indices of qu for each BFS, qts; and the
% number of connected components, ncc.  Specifically,  qu(qts(i) :
% qts(i+1)-1) are linear indices to the i'th connected component of G,
% ordered by their discovery time in the i'th BFS, for i = 1 : ncc.  Note,
% qts(1) = 1.
%
% [...] = bfs_conncomp(..., qu)
% works on qu using memory in-place, where the fourth argument qu is a
% vector of length N and class double.  Pass r as [] to use this with the
% second form above.
%
% [qu, qts, ncc, G] = bfs_conncomp(...)
% also returns a modified G that is true only where the BFS has walked.
% When r is not provided, the output G should match the input G.
%
% [qu, qts, ncc, G, L] = bfs_conncomp(...)
% also returns a label array L, the same size as G. Specifically, L(n) = i
% for all n in the i'th connected component, for i = 1 : ncc; otherwise
% L(n) = 0.
%
%
% --- Input:
% G [array with N elements]: true where there are nodes, false elsewhere
% A [D,N]: adjacency, where D is the most neighbours possible
% r [1,1]: perform one BFS from this root node (optional)
% qu [N,1]: vector to work in-place (optional)
%
%
% --- Output:
% qu [N,1]: the nodes visited by the BFS's in order from 1 to qts(end)
% qts [ncc+1,1]: the queue tail indices for each connected component
% ncc [1,1]: the number of connected components of G
% G [array with N elements]: true where there are nodes that have been discovered
% L [array with N elements]: label array
%
%
% See also GRID_ADJACENCY

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


N = numel(G);

if nargin < 4 || isempty(qu)
  qu = zeros(N, 1); % pre-allocate queue storing linear indices to nodes
end


qt = 0; % Queue Tail
qh = 0; % Queue Head

if nargin < 3 || isempty(r)
  % No root provided.  Sweep over all nodes, starting a new BFS rooted at
  % any previously undiscovered node.
  ncc = 0;             % Number of Connected Components
  len = ceil(sqrt(N)); % A # hopefully greater than the final ncc
  qts = zeros(len,1);  % Queue TailS
  
  for r = 1 : N
    if G(r)
      
      ncc = ncc + 1; % Found another connected component
      if ncc > len   % Too many components.  Double length of qts. Crude.
        qts = [qts; zeros(len,1)]; %#ok<AGROW>
        len = len * 2;
      end
      qts(ncc) = qt + 1; % Record search queue tail from this connected component
      
      [G, qu, qh, qt] = bfs(r, G, A, qu, qh, qt);
      
    end
  end
  
  qts(ncc+1) = qt + 1;  % Record search queue tail for final connected component
  qts = qts(1 : ncc+1); % Trim
  
else
  % Root location specified.  Do one BFS from that root node.
  
  [G, qu, ~, qt] = bfs(r, G, A, qu, qh, qt);
  
  ncc = 1;
  qts = [1; qt+1];
end

if nargout >= 4
  % Create the Label
  G(:) = false;
  G(qu(1 : qts(end)-1)) = true;
  
  if nargout == 5
    % Create the Label
    L = zeros(size(G));
    for i = 1 : ncc
      L(qu(qts(i) : qts(i+1)-1)) = i;
    end
  end
  
end


end

function [G, qu, qh, qt] = bfs(r, G, A, qu, qh, qt)
D = size(A,1); % maximal degree
N = numel(G);

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

end


