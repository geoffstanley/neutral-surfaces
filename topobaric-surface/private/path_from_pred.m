function path = path_from_pred(pred,d,t)
%PATH_FROM_PRED  Convert a predecessor array into a path to a vertex.
%
%
% path = path_from_pred(pred,d,t) 
% returns the list of nodes on a walk from a source node to a target node
% t. The predecessor array pred is a row vector such that
%   pred(i) = 0 if node i is the source node,
%   pred(i) = 0 if node i has no predecessor associated with it, or
%   pred(i) = j where node j preceeds node i on the shortest path from the
%               source node to i.
% The discovery time d is a vector such that
%   d(i) = -1 where node i is unreachable from the source node; otherwise,
%   d(i) = the number of nodes between the source node and node i.
% The target node, t, is an integer from 1 to n, where n is the length of
% the predecessor array.
% path(i) is the i'th node in the walk from the source node to node t.
% Typically, pred and d are output from a breadth-first search such as bfs
% in the GAIMC toolbox.
% 
% In the case that the node is the source or unreachable from the source,
% then path is a scalar, the source node. These input cases must be handled
% with care, for the output of path_from_pred in these cases are
% indistinguishable.
%
%
% --- Input:
% pred [1, N]: predecessor array
% d [1, N]: discovery time
% t [1, 1]: target node
% 
%
% --- Output:
% path [*, 1]: list of nodes in the walk from source node to target node.

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

if d(t) < 0
    % target unreachable from source vertex
    path = find(d == 0);
else
    % target reachable from source vertex.
    % Build list of nodes visited from source to target.
    path = zeros(d(t) + 1, 1);
    path(end) = t;
    for k = length(path)-1:-1:1
        path(k) = pred(path(k+1));
    end
end