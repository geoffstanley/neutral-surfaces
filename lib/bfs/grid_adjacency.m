function A = grid_adjacency(sz, wrap, D)
% GRID_ADJACENCY  Linear indices to each neighbour of each grid point on a grid
%
%
% A = grid_adjacency(sz, wrap)
% builds a matrix A giving linear indices to each neighbour of each grid
% point in a grid of size specified by sz and periodicity specified by
% wrap.  The grid is periodic in the i'th dimension if and only if wrap(i)
% is true. The n'th column of A is a vector of linear indices to grid
% points that are adjacent to the grid point whose linear index is n.  The
% maximum number of neighbours of any grid point is D, which is the number
% of rows of A.  If grid point n has fewer neighbours than D, A(:,n) will
% contain some 0's.
%
% The connectivity D specifies the number of neighbours to a central grid
% point, possibly including itself.  For i between 1 and D, the i'th
% neighbour is located relative to the central grid point according to the
% following diagram:
%    D == 4       D == 5    D == 8    D == 9
% +---------> j
% |   . 2 .       . 2 .     5 2 6     1 4 7
% |   1 . 4       1 . 4     1 . 4     2 5 8
% v   . 3 .       . 3 .     7 3 8     3 6 9
% i
% Here, i increases downward and j increases right.  For example, if D ==
% 4, the 2'nd neighbour of the central grid point at (i,j) is located at
% (i-1,j) .
%
%
% --- Input:
% sz: vector specifying dimensions of the grid, i.e. the output of size()
% wrap: vector specifying periodicity of the grid.
% D: the connectivity: 4 (default), 5, 8, or 9.
%
%
% --- Output:
% A [D,N]: the adjacency matrix

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2019  Geoff Stanley
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
% Version   : 2.1.1
%
% Modified by : --
% Date        : --
% Changes     : --

dim = length(sz); % Number of dimensions in the grid
assert(dim == 2, 'grid_adjacency currently only works for 2D grids.');
assert(length(wrap) == dim, 'wrap must be a vector the same length as sz.');

if nargin < 3 || isempty(D)
  D = 4;
end

N = prod(sz);
ni = sz(1);
nj = sz(2);


% Build adjacency matrix and handle periodicity
if D == 4
  % . 2 .
  % 1 . 4
  % . 3 .
  DIR = [2 4 6 8];
  A = helper(ni, nj, N, D, DIR);
  
  if ~wrap(1)
    A(2, 1 , :) = 0; % i-1 hits a wall when i = 1
    A(3, ni, :) = 0; % i+1 hits a wall when i = ni
  end
  if ~wrap(2)
    A(1, :, 1 ) = 0; % j-1 hits a wall when j = 1
    A(4, :, nj) = 0; % j+1 hits a wall when j = nj
  end
  
elseif D == 5
  % . 2 .
  % 1 3 5
  % . 4 .
  DIR = [2 4 5 6 8];
  A = helper(ni, nj, N, D, DIR);
  
  if ~wrap(1)
    A(2, 1 , :) = 0; % i-1 hits a wall when i = 1
    A(3, ni, :) = 0; % i+1 hits a wall when i = ni
  end
  if ~wrap(2)
    A(1, :, 1 ) = 0; % j-1 hits a wall when j = 1
    A(4, :, nj) = 0; % j+1 hits a wall when j = nj
  end
  
elseif D == 8
  % 5 2 6
  % 1 . 4
  % 7 3 8
  DIR = [2 4 6 8 1 7 3 9];
  A = helper(ni, nj, N, D, DIR);
  
  if ~wrap(1)
    A([2,5,6], 1 , :) = 0; % i-1 hits a wall when i = 1
    A([3,7,8], ni, :) = 0; % i+1 hits a wall when i = ni
  end
  if ~wrap(2)
    A([1,5,7], :, 1 ) = 0; % j-1 hits a wall when j = 1
    A([4,6,8], :, nj) = 0; % j+1 hits a wall when j = nj
  end
  
elseif D == 9
  % 1 4 7
  % 2 5 8
  % 3 6 9
  DIR = 1:9;
  A = helper(ni, nj, N, D, DIR);
  
  if ~WRAP(1)
    A([1 4 7], 1 , :) = 0; % i-1 hits a wall when i = 1
    A([3 6 9], ni, :) = 0; % i+1 hits a wall when i = ni
  end
  if ~WRAP(2)
    A([1 2 3], :, 1 ) = 0; % j-1 hits a wall when j = 1
    A([7 8 9], :, nj) = 0; % j+1 hits a wall when j = nj
  end
  
else
  error('Unknown number of neighbours.  D should be one of 4, 5, 8, or 9.');
end


% Reshape A to a matrix of dimensions [D, N]
A = reshape(A, D, N);

end

function A = helper(ni, nj, N, D, DIR)

% Prepare to circshift linear indices to some subset of its neighbours,
% generally ordered as follows:
% 1 4 7
% 2 5 8
% 3 6 9
spin = zeros(3, 9);
spin(2:3,:) = -[ ... % negative prepares for circshift
  1,     0,    -1,     1,     0,    -1,     1,     0,    -1;
  1,     1,     1,     0,     0,     0,    -1,    -1,    -1];

% Build linear index to each grid point, and repeat them D times
A = repmat(reshape(1:N, [1, ni, nj]), [D 1 1]);  % D x ni x nj

% Shift these linear indices so they refer to their neighbours.
for d = 1 : D
  A(d,:,:) = circshift(A(d,:,:), -spin(:, DIR(d)));
end
end