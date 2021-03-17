function ADJ = grid_adjacency(SZ, CONN, WRAP)
% GRID_ADJACENCY  Linear indices to each neighbour of each grid point on a grid
%
%
% ADJ = grid_adjacency(SZ, CONN, WRAP)
% builds a matrix ADJ giving linear indices to each of CONN neighbours of
% each grid point in a grid of size specified by SZ.  The grid is periodic
% in the i'th dimension if and only if WRAP(i) is true. The n'th column of
% ADJ is a vector of linear indices to grid points that are adjacent to the
% grid point whose linear index is n.  The maximum number of neighbours of
% any grid point is CONN, which is the number of rows of ADJ.  If grid
% point n has fewer neighbours than CONN, due to being near a non-periodic
% boundary, ADJ(:,n) will contain some flag values.  The flag value is
% ni*nj+1, which can be used to index some special value.
%
% The connectivity CONN specifies the number of neighbours to a central grid
% point, possibly including itself.  For i between 1 and CONN, the i'th
% neighbour is located relative to the central grid point according to the
% following diagram:
%  CONN==4        CONN==5   CONN==8   CONN==9
% +---------> j
% |   . 2 .       . 2 .     5 2 6     1 4 7
% |   1 . 4       1 3 5     1 . 4     2 5 8
% v   . 3 .       . 4 .     7 3 8     3 6 9
% i
% Here, i increases downward and j increases right.  For example, if CONN
% == 4, the 2'nd neighbour of the central grid point at (i,j) is located at
% (i-1,j) .
%
%
% --- Input:
% SZ: vector specifying dimensions of the grid, i.e. the output of size()
% CONN: the connectivity: 4, 5, 8, or 9.
% WRAP: vector specifying periodicity of the grid.
%
%
% --- Output:
% ADJ [CONN, SZ]: the adjacency matrix

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
% MERCHANTABILITY or FITNESS FOR ADJ PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


dim = length(SZ); % Number of dimensions in the grid
assert(dim == 2, 'grid_adjacency currently only works for 2D grids.');
assert(length(WRAP) == dim, 'WRAP must be a vector the same length as SZ.');

ni = SZ(1);
nj = SZ(2);

WALLVAL = ni * nj + 1;


% Build adjacency matrix and handle periodicity
if CONN == 4
  % . 2 .
  % 1 . 4
  % . 3 .
  DIR = [2 4 6 8];
  ADJ = helper(ni, nj, DIR);
  
  if ~WRAP(1)
    ADJ(2, 1 , :) = WALLVAL; % i-1 hits a wall when i = 1
    ADJ(3, ni, :) = WALLVAL; % i+1 hits a wall when i = ni
  end
  if ~WRAP(2)
    ADJ(1, :, 1 ) = WALLVAL; % j-1 hits a wall when j = 1
    ADJ(4, :, nj) = WALLVAL; % j+1 hits a wall when j = nj
  end
  
elseif CONN == 5
  % . 2 .
  % 1 3 5
  % . 4 .
  DIR = [2 4 5 6 8];
  ADJ = helper(ni, nj, DIR);
  
  if ~WRAP(1)
    ADJ(2, 1 , :) = WALLVAL; % i-1 hits a wall when i = 1
    ADJ(4, ni, :) = WALLVAL; % i+1 hits a wall when i = ni
  end
  if ~WRAP(2)
    ADJ(1, :, 1 ) = WALLVAL; % j-1 hits a wall when j = 1
    ADJ(5, :, nj) = WALLVAL; % j+1 hits a wall when j = nj
  end
  
elseif CONN == 8
  % 5 2 6
  % 1 . 4
  % 7 3 8
  DIR = [2 4 6 8 1 7 3 9];
  ADJ = helper(ni, nj, DIR);
  
  if ~WRAP(1)
    ADJ([2,5,6], 1 , :) = WALLVAL; % i-1 hits a wall when i = 1
    ADJ([3,7,8], ni, :) = WALLVAL; % i+1 hits a wall when i = ni
  end
  if ~WRAP(2)
    ADJ([1,5,7], :, 1 ) = WALLVAL; % j-1 hits a wall when j = 1
    ADJ([4,6,8], :, nj) = WALLVAL; % j+1 hits a wall when j = nj
  end
  
elseif CONN == 9
  % 1 4 7
  % 2 5 8
  % 3 6 9
  DIR = 1:9;
  ADJ = helper(ni, nj, DIR);
  
  if ~WRAP(1)
    ADJ([1 4 7], 1 , :) = WALLVAL; % i-1 hits a wall when i = 1
    ADJ([3 6 9], ni, :) = WALLVAL; % i+1 hits a wall when i = ni
  end
  if ~WRAP(2)
    ADJ([1 2 3], :, 1 ) = WALLVAL; % j-1 hits a wall when j = 1
    ADJ([7 8 9], :, nj) = WALLVAL; % j+1 hits a wall when j = nj
  end
  
else
  error('Unknown number of neighbours.  CONN must be one of 4, 5, 8, or 9.');
end


% Reshape ADJ to a matrix of dimensions [CONN, N]
ADJ = reshape(ADJ, CONN, []);

end

function ADJ = helper(ni, nj, DIR)

D = length(DIR);

% Prepare to circshift linear indices to some subset of its neighbours,
% generally ordered as follows:
% 1 4 7
% 2 5 8
% 3 6 9
spin = [ ...
  0 , 0, 0,   0, 0, 0,   0, 0, 0; ...
  -1, 0, 1,  -1, 0, 1,  -1, 0, 1; ...
  -1,-1,-1,   0, 0, 0,   1, 1, 1];

% Build linear index to each grid point, and repeat them D times
ADJ = repmat(reshape(1:ni*nj, [1, ni, nj]), [D 1 1]);  % D x ni x nj

% Shift these linear indices so they refer to their neighbours.
for d = 1 : D
  ADJ(d,:,:) = circshift(ADJ(d,:,:), -spin(:, DIR(d)));
end
end