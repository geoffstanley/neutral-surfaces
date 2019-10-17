function A = grid_adjacency(sz, wrap, D)
% GRID_ADJACENCY  Linear indices to each neighbour of each grid point on a grid
%
%
% A = grid_adjacency(sz, wrap)
% builds a matrix A giving linear indices to each neighbour of each grid
% point in a grid of size specified by sz and periodicity specified by
% wrap.  The grid is periodic in the i'th dimension if and only if wrap(i)
% is true. The n'th row of A is a vector of linear indices to grid points
% that are adjacent to the grid point whose linear index is n.  The maximum
% number of neighbours of any grid point is D, which is the number of
% columns of A.  If grid point n has fewer neighbours than D, A(n,:) will
% contain some 0's.
%
% Currently, this function only handles 2D grids with 4-connectivity.
% Specifically, the grid is laid out as a matrix, and, in order, A(n,:)
% gives linear indices to the grid points left (<), above (^), below (v),
% and right (>) of grid point n.
%
%
% --- Input:
% sz: vector specifying dimensions of the grid, i.e. the output of size()
% wrap: vector specifying periodicity of the grid.
%
%
% --- Output:
% A [N,D]: the adjacency matrix

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
% Version   : 2.0.0
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

% Build linear index to each grid point
idx = reshape(1:N, ni, nj);

% Build adjacency matrix, currently as a multidimensional array of dimensions [sz, 4]
shift_im1 = @(F) circshift(F, [+1, 0]);
shift_jm1 = @(F) circshift(F, [0, +1]);
shift_ip1 = @(F) circshift(F, [-1, 0]);
shift_jp1 = @(F) circshift(F, [0, -1]);
if D == 4
    A = cat(dim+1, ...
        shift_jm1(idx), ... % Linear index to the grid point left  (<) of the local grid point
        shift_im1(idx), ... % Linear index to the grid point above (^)    the local grid point
        shift_ip1(idx), ... % Linear index to the grid point below (v)    the local grid point
        shift_jp1(idx));    % Linear index to the grid point right (>) of the local grid point
else % D == 8
    shift_1mm = @(F) circshift(F, [+1 +1]);
    shift_1mp = @(F) circshift(F, [+1 -1]);
    shift_1pm = @(F) circshift(F, [-1 +1]);
    shift_1pp = @(F) circshift(F, [-1 -1]);
    A = cat(dim+1, ...
        shift_jm1(idx), ... % Linear index to the grid point left  (<) of the local grid point
        shift_im1(idx), ... % Linear index to the grid point above (^)    the local grid point
        shift_ip1(idx), ... % Linear index to the grid point below (v)    the local grid point
        shift_jp1(idx), ... % Linear index to the grid point right (>) of the local grid point
        shift_1mm(idx), ... % Linear index to the grid point left  & above (<^) the local grid point
        shift_1mp(idx), ... % Linear index to the grid point right & above (^>) the local grid point
        shift_1pm(idx), ... % Linear index to the grid point left  & below (<v) the local grid point
        shift_1pp(idx));    % Linear index to the grid point right & below (v>) the local grid point
end

% Handle periodicity
if ~wrap(1)
    if D == 4
        A(1 ,:,2) = 0; % i-1 hits a wall when i = 1
        A(ni,:,3) = 0; % i+1 hits a wall when i = ni
    else
        A(1 ,:,[2,5,6]) = 0; % i-1 hits a wall when i = 1
        A(ni,:,[3,7,8]) = 0; % i+1 hits a wall when i = ni
    end
end
if ~wrap(2)
    if D == 4
        A(:,1 ,1) = 0; % j-1 hits a wall when j = 1
        A(:,nj,4) = 0; % j+1 hits a wall when j = nj
    else
        A(:,1 ,[1,5,7]) = 0; % j-1 hits a wall when j = 1
        A(:,nj,[4,6,8]) = 0; % j+1 hits a wall when j = nj
    end
end

% Reshape A to a matrix of dimensions [N, D]
A = reshape(A, N, D);
