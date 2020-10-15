function [s, t, x, freshly_wet, qu] = bfs_wet(SppX, TppX, X, s, t, x, X_TOL, A, BotK, qu) %#codegen
%BFS_WET  Test neutral tangent plane connections from the edges of a
%         surface, using breadth first search
%
%
% [s, t, x, freshly_wet] = bfs_wet(SppX, TppX, X, s, t, x, X_TOL, A, BotK)
% expands the valid parts of an approximately neutral surface --- of depth
% or pressure x on which the salinity is s and the temperature is t --- by
% testing for neutral tangent plane (NTP) connections from its perimeter to
% adjacent water columns.  The ocean water columns are given as piecewise
% polynomials SppX and TppX for the salinity and temperature, with knots at
% pressure or depth values X.  NTP connections are accurate to a pressure
% or depth of X_TOL. The adjacency information of the water columns is
% encoded in A, the output of grid_adjacency(). See also BFS_CONNCOMP for
% more information on A. The number of water columns added to the surface
% is given by freshly_wet.  The method begins with a breadth first search
% (BFS) on the valid pixels (where x,s,t are not NaN), discovering all
% invalid pixels (where x,s,t are NaN) that are adjacent to valid pixels.
% Then, a second BFS is initialized with all of these pixels. This second
% BFS tests NTP connections from all of its valid adjacent pixels to the
% current water column. If no such NTP connection is successful, the BFS
% continues.  Otherwise, the depths or pressures x of all the successful
% NTP connections is averaged and assigned to be the new depth or pressure
% of this water column in the surface, and then the neighbouring invalid
% water columns are added to the BFS queue.
%
% [s, t, x, freshly_wet, qu] = bfs_wet(SppX, TppX, X, s, t, x, X_TOL, A, BotK, qu)
% also provides a vector qu for working in-place, saving a small amount of
% memory allocation.  Note the contents of qu will be overwritten.
%
%
% --- Input:
% S [nk,ni,nj]: practical/Absolute salinity
% T [nk,ni,nj]: potential/Conservative temperature
% X [nk,ni,nj] or [nk,1]: pressure [dbar] or depth [m, positive]
% s [ni,nj]: practical/Absolute salinity on the surface
% t [ni,nj]: potential/Conservative temperature on the surface
% x [ni,nj]: pressure [dbar] or depth [m] of the surface
% X_TOL [1,1]: tolerance in x for finding neutral connections
% A [D, ni*nj]: adjacency, where D is the most neighbours possible
% BotK [ni,nj]: number of valid bottles on each cast
% qu [N,1]: vector to work in-place (optional)
%
%
% --- Output:
% s [ni,nj]: updated practical/Absolute salinity on the surface
% t [ni,nj]: updated potential/Conservative temperature on the surface
% x [ni,nj]: updated pressure [dbar] or depth [m] of the surface
% freshly_wet [1,1]: number of casts freshly wet
% qu [N,1]: vector serving no purpose, but for qu to be worked in-place
%
%
% See also BFS_CONNCOMP, GRID_ADJACENCY
%
%
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


[ni,nj] = size(x);
nij = ni * nj;
Xmat = double(~isvector(X)); % 0 if X is a vector, 1 otherwise
freshly_wet = 0; % count number of casts that wetting adds to the surface

if nargin < 10 || isempty(qu)
  qu = zeros(nij, 1); % pre-allocate queue storing linear indices to pixels
end

D = size(A,1); % maximal degree

G = isfinite(x); % Good nodes

test = (BotK > 1) & ~G; % Try wetting only these locations: ocean and not currently in the surface

qt = 0; % Queue Tail
qh = 0; % Queue Head

% Scan for Bad nodes that are not land and are neighbours of Good nodes
for m = 1 : nij % r = root node
  if ~G(m) && BotK(m) > 1
    for d = 1 : D
      n = A(d,m); % neighbour node
      if n && G(n) % if n is good
        qt = qt + 1; % advance the tail of the queue
        qu(qt) = m;  % push n node onto queue
      end
    end
  end
end

% --- BEGIN BFS
Xm = X(:,1);
while qt > qh
  qh = qh + 1; % advance head of the queue
  m = qu(qh); % me node; pop from head of queue
  
  km = BotK(m);
  if Xmat
    Xm = X(:,m);
  end
  
  xm = 0;
  nm = 0;
  
  % Check all possible neutral connections from m cast to its neighbours
  for d = 1 : D
    n = A(d,m); % neighbour node
    if n && G(n)
      x_ = ntp_bottle_to_cast(SppX(:,:,m), TppX(:,:,m), Xm, km, s(n), t(n), x(n), X_TOL);
      if isfinite(x_)
        xm = xm + x_;
        nm = nm + 1;
      end
    end
  end
  
  if nm > 0
    % Cast m is NTP-connected to at least one of its neighbours, so wet it.
    x(m) = xm / nm;
    [s(m), t(m)] = ppc_val2(Xm, SppX(:,:,m), TppX(:,:,m), x(m));
    
    % Add Bad neighbours of this cast to BFS
    for d = 1 : D
      n = A(d,m); % neighbour node
      if n && test(n)
        qt = qt + 1; % advance tail of the queue
        qu(qt) = n;  % push n onto queue
        test(n) = false; % ensure n is only added at most once to the queue
      end
    end
    
    G(m) = true;  % set m node as Good
        
    freshly_wet = freshly_wet + 1;  % augment counter of freshly wet casts
  end
end

end
