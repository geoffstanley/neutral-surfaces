function [s, t, x, freshly_wet, qu, qt] = bfs_conncomp1_wet(SppX, TppX, X, s, t, x, X_TOL, A, BotK, r, qu) %#codegen
%BFS_CONNCOMP1_WET  Find one connected component using Breadth First Search,
%                   and test neutral tangent plane connections from the perimeter
%
%
% [s, t, x, freshly_wet] = bfs_conncomp1_wet(SppX, TppX, X, s, t, x, X_TOL, A, BotK, r, qu)
% works as bfs_conncomp1.m, with G in that function given by isfinite(x).
% That is, nodes in the surface are walked from the root node r.  Where an
% invalid node is reached, that is not part of the surface but nonetheless
% ocean (BotK > 1 here), a neutral tangent plane calculation is performed
% from the surface to the cast at this node.  If successful, this node is
% added to the surface, the salinity and temperature are interpolated using
% the piecewise polynomials SppX and TppX with knots at X, and the BFS
% continues.  NTP connections are accurate to a pressure or depth of X_TOL.
% The input qu is optional, simply to save memory by working in-place.  On
% output, qu(1:qt) are linear indices to the nodes on the surface (whether
% or not they are freshly wet), in the order that they were added to the
% queue.  Note qu(1) == r.
%
%
% --- Input:
% SppX [O,K-1,ni,nj]: coefficients for piecewise polynomial for practical
%                   / Absolute Salinity in terms of X
% TppX [O,K-1,ni,nj]: coefficients for piecewise polynomial for potential
%                   / Conservative Temperature in terms of X
% X [K,ni,nj]: knots for the pressure or depth of the casts
% s [ni,nj]: practical / Absolute salinity on the surface
% t [ni,nj]: potential / Conservative temperature on the surface
% x [ni,nj]: pressure or depth on the surface
% X_TOL [1,1]: tolerance in x for finding neutral connections
% A [D, ni*nj]: adjacency, where D is the most neighbours possible
% qu [N,1]: vector to work in-place (optional)
% BotK [ni,nj]: number of valid data points on each cast
% r [1, 1]: linear index to the reference cast
% qu [N,1]: vector to work in-place (optional)
%
%
% --- Output:
% s [ni,nj]: updated practical/Absolute salinity on the surface
% t [ni,nj]: updated potential/Conservative temperature on the surface
% x [ni,nj]: updated pressure [dbar] or depth [m] of the surface
% freshly_wet [1,1]: number of casts freshly wet
% qu [N,1]: the nodes visited by the BFS's in order from 1 to qt
% qt [1,1]: the queue tail index
%
%
% See also BFS_CONNCOMP1, GRID_ADJACENCY
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

dry = (BotK > 1) & ~G; % Try wetting only these locations: ocean and not currently in the surface

qt = 0; % Queue Tail
qh = 0; % Queue Head

qt = qt + 1; % Add r to queue
qu(qt) = r;
G(r) = false; % mark r as discovered

Xn = X(:,1);
while qt > qh
  qh = qh + 1; % advance head of the queue
  m = qu(qh); % me node; pop from head of queue
  
  for d = 1 : D
    n = A(d,m); % neighbour node
    if n <= N && G(n) % First condition checks n is not a neighbour across a non-periodic boundary
      % n is good, and undiscovered
      
      qt = qt + 1;  % Add n to queue
      qu(qt) = n;
      G(n) = false; % mark n as discovered
      
    elseif dry(n)
      % n is "dry".  Try wetting.
      
      if Xmat
        Xn = X(:,n);
      end
      x(n) = ntp_bottle_to_cast(SppX(:,:,n), TppX(:,:,n), Xn, BotK(n), s(m), t(m), x(m), X_TOL);
      if isfinite(x(n))
        % The NTP connection was successful
        
        [s(n), t(n)] = ppc_val2(Xn, SppX(:,:,n), TppX(:,:,n), x(n));
        
        qt = qt + 1;  % Add n to queue
        qu(qt) = n;
        G(n) = false; % mark n as discovered
        dry(n) = false;
        
        freshly_wet = freshly_wet + 1;  % augment counter of freshly wet casts
        
      end
    end
  end
  
end
end
