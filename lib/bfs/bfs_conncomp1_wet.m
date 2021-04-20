function [s, t, p, freshly_wet, qu, qt] = bfs_conncomp1_wet(Sppc, Tppc, P, s, t, p, TOL_P, A, BotK, r, qu) %#codegen
%BFS_CONNCOMP1_WET  Find one connected component using Breadth First Search,
%                   and test neutral tangent plane connections from the perimeter
%
%
% [s, t, p, freshly_wet] = bfs_conncomp1_wet(Sppc, Tppc, P, s, t, p, TOL_P, A, BotK, r, qu)
% works as bfs_conncomp1.m, with G in that function given by isfinite(p).
% That is, nodes in the surface are walked from the root node r.  Where an
% invalid node is reached, that is not part of the surface but nonetheless
% ocean (BotK > 1 here), a neutral tangent plane calculation is performed
% from the surface to the cast at this node.  If successful, this node is
% added to the surface, the salinity and temperature are interpolated using
% the piecewise polynomials Sppc and Tppc with knots at P, and the BFS
% continues.  NTP connections are accurate to a pressure or depth of TOL_P.
% The input qu is optional, simply to save memory by working in-place.  On
% output, qu(1:qt) are linear indices to the nodes on the surface (whether
% or not they are freshly wet), in the order that they were added to the
% queue.  Note qu(1) == r.
%
%
% --- Input:
% Sppc [O,K-1,ni,nj]: coefficients for piecewise polynomial for practical
%                   / Absolute Salinity in terms of P
% Tppc [O,K-1,ni,nj]: coefficients for piecewise polynomial for potential
%                   / Conservative Temperature in terms of P
% P [K,ni,nj]: knots for the pressure or depth of the casts
% s [ni,nj]: practical / Absolute salinity on the surface
% t [ni,nj]: potential / Conservative temperature on the surface
% p [ni,nj]: pressure or depth on the surface
% TOL_P [1,1]: tolerance in p for finding neutral connections
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
% p [ni,nj]: updated pressure [dbar] or depth [m] of the surface
% freshly_wet [1,1]: number of casts freshly wet
% qu [N,1]: the nodes visited by the BFS's in order from 1 to qt
% qt [1,1]: the queue tail index
%
%
% See also BFS_CONNCOMP1, GRID_ADJACENCY
%
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


[ni,nj] = size(p);
nij = ni * nj;
Pmat = double(~isvector(P)); % 0 if P is a vector, 1 otherwise
freshly_wet = 0; % count number of casts that wetting adds to the surface

if nargin < 10 || isempty(qu)
  qu = zeros(nij, 1); % pre-allocate queue storing linear indices to pixels
end

D = size(A,1); % maximal degree

G = isfinite(p); % Good nodes

dry = (BotK > 1) & ~G; % Try wetting only these locations: ocean and not currently in the surface

qt = 0; % Queue Tail
qh = 0; % Queue Head

qt = qt + 1; % Add r to queue
qu(qt) = r;
G(r) = false; % mark r as discovered

Pn = P(:,1);
while qt > qh
  qh = qh + 1; % advance head of the queue
  m = qu(qh); % me node; pop from head of queue
  
  for d = 1 : D
    n = A(d,m); % neighbour node
    if n <= nij  % check n is not a neighbour across a non-periodic boundary
      if G(n) 
        % n is good, and undiscovered

        qt = qt + 1;  % Add n to queue
        qu(qt) = n;
        G(n) = false; % mark n as discovered

      elseif dry(n)
        % n is "dry".  Try wetting.

        if Pmat
          Pn = P(:,n);
        end
        p(n) = ntp_bottle_to_cast(Sppc(:,:,n), Tppc(:,:,n), Pn, BotK(n), s(m), t(m), p(m), TOL_P);
        if isfinite(p(n))
          % The NTP connection was successful

          [s(n), t(n)] = ppc_val2(Pn, Sppc(:,:,n), Tppc(:,:,n), p(n));

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
end
