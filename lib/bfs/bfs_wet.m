function [qu, qt, s, t, x, freshly_wet] = bfs_wet(SppX, TppX, X, s, t, x, X_TOL, A, BotK, r, qu) %#codegen
%BFS_WET  Test neutral tangent plane connections from the edges of a
%         surface, using breadth first search
%
%
% [qu, qt, s, t, x] = bfs_wet(S, T, X, BotK, s, t, x, X_TOL, r, qu)
% performs a breadth first search (BFS), as in bfs_conncomp, on the valid
% parts of a surface, except steps are also made from a point on the
% (current) perimeter of the surface to a point in an adjacent water column
% if these points are neutrally connected.  The ocean has practical /
% Absolute salinity, potential / Conservative temperature, and pressure or
% depth x given by (S, T, X).  For the surface, these are (s, t, x).
% Neutral connections are given by the neutral tangent plane discretization
% where two points have the same potential density referenced to their
% average pressure.  This connection is numerically tested with an error
% tolerance, in pressure or depth, of X_TOL. The number of valid data
% points in each water column is BotK, which should be BotK = sum(isfinite(S),1).
% Inputs A, r, and qu are as in bfs_conncomp.  The outputs are: the search
% queue for the BFS, qu; the tail index of qu for the BFS, qt; the S, T,
% and X values on the updated wet surface, (s, t, x); the number of water
% columns added to the surface, freshly_wet.  Note, the points where the
% updated surface exists are given by linear indices qu(1:qt).  Note, this
% function changes the connected components of the surface, and hence the
% tail indices of the BFS search queue are not given as output here; that
% is, only qt is provided, the total number of valid points discovered or
% wet, rather than qts as provided by bfs_conncomp.  To obtain indices for
% the different connected components, run bfs_conncomp on the output
% surface, i.e. ``wet = false(size(x)); wet(qu(1:qt)) = true; [qu, qts, ncc]
% = bfs_conncomp(wet, A);``.  
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
% r [1,1]: perform one BFS from this root node (optional)
% qu [N,1]: vector to work in-place (optional)
%
%
% --- Output:
% qu [N,1]: the nodes visited by the BFS's in order from 1 to qt
% qt [1,1]: index for the tail element of qu
% s [ni,nj]: updated practical/Absolute salinity on the surface
% t [ni,nj]: updated potential/Conservative temperature on the surface
% x [ni,nj]: updated pressure [dbar] or depth [m] of the surface
% freshly_wet [1,1]: number of casts freshly wet
%
%
% See also BFS_CONNCOMP, GRID_ADJACENCY
%
%
% Note:  bfs_wet does not return qts or ncc like bfs_conncomp does, because
% wetting can change the connectivity structure of the surface, unaware of
% this algorithm.  If the connected components are required, simply use
% bfs_conncomp after using bfs_wet.

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
% Version   : 2.1.0
%
% Modified by : --
% Date        : --
% Changes     : --

[ni,nj] = size(x);
nij = ni * nj;
Xmat = double(~isvector(X)); % 0 if X is a vector, 1 otherwise
freshly_wet = 0; % count number of casts that wetting adds to the surface

if nargin < 11 || isempty(qu)
    qu = zeros(nij, 1); % pre-allocate queue storing linear indices to pixels
end

D = size(A,1); % maximal degree

G = isfinite(x); % good pixels

test = (BotK > 1) & ~G; % Try wetting only these locations: ocean and not currently in the surface

qt = 0; % Queue Tail
qh = 0; % Queue Head


if nargin < 10 || isempty(r)
    % No root provided.  Sweep over all nodes, starting a new BFS rooted at
    % any previously undiscovered node.
    for r = 1 : nij % r = root node
        if G(r)
            
            % --- BEGIN BFS
            % Push r into queue
            qt = qt + 1; % advance tail of the queue
            qu(qt) = r;
            G(r) = false; % set root node as discovered by the BFS
            
            while qt > qh
                qh = qh + 1; % advance head of the queue
                m = qu(qh); % me node; pop from head of queue
                for d = 1 : D
                    n = A(d,m); % neighbour node
                    if n  % check that n is not 0, i.e. a non-periodic boundary
                        if G(n)
                            % n is on the surface, and undiscovered
                            qt = qt + 1;  % Add n to queue
                            qu(qt) = n;
                            G(n) = false; % mark n as discovered
                        elseif test(n)
                            k = BotK(n);
                            % n is off the surface but in the ocean. Check for neutral connection
                            nX = (n-1) * Xmat + 1; % = n if Xmat, 1 if not Xmat
                            [x(n), s(n), t(n), success] = ntp_bottle_to_cast(SppX(:,:,n), TppX(:,:,n), X(:,nX), s(m), t(m), x(m), X_TOL, k);
                            if success
                                qt = qt + 1;     % Add n to queue
                                qu(qt) = n;
                                test(n) = false; % do not try wetting here again
                                freshly_wet = freshly_wet + 1;
                            end
                        end
                    end
                end
            end
            % --- END BFS
            
        end
    end
else
    % Root location specified.  Do one BFS from that root node.
    
    % --- BEGIN BFS
    % Push r into queue
    qt = qt + 1; % advance tail of the queue
    qu(qt) = r;
    G(r) = false; % set root node as discovered by the BFS
    
    while qt > qh
        qh = qh + 1; % advance head of the queue
        m = qu(qh); % me node; pop from head of queue
        for d = 1 : D
            n = A(d,m); % neighbour node
            if n  % check that n is not 0, i.e. a non-periodic boundary
                if G(n)
                    % n is on the surface, and undiscovered
                    qt = qt + 1;  % Add n to queue
                    qu(qt) = n;
                    G(n) = false; % mark n as discovered
                elseif test(n)
                    k = BotK(n);
                    % n is off the surface but in the ocean. Check for neutral connection
                    nX = (n-1) * Xmat + 1; % = n if Xmat, 1 if not Xmat
                    [x(n), s(n), t(n), success] = ntp_bottle_to_cast(SppX(:,:,n), TppX(:,:,n), X(:,nX), s(m), t(m), x(m), X_TOL, k);
                    if success
                        qt = qt + 1;     % Add n to queue
                        qu(qt) = n;
                        test(n) = false; % do not try wetting here again
                        freshly_wet = freshly_wet + 1;
                    end
                end
            end
        end
    end
    % --- END BFS
end

