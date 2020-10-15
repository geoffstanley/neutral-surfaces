function [v, f, Qmap, nV] = latlon2vertsfaces(Q, WRAP, X, Y, method, select)
%LATLON2VERTSFACES  Make simplical mesh from rectilinear data.
%
%
% [v, f] = latlon2vertsfaces(Q)
% creates triangular mesh, specified by vertices v and triangles f, from
% the rectilinear data Q: each rectangle is made into two triangles, except
% where precisely one of the four Q data points on the rectangle is NaN, in
% which case only one triangle is produced.  When two valid (non-NaN)
% points are adjacent and have four NaN's on either side of them (an
% "alleyway"), an additional vertex is added between them and an additional
% face connecting these three vertices.
%
% [v, f, Qmap, nV] = latlon2vertsfaces(Q)
% also returns Qmap that maps an (i,j) coordinate for accessing Q onto a
% vertex coordinate for accessing v. Specifically, v(I,3) == Q(i,j) where I
% = Qmap(i,j). Also returns nV, the number of vertices from the original Q
% mesh.
%
% ... = latlon2vertsfaces(Q, WRAP)
% as above if WRAP is [false false]. This specifies that Q data is periodic
% in its first dimension iff WRAP(1) is true, and periodic in its second
% dimension iff WRAP(2) is true.
%
% ... = latlon2vertsfaces(Q, WRAP, X, Y)
% as above if X is 1:size(Q,1) and Y is 1:size(Q,2). This uses vectors X and
% Y to specify the physical coordinates of the Q data, and hence also the
% vertices v.
%
% ... = latlon2vertsfaces(Q, WRAP, X, Y, method)
% as above if method(1) == 'd' (for "diagonal"). If method(1) == 'c' (for
% "cross"), each rectangle of Q data is split into four triangles by
% introducing a new vertex at the centre of the rectangle, with Q value
% obtained by averaging the four surrounding Q points. If precisely one of
% the four Q data points on the rectangle is NaN, only one triangle is
% produced, and no extra vertices.
%
% ... = latlon2vertsfaces(Q, WRAP, X, Y, method, select)
% as above if select == false. If select == true, only the largest
% connected component of the triangular mesh is kept. If select == [x, y],
% only the connected component of the triangular mesh containing point
% (x,y) is kept.
%
%
% --- Input:
% Q [ni, nj]: rectilinear data
% WRAP [1, 2]: specification for which dimensions of Q are periodic [logical]
% X [ni, 1]: spatial positions for the first dimension of Q
% Y [ni, 1]: spatial positions for the second dimension of Q
% method [string]:
%   method(1) == 'd' (for diagonal) splits a rectangular face of Q into two
%       triangles.
%   method(1) == 'c' (for cross) splits a rectangular face of Q into four
%       triangles, with the mean of Q on the four corners taken as the
%       value in the centre.
% select [1, 1] or [1, 2]:
%   true, to select the largest connected component;
%   false, to keep all components, even if they are disconnected;
%   [x y], to select the connected component containing vertex nearest to (x,y).
%
%
% --- Output:
% v [N, 3]: vertices of the triangular mesh. v(:,1:2) gives the X and Y
%   coordinates of the vertices. v(:,3) gives the Q values.
% f [M, 3]: three integers in each row index v to form a triangle.
% Qmap [ni, nj]: mapping from Q indices to v index.
% nV [1, 1]: the number of vertices in the mesh originating from the
%   original Q grid.
%
%
% Note: to view the mesh, uncomment and comment out lines as specified in
% the code, run the function, then run the following:
% figure; patch('Faces',faces, 'Vertices', [v(:,1:2), zeros(size(v,1),1)], ...
%     'edgecolor', 'k', 'facecolor', 'interp', 'FaceVertexCData', v(:,3));

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


% --- Setup and check arguments
[ni,nj] = size(Q);
nij = ni * nj;

if nargin < 2 || isempty(WRAP)
    WRAP = [false false];
end
if nargin < 3 || isempty(X)
    X = 1:ni;
end
if nargin < 4 || isempty(Y)
    Y = 1:nj;
end
if nargin < 5 || isempty(method)
    method = 'd';
end
if nargin < 6 || isempty(select)
    select = false;
end

method = lower(method(1));
assert(method == 'c' || method == 'd', 'method should be either ''diagonal'' or ''cross''');


% --- Build Vertex list and a 2D map of the vertex indices
goodQ = isfinite(Q);

% Qmap(i,j) gives the index of vertex (i,j) on the Q grid. Same size as Q.
% i.e. Qvec(Qmap(i,j)) == Q(i,j), where Qvec = Q(goodQ).
Qmap = reshape(cumsum(goodQ(:)), ni, nj);
nV = Qmap(end);
Qmap(~goodQ) = 0; % Within this function, these will never be accessed

% --- Count number of vertices and number of faces (for pre-allocation of v and f)
% In the following diagrams:
%   + is a valid grid point (ocean),
%   % is invalid (ground), and
%   -- and | represent connections between grid points.

% Use bitwise-and to determine, for each (i,j), which of the neighbouring 4
% [(i,j), (i+1,j), (i,j+1), (i+1,j+1)] are valid points. Except, in
% MATLAB, it's faster to just do arithmetic. If (i,j) is the upper-left
% point in the following diagram, then b(i,j) is 1 + the sum of a subset of
% {1, 2, 4, 8} depending on which points are valid.
%   1 -- 4
%   |    |
%   2 -- 8
% The additional +1 to b simply makes indexing easier later, since
% MATLAB starts indexing at 1 (not 0)
goodQ8 = uint8(goodQ);
b = 1 + goodQ8 + 2*ip1(goodQ8) + 4*jp1(goodQ8) + 8*ip1jp1(goodQ8);

% Group indices that share a value of b.
bgroup = accumarray(b(:), 1:numel(b), [16, 1], @(x) {x});

% Alleyways: two valid points with invalid points on either side.
% Find valid points (+) that match the following arrangement:
%   % --(+)-- %
%   |    |    |
%   % -- + -- %
% The left square is 1+4+8=13 and the right square is 1+1+2=4
indsUD = find(b == 4 & jm1(b) == 13); % up-down
nUD = length(indsUD);

% Find grid points (+) that match the following arrangement:
%   % -- %
%   |    |
%  (+)-- +
%   |    |
%   % -- %
% The top square is 1+2+8 = 11 and the bottom square is 1+1+4=6
indsLR = find(b == 6 & im1(b) == 11); % left-right
nLR = length(indsLR);

if method == 'c'
    % 4 way average the data values
    M = Q + jp1(Q);
    M = (M + ip1(M)) / 4;
    goodM = isfinite(M);
    
    Mmap = nV + reshape(cumsum(goodM(:)), ni, nj);
    nVerts = Mmap(end);
    Mmap(~goodM) = 0;
    
else % method == 'd'
    
    nVerts = nV;
    
end

len = cellfun('length', bgroup);
%nFaces = len(8) + len(12) + len(14) + len(15) + len(16) * (2+2*(method == 'c'));
nFaces = len(8) + len(12) + len(14) + len(15) + len(16) * (2+2*(method == 'c')) + nUD + nLR;
nVerts = nVerts + nUD + nLR;

% --- Build Vertices and Faces
v = zeros(nVerts, 3, 'like', Q); % pre-allocate
%[xq, yq] = ndgrid(X,Y);                                                   % <-- Uncomment to give vertices spatial coordinates.
%v(1:nV,1) = xq(goodQ);                                                    % <-- Uncomment to give vertices spatial coordinates.
%v(1:nV,2) = yq(goodQ);                                                    % <-- Uncomment to give vertices spatial coordinates.
v(1:nV,3) = Q(goodQ);

if method == 'c'
    %XM = [(X(1:end-1) + X(2:end))/2, X(end) + (X(end)-X(end-1))/2];       % <-- Uncomment to give vertices spatial coordinates.
    %YM = [(Y(1:end-1) + Y(2:end))/2, Y(end) + (Y(end)-Y(end-1))/2];       % <-- Uncomment to give vertices spatial coordinates.
    %[xm, ym] = ndgrid(XM,YM);                                             % <-- Uncomment to give vertices spatial coordinates.
    %v(nV+1:end,1) = xm(goodM);                                            % <-- Uncomment to give vertices spatial coordinates.
    %v(nV+1:end,2) = ym(goodM);                                            % <-- Uncomment to give vertices spatial coordinates.
    v(nV+1:end,3) = M(goodM);
end

iF = 0;  % iterate for f
iV = nV; % iterate for v
f = zeros(nFaces,3); % pre-allocate

Qmap_ip1 = ip1(Qmap);
Qmap_jp1 = jp1(Qmap);
Qmap_ip1jp1 = jp1(Qmap_ip1);

% 16: all four corners valid:
I = bgroup{16}; nF = len(16);
if method == 'd'
    % Diagonal simplical mesh.
    f(1+iF : nF+iF,:) = [Qmap(I), Qmap_ip1(I), Qmap_jp1(I)];
    iF = iF + nF;
    f(1+iF : nF+iF,:) = [Qmap_ip1(I), Qmap_ip1jp1(I), Qmap_jp1(I)];
    iF = iF + nF;
else % method == 'c'
    % Cross simplical mesh.
    f(1+iF : nF+iF,:) = [Qmap(I), Qmap_ip1(I), Mmap(I)];
    iF = iF + nF;
    f(1+iF : nF+iF,:) = [Qmap_ip1(I), Qmap_ip1jp1(I), Mmap(I)];
    iF = iF + nF;
    f(1+iF : nF+iF,:) = [Qmap_ip1jp1(I), Qmap_jp1(I), Mmap(I)];
    iF = iF + nF;
    f(1+iF : nF+iF,:) = [Qmap_jp1(I), Qmap(I), Mmap(I)];
    iF = iF + nF;
end

% 15: (i,j) is NaN:
I = bgroup{15}; nF = len(15);
f(1+iF : nF+iF,:) = [Qmap_ip1(I), Qmap_ip1jp1(I), Qmap_jp1(I)];
iF = iF + nF;

% 14: (ip1,j) is NaN:
I = bgroup{14}; nF = len(14);
f(1+iF : nF+iF,:) = [Qmap(I), Qmap_jp1(I), Qmap_ip1jp1(I)];
iF = iF + nF;

% 12: (i,jp1) is NaN:
I = bgroup{12}; nF = len(12);
f(1+iF : nF+iF,:) = [Qmap(I), Qmap_ip1(I), Qmap_ip1jp1(I)];
iF = iF + nF;

% 8: (ip1,jp1) is NaN:
I = bgroup{8}; nF = len(8);
f(1+iF : nF+iF,:) = [Qmap(I), Qmap_jp1(I), Qmap_ip1(I)];
iF = iF + nF;


% Add a vertex between two vertices taking the average value of the two,
% and add a face using these two vertices and the new vertex.
I = indsUD;
J = I + 1 - ni * (mod(I,ni)==0); % Index to (i+1,j) from (+), with i modulo ni; note i will never equal ni if WRAP(1) is false
%v(nV+1 : nV+nUD, :) = 0.5 * [xq(I) + xq(J), yq(I) + yq(J), Q(I) + Q(J)];  % <-- Uncomment to give vertices spatial coordinates.
v(iV+1 : iV+nUD, :) = [zeros(nUD, 2), 0.5 * (Q(I) + Q(J))];                   % <-- Comment   to give vertices spatial coordinates.
f(iF+1 : iF+nUD, :) = [Qmap(I), Qmap(J), (iV+1 : iV+nUD).'];
iV = iV + nUD;
iF = iF + nUD;

I = indsLR;
J = I + ni - nij * (I > nij-ni); % Index to (i,j+1) from (+), with j modulo nj; note j will never equal nj if WRAP(2) is false
%v(nV+1 : nV+nLR, :) = 0.5 * [xq(I) + xq(J), yq(I) + yq(J), Q(I) + Q(J)];  % <-- Uncomment to give vertices spatial coordinates.
v(iV+1 : iV+nLR, :) = [zeros(nLR, 2), 0.5 * (Q(I) + Q(J))];                   % <-- Comment   to give vertices spatial coordinates.
f(iF+1 : iF+nLR, :) = [Qmap(I), Qmap(J), (iV+1 : iV+nLR).'];
%iV = iV + nLR;
%iF = iF + nLR;


if ~select
    return
end


% --- Select one component
G = sparse(f(:), reshape(f(:,[2 3 1]), [], 1), ones(numel(f),1), nVerts, nVerts);
G = G + G'; % Make symmetric, ie. an undirected graph.
% G has a few 1's, at the perimeter of simplical mesh. Most edges
% (between two vertices v_i and v_j, say) are part of two faces, and
% G(v_i,v_j) == 2. But at the perimeter, the edge is only part of one
% face, where G(v_i,v_j) == 1.

% Find disjoint regions.
[C,sizeC] = scomponents(G); % from the gaimc toolbox.
nCpts = length(sizeC);
C = C(:);
if nCpts == 1
    return
end

if isvector(select) && length(select) == 2 && isnumeric(select)
    % (x,y) given:  Select the component containing the specified point
    x = select(1);
    y = select(2);
    i = interp1(X, 1:ni, x, 'nearest');
    j = interp1(Y, 1:nj, y, 'nearest');
    c = C(Qmap(i,j)); % c is the label for that component, in C.
else
    % Select the largest component:
    [~,c] = max(sizeC); % c is the label for that component, in C.
end

% Remove all vertices not in the chosen connected region
% Remove all faces connected to a vertex not in the chosen connected region
keep_v = (C == c);
keep_f = keep_v(f(:,1));
v = v(keep_v,:);
f = f(keep_f,:);

% Update f to index the new, reduced set of vertices
vertmap = cumsum(keep_v);
f = vertmap(f);

if nargout < 3
    return
end

% Finally, update Qmap to map from the reduced set of vertices back to
% physical space.

% indsQ is the inverse to Qmap.  That is,
% Qmap(indsQ(v)) = v  for all I, and
% indsQ(Qmap(I)) = I  for all I such that Qmap(I) is non-NaN.
indsQ = find(Qmap > 0);

% Recreate Qmap's values and nV:
goodQ = false(size(Q));
goodQ(indsQ(keep_v(1:nV))) = true;
Qmap = reshape(cumsum(goodQ(:)), ni, nj);
nV = Qmap(end);
Qmap(~goodQ) = 0;


    function out = im1(in)
        out = circshift(in, [+1 0]);
        if ~WRAP(1)
            out(1,:) = nan;
        end
    end

    function out = ip1(in)
        out = circshift(in, [-1 0]);
        if ~WRAP(1)
            out(end,:) = nan;
        end
    end

    function out = jm1(in)
        out = circshift(in, [0 +1]);
        if ~WRAP(2)
            out(:,1) = nan;
        end
    end

    function out = jp1(in)
        out = circshift(in, [0 -1]);
        if ~WRAP(2)
            out(:,end) = nan;
        end
    end

    function out = ip1jp1(in)
        out = circshift(in, [-1 -1]);
        if ~WRAP(1)
            out(end,:) = nan;
        end
        if ~WRAP(2)
            out(:,end) = nan;
        end
    end

end