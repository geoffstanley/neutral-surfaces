function [v, f, Qmap, nQ] = latlon2vertsfaces(Q, WRAP, X, Y, method, select)
%LATLON2VERTSFACES  Make simplical mesh from rectilinear data.
%
%
% [v, f] = latlon2vertsfaces(Q)
% creates triangular mesh, specified by vertices v and triangles f, from
% the rectilinear data Q: each rectangle is made into two triangles, except
% where precisely one of the four Q data points on the rectangle is NaN, in
% which case only one triangle is produced.
%
% [v, f, Qmap, nQ] = latlon2vertsfaces(Q)
% also returns Qmap that maps an (i,j) coordinate for accessing Q onto a
% vertex coordinate for accessing v. Specifically, v(I,3) == Q(i,j) where 
% I = Qmap(i,j). Also returns nQ == size(v,1).
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
% produced, and no extra vertices. Also, nQ specifies that v(1:nQ,:) are
% the vertices corresponding to the original Q points; the extra vertices
% are v(nQ+1,:).
% 
% ... = latlon2vertsfaces(Q, WRAP, X, Y, method, select)
% as above if select == false. If select == true, only the largest
% connected component of the triangular mesh is kept. If select == [x, y],
% only the connected component of the triangular mesh containing point
% (x,y) is kept.
% 
%
% --- Input:
% Q [nx, ny]: rectilinear data
% WRAP [1, 2]: specification for which dimensions of Q are periodic [logical]
% X [nx, 1]: spatial positions for the first dimension of Q
% Y [nx, 1]: spatial positions for the second dimension of Q
% method [string]:
%   method(1) == 'd' (for diagonal) splits a rectangular face of Q into two
%       triangles.
%   method(1) == 'c' (for cross) splits a rectangular face of Q into four
%       triangles, with the mean of Q on the four corners taken as the
%       value in the centre.
%   Where Q is NaN on 1 corner of a rectangle, both methods both produce a
%       single triangle.
%   Where Q is NaN on 2 or more corners, both methods produce no triangle.
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
% Qmap [nx, ny]: mapping from Q indices to v index. 
% nQ [1, 1]: the number of vertices in the mesh originating from the
%   original Q grid.
%
%
% Note: to view the mesh, run the following:
% figure; patch('Faces',faces, 'Vertices', [v(:,1:2), zeros(size(v,1),1)], ...
%     'edgecolor', 'k', 'facecolor', 'interp', 'FaceVertexCData', v(:,3));
%
%
% --- Requirements:
% bfs, scomponents - https://www.mathworks.com/matlabcentral/fileexchange/24134

% --- Copyright:
% Copyright 2019 Geoff Stanley
%
% This file is part of Topobaric Surface.
% 
% Topobaric Surface is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
% 
% Topobaric Surface is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with Topobaric Surface.  If not, see
% <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au 
% Email     : geoffstanley@gmail.com
% Version   : 1.0
%
% Modified by : --
% Date        : --
% Changes     : --

% --- Setup and check arguments
[ni,nj] = size(Q);
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


% --- Build Vertex list and a 2D map of the vertex indices
goodQ = isfinite(Q);

% Qmap(i,j) gives the vertex number of vertex (i,j) on the Q grid. Same
% size as Q. Similarly for Mmap. ie. Qvec(Qmap(i,j)) == Q(i,j), where Qvec
% == Q(goodQ).
Qmap = reshape(cumsum(goodQ(:)), ni, nj);
nQ = Qmap(end);
Qmap(~goodQ) = nan; % These ought never to be accessed anyway.


% Do a 4 way average
if method == 'c'
    XM = [(X(1:end-1) + X(2:end))/2, X(end) + (X(end)-X(end-1))/2];
    YM = [(Y(1:end-1) + Y(2:end))/2, Y(end) + (Y(end)-Y(end-1))/2];
    M = Q + jp1(Q,WRAP);
    M = (M + ip1(M,WRAP)) / 4;
    goodM = isfinite(M);
    
    Mmap = nQ + reshape(cumsum(goodM(:)), ni, nj);
    nVerts = Mmap(end);
    Mmap(~goodM) = 0;
    
    nFaces = 4*ni*nj; % Upper bound
    
elseif method == 'd'
    
    nVerts = nQ;
    nFaces = 2*ni*nj; % Upper bound
    
end

v = nan(nVerts, 3, 'like', Q);
[icm, jcm] = ndgrid(X,Y);
v(1:nQ,1) = icm(goodQ);
v(1:nQ,2) = jcm(goodQ);
v(1:nQ,3) = Q(goodQ);

if method == 'c'
    [icm, jcm] = ndgrid(XM,YM);
    v(nQ+1:end,1) = icm(goodM);
    v(nQ+1:end,2) = jcm(goodM);
    v(nQ+1:end,3) = M(goodM);
end
clear icm jcm M


% --- Build Faces
goodQ8 = uint8(goodQ);
iF = 0; % iterate for Face
f = nan(nFaces,3);

Qmap_ip1 = ip1(Qmap,WRAP);
Qmap_jp1 = jp1(Qmap,WRAP);
Qmap_ip1jp1 = jp1(Qmap_ip1,WRAP);

% Use bitwise-and to determine, for each (i,j), which of the neighbouring 4
% [(i,j), (i+1,j), (i,j+1), (i+1,j+1)] are good data. Except, in MATLAB,
% it's faster to just do arithmetic. 
b = 1 + goodQ8 + 2*ip1(goodQ8,WRAP) + 4*jp1(goodQ8,WRAP) + 8*ip1jp1(goodQ8,WRAP);

% Group indices that share a value of b.
bgroup = accumarray(b(:), 1:numel(b), [16, 1], @(x) {x});


% 16: all four corners valid:
inds = bgroup{16}; nF = length(inds);
if method == 'd'
    % Diagonal simplical mesh.
    f((1+iF):(nF+iF),:) = [Qmap(inds), Qmap_ip1(inds), Qmap_jp1(inds)];
    iF = iF + nF;
    f((1+iF):(nF+iF),:) = [Qmap_ip1(inds), Qmap_ip1jp1(inds), Qmap_jp1(inds)];
    iF = iF + nF;
elseif method == 'c'
    % Cross simplical mesh.
    f((1+iF):(nF+iF),:) = [Qmap(inds), Qmap_ip1(inds), Mmap(inds)];
    iF = iF + nF;
    f((1+iF):(nF+iF),:) = [Qmap_ip1(inds), Qmap_ip1jp1(inds), Mmap(inds)];
    iF = iF + nF;
    f((1+iF):(nF+iF),:) = [Qmap_ip1jp1(inds), Qmap_jp1(inds), Mmap(inds)];
    iF = iF + nF;
    f((1+iF):(nF+iF),:) = [Qmap_jp1(inds), Qmap(inds), Mmap(inds)];
    iF = iF + nF;
end

% 15: (i,j) is NaN:
inds = bgroup{15}; nF = length(inds);
f((1+iF):(nF+iF),:) = [Qmap_ip1(inds), Qmap_ip1jp1(inds), Qmap_jp1(inds)];
iF = iF + nF;

% 14: (ip1,j) is NaN:
inds = bgroup{14}; nF = length(inds);
f((1+iF):(nF+iF),:) = [Qmap(inds), Qmap_jp1(inds), Qmap_ip1jp1(inds)];
iF = iF + nF;

% 12: (i,jp1) is NaN:
inds = bgroup{12}; nF = length(inds);
f((1+iF):(nF+iF),:) = [Qmap(inds), Qmap_ip1(inds), Qmap_ip1jp1(inds)];
iF = iF + nF;

% 8: (ip1,jp1) is NaN:
inds = bgroup{8}; nF = length(inds);
f((1+iF):(nF+iF),:) = [Qmap(inds), Qmap_jp1(inds), Qmap_ip1(inds)];
iF = iF + nF;

clear goodQ8
clear Qmap_ip1 Qmap_jp1 Qmap_ip1jp1

% Trim the faces matrix:
f(iF+1:end,:) = [];

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

% indsQ is the inverse to Qmap. [i,j] = sub2ind([ni nj], indsQ(v)), then Qmap(i,j) = v.
indsQ = find(Qmap > 0);

% Recreate Qmap's values and nQ:
goodQ = false(size(Q));
if method == 'd'
    goodQ(indsQ(keep_v)) = true;
else
    goodQ(indsQ(keep_v(1:nQ))) = true;
end
Qmap = reshape(cumsum(goodQ(:)), ni, nj);
nQ = Qmap(end);
Qmap(~goodQ) = 0;


end


function out = ip1(in,WRAP)
out = circshift(in, [-1 0]);
if ~WRAP(1)
    out(end,:) = nan;
end
end

function out = jp1(in,WRAP)
out = circshift(in, [0 -1]);
if ~WRAP(2)
    out(:,end) = nan;
end
end

function out = ip1jp1(in,WRAP)
out = circshift(in, [-1 -1]);
if ~WRAP(1)
    out(end,:) = nan;
end
if ~WRAP(2)
    out(:,end) = nan;
end
end