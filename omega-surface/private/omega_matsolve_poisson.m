function phi = omega_matsolve_poisson(s, t, x, DIST2on1_iJ, DIST1on2_Ij, WRAP, A4, qu, qt, mr)
% OMEGA_MATSOLVE_POISSON  Build & solve the sparse matrix Poisson problem for omega surfaces
%
%
% phi = omega_matsolve_poisson(s, t, x, DIST2on1_iJ, DIST1on2_Ij, WRAP, A4, qu, qt, m_ref)
% builds and solves the sparse matrix problem for omega surfaces in Poisson
% form. 
%
%
% --- Input:
%  s [ni, nj]: practical / Absolute Salinity on the surface
%  t [ni, nj]: potential / Conservative Temperature on the surface
%  x [ni, nj]: pressure or depth on the surface
% DIST2on1_iJ [ni, nj]: the area of a grid cell centred at (I-1/2, J),
%   divided by the distance, squared, from (I-1,J) to (I,J).  Equivalently,
%   this is the grid distance in second dimension divided by grid distance
%   in first dimension, both centred at (I-1/2, J).
% DIST1on2_Ij [ni, nj]: the area of a grid cell centred at (I, J-1/2),
%   divided by the distance, squared, from (I,J-1) to (I,J).  Equivalently,
%   this is the grid distance in first dimension divided by grid distance
%   in second dimension, both centred at (I, J-1/2).
% WRAP [2 element array]: WRAP(i) is true iff the domain is periodic in the
%                         i'th lateral dimension.  
% A4 [4, ni*nj]: adjacency matrix  (see grid_adjacency.m)
% qu [ni*nj,1]: the nodes visited by the BFS's in order from 1 to qt (see bfs_conncomp1.m)
% qt [1,1]: the last valid index of qu (see bfs_conncomp1.m)
% m_ref [1,1]  : linear index to a reference cast at which phi will be zero.
%
%
% --- Output:
% phi [ni, nj]: density perturbation satisfying the discrete version of
%               div grad phi = - div epsilon,
%               where epsilon is the neutrality error (see ntp_errors.m).

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
WALL = ni * nj + 1; % a flag values used by A4 to index neighbours that would go across a non-periodic boundary

autoexp = @(x) repmat(x, ni / size(x,1), nj / size(x,2)); % automatic expansion to [ni,nj]

% If both gridding variables are 1, then grid is uniform
UNIFORM_GRID = ...
  isscalar(DIST2on1_iJ) && DIST2on1_iJ == 1 && ...
  isscalar(DIST1on2_Ij) && DIST1on2_Ij == 1;

%% Begin building D = divergence of epsilon, and L = Laplacian (compact representation)

% L refers to neighbours in this order (so does A4, except without the 5'th entry):
% . 2 .
% 1 5 4
% . 3 .
IM = 1; % (I  ,J-1)
MJ = 2; % (I-1,J  )
PJ = 3; % (I+1,J  )
IP = 4; % (I  ,J+1)
IJ = 5; % (I  ,J  )
L = zeros(5, ni, nj); % pre-alloc space


% Anonymous functions for shifting data.  (Faster than indexing with A4.)
im1 = @(D) circshift(D, [+1 0]);
ip1 = @(D) circshift(D, [-1 0]);
jm1 = @(D) circshift(D, [0 +1]);
jp1 = @(D) circshift(D, [0 -1]);

% Aliases
sm = s;
tm = t;
xm = x;

%% m = (i, j) and n = (i-1, j),  then also n = (i+1, j) by symmetry
sn = im1(sm);
tn = im1(tm);
xn = im1(xm);

if ~WRAP(1)
  sn(1,:) = nan;
end

% A stripped down version of ntp_errors(s,t,x,1,1,true,false,true);
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (xm + xn) );
% [vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 1500 );  % DEV: testing omega software to find potential density surface
eps = vs .* (sm - sn) + vt .* (tm - tn);

bad = isnan(eps);
eps(bad) = 0;

if UNIFORM_GRID
  fac = double(~bad); % 0 and 1
else
  fac = autoexp(DIST2on1_iJ);
  fac(bad) = 0;
  eps = eps .* fac; % scale eps
end

D = eps - ip1(eps);

L(IJ,:,:) = fac + ip1(fac);                     

L(MJ,:,:) = -fac;

L(PJ,:,:) = -ip1(fac); 

%% m = (i, j) and n = (i, j-1),  then also n = (i, j+1) by symmetry
sn = jm1(sm);
tn = jm1(tm);
xn = jm1(xm);

if ~WRAP(2)
  sn(:,1) = nan;
end

% A stripped down version of ntp_errors(s,t,x,1,1,true,false,true);
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (xm + xn) );
% [vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 1500 );  % DEV: testing omega software to find potential density surface

eps = vs .* (sm - sn) + vt .* (tm - tn);
bad = isnan(eps);
eps(bad) = 0;

if UNIFORM_GRID
  fac = double(~bad); % 0 and 1
else
  fac = autoexp(DIST1on2_Ij);
  fac(bad) = 0;
  eps = eps .* fac; % scale eps
end

D = D + eps - jp1(eps);

L(IJ,:,:) = squeeze(L(IJ,:,:)) + fac + jp1(fac);       

L(IM,:,:) = -fac;

L(IP,:,:) = -jp1(fac); 

%% Finish building L
% For any m where all neighbours are NaN, set L(IJ,m) to 1 so that this
% equation amounts to:  1 * phi(m) = 0. This keeps phi(m) = 0, rather than
% becoming NaN and infecting its neighbours.
L(IJ, L(IJ,:,:) == 0) = 1;

%% Build and solve sparse matrix problem
phi = nan(ni, nj);     % solution to matrix problem

% Collect and sort linear indices to all pixels in this region
m = sort(qu(1:qt));  % sorting here makes matrix better structured; overall speedup.

N = length(m);  % Number of water columns
if N <= 1  % There are definitely no equations to solve
  phi(m) = 0; % Leave this isolated pixel at current pressure
  return
end

% remap changes from linear indices for the entire 2D space into linear
% indices for the current connected component.  
% If the domain were doubly periodic, we would want remap to be a 2D array
% of size [ni,nj]. However, with a potentially non-periodic domain, we need
% one more value for A4 to index into, hence use a vector with ni*nj+1
% elements.
% All other water columns are left with remap values of 0,
% including the "water column" corresponding to non-periodic boundaries.
remap = zeros(1, WALL);  % pre-alloc space and leave land as 0.
remap(m) = 1:N;          % Label the water columns in this region alone by 1, 2, ... N.

% Pin surface at mr by changing the mr'th equation to be 1 * phi[mr] = 0.
D(mr) = 0;
L(:,mr) = 0;
L(IJ,mr) = 1;

% The above change renders the mr'th column on all rows irrelevant,
% since phi[mr] will be zero.  So, we may also set this column to 0,
% which we do here by setting the appropriate links in L to 0. This
% maintains symmetry of the matrix, and speeds up solution by a
% factor of about 2.
if A4(IP,mr) ~= WALL
  L(IM,A4(IP,mr)) = 0;
end
if A4(PJ,mr) ~= WALL
  L(MJ,A4(PJ,mr)) = 0;
end
if A4(MJ,mr) ~= WALL
  L(PJ,A4(MJ,mr)) = 0;
end
if A4(IM,mr) ~= WALL
  L(IP,A4(IM,mr)) = 0;
end

% Build the RHS of the matrix problem
rhs = D(m);

% Build indices for the rows of the sparse matrix
r = repmat(1:N, 5, 1);

% Build indices for the columns of the sparse matrix
% remap() changes global indices to local indices for this region, numbered 1, 2, ... N
c = [remap(A4(IM,m)); remap(A4(MJ,m)); remap(A4(PJ,m)); remap(A4(IP,m)); 1:N];

% Build the values of the sparse matrix
v = L(:,m);

% Build the sparse matrix, with N rows and N columns
good = c > 0; % Ignore connections to dry pixels (though they should have zero J anyway, this is faster)
mat = sparse( r(good), c(good), v(good), N, N );


% Solve the matrix problem
sol = mat \ rhs;  % Cholesky direct solver

% Save solution
phi(m) = sol;

