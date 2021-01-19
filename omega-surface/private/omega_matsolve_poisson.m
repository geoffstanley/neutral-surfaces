function phi = omega_matsolve_poisson(s, t, x, DIST2on1_iJ, DIST1on2_Ij, A4, qu, qt, mr)
% OMEGA_MATSOLVE_POISSON  Build & solve the sparse matrix Poisson problem for omega surfaces
%
%
% phi = omega_matsolve_poisson(s, t, x, DIST2on1_iJ, DIST1on2_Ij, A4, qu, qts)
% builds and solves the sparse matrix problem for omega surfaces in Poisson
% form. In each connected component of the surface, phi has zero
% arithmetic, unweighted mean.
%
% phi = omega_matsolve_poisson(s, t, x, DIST2on1_iJ, DIST1on2_Ij, A4, qu, qts, Lcc, m_ref)
% as above, but makes phi be zero at m_ref.  Other connected components not
% containing m_ref are as above, with mean phi = 0 in these regions.
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
% A4 [4, ni*nj]: adjacency matrix  (see grid_adjacency.m)
% qu [ni*nj,1]: the nodes visited by the BFS's in order from 1 to qts(end) (see bfs_conncomp.m)
% qts [ncc+1,1]: the queue tail indices for each connected component (see bfs_conncomp.m)
% Lcc: [array with ni*nj elements]: label array, giving unique integer 
%   labels to each connected component of isfinite(x) (see bfs_conncomp.m)
% m_ref [1,1]:  linear index to a reference cast at which phi will be zero.
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

% A stripped down version of ntp_errors(s,t,x,1,1,true,false,true);
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (xm + xn) );
% [vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 1500 );  % DEV: testing omega software to find potential density surface
eps = vs .* (sm - sn) + vt .* (tm - tn);

bad = isnan(eps);
eps(bad) = 0;

if UNIFORM_GRID
  fac = ~double(bad); % 0 and 1
else
  fac = DIST2on1_iJ;
  fac = repmat(fac, ni / size(fac,1), nj / size(fac,2)); % automatic expansion to [ni,nj]
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

% A stripped down version of ntp_errors(s,t,x,1,1,true,false,true);
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (xm + xn) );
% [vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 1500 );  % DEV: testing omega software to find potential density surface

eps = vs .* (sm - sn) + vt .* (tm - tn);
bad = isnan(eps);
eps(bad) = 0;

if UNIFORM_GRID
  fac = ~double(bad); % 0 and 1
else
  fac = DIST1on2_Ij;
  fac = repmat(fac, ni / size(fac,1), nj / size(fac,2)); % automatic expansion to [ni,nj]
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
remap = zeros(ni,nj);  % remap indices from 2D into linear indices for the current connected component

% Collect and sort linear indices to all pixels in this region
m = sort(qu(1:qt));  % sorting here makes matrix better structured; overall speedup.

N = length(m);  % Number of water columns
if N <= 1  % There are definitely no equations to solve
  phi(m) = 0; % Leave this isolated pixel at current pressure
  return
end

% Label the water columns in this region alone by 1, 2, ... N
% No need to reset other elements of remap to 0.
remap(m) = 1:N;

% Pin surface at mr by changing the mr'th equation to be 1 * phi[mr] = 0.
D(mr) = 0;
L(:,mr) = 0;
L(IJ,mr) = 1;

% The above change renders the mr'th column on all rows irrelevant,
% since phi[mr] will be zero.  So, we may also set this column to 0,
% which we do here by setting the appropriate links in L to 0. This
% maintains symmetry of the matrix, and speeds up solution by a
% factor of about 2.
L(IM,A4(IP,mr)) = 0;
L(MJ,A4(PJ,mr)) = 0;
L(PJ,A4(MJ,mr)) = 0;
L(IP,A4(IM,mr)) = 0;

% Build the RHS of the matrix problem
rhs = D(m);

% Build indices for the rows of the sparse matrix
r = repmat(1:N, 5, 1);

% Build indices for the columns of the sparse matrix
% remap() changes global indices to local indices for this region, numbered 1, 2, ... N
c = [remap(A4(IM,m)); remap(A4(MJ,m)); remap(A4(PJ,m)); remap(A4(IP,m)); 1:N];

% Build the values of the sparse matrix
v = L(:,m);

% Pin surface at mr = m_ref or mr = m(1), by adding a 1 to (mr, mr) entry
% v(IJ, remap(mr)) = v(IJ, remap(mr)) + 1;

% Build the sparse matrix, with N rows and N columns
good = c > 0; % Ignore connections to dry pixels (though they should have zero J anyway, this is faster)
mat = sparse( r(good), c(good), v(good), N, N );


% Solve the matrix problem
sol = mat \ rhs;  % Cholesky direct solver

% Save solution
phi(m) = sol;

