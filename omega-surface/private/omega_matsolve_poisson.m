function phi = omega_matsolve_poisson(s, t, p, DIST2on1_iJ, DIST1on2_Ij, WRAP, A4, qu, N, mr)
% OMEGA_MATSOLVE_POISSON  Build & solve the sparse matrix Poisson problem for omega surfaces
%
%
% phi = omega_matsolve_poisson(s, t, p, DIST2on1_iJ, DIST1on2_Ij, WRAP, A4, qu, N, m_ref)
% builds and solves the sparse matrix problem for omega surfaces in Poisson
% form. 
%
%
% --- Input:
%  s [ni, nj]: practical / Absolute Salinity on the surface
%  t [ni, nj]: potential / Conservative Temperature on the surface
%  p [ni, nj]: pressure (or depth) on the surface
% DIST2on1_iJ [ni, nj]: The grid distance in the second dimension divided
%   by the grid distance in the first dimension, both centred at (I-1/2,J).
%   Equivalently, the square root of the area of a grid cell centred at
%   (I-1/2,J), divided by the distance from (I-1,J) to (I,J).   
% DIST1on2_Ij [ni, nj]: The grid distance in the first dimension divided by
%   the grid distance in the second dimension, both centred at (I-1/2,J).
%   Equivalently, the square root of the area of a grid cell centred at
%   (I,J-1/2), divided by the distance from (I,J-1) to (I,J).   
% WRAP [2 element array]: WRAP(i) is true iff the domain is periodic in the
%                         i'th lateral dimension.  
% A4 [4, ni*nj]: adjacency matrix  (see grid_adjacency.m)
% qu [ni*nj,1]: the nodes visited by the BFS's in order from 1 to N (see bfs_conncomp1.m)
% N [1,1]: the number of water columns in the region; 
%          also, the last valid index of qu, i.e. the queue tail (see bfs_conncomp1.m)
% mr [1,1]  : linear index to a reference cast at which phi will be zero.
%
%
% --- Output:
% phi [ni, nj]: density perturbation satisfying the discrete version of
%               div grad phi = - div epsilon,
%               where epsilon is the neutrality error (see ntp_errors.m).

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com

[ni,nj] = size(p);

phi = nan(ni, nj);

% If there is only one water column, there are no equations to solve,
% and the solution is simply phi = 0 at that water column, and nan elsewhere.
% Note, N > 0 should be guaranteed by omega_surf(), so N <= 1 should
% imply N == 1.  If N > 0 weren't guaranteed, this could throw an error.
if N <= 1
  phi(qu(1)) = 0; % Leave this isolated pixel at current pressure
  return
end

WALL = ni * nj + 1; % a flag values used by A4 to index neighbours that would go across a non-periodic boundary

% If both gridding variables are 1, then grid is uniform
UNIFORM_GRID = ...
  isscalar(DIST2on1_iJ) && DIST2on1_iJ == 1 && ...
  isscalar(DIST1on2_Ij) && DIST1on2_Ij == 1;

%% Begin building D = divergence of epsilon, and L = negative Laplacian (compact representation)

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
pm = p;

%% m = (i, j) and n = (i-1, j),  then also n = (i+1, j) by symmetry
sn = im1(sm);
tn = im1(tm);
pn = im1(pm);

if ~WRAP(1)
  sn(1,:) = nan;
end

% A stripped down version of ntp_errors(s,t,p,1,1,true,false,true);
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (pm + pn) );
% [vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 1500 );  % DEV: testing omega software to find potential density surface
eps = vs .* (sm - sn) + vt .* (tm - tn);

bad = isnan(eps);
eps(bad) = 0;

if UNIFORM_GRID
  fac = double(~bad); % 0 and 1
else
  fac = DIST2on1_iJ;
  fac(bad) = 0;
  eps = eps .* fac; % scale eps
end

D = -eps + ip1(eps);

L(IJ,:,:) = fac + ip1(fac);                     

L(MJ,:,:) = -fac;

L(PJ,:,:) = -ip1(fac); 

%% m = (i, j) and n = (i, j-1),  then also n = (i, j+1) by symmetry
sn = jm1(sm);
tn = jm1(tm);
pn = jm1(pm);

if ~WRAP(2)
  sn(:,1) = nan;
end

% A stripped down version of ntp_errors(s,t,p,1,1,true,false,true);
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (pm + pn) );
% [vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 1500 );  % DEV: testing omega software to find potential density surface

eps = vs .* (sm - sn) + vt .* (tm - tn);
bad = isnan(eps);
eps(bad) = 0;

if UNIFORM_GRID
  fac = double(~bad); % 0 and 1
else
  fac = DIST1on2_Ij;
  fac(bad) = 0;
  eps = eps .* fac; % scale eps
end

D = D - eps + jp1(eps);

L(IJ,:,:) = squeeze(L(IJ,:,:)) + fac + jp1(fac);       

L(IM,:,:) = -fac;

L(IP,:,:) = -jp1(fac); 

%% Build and solve sparse matrix problem

% Collect and sort linear indices to all pixels in this region
m = sort(qu(1:N));  % sorting here makes matrix better structured; overall speedup.

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

