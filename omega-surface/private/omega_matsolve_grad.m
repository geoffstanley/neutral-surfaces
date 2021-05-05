function phi = omega_matsolve_grad(s, t, p, sqrtDIST2on1_iJ, sqrtDIST1on2_Ij, WRAP, A4, qu, N, mr)
% OMEGA_MATSOLVE_GRAD  Build & solve the sparse matrix gradient equations for omega surfaces
%
%
% phi = omega_matsolve_grad(s, t, p, DIST2on1_iJ, DIST1on2_Ij, WRAP, A4, qu, N, mr)
% builds and solves the sparse matrix problem for omega surfaces in
% gradient equation form, ensuring phi is zero at mr.
%
%
% --- Input:
%  s [ni, nj]: practical / Absolute Salinity on the surface
%  t [ni, nj]: potential / Conservative Temperature on the surface
%  p [ni, nj]: pressure (or depth) on the surface
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
% qu [ni*nj,1]: the nodes visited by the BFS's in order from 1 to N (see bfs_conncomp1.m)
% N [1,1]: the number of water columns in the region; 
%          also, the last valid index of qu, i.e. the queue tail (see bfs_conncomp1.m)
% mr [1,1]:  linear index to a reference cast at which phi will be zero.
%
%
% --- Output:
% phi [ni, nj]: density perturbation, providing the least-squares solution
%               to the discrete version of
%               grad phi = - epsilon,
%               where epsilon is the neutrality error (see ntp_errors.m).

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


%#ok<*UNRCH>


%% --- Hard-coded parameters

% Relative tolerance for LSQR. Since the matrix problem is overdetermined,
% the relative residual will, in general, exceed this tolerance bound. As
% such, it is challenging to relate LSQR's relative tolerance to physical
% tolerances on density or pressure.  From several numerical tests,
% the default relative tolerance of 1e-6 seems to work well.
TOL_LSQR_REL = 1e-6;

% ZERO_MEAN = true; % Add extra row specifying mean(phi) = 0
ZERO_MEAN = false; % PINNING.  Add extra row specifying phi = 0 for the first cast in each region

% The matrix has 1's and 0's everywhere, except this value in the final row
F = 1e-2; % chosen empirically from tests on 1x1deg OCCA data

% A4 refers to neighbours in this order:
% . 2 .
% 1 . 4
% . 3 .
IM = 1; % (I  ,J-1)
MJ = 2; % (I-1,J  )


%% Process Inputs
[ni,nj] = size(p);
WALL = ni * nj + 1; % a flag values used by A4 to index neighbours that would go across a non-periodic boundary

% If both gridding variables are 1, then grid is uniform
UNIFORM_GRID = ...
  isscalar(sqrtDIST2on1_iJ) && sqrtDIST2on1_iJ == 1 && ...
  isscalar(sqrtDIST1on2_Ij) && sqrtDIST1on2_Ij == 1;

phi = nan(ni, nj);     % solution to matrix problem

%% Check for trivial solution

% If there is only one water column, there are no equations to solve,
% and the solution is simply phi = 0 at that water column, and nan elsewhere.
% Note, N > 0 should be guaranteed by omega_surf(), so N <= 1 should
% imply N == 1.  If N > 0 weren't guaranteed, this could throw an error.
if N <= 1
  phi = nan(ni, nj);
  phi(qu(1)) = 0; % Leave this isolated pixel at current pressure
  return
end

%% --- Compute epsilon neutrality errors
% Aliases
sm = s;
tm = t;
pm = p;

% Anonymous functions for shifting data.  (Faster than indexing with A4.)
im1 = @(D) circshift(D, [+1 0]);
jm1 = @(D) circshift(D, [0 +1]);
flat = @(F) F(:);

%% m = (i, j) and n = (i-1, j)
sn = im1(sm);
tn = im1(tm);
pn = im1(pm);
if ~WRAP(1)
  sn(1,:) = nan;
end
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (pm + pn) );
eps_iJ = vs .* (sm - sn) + vt .* (tm - tn);
if ~UNIFORM_GRID
  eps_iJ = eps_iJ .* sqrtDIST2on1_iJ;
end

%% m = (i, j) and n = (i, j-1)
sn = jm1(sm);
tn = jm1(tm);
pn = jm1(pm);
if ~WRAP(2)
  sn(:,1) = nan;
end
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (pm + pn) );
eps_Ij = vs .* (sm - sn) + vt .* (tm - tn);
if ~UNIFORM_GRID
  eps_Ij = eps_Ij .* sqrtDIST1on2_Ij;
end

%% --- Build and solve sparse matrix problem

% Collect and sort linear indices to all pixels in this region
m = sort(qu(1 : N));  % sorting here makes mat better structured; overall speedup.

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


% Build the RHS of the matrix problem, containing neutrality errors
rhs = -[eps_iJ(m); eps_Ij(m)];
good_eq = isfinite(rhs);
rhs = [rhs(good_eq); 0]; % Ignore equations where the east or the north water column is invalid.
E = length(rhs) - 1; % number of Equations, excluding density conserving equation.  Note E > 0 is guaranteed, because bfs_conncomp used 4-connectivity

% Begin building the sparse matrix.  The k'th value is v(k) in row r(k)
% and column c(k).  For each k from 1 to E, there is a central grid point
% at (i, j), say. The first E entries of r,c,v are for the points (i-1,j)
% and then (i,j-1).  The second set of E entries are for the points (i,j)
% and (i,j) again.

% Build indices for the rows of the sparse matrix.
if ZERO_MEAN
  r = [flat(repelem((1:E), 1, 2)); repelem(E+1, N, 1)]; % [1; 1; 2; 2; ...; E; E; E+1; ...; E+1]
else % PINNING
  r = [flat(repelem((1:E), 1, 2)); E+1];                % [1; 1; 2; 2; ...; E; E; E+1;]
end

% Build indices for the columns of the sparse matrix
c = [remap(A4(MJ,m)), remap(A4(IM,m)); 1:N, 1:N];
if ZERO_MEAN
  c = [flat(c(:,good_eq)); (1:N).']; % Augment with equation for the specific volume mean = 0
else % PINNING
  c = [flat(c(:,good_eq)); 1];       % Augment with equation for the specific volume at cast 1 = 0
end


% Build the values of the sparse matrix
if UNIFORM_GRID
  if ZERO_MEAN
    v = [flat(repelem([-1; 1], 1, E)); repelem(F, N, 1)];  % [-1; 1; -1; 1; ...; -1; 1; F; ...; F];
  else % PINNING
    v = [flat(repelem([-1; 1], 1, E)); F ];                % [-1; 1; -1; 1; ...; -1; 1; F];
  end
else
  d = [sqrtDIST2on1_iJ(m).', sqrtDIST1on2_Ij(m).'];
  d = d(good_eq);
  if ZERO_MEAN
    v = [flat([-d; d]); repelem(F, N, 1)];
  else % PINNING
    v = [flat([-d; d]); F];
  end
end

% Build the sparse matrix, with E+1 rows and N columns
mat = sparse( r, c, v, E+1, N );

% Solve the LSQR problem
[sol,flag] = omega_lsqr(mat, rhs, TOL_LSQR_REL);
if flag > 0
  warning('omega_surface:lsqr did not converge.');
end

% Force phi(mr) = 0 exactly.  This keeps the surface pinned at its
% initial depth in water column mr.  This does not affect other
% connected components.
sol = sol - sol(remap(mr));

% Save solution
phi(m) = sol;