function phi = omega_matsolve_grad(s, t, x, sqrtDIST2on1_iJ, sqrtDIST1on2_Ij, A4, qu, qts, Lcc, m_ref)
% OMEGA_MATSOLVE_GRAD  Build & solve the sparse matrix gradient equations for omega surfaces
%
%
% phi = omega_matsolve_grad(s, t, x, DIST2on1_iJ, DIST1on2_Ij, A4, qu, qts)
% builds and solves the sparse matrix problem for omega surfaces in
% gradient equation form. In each connected component of the surface, phi
% has zero arithmetic, unweighted mean.
%
% phi = omega_matsolve_grad(s, t, x, DIST2on1_iJ, DIST1on2_Ij, A4, qu, qts, Lcc, m_ref)
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
% phi [ni, nj]: density perturbation, providing the least-squares solution
%               to the discrete version of
%               grad phi = - epsilon,
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
% Version   : 2.2.0
%
% Modified by : --
% Date        : --
% Changes     : --

%#ok<*UNRCH>

%% --- Hard-coded parameters

% Relative tolerance for LSQR. Since the matrix problem is overdetermined,
% the relative residual will, in general, exceed this tolerance bound. As
% such, it is challenging to relate LSQR's relative tolerance to physical
% tolerances on density or pressure.  From several numerical tests,
% the default relative tolerance of 1e-6 seems to work well.
TOL_LSQR_REL = 1e-6;
%TOL_LSQR_REL = 1e-10; % DEV

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
[ni,nj] = size(x);

% If both gridding variables are 1, then grid is uniform
UNIFORM_GRID = ...
  isscalar(sqrtDIST2on1_iJ) && sqrtDIST2on1_iJ == 1 && ...
  isscalar(sqrtDIST1on2_Ij) && sqrtDIST1on2_Ij == 1;

%% --- Compute epsilon neutrality errors
% Aliases
sm = s;
tm = t;
xm = x;

% Anonymous functions for shifting data.  (Faster than indexing with A5.)
im1 = @(D) circshift(D, [+1 0]);
jm1 = @(D) circshift(D, [0 +1]);
flat = @(F) F(:);

%% m = (i, j) and n = (i-1, j)
sn = im1(sm);
tn = im1(tm);
xn = im1(xm);
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (xm + xn) );
eps_iJ = vs .* (sm - sn) + vt .* (tm - tn);
if ~UNIFORM_GRID
  eps_iJ = eps_iJ .* sqrtDIST2on1_iJ;
end

%% m = (i, j) and n = (i, j-1)
sn = jm1(sm);
tn = jm1(tm);
xn = jm1(xm);
[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (xm + xn) );
eps_Ij = vs .* (sm - sn) + vt .* (tm - tn);
if ~UNIFORM_GRID
  eps_Ij = eps_Ij .* sqrtDIST1on2_Ij;
end

%% --- Build and solve sparse matrix problem
phi = nan(ni, nj);     % solution to matrix problem
remap = zeros(ni,nj);  % remap indices from 2D into linear indices for the current connected component
PIN = ~isempty(m_ref); % whether to phi(m_ref) = 0 or not
ncc = length(qts) - 1; % number connected components
for icc = 1 : ncc % icc = index of region
    
  % Collect and sort linear indices to all pixels in this region
  m = sort(qu(qts(icc) : qts(icc+1)-1));  % sorting here makes mat better structured; overall speedup.
  
  N = length(m);  % Number of water columns
  if N <= 1  % There are definitely no equations to solve
    phi(m) = 0; % Leave this isolated pixel at current pressure
    continue
  end
  
  % Label the water columns in this region alone by 1, 2, ...
  % No need to reset other elements of remap to 0.
  remap(m) = 1:N;
  
  
  % Build the RHS of the matrix problem, containing neutrality errors
  rhs = [eps_iJ(m); eps_Ij(m)];
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
    r = [flat(repelem((1:E), 1, 2)); repelem(E+1, N, 1)]; % [1; 1; 2; 2;, ...; E; E; ..., E; E+1; ...; E+1]
  else % PINNING
    r = [flat(repelem((1:E), 1, 2)); E+1];                % [1; 1; 2; 2;, ...; E; E; ..., E; E+1;]
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
    warning('omega_surface:lsqr did not converge in region %d', icc);
  end
  
  
  if PIN && icc == Lcc(m_ref)
    % Force phi(m_ref) = 0 exactly.  This keeps the surface pinned at its
    % initial depth in water column m_ref.  This does not affect other
    % connected components.
    sol = sol - sol(remap(m_ref));
  else
    % Force mean(sol) = 0 exactly.  This can be used even if using pinning
    % in the matrix equation. 
    sol = sol - mean(sol);
  end
  
  
  % Save solution
  phi(m) = sol;
  
end % icc