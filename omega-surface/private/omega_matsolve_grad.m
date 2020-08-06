function phi = omega_matsolve_grad(eps_i, eps_j, idx_im1, idx_jm1, L, I0, qu, qts)
%
%
% phi = omega_matsolve_grad(eps_i, eps_j, idx_im1, idx_jm1, L, I0, qu, qts)
% uses the grid adjacency information idx_im1, idx_jm1 to convert the
% epsilon errors eps_i, eps_j into a sparse matrix, doing so for each
% connected component as provided by qu and qts. Pinning is determined by
% I0, with help from L.  See documentation for ntp_errors and bfs_conncomp
% for full details.  The matrix problem solved here is rectangular, with
% approximately twice as many rows as columns.

% --- Hard-coded parameters

% Relative tolerance for LSQR. Since the matrix problem is overdetermined,
% the relative residual will, in general, exceed this tolerance bound. As
% such, it is challenging to relate LSQR's relative tolerance to physical
% tolerances on density or pressure.  From several numerical tests,
% the default relative tolerance of 1e-6 seems to work well.
TOL_LSQR_REL = 1e-6;

% The matrix has 1's and 0's everywhere, except this value in the final row
FINAL_ROW_VALUES = 1e-2; % chosen empirically from tests on 1x1deg OCCA data


% --- Build and solve sparse matrix problem
[ni,nj] = size(eps_i);
phi = nan(ni, nj);     % solution to matrix problem
idx = zeros(ni,nj);    % index from 2D map into linear indices for the current connected component
PIN = ~isempty(I0);    % whether to phi(I0) = 0 or not
ncc = length(qts) - 1; % number connected components
for icc = 1 : ncc % icc = index of region
    
  % Collect and sort linear indices to all pixels in this region
  I = sort(qu(qts(icc) : qts(icc+1)-1));  % sorting here makes mat better structured; overall speedup.
  
  nwc = length(I);  % Number of water columns
  if nwc <= 1  % There are definitely no equations to solve
    phi(I) = 0; % Leave this isolated pixel at current pressure
    continue
  end
  
  % Label the water columns in this region alone by 1, 2, ...
  % No need to reset other elements of idx to 0.
  idx(I) = 1:nwc;
  
  
  % Build the RHS of the matrix problem, containing neutrality errors
  rhs = [eps_i(I); eps_j(I)];
  good_eq = ~isnan(rhs);
  rhs = [rhs(good_eq); 0]; % Ignore equations where the east or the north water column is invalid.
  neq = length(rhs) - 1; % Number of equations, excluding density conserving equation.  Note neq > 0 is guaranteed, because bfs_conncomp used 4-connectivity
  
  % Build indices for the columns of the sparse matrix
  jj = [[idx_im1(I); idx_jm1(I)], repmat(I, 2, 1)];
  jj = idx(jj);  % Remap global indices to local indices for this region, numbered 1, 2, ...
  jj = [flat(jj(good_eq,:)); (1:nwc).']; % Augment with the specific volume conserving equation
  
  % Build indices for the rows of the sparse matrix
  ii = [flat(repelem((1:neq).', 1, 2)); repelem(neq+1, nwc, 1)];
  
  % Build the values of the sparse matrix
  vv = [flat(repelem([-1, 1], neq, 1)); repelem(FINAL_ROW_VALUES, nwc, 1)];
  
  
  % Build the sparse matrix, with neq+1 rows and nwc columns
  mat = sparse( ii, jj, vv, neq+1, nwc );
  
  % Solve the LSQR problem
  [sol,flag] = omega_lsqr(mat, rhs, TOL_LSQR_REL);
  if flag > 0
    warning('omega_surface:lsqr did not converge on iter %d, region %d', iter, icc);
  end
  
  
  if PIN && icc == L(I0)
    % Force phi(I0) = 0 exactly.  This keeps the surface pinned at
    % its initial depth in water column I0.  Note, there could be
    % multiple connected components, and the reference column exists in
    % only one of them.  The other connected components are heaved up
    % or down by this density (phi) amount.
    sol = sol - sol(idx(I0));
  else
    % Force mean(phi) = 0 exactly. The LSQR solution only tries to keep
    % mean(phi) near zero in a least-squares sense, while also trying to
    % make the neutrality errors zero. When FINAL_ROW_VALUES is badly
    % chosen, |phi|_1 can wobble as iterations proceed, but forcing
    % mean(phi) = 0 seems to help |phi|_1 to decrease monotonically.
    sol = sol - mean(sol);
  end
  
  
  % Save solution
  phi(I) = sol;
  
end % icc