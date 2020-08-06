function [phi, G, H] = omega_matsolve_poisson(G, H, A5, L, I0, qu, qts)
%
%
% [phi, G, H] = omega_matsolve_poisson(G, H, A5, L, I0, qu, qts)
% uses the grid adjacency information A5 to convert the compact
% representation of the matrix problem from G and H into a sparse matrix,
% doing so for each connected component as provided by qu and qts. Pinning
% is determined by I0, with help from L.  See documentation for
% linprob_dens and bfs_conncomp for full details.  The matrix problem
% solved here is square, sparse, and symmetric.  This function modifies G
% and H internally, so they are returned as outputs so that MATLAB can
% effectively pass them by reference rather than by value, to avoid making
% duplicate copies inside this function.

% --- Build and solve sparse matrix problem
[ni,nj] = size(G);
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
  
  
  %PIN_HERE = PIN && I(binsrchleft(I0, I)) == I0;  % Alternative test option
  PIN_HERE = PIN && icc == L(I0);
  if PIN_HERE
    % I0 is in this region.
    I1 = I0;
  else
    % I0 is not in this region.  Pin the first cast we see, and adjust later
    I1 = I(1);
  end
  
  % Pin surface at I1 = I0 or I1 = I(1), by changing the I1'th equation
  % to be phi[I1], to be 1 * phi[I1] = 0.
  G(I1) = 0;
  H(:,I1) = 0;
  H(3,I1) = 1;
  
  % The above change renders the I1'th column on all rows irrelevant,
  % since phi[I1] will be zero.  So, we may also set this column to 0,
  % which we do here by setting the appropriate links in H to 0. This
  % maintains symmetry of the matrix, and speeds up solution by a
  % factor of about 2.
  H(1,A5(5,I1)) = 0;
  H(2,A5(4,I1)) = 0;
  H(4,A5(2,I1)) = 0;
  H(5,A5(1,I1)) = 0;
  
  % Build the RHS of the matrix problem
  rhs = G(I);
  
  % Build indices for the rows of the sparse matrix
  ii = repmat(1:nwc, 5, 1);
  
  % Build indices for the columns of the sparse matrix
  % idx() remaps global indices to local indices for this region, numbered 1, 2, ...
  jj = idx(A5(:,I));
  
  % Build the values of the sparse matrix
  vv = H(:,I);
  
  % Build the sparse matrix, with nwc rows and nwc columns
  good = jj > 0; % Ignore connections to dry pixels (though they should have zero J anyway, this is faster)
  mat = sparse( ii(good), jj(good), vv(good), nwc, nwc );
  
  % Solve the matrix problem
  sol = mat \ rhs;
  
  if ~PIN_HERE
    % Adjust solution, not in the region containing the reference cast, to have zero mean.
    sol = sol - mean(sol);
  end
  
  % Save solution
  phi(I) = sol;
  
end % icc