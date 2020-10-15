cd ~/work/projects-gfd/neutral-surfaces/omega-surface/private/
format long
ni = 5; 
nj = 5;
nij = ni * nj;
WRAP = [1; 1];
% Set up indices used when building the matrix mat
shift_im1 = @(F) circshift(F, [+1, 0]); % if x is [lat x lon], this shifts south
shift_jm1 = @(F) circshift(F, [0, +1]); % if x is [lat x lon], this shifts west
idx = reshape(1:nij, ni, nj); % Linear index to each pixel on the full grid (ocean or not)
idx_im1 = shift_im1(idx);     % Linear index to the pixel above (^)   the local pixel
idx_jm1 = shift_jm1(idx);     % Linear index to the pixel left (<) of the local pixel


good = logical([ ...
  0 0 0 0 0;
  0 0 1 1 0;
  0 1 1 1 0;
  0 0 1 1 0;
  0 0 0 0 0]);

good = flipud(good)'; % convert from map thinking to matrix layout

i0 = 3;
j0 = 3;
I0 = sub2ind([ni nj], i0, j0);

x = 2000 + 50 * randn(ni, nj);
s = 35 + .1 * randn(ni,nj);
t = 3 + .6 * randn(ni, nj);
x(~good) = nan;
s(~good) = nan;
t(~good) = nan;


%%
TOL_X_UPDATE = 1e-12;

% Pre-calculate things for Breadth First Search
qu = zeros(nij, 1); % queue storing linear indices to pixels
A5 = grid_adjacency([ni,nj], 5, WRAP); % all grid points that are adjacent to all grid points, using 5-connectivity
A4 = A5([1,2,4,5],:); % all grid points that are adjacent to all grid points, using 4-connectivity
idx = nan(ni,nj);

[eps_i, eps_j] = ntp_errors(s, t, x, 1, 1, true, false, WRAP);

[qu, qts, ncc, ~, Lcc] = bfs_conncomp(good, A4, [], qu);
I = sort(qu(1:qts(2)-1));
N = length(I);

% Label the water columns in this region alone by 1, 2, ...
% No need to reset other elements of idx to 0.
idx(I) = 1:N;


rhs = [eps_i(I); eps_j(I)];
good_eq = ~isnan(rhs);
rhs = [rhs(good_eq); 0];

E = length(rhs) - 1; % number of Equations, excluding density conserving equation.  Note E > 0 is guaranteed, because bfs_conncomp used 4-connectivity

% Begin building the sparse matrix.  The k'th value is v(k) in row r(k)
% and column c(k).  For each k from 1 to E, there is a central grid point
% at (i, j), say. The first E entries of r,c,v are for the points (i-1,j)
% and then (i,j-1).  The second set of E entries are for the points (i,j)
% and (i,j) again.

% Build indices for the rows of the sparse matrix.
% r is [1; 2; ...; E; 1; 2; ..., E; E+1; ...; E+1]
r = [flat(repelem((1:E).', 1, 2)); repelem(E+1, N, 1)];

% Build indices for the columns of the sparse matrix
c = [[idx_im1(I); idx_jm1(I)], repmat(I, 2, 1)];
c = idx(c);  % Remap global indices to local indices for this region, numbered 1, 2, ...
c = [flat(c(good_eq,:)); (1:N).']; % Augment with the specific volume conserving equation

% Build the values of the sparse matrix
FINAL_ROW_VALUES = 1e-2;
%FINAL_ROW_VALUES = 1;
v = [flat(repelem([-1, 1], E, 1)); repelem(FINAL_ROW_VALUES, N, 1)];

% Build the sparse matrix, with E+1 rows and N columns
matG = sparse( r, c, v, E+1, N );

solG = omega_lsqr(matG, rhs, 1e-6)

matG = full(matG);

% Change last row to be pinning
%matG(end, :) = 0;
%matG(end, 4) = 1;

solG = matG \ rhs;
phiG = nan(ni,nj); 
phiG(good) = solG;
phiG = phiG - phiG(i0,j0);

solG

%%
[D, L] = omega_build_poisson(x, s, t);

%[phiP, D, L, matP] = omega_matsolve_poisson(D, L, A5, Lcc, I0, qu, qts);
[phiP, D, L, matP] = omega_matsolve_poisson(D, L, A5, Lcc, [], qu, qts);
%[xP, sP, tP] = omega_vertsolve_mex(SppX, TppX, X, BotK, s, t, x, TOL_X_UPDATE, phiP);
matP = full(matP);

phiP2 = (matG' * matG) \ D(good); 

phiP(good) - phiP2 %  they match!  
%% WRONG
% remove mean phi conserving row
% replace a row with pinning of cast 4
matG = matG(1:end-1,:); 
matG(end, :) = 0;
matG(end, 4) = 1;

matG' * matG