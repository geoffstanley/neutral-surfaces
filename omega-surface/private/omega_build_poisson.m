function [D, L, epsL2, eps2] = omega_build_poisson(pm, sm, tm)
% OMEGA_BUILD_POISSON  Build the Poisson matrix problem for omega surfaces
% 
% 
% [D, L] = omega_build_poisson(p, s, t)
% evaluates D, the divergence of epsilon, and L, the matrix representation
% of the discrete Laplacian, for a surface of pressure p, salinity s, and
% temperature t.  L is a compact representation, to be built into a sparse
% matrix by omega_matsolve_poisson.
%
% [D, L, epsL2, eps2] = omega_build_poisson(p, s, t)
% also calculates the sum of the squares of epsilon's components values,
% eps2, and the L2 norm of epsilon, epsL2.
%
% Assumes domain is periodic in both dimensions.

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2019  Geoff Stanley
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


[ni,nj] = size(pm);

eps2 = 0;             % For accumulating the square of epsilon
n_edges = 0;          % The total number of epsilon links
L = zeros(5, ni, nj); % pre-alloc space

CALC_EPSL2 = (nargout >= 3);

% . 2 .
% 1 3 5
% . 4 .
sL = 1; % Left
sU = 2; % Up
sC = 3; % Centre
sD = 4; % Down
sR = 5; % Right

% Anonymous functions for shifting data
im1 = @(D) circshift(D, [+1 0]);
ip1 = @(D) circshift(D, [-1 0]);
jm1 = @(D) circshift(D, [0 +1]);
jp1 = @(D) circshift(D, [0 -1]);


%% m = (i, j) and n = (i, j-1),  then also n = (i, j+1) by symmetry
sn = jm1(sm);
tn = jm1(tm);
pn = jm1(pm);

[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (pm + pn) );

eps = vs .* (sm - sn) + vt .* (tm - tn);

good = ~isnan(eps);
eps(~good) = 0;
if CALC_EPSL2
  n_edges = n_edges + sum(good(:));
  eps2 = eps2 + sum(eps(:) .* eps(:));
end

D =     eps - jp1(eps);   


L(sC,:,:) = good + jp1(good);       

L(sL,:,:) = -good;

L(sR,:,:) = -jp1(good); 


%% m = (i, j) and n = (i-1, j),  then also n = (i+1, j) by symmetry
sn = im1(sm);
tn = im1(tm);
pn = im1(pm);

[vs, vt] = eos_s_t( 0.5 * (sm + sn), 0.5 * (tm + tn), 0.5 * (pm + pn) );

eps = vs .* (sm - sn) + vt .* (tm - tn);

bad = isnan(eps);
good = double(~bad);
eps(bad) = 0;
if CALC_EPSL2
  n_edges = n_edges + sum(good(:));
  eps2 = eps2 + sum(eps(:) .* eps(:));
end

D = D + eps - ip1(eps);   

L(sC,:,:) = squeeze(L(sC,:,:)) + good + ip1(good);                     

L(sU,:,:) = -good;    

L(sD,:,:) = -ip1(good); 


%% Finish up

if CALC_EPSL2
  epsL2 = sqrt(eps2 / n_edges); % This is sqrt(\mathcal{E}), where E = 1/2 sum_m sum_{n next to m} \epsilon_{m,n}^2)
end

% For any m where all neighbours are NaN, set L(sC,m) to 1 so that this
% equation amounts to:  1 * phi(m) = 0. This keeps phi(m) = 0, rather than
% becoming NaN and infecting its neighbours.
L(sC, L(sC,:,:) == 0) = 1;

end
