function [x,s,t] = neutral_trajectory(S, T, X, x0, s0, t0, interpfn, tolx)
% NEUTRAL_TRAJECTORY  Calculate a neutral trajectory through a sequence of casts.
%
%
% [x,s,t] = neutral_trajectory(S, T, X, x0)
% calculates a discrete neutral trajectory through the consecutive casts
% (S(:,c), T(:,c), X(:,c)) for increasing c, beginning at depth or pressure
% x0 on cast c=1.  The output are 1D arrays x, s, and t, whose c'th
% elements provide the depth / pressure, salinity, and temperature values
% on the c'th cast along the neutral trajectory. The equation of state for
% density is given by eos.m in the PATH.
%
% [x,s,t] = neutral_trajectory(S, T, X, x0, s0, t0)
% as above, but the first step is a discrete neutral trajectory from the
% bottle (s0, t0, x0) to the cast (S(:,1), T(:,1), X(:,1)).
%
% ... = neutral_trajectory(..., interpfn)
% uses interpfn (a function handle) to interpolate S and T as piecewise
% polynomials of X. By default, interpfn = @ppc_linterp to use linear
% interpolation. Other functions from the PPC toolbox can be used, e.g.
% @ppc_pchip and @ppc_makima.
%
% ... = neutral_trajectory(..., tolx)
% evaluates the discrete neutral trajectory with an accuracy of tolx [m or dbar].
% By default, tolx = 1e-6 [m or dbar]. 
%
% Provide [] for any optional arguments that are required only to provide 
% a value for an argument later in the list. 
%
%
% --- Input:
%  S [nk, nc]: practical / Absolute Salinity values on a cast
%  T [nk, nc]: potential / Conservative Temperature values on a cast
%  X [nk, nc]: pressure / depth values on a cast
%  x0 [1, 1]: pressure / depth of starting bottle
%  s0 [1, 1]: practical / Absolute Salinity of starting bottle
%  t0 [1, 1]: potential / Conservative Temperature of starting bottle
%  interpfn [function handle]: function to calcualte piecewise polynomial coefficients
%  tolx [1, 1]: error tolerance in vertical for neutral trajectory calculations
%
%
% --- Output:
%  x [1, nc]: pressure / depth along the neutral trajectory
%  s [1, nc]: practical / Absolute Salinity along the neutral trajectory
%  t [1, nc]: potential / Conservative Temperature along the neutral trajectory

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


[nk,nc] = size(S);
assert(all(size(T) == size(S)), 'T must be same size as S')
assert(all(size(X) == size(S)) || all(size(X) == [nk, 1]), 'X must be [nk,nc] or [nk,1]');
X_2D = ~isvector(X);
Xc = X(:,1);
c1 = 1;

s = nan(1,nc);
t = nan(1,nc);
x = nan(1,nc);

if nargin < 7 || isempty(interpfn)
  interpfn = @ppc_linterp; % default: piecewise linear interpolation
end

if nargin < 8 || isempty(tolx)
  tolx = 1e-6;
end

if nargin < 6 || isempty(s0) || isempty(t0)
  % Only a depth x0 was given, not a bottle (s0,t0,x0).
  % Get (s0,t0) by evaluating the first cast at x0. 
  Sc = S(:,1);
  Tc = T(:,1);
  SppX = interpfn(Xc, Sc);
  TppX = interpfn(Xc, Tc);
  [s0, t0] = ppc_val2(Xc, SppX, TppX, x0);
  
  % Skip the first cast, because the ntp from (s0,t0,x0) to the first cast
  % will just return (s0,t0,x0).  Also record (s0,t0,x0) as the answers for
  % the first cast.
  c1 = 2;
  s(1) = s0;
  t(1) = t0;
  x(1) = x0;
  
end

% Use precompiled mex version if it exists. Beware with this, however!  It
% could be using an equation of state other than the current eos.m on your
% path.  If you change eos.m, be sure to recompile ntp_bottle_to_cast_mex.
if which('ntp_bottle_to_cast_mex')
  b2c = @ntp_bottle_to_cast_mex;
else
  b2c = @ntp_bottle_to_cast;
end

% Loop over casts
for c = c1 : nc
  
  Sc = S(:,c);
  Tc = T(:,c);
  if X_2D
    Xc = X(:,c);
  end
  
  % Interpolate Sc and Tc as piecewise polynomials of X
  SppX = interpfn(Xc, Sc);
  TppX = interpfn(Xc, Tc);
  
  % Make a neutral connection from (s0,t0,x0) to the cast (S(:,c), T(:,c), X(:,c))
  K = sum(isfinite(Sc));
  [x(c), s(c), t(c)] = b2c(SppX, TppX, Xc, K, s0, t0, x0, tolx); 
  
  if isnan(x(c))
    % The neutral trajectory incropped or outcropped
    break
  end
  
  % Prepare for the next cast (mix with environment)
  s0 = s(c);
  t0 = t(c);
  x0 = x(c);

end

end
