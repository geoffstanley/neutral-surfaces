function [p,s,t] = neutral_trajectory(S, T, P, p0, s0, t0, interpfn, tolp)
% NEUTRAL_TRAJECTORY  Calculate a neutral trajectory through a sequence of casts.
%
%
% [p,s,t] = neutral_trajectory(S, T, P, p0)
% calculates a discrete neutral trajectory through the consecutive casts
% (S(:,c), T(:,c), P(:,c)) for increasing c, beginning at depth or pressure
% p0 on cast c=1.  The output are 1D arrays p, s, and t, whose c'th
% elements provide the depth / pressure, salinity, and temperature values
% on the c'th cast along the neutral trajectory. The equation of state for
% density is given by eos.m in the PATH.
%
% [p,s,t] = neutral_trajectory(S, T, P, p0, s0, t0)
% as above, but the first step is a discrete neutral trajectory from the
% bottle (s0, t0, p0) to the cast (S(:,1), T(:,1), P(:,1)).
%
% ... = neutral_trajectory(..., interpfn)
% uses interpfn (a function handle) to interpolate S and T as piecewise
% polynomials of P. By default, interpfn = @ppc_linterp to use linear
% interpolation. Other functions from the PPC toolbox can be used, e.g.
% @ppc_pchip and @ppc_makima.
%
% ... = neutral_trajectory(..., tolp)
% evaluates the discrete neutral trajectory with an accuracy of tolp [m or dbar].
% By default, tolp = 1e-6 [m or dbar]. 
%
% Provide [] for any optional arguments that are required only to provide 
% a value for an argument later in the list. 
%
%
% --- Input:
%  S [nk, nc]: practical / Absolute Salinity values on a cast
%  T [nk, nc]: potential / Conservative Temperature values on a cast
%  P [nk, nc]: pressure / depth values on a cast
%  p0 [1, 1]: pressure / depth of starting bottle
%  s0 [1, 1]: practical / Absolute Salinity of starting bottle
%  t0 [1, 1]: potential / Conservative Temperature of starting bottle
%  interpfn [function handle]: function to calcualte piecewise polynomial coefficients
%  tolp [1, 1]: error tolerance in vertical for neutral trajectory calculations
%
%
% --- Output:
%  p [1, nc]: pressure / depth along the neutral trajectory
%  s [1, nc]: practical / Absolute Salinity along the neutral trajectory
%  t [1, nc]: potential / Conservative Temperature along the neutral trajectory

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


[nk,nc] = size(S);
assert(all(size(T) == size(S)), 'T must be same size as S')
assert(all(size(P) == size(S)) || all(size(P) == [nk, 1]), 'P must be [nk,nc] or [nk,1]');
P_2D = ~isvector(P);
Pc = P(:,1);
c1 = 1;

s = nan(1,nc);
t = nan(1,nc);
p = nan(1,nc);

if nargin < 7 || isempty(interpfn)
  interpfn = @ppc_linterp; % default: piecewise linear interpolation
end

if nargin < 8 || isempty(tolp)
  tolp = 1e-6;
end

if nargin < 6 || isempty(s0) || isempty(t0)
  % Only a depth p0 was given, not a bottle (s0,t0,p0).
  % Get (s0,t0) by evaluating the first cast at p0. 
  Sc = S(:,1);
  Tc = T(:,1);
  Sppc = interpfn(Pc, Sc);
  Tppc = interpfn(Pc, Tc);
  [s0, t0] = ppc_val2(Pc, Sppc, Tppc, p0);
  
  % Skip the first cast, because the ntp from (s0,t0,p0) to the first cast
  % will just return (s0,t0,p0).  Also record (s0,t0,p0) as the answers for
  % the first cast.
  c1 = 2;
  s(1) = s0;
  t(1) = t0;
  p(1) = p0;
  
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
  if P_2D
    Pc = P(:,c);
  end
  
  % Interpolate Sc and Tc as piecewise polynomials of P
  Sppc = interpfn(Pc, Sc);
  Tppc = interpfn(Pc, Tc);
  
  % Make a neutral connection from (s0,t0,p0) to the cast (S(:,c), T(:,c), P(:,c))
  K = sum(isfinite(Sc));
  [p(c), s(c), t(c)] = b2c(Sppc, Tppc, Pc, K, s0, t0, p0, tolp); 
  
  if isnan(p(c))
    % The neutral trajectory incropped or outcropped
    break
  end
  
  % Prepare for the next cast (mix with environment)
  s0 = s(c);
  t0 = t(c);
  p0 = p(c);

end

end
