function d1 = veronis_density(x_ref, S, T, X, x0, x1, s1, t1, dx, interpfn)
% VERONIS_DENSITY  The surface density plus the integrated vertical 
%                  gradient of Locally Referenced Potential Density.
%
%
% d1 = veronis_density(x_ref, S, T, X, x0, x1)
% determines the Veronis density d1 at vertical position x1 on a cast with
% practical / Absolute salinity S and potential / Conservative temperature
% T values at depth or pressure values X.  The Veronis density is given by
% the potential density (with reference pressure / depth x_ref) evaluated
% at x0 on the cast, plus the integral of the vertical (d/dX) derivative of
% Locally Referenced Potential Density (LRPD) from X = x0 to X = x1. The
% vertical (d/dX) derivative of LRPD is rho_S dS/dX + rho_T dT/dX where
% rho_S and rho_T are the partial derivatives of density with respect to S
% and T, and dS/dX and dT/dX are the derivatives of S and T with respect to
% X.  The equation of state for density is given by eos.m in the PATH, and
% its partial derivatives with respect to S and T are given by eos_s_t.m in
% the PATH. It is assumed here that S and T are piecewise linear
% interpolants in X.
% 
% d1 = veronis_density(x_ref, S, T, X, x0, x1, s1, t1)
% first finds the pressure or depth x2 at which a (discrete) neutral
% trajectory from the bottle (s1,t1,x1) intersects the cast (S,T,X), then
% returns veronis_density(x_ref, S, T, X, x0, x2).  This form allows the cast
% (S,T,X) to not necessarily be from the same dataset as the bottle (s1,
% t1, x1).  If it is, i.e. if (x1,s1,t1) is a point on the (S,T,X) cast,
% then x2 = x1.
%
% d1 = veronis_density(..., dx)
% specifies the maximum interval size used in the trapezoidal numerical
% integration.  If omitted, the default size is 1 unit of X (1m or 1 dbar).

% d1 = veronis_density(..., interpfn)
% uses interpfn (a function handle) to interpolate S and T as piecewise polynomials of X.
% If interpfn = @ppc_linterp, the result is the same as if interpfn were omitted
% and linear interpolation were performed native to this code.  Other functions
% from the PPC toolbox can be used, e.g. ppc_pchip and ppc_makima.
%
% Provide [] for any optional arguments that are required only to provide 
% a value for an argument later in the list. 
%
%
% --- Input:
% x_ref [1,1]: reference pressure / depth to evaluate potential density at x0
%  S [nk, 1]: practical / Absolute Salinity values on a cast
%  T [nk, 1]: potential / Conservative Temperature values on a cast
%  X [nk, 1]: pressure / depth values on a cast
%  x0 [1, 1]: pressure / depth that starts the integral
%  x1 [1, 1]: pressure / depth of point on surface
%  s1 [1, 1]: practical / Absolute Salinity value of point on surface
%  t1 [1, 1]: potential / Conservative Temperature value of point on surface
%  dx [1, 1]: maximum interval of pressure / depth in numerical integration
%  interpfn [function handle]: function to calcualte piecewise polynomial coefficients
%
%
% --- Output:
%  d1 [1, 1]: Veronis density
%
%
% --- Discussion:
% The result of this function can serve as a density label for an
% approximately neutral surface. However, this is NOT the same as a value
% of the Jackett and McDougall (1997) Neutral Density variable. This is
% true even if you were to provide this function with the same cast that
% Jackett and McDougall (1997) used to initially label their Neutral
% Density variable, namely the cast at 188 deg E, 4 deg S, from the Levitus
% (1982) ocean atlas. Some difference would remain, because of differences
% in numerics, and because of a subsequent smoothing step in the Jackett
% and McDougall (1997) algorithm. This function merely allows one to label
% an approximately neutral surface with a density value that is INTERNALLY
% consistent within the dataset where one's surface lives. This function is
% NOT to compare density values against those from any other dataset, such
% as 1997 Neutral Density.
%
%
% --- References:
% Veronis, G. (1972). On properties of seawater defined by temperature,
% salinity, and pressure. Journal of Marine Research, 30(2), 227.
%
% Stanley, G. J., McDougall, T. J., & Barker, P. M. (2020). Algorithmic
% improvements to finding approximately neutral surfaces. Journal of
% Advances in Modeling Earth Systems, submitted.

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



k0 = discretize(x0, X);  % X(k0) <= x0 < X(k0+1);  but if x0 == X(end), then k0 = length(X) - 1;  else, k0 = NaN.
k1 = discretize(x1, X);  % X(k1) <= x1 < X(k1+1);  but if x1 == X(end), then k1 = length(X) - 1;  else, k1 = NaN.

assert(~isnan(k0), 'x0 must be in the range of X');
assert(~isnan(k1), 'x1 must be in the range of X');

if nargin < 10 
  % Piecewise linear interpolation, by default
  interpfn = @ppc_linterp;
end

% Interpolate S and T as piecewise polynomials of X
SppX = interpfn(X, S);
TppX = interpfn(X, T);
SxppX = ppc_deriv(X, SppX);
TxppX = ppc_deriv(X, TppX);


if nargin >= 8 && isscalar(s1) && isscalar(t1)
  % Make a neutral connection from (s1,t1,x1) on the approximately neutral
  % surface to the cast (S,T,X), which need not be from the same oceanic
  % dataset as the approximately neutral surface lives in.
  K = sum(isfinite(S));

  x1 = ntp_bottle_to_cast(SppX, TppX, X, K, s1, t1, x1, 1e-6);
end


if nargin < 9 || isempty(dx)
  dx = 1;  % Upper limit for numerical integration step size
end

% Calculate ipotential density, referenced to x_ref, at x0
[s0, t0] = ppc_val2(X, SppX, TppX, x0);
d0 = eos(s0, t0, x_ref);

% Integrate from x0 to X(k0+1)
d1 =   d0 + int_x_kp1(X, x0, k0, dx, SppX, TppX, SxppX, TxppX);

% Integrate from X(k0+1) to X(k1+1)
for k = k0+1 : k1
  d1 = d1 + int_k_kp1(X, k, dx, SppX, TppX, SxppX, TxppX);
end

% Integrate from x1 to X(k1+1), and subtract this
d1 =   d1 - int_x_kp1(X, x1, k1, dx, SppX, TppX, SxppX, TxppX);

end


function d = int_x_kp1(X, x, k, dx, SppX, TppX, SxppX, TxppX)
% Integrate from x to X(k+1) using trapezoidal integration with spacing dx

n = ceil((X(k+1) - x) / dx) + 1; % # points between x and X(k+1), inclusive
xx = linspace(x, X(k+1), n);     % intervals are not larger than dx

% Use piecewise polynomial coefficients as provided
[ss, tt] = ppc_val2_(X, SppX, TppX, xx, k+1);
[dsdx, dtdx] = ppc_val2_(X, SxppX, TxppX, xx, k+1);

% To use linear interpolation internally, replace the above two lines with the following 7 lines (and also pass S and T as arguments):
%   dX = X(k+1) - X(k);
%   dsdx = (S(k+1) - S(k)) / dX;
%   dtdx = (T(k+1) - T(k)) / dX;
%   s0 = ( S(k) * (X(k+1) - x) + S(k+1) * (x - X(k)) ) / dX;
%   ss = linspace(s0, S(k+1), n);
%   t0 = ( T(k) * (X(k+1) - x) + T(k+1) * (x - X(k)) ) / dX;
%   tt = linspace(t0, T(k+1), n);

[rs, rt] = eos_s_t(ss, tt, xx);
yy = rs .* dsdx + rt .* dtdx;
d = trapz(xx, yy);

end

function d = int_k_kp1(X, k, dx, SppX, TppX, SxppX, TxppX)
% Integrate from X(k) to X(k+1) using trapezoidal integration with spacing dx

n = ceil((X(k+1) - X(k)) / dx) + 1; % # points between x and X(k+1), inclusive
xx = linspace(X(k), X(k+1), n);     % intervals are not larger than dx

% Use piecewise polynomial coefficients as provided
[ss, tt] = ppc_val2_(X, SppX, TppX, xx, k+1);
[dsdx, dtdx] = ppc_val2_(X, SxppX, TxppX, xx, k+1);

% To use linear interpolation internally, replace the above two lines with the following 5 lines (and also pass S and T as arguments):
%   dX = X(k+1) - X(k);
%   dsdx = (S(k+1) - S(k)) / dX;
%   dtdx = (T(k+1) - T(k)) / dX;
%   ss = linspace(S(k), S(k+1), n);
%   tt = linspace(T(k), T(k+1), n);

[rs, rt] = eos_s_t(ss, tt, xx);
yy = rs .* dsdx + rt .* dtdx;
d = trapz(xx, yy);

end


function [y,z] = ppc_val2_(X, C, D, x, i)
% like pcc_val2, but supplies i such that
% i = 1                      if x <= X(1), or
% i = M                      if X(M) < x
% X(i-1) < x <= X(i)         otherwise
%
% where M = length(X)
% Also assumes X is 1D, and that C and D are coefficients for 1D piecewise
% polynomials of X.


O = size(C,1); % order

y = nan(size(x));
z = nan(size(x));

if i == 1
  % Note: X(1) == x   is guaranteed
  y(:) = C(O, i);
  z(:) = D(O, i);
  
else
  
  % Evaluate this piece of the polynomial
  t = x - X(i-1); % Switch to local coordinates
  
  y(:) = C(1, i-1);
  z(:) = D(1, i-1);
  
  for o = 2:O
    y = t .* y + C(o, i-1);
    z = t .* z + D(o, i-1);
  end
end
end