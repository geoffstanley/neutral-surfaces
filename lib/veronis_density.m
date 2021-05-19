function d1 = veronis_density(p_ref, S, T, P, p0, p1, dp, interpfn)
% VERONIS_DENSITY  The surface density plus the integrated vertical 
%                  gradient of Locally Referenced Potential Density.
%
%
% d1 = veronis_density(p_ref, S, T, P, p0, p1)
% determines the Veronis density d1 at vertical position p1 on a cast with
% practical / Absolute salinity S and potential / Conservative temperature
% T values at depth or pressure values P.  The Veronis density is given by
% the potential density (with reference pressure / depth p_ref) evaluated
% at p0 on the cast, plus the integral of the vertical (d/dX) derivative of
% Locally Referenced Potential Density (LRPD) from P = p0 to P = p1. The
% vertical (d/dP) derivative of LRPD is rho_S dS/dP + rho_T dT/dP where
% rho_S and rho_T are the partial derivatives of density with respect to S
% and T, and dS/dP and dT/dP are the derivatives of S and T with respect to
% P.  The equation of state for density is given by eos.m in the PATH, and
% its partial derivatives with respect to S and T are given by eos_s_t.m in
% the PATH.  If p0 or p1 are outside the range of P, d1 is returned as NaN.
% 
% d1 = veronis_density(..., dp)
% specifies the maximum interval size used in the trapezoidal numerical
% integration.  If omitted, the default size is 1 unit of P (1m or 1 dbar).

% d1 = veronis_density(..., interpfn)
% uses interpfn (a function handle) to interpolate S and T as piecewise polynomials of P.
% If interpfn = @ppc_linterp, the result is the same as if interpfn were omitted
% and linear interpolation were performed native to this code.  Other functions
% from the PPC toolbox can be used, e.g. ppc_pchip and ppc_makima.
%
% Provide [] for any optional arguments that are required only to provide 
% a value for an argument later in the list. 
%
%
% --- Input:
% p_ref [1,1]: reference pressure / depth to evaluate potential density at p0
% S [nk, nt] : practical / Absolute Salinity values on a cast
% T [nk, nt] : potential / Conservative Temperature values on a cast
% P [nk, nt] : pressure / depth values on a cast
% p0 [1, 1]  : pressure / depth that starts the integral
% p1 [1, 1]  : pressure / depth that ends the integral
% dp [1, 1]  : maximum interval of pressure / depth in numerical integration
% interpfn [function handle]: function to calcualte piecewise polynomial coefficients
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
% Stanley, G. J., McDougall, T. J., & Barker, P. M. (2021). Algorithmic
% improvements to finding approximately neutral surfaces. Journal of
% Advances in Modeling Earth Systems, submitted.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com

assert(all(size(T) == size(S)), 'T must be same size as S')
assert(all(size(P) == size(S)), 'P must be same size as S')
assert(isvector(S), 'S, T, P must be 1D. (Veronis density is only useful for one water column at a time!)');
assert(isscalar(p0), 'p0 must be a scalar');
assert(isscalar(p1), 'p1 must be a scalar');

if nargin < 7 || isempty(dp)
  dp = 1; % default: 1m or 1dbar step size for trapezoidal integration
end

if nargin < 8 || isempty(interpfn)
  interpfn = @ppc_linterp; % default: piecewise linear interpolation
end


% Interpolate S and T as piecewise polynomials of P
Sppc = interpfn(P, S);
Tppc = interpfn(P, T);

% Evaluate the derivatives of these piecewise polynomials as functions of P
SPppc = ppc_deriv(P, Sppc);
TPppc = ppc_deriv(P, Tppc);

% Calculate potential density, referenced to p_ref, at p0
[s0, t0] = ppc_val2(P, Sppc, Tppc, p0);
d0 = eos(s0, t0, p_ref);

% Find indices at which p0 and p1 sit within P
k0 = discretize(p0, P);  % P(k0) <= p0 < P(k0+1);  but if p0 == P(end), then k0 = length(P) - 1;  else, k0 = NaN.
k1 = discretize(p1, P);  % P(k1) <= p1 < P(k1+1);  but if p1 == P(end), then k1 = length(P) - 1;  else, k1 = NaN.
if isnan(k0) || isnan(k1)
  d1 = nan;
  return
end

% Integrate from p0 to P(k0+1)
d1 =   d0 + int_x_kp1(P, p0, k0, dp, Sppc, Tppc, SPppc, TPppc);

% Integrate from P(k0+1) to P(k1+1)
for k = k0+1 : k1
  d1 = d1 + int_x_kp1(P, P(k), k, dp, Sppc, Tppc, SPppc, TPppc);
end

% Integrate from p1 to P(k1+1), and subtract this
d1 =   d1 - int_x_kp1(P, p1, k1, dp, Sppc, Tppc, SPppc, TPppc);

end


function d = int_x_kp1(P, p, k, dp, Sppc, Tppc, SPppc, TPppc)
% Integrate from p to P(k+1) using trapezoidal integration with spacing dp

n = ceil((P(k+1) - p) / dp) + 1; % # points between p and P(k+1), inclusive
p_ = linspace(p, P(k+1), n);     % intervals are not larger than dp

% Use piecewise polynomial coefficients as provided
[s_, t_] = ppc_val2_(P, Sppc, Tppc, p_, k+1);
[dsdp_, dtdp_] = ppc_val2_(P, SPppc, TPppc, p_, k+1);

% To use linear interpolation internally, replace the above 2 lines with the following 9 lines:
% S = Sppc(end,:);
% T = Tppc(end,:);
% dP = P(k+1) - P(k);
% dsdp_ = (S(k+1) - S(k)) / dP;
% dtdp_ = (T(k+1) - T(k)) / dP;
% s0 = ( S(k) * (P(k+1) - p) + S(k+1) * (p - P(k)) ) / dP;
% t0 = ( T(k) * (P(k+1) - p) + T(k+1) * (p - P(k)) ) / dP;
% s_ = linspace(s0, S(k+1), n);
% t_ = linspace(t0, T(k+1), n);

[rs_, rt_] = eos_s_t(s_, t_, p_);
y_ = rs_ .* dsdp_ + rt_ .* dtdp_;
d = trapz(p_, y_);

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