function [x, a, b] = fzero_bisect(f, a, b, t, varargin) %#codegen
%FZERO_BISECT  Root-finding by bisection
%
%
% [x,a,b] = fzero_bisect(f,lb,ub,t) finds a and b such that f(a) and f(b)
% have opposite signs, and lb <= a < b <= ub, and b - a <= t. Also
% returns x = (a+b)/2, satisfying |x-y| <= t/2 where f(y) = 0.
%
% ... = fzero_bisect(f,lb,ub,t,x,...) passes additional inputs directly
% to f.
%
%
% --- Notes:
% f is only evaluated within the interval [lb, ub].
%
% Any NaN or 0 returned by f is treated as a negative, finite number.
%
% If no sign change is found between lb and ub, then NaN is returned for x,
% a, and b.
%
% This function does not test for f(x) == 0 identically. That is
% "infinitely unlikely" to happen in real-world problems, and testing for
% that would only slow down execution. If f(x) actually is zero, bisection
% will continue to narrow in on that root. Similarly, if lb == ub, then NaN
% is returned for x, a, and b, even if f(x) == 0 exactly.
%
% This function is compatible with MATLAB's code generation -- so long f is
% similarly compatible. Many root-finding problems can be solved with
% fzero_bisect by writing another function which calls fzero_bisect inside a
% for loop, and this can be made fast by using code generation on that
% wrapper function. Note that f will only be called with a scalar input as
% its first argument; codegen knows this, and might strip out unnecessary
% code from the function definition underlying f.
%
%
% --- Input:
%   f       : handle to a function that accepts a real scalar as its first
%             input and returns a real scalar
%   lb      : scalar lower bound
%   ub      : scalar upper bound
%   t       : scalar error tolerance in x
%   x       : initial scalar guess for a root of f
%   varargin: All additional inputs are passed directly to f
%
%
% --- Output:
%   x : estimated root, scalar
%   a : lower bound for interval containing a root, scalar
%   b : upper bound for interval containing a root, scalar
%
%
% --- Example:
% Simple bisection between bounds of -0.5 and 1.5 would fail to find the
% root. Starting with a guess of .85 and expanding outwards finds a root.
% This example is shown graphically on the MATLAB File Exchange.
%
% c = poly([0 1]); % a polynomial with roots at 0 and 1
% f = @(x) polyval(c,x); % encapsulate extra parameters
% [a,b] = fzero_guess_to_bounds(f, .85, -.5, 1.5);
% root = fzero_bisect(f, a, b, .05);
%
%
% --- Acknowledgements:
% Bisection code adapted from Sky Sartorius' bisection.m, available at
% https://www.mathworks.com/matlabcentral/fileexchange/28150
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 13/12/2018 - initial release  as bisectguess.m
%             06/08/2020 - conversion into fzero_bisect.m, for use with
%             fzero_guess_to_bounds.m

% Bump up tolerance to ensure input isn't smaller than machine precision
tol = 2.0 * eps * abs ( b ) + t;

% Note: comparing to 0 is much faster than calling sign()
fbpos = f(b, varargin{:}) > 0;
fapos = f(a, varargin{:}) > 0;

if xor(fapos, fbpos)
  
  % Bisect to narrow in on the root
  while (b - a) > tol
    
    x = (a + b) / 2;
    fxpos = f(x, varargin{:}) > 0;
    
    if xor(fxpos, fbpos)
      a = x;
    else
      b = x;
      fbpos = fxpos;
    end
    
  end
  
  x = (a + b) / 2; % For good measure
  
else
  
  a = nan;
  b = nan;
  x = nan;
end

end