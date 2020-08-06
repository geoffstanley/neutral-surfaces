function [a, b] = fzero_guess_to_bounds(f, x, lb, ub, varargin) %#codegen
%FZERO_GUESS_TO_BOUNDS  Search for a sign change bounding a zero of a
%                       univariate function, expanding geometrically
%                       outward from an initial guess.
%
%
% [a, b] = fzero_guess_to_bounds(f, x)
% finds a < b such that f(a) and f(b) have different sign*, meaning a
% solution exists within the interval [a,b].  The bounds a,b are expanded
% outward in geometric progression from an initial guess for the root of f
% at x. If f evaluates to NaN at any point during the search, then a = nan
% and b = nan are immediately returned.  If the function is genuinely
% single-signed, or even if it is not but its values of opposite sign are
% skipped over, it is possible to enter an infinite loop.  Calling the
% function in this form is therefore not recommended unless you know the
% function will not result in such an infinite loop.
%
% [a, b] = fzero_guess_to_bounds(f, x, lb, ub)
% as above, but limits [a,b] to remain inside the subset [lb, ub].  If x is
% outside of [lb, ub], it is immediately moved into this range. If no
% sign-change is found within [lb, ub], then a = nan and b = nan are
% returned.  Note, as above, it is possible that a sign-change is skipped
% over as f is only evaluated at finitely many x values.
%
% [a,b] = fzero_guess_to_bounds(f, x, lb, ub, ...)
% passes all additional arguments to the function f. 
%
% * Note: for computational speed, herein the "sign" of 0 is considered the
% same as the sign of a negative number.
%
% This function is compatible with MATLAB's code generation.
%
%
% --- Input:
%   f       : handle to a function that accepts a real scalar as its first
%             input and returns a real scalar
%   x       : initial scalar guess for a root of f
%   lb      : scalar lower bound
%   ub      : scalar upper bound
%   varargin: All additional inputs are passed directly to f
%
%
% --- Output:
%   a : lower bound for interval containing a root, scalar
%   b : upper bound for interval containing a root, scalar
%
%
% --- Acknowledgements:
% Expansion from initial guess inspired by MATLAB's fzero.m.
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 01/07/2020 - initial release
%           : 22/07/2020 - fix infinite loop in bounded case, arising from machine precision rounding

% Geometrically expand from the guess x, until a sign change is found
sqrttwo = 1.414213562373095;


if nargin >= 4 && isscalar(lb) && isscalar(ub)

  x = min(max(x, lb), ub);
  
  % bounds are given
  dxp = (ub - x) / 50;
  dxm = (x - lb) / 50;
  
  % Set a = x, except when x is so close to lb that machine roundoff makes dxm identically 0,
  % which would lead to an infinite loop below.  In this case, set a = lb.
  if dxm == 0
    a = lb;
  else
    a = x;
  end
  
  % Similarly, set b = x, except for machine precision problems.
  if dxp == 0
    b = ub;
  else
    b = x;
  end
  
  fbpos = f(b, varargin{:}) > 0;
  fapos = fbpos; % needed for codegen
  while true
    if a > lb
      % Move a left, and test for a sign change
      dxm = sqrttwo * dxm;
      a = max(x - dxm, lb);
      fapos = f(a, varargin{:}) > 0;
      if xor(fapos, fbpos) % fa and fb have different signs
        break
      end
    elseif b == ub % also a == lb, so cannot expand anymore
      if xor(fapos, fbpos) % one last test for sign change
        break
      else % no sign change found
        a = nan; b = nan;
        return
      end
    end
    
    if b < ub
      % Move b right, and test for a sign change
      dxp = sqrttwo * dxp;
      b = min(x + dxp, ub);
      fbpos = f(b, varargin{:}) > 0;
      if xor(fapos, fbpos) % fa and fb have different signs
        break
      end
    elseif a == lb % also b == ub, so cannot expand anymore
      if xor(fapos, fbpos) % one last test for sign change
        break
      else % no sign change found
        a = nan; b = nan;
        return
      end
    end
  end
  
else
  
  % no bounds given
  if x ~= 0
    dx = abs(x) / 50;
  else
    dx = 1/50;
  end
  
  a = x;
  b = x;
  fb = f(b, varargin{:});
  fa = fb; %#ok<NASGU> % needed for codegen
  while true
    
    dx = sqrttwo * dx;
    
    % Move a left, and test for a sign change
    a = x - dx;
    fa = f(a, varargin{:});
    if isnan(fa)  % Outside the valid bounds of the function.  Give up.
      a = nan; 
      b = nan;
      break
    elseif xor(fa > 0, fb > 0) % fa and fb have different signs
      break
    end
    
    % Move b right, and test for a sign change
    b = x + dx;
    fb = f(b, varargin{:});
    if isnan(fb)  % Outside the valid bounds of the function.  Give up.
      a = nan; 
      b = nan;
      break
    elseif xor(fa > 0, fb > 0) % fa and fb have different signs
      break
    end
    
  end
  
end

end