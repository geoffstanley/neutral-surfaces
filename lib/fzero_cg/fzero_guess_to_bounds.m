function [a, b] = fzero_guess_to_bounds(f, x, A, B, varargin) %#codegen
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
% [a, b] = fzero_guess_to_bounds(f, x, A, B)
% as above, but limits [a,b] to remain inside the subset [A, B].  If x is
% outside of [A, B], it is immediately moved into this range. If no
% sign-change is found within [A, B], then a = nan and b = nan are
% returned.  Note, as above, it is possible that a sign-change is skipped
% over as f is only evaluated at finitely many x values.
%
% [a,b] = fzero_guess_to_bounds(f, x, A, B, ...)
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
%   A       : scalar lower bound
%   B       : scalar upper bound
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


if nargin >= 4 && isscalar(A) && isscalar(B)

  fa = f(A, varargin{:});
  if fa == 0
    a = A;
    b = A;
    return
  end
  
  fb = f(B, varargin{:});
  if fb == 0
    a = B;
    b = B;
    return
  end
  
  x = min(max(x, A), B);
  
  % bounds are given
  dxp = (B - x) / 50;
  dxm = (x - A) / 50;
  
  % Set a = x, except when x is so close to A that machine roundoff makes dxm identically 0,
  % which would lead to an infinite loop below.  In this case, set a = A.
  if dxm == 0
    a = A;
  else
    a = x;
  end
  fapos = f(a, varargin{:}) > 0;
  
  % Similarly, set b = x, except for machine precision problems.
  if dxp == 0
    b = B;
    fbpos = f(b, varargin{:}) > 0;
  else
    b = x;
    if dxm == 0
      fbpos = fapos; % since a = b = x
    else
      fbpos = f(b, varargin{:}) > 0;
    end
  end
  
  while true
    if a > A
      % Move a left, and test for a sign change
      dxm = sqrttwo * dxm;
      a = max(x - dxm, A);
      fapos = f(a, varargin{:}) > 0;
      if xor(fapos, fbpos) % fa and fb have different signs
        break
      end
    elseif b == B % also a == A, so cannot expand anymore
      if xor(fapos, fbpos) % one last test for sign change
        break
      else % no sign change found
        a = nan; b = nan;
        return
      end
    end
    
    if b < B
      % Move b right, and test for a sign change
      dxp = sqrttwo * dxp;
      b = min(x + dxp, B);
      fbpos = f(b, varargin{:}) > 0;
      if xor(fapos, fbpos) % fa and fb have different signs
        break
      end
    elseif a == A % also b == B, so cannot expand anymore
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