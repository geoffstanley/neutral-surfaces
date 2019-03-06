function [x,a,b] = bisectguess(f,lb,ub,tolx,x,varargin) %#codegen
%BISECTGUESS  Root-finding by bisection, expanding from an initial guess.
%
%
% [x,a,b] = bisectguess(f,lb,ub,tolx) finds a and b such that f(a) and f(b)
% have opposite signs, and lb <= a < b <= ub, and b - a <= tolx. Also
% returns x = (a+b)/2, satisfying |x-y| <= tolx/2 where f(y) = 0.
%
% [x,a,b] = bisectguess(f,lb,ub,tolx,x) uses an initial guess x for a root
% of f, expanding an interval outward from x until a sign change is found
% for f evaluated at the limits of this interval, then bisection proceeds
% within this interval.
%
% ... = bisectguess(f,lb,ub,tolx,x,...) passes additional inputs directly
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
% bisectguess by writing another function which calls bisectguess inside a
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
%   tolx    : scalar error tolerance in x
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
% --- Examples:
% 1.    A very simple example. Simple bisection between bounds of -0.5 and
%       1.5 would fail to find the root. Starting with a guess of .85 and
%       expanding outwards finds a root. This example is shown graphically
%       on the MATLAB File Exchange.
%
% c = poly([0 1]); % a polynomial with roots at 0 and 1
% f = @(x) polyval(c,x); % encapsulate extra parameters 
% root = bisectguess(f, -.5, 1.5, .05, .85);
%
% 2.    Use of codegen to solve many similar problems quickly inside a for
%       loop. Create a file bisectionguessloop.m containing the following
%       function definition:
% function found_roots = bisectionguessloop(actual_roots, lb, ub, tolx, x)
% N = size(actual_roots,2);
% found_roots = nan(1, N);
% for i = 1 : N
%     coefs = poly(actual_roots(:,i));
%     found_roots(i) = bisectguess(@myfun, lb(i), ub(i), tolx, x(i), coefs);
% end
%     function out = myfun(x, c)
%         out = polyval(c, x);
%     end
% end
%
%       Then, run the following script:
% type_r = coder.typeof(0, [2, 2^16], [false, true]); % We'll always have 2 roots, but might solve many such problems
% type_x = coder.typeof(0, [1, 2^16], [false, true]);
% codegen('bisectionguessloop', '-args', {type_r, type_x, type_x, 0, type_x}, '-o', 'bisectionguessloop_mex');
% N = 10; % Number of similar problems to solve. 
% actual_roots = randn(2, N);      % Each problem is to find a root of a quadratic having 2 random roots
% lower_bounds = repmat(-3, 1, N); % Assume each problem has a root greater than -3, say
% upper_bounds = repmat( 3, 1, N); % ... and less than +3, say.
% tolx = 1e-6;                     % get this close to the true root.
% guess_roots = zeros(1, N);       % Initial guess for each problem -- try 0. 
% found_roots = bisectionguessloop_mex(actual_roots, lower_bounds, upper_bounds, tolx, guess_roots);
% disp(actual_roots' - found_roots') % sometimes bisection finds the left root, sometimes the right
%
%
% --- Acknowledgements:
% Expansion from initial guess inspired by MATLAB's fzero.m.
% Bisection code adapted from Sky Sartorius' bisection.m, available at
% https://www.mathworks.com/matlabcentral/fileexchange/28150
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 13/12/2018 - initial release

if nargin >= 5 && isscalar(x) && isfinite(x) && lb < x && x < ub
    % Geometrically expand from the guess x, until a sign change is found
    sqrttwo = 1.414213562373095;
    dxp = (ub - x) / 50;
    dxm = (x - lb) / 50;
    
    a = x;
    b = x;
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
                x = nan; a = nan; b = nan;
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
                x = nan; a = nan; b = nan;
                return
            end
        end
    end
    
else
    
    % No guess provided. Bisection will proceed between the lower and upper
    % bounds, but only if f has opposite signs at these boundaries.
    a = lb;
    b = ub;
    fa = f(a, varargin{:});
    fapos = fa > 0;
    fbpos = f(b, varargin{:}) > 0;
    
    if fapos == fbpos
        x = nan; a = nan; b = nan;
        return
    end
    
end

% Bisect to narrow in on the root
% Note: comparing to 0 is much faster than calling sign()
while (b - a) > tolx
    
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
end