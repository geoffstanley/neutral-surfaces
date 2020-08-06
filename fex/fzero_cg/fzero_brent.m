function b = fzero_brent(f, a, b, t, varargin) %#codegen
%FZERO_BRENT  Find a root of a univariate function within a given interval
%             using Brent's method
%
% x = fzero_brent(f,a,b,t) 
% finds x in the interval [a,b] satisfying |x-y| <= t/2 where f(y) = 0.
% f(a) and f(b) must have opposite signs.
%
% ... = fzero_brent(f,lb,ub,t,x,...) passes additional inputs directly
% to f.
%
% This function is compatible with MATLAB's code generation -- so long f is
% similarly compatible. Many root-finding problems can be solved with
% fzero_brent by writing another function which calls fzero_brent inside a
% for loop, and this can be made fast by using code generation on that
% wrapper function. Note that f will only be called with a scalar input as
% its first argument; codegen knows this, and might strip out unnecessary
% code from the function definition underlying f.
% --- Example:
% Simple bisection between bounds of -0.5 and 1.5 would fail to find the
% root. Starting with a guess of .85 and expanding outwards finds a root.
% This example is shown graphically on the MATLAB File Exchange.
%
% c = poly([0 1]); % a polynomial with roots at 0 and 1
% f = @(x) polyval(c,x); % encapsulate extra parameters
% [a,b] = fzero_guess_to_bounds(f, .85, -.5, 1.5);
% root = fzero_brent(f, a, b, .05);
%
%  Discussion:
%
%    The interval [A,B] must be a change of sign interval for F.
%    That is, F(A) and F(B) must be of opposite signs.  Then
%    assuming that F is continuous implies the existence of at least
%    one value C between A and B for which F(C) = 0.
%
%    The location of the zero is determined to within an accuracy
%    of 6 * EPS * abs ( C ) + 2 * T, where EPS is the machine epsilon.
%
%    Thanks to Thomas Secretin for pointing out a transcription error in the
%    setting of the value of P, 11 February 2013.
%
%    Additional parameters given by varargin are passed directly to F.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 September 2019
%    30 June 2020 - Geoff Stanley
%
%  Author:
%
%    Original FORTRAN77 version by Richard Brent.
%    MATLAB version by John Burkardt.
%    Minor changes (passing extra arguments) by Geoff Stanley.
%
%  Reference:
%
%    Richard Brent,
%    Algorithms for Minimization Without Derivatives,
%    Dover, 2002,
%    ISBN: 0-486-41998-3,
%    LC: QA402.5.B74.
%
%  Parameters:
%
%    Input, real A, B, the endpoints of the change of sign interval.
%
%    Input, real T, a positive error tolerance.
%
%    Input, real value = F ( x ), the name of a user-supplied
%    function which evaluates the function whose zero is being sought.
%
%    Output, real VALUE, the estimated value of a zero of
%    the function F.
%

fa = f(a, varargin{:});
fb = f(b, varargin{:});

c = a;
fc = fa;
e = b - a;
d = e;

while ( true )
  
  if ( abs ( fc ) < abs ( fb ) )
    
    a = b;
    b = c;
    c = a;
    fa = fb;
    fb = fc;
    fc = fa;
    
  end
  
  tol = 2.0 * eps * abs ( b ) + t;
  m = 0.5 * ( c - b );
  
  if ( abs ( m ) <= tol || fb == 0.0 )
    break
  end
  
  if ( abs ( e ) < tol || abs ( fa ) <= abs ( fb ) )
    
    e = m;
    d = e;
    
  else
    
    s = fb / fa;
    
    if ( a == c )
      
      p = 2.0 * m * s;
      q = 1.0 - s;
      
    else
      
      q = fa / fc;
      r = fb / fc;
      p = s * ( 2.0 * m * q * ( q - r ) - ( b - a ) * ( r - 1.0 ) );
      q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
      
    end
    
    if ( 0.0 < p )
      q = - q;
    else
      p = - p;
    end
    
    s = e;
    e = d;
    
    if ( 2.0 * p < 3.0 * m * q - abs ( tol * q ) && p < abs ( 0.5 * s * q ) )
      d = p / q;
    else
      e = m;
      d = e;
    end
    
  end
  
  a = b;
  fa = fb;
  
  if ( tol < abs ( d ) )
    b = b + d;
  elseif ( 0.0 < m )
    b = b + tol;
  else
    b = b - tol;
  end
  
  fb = f(b, varargin{:});
  
  if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
    c = a;
    fc = fa;
    e = b - a;
    d = e;
  end
  
end

end