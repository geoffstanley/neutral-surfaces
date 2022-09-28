function C = ppc_spline(X, Y, x) %#codegen
%PPC_SPLINE Cubic Spline Interpolant
%
%
% C = ppc_spline(X,Y)
% builds the coefficients C of a Spline that interpolates each column of Y
% as a function of each column of X.
%
% y = ppc_spline(X,Y,x)
% evaluates the Spline formed from C = ppc_spline(X,Y) at each column of x. That
% is, y(l,n) interpolates Y(:,n), as a function of X(:,n), at x(l,n). No
% extrapolation is done: when x(l,n) is out of bounds of X(:,n), y(l,n) is
% NaN.
%
%
% --- Input:
% X [K x N], the data sites
% Y [K x N], the data values at the the data sites
%
%
% --- Output:
% C [4 x K-1 x N], the piecewise polynomial coefficients of order 4
%
%
% --- Notes:
% This uses not-a-knot conditions for the end slopes.
%
% X(:,n) must be monotonically increasing for all n.
%
% NaN's in X are treated as +Inf (and as such must come at the end of each
% column).
%
% Any input can have a singleton second dimension (N = 1), in which case
% that single value is used for each interpolation problem.
%
% Any dimension N can actually be higher-dimensional, so long as it has N
% elements.
%
% Even if L == 1, x does need a leading singleton dimension.
%
% If L == 1, this dimension is squeeze()'d out of y.
%
%
% --- Acknowledgements:
% This code is adapted from MATLAB's spline.m and pwch.m

% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 14/09/2022- initial release

szX = size(X);
szY = size(Y);
K = szX(1);
XN = prod(szX(2:end));
YN = prod(szY(2:end));
if XN > YN
  N = XN;
  szC = [4, K-1, szX(2:end)];
else
  N = YN;
  szC = [4, K-1, szY(2:end)];
end
assert(szY(1) == K, 'X and Y must have the same number of rows.');
assert(XN == 1 || XN == N, 'X must have N=%d columns or 1 column', N);
assert(YN == 1 || YN == N, 'Y must have N=%d columns or 1 column', N);

C = nan([4, K-1, N]);   % Build with 3 dimensions only, to speed execution
Y = reshape(Y, K, YN);  % Just to speed execution

mmdflag = spparms('autommd');
spparms('autommd',0);

% Calculations on all X and Y
h = reshape(diff(X, 1, 1), K-1, XN);            % distance between data sites
delta = diff(Y, 1, 1) ./ h;                     % slope of linear interpolant between data sites
BotK = min(sum(isfinite(X)), sum(isfinite(Y))); % Index of bottom valid data point in each column

% Prepare things for linear indexing
OKm1 = 4 * (K-1);           
Kh = (K-1) * double(XN > 1);
KX =  K    * double(XN > 1);
KY =  K    * double(YN > 1);
nh = 0 ;
nX = 0 ;
nY = 0 ;

% Loop over columns, each of which is an interpolation problem
for n = 1:N
  k = BotK(n);
  
  if k >= 3
    
    dx = h((nh+1 : nh+k-1)');  % column vector
    divdif = delta(:,n);  % column vector
    x31 = X(3 + nX) - X(1 + nX);
    
    if k > 3
      
      % set up the sparse, tridiagonal, linear system b = ?*c for the slopes
      b = zeros(1,k);
      b(2:k-1) = 3*(dx(2:k-1).*divdif(1:k-2)+dx(1:k-2).*divdif(2:k-1));
      
      % Do not-a-knot conditions -- isempty(endslopes) is true
      xk = X(k + nX) - X(k-2 + nX);
      b(1,1) = ((dx(1)+2*x31)*dx(2)*divdif(1)+dx(1)^2*divdif(2))/x31;
      b(1,k) = (dx(k-1)^2*divdif(k-2)+(2*xk+dx(k-1))*dx(k-2)*divdif(k-1))/xk;
      
      c = spdiags([ [x31;dx(1:k-2);0] ...
        [dx(2);2*(dx(2:k-1)+dx(1:k-2));dx(k-2)] ...
        [0;dx(2:k-1);xk] ],[-1 0 1],k,k);
      
      % sparse linear equation solution for the slopes
      s = b/c;
      
      % Build piecewise cubic Hermite polynomial (see pwch.m)
      nn = (n-1)*OKm1;
      for i = 1 : k-1
        dzzdx = (divdif(i) - s(i)) / dx(i);
        dzdxdx = (s(i+1) - divdif(i)) / dx(i);
        C(nn + i*4 - 3) = (dzdxdx - dzzdx) / dx(i);   % C(1, i, n)
        C(nn + i*4 - 2) = 2*dzzdx - dzdxdx;           % C(2, i, n)
        C(nn + i*4 - 1) = s(i);                       % C(3, i, n)
        C(nn + i*4    ) = Y(i + nY);                  % C(4, i, n) = Y(i,n)
      end
      
    else % k == 3
      % Special case k = 3: interpolant is a parabola
      
      % Build coeffs for polynomial, f(x) = c2 (x - x1)^2 + c1 (x - x1) + c0
      c0 = Y(1 + nY); % constant coeff
      c2 = (divdif(2) - divdif(1)) / x31;  % quadratic coeff
      c1 = divdif(1) - c2 * dx(1); % linear coeff
      
      nn = (n-1)*OKm1;
      C(nn + 1) = 0;              % C(1, 1, n)
      C(nn + 2) = c2;             % C(2, 1, n)
      C(nn + 3) = c1;             % C(3, 1, n)
      C(nn + 4) = c0;             % C(4, 1, n)
      
      % Convert coeffs for polynomial, f(x) = c2 (x - x1)^2 + c1 (x - x1) + c0
      c1 = c1 + 2 * c2 * dx(1);
      c0 = c0 + c1 * dx(1) - c2 * dx(1)^2;
      
      C(nn + 4 + 1) = 0;              % C(1, 2, n)
      C(nn + 4 + 2) = c2;             % C(2, 2, n)
      C(nn + 4 + 3) = c1;             % C(3, 2, n)
      C(nn + 4 + 4) = c0;             % C(4, 2, n)
    end
    
  elseif k == 2
    % Special case k = 2, use linear interpolation.
    nn = (n-1)*OKm1;
    C(nn + 1) = 0;              % C(1, 1, n)
    C(nn + 2) = 0;              % C(2, 1, n)
    C(nn + 3) = delta(1,n);     % C(3, 1, n)
    C(nn + 4) = Y(1 + nY);      % C(4, 1, n),  Y(i,n)
    
  end
  
  nh = nh + Kh;
  nX = nX + KX;
  nY = nY + KY;
end % for n

spparms('autommd',mmdflag);

% Reshape trailing dimensions of C to be like the input
C = reshape(C, szC);

% Evaluate the piecewise polynomial, if requested
if nargin == 3
  C = ppc_val(X, C, x);
end


