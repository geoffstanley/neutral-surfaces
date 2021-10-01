function C = ppc_linterp(X, Y, x) %#codegen
%PPC_LINTERP  Linear Interpolant
%
%
% C = ppc_linterp(X,Y)
% builds the coefficients C of a piecewise linear function that
% interpolates each column of Y as a function of each column of X.
%
% y = ppc_linterp(X,Y,x)
% evaluates the piecewise linear function formed from C = ppc_linterp(X,Y)
% at each column of x. That is, y(l,n) interpolates Y(:,n), as a function
% of X(:,n), at x(l,n). No extrapolation is done: when x(l,n) is out of
% bounds of X(:,n), y(l,n) is NaN.
%
%
% --- Input:
% X [K x N], the data sites
% Y [K x N], the data values at the the data sites
%
%
% --- Output:
% C [O x K-1 x N], the piecewise polynomial coefficients of order O
%
%
% --- Notes:
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

% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 24/10/2019 - initial release

% O = 2;  % Order of the polynomial pieces

szX = size(X);
szY = size(Y);
K = szX(1);
XN = prod(szX(2:end));
YN = prod(szY(2:end));
if XN > YN
    N = XN;
    szC = [2, K-1, szX(2:end)];
else
    N = YN;
    szC = [2, K-1, szY(2:end)];
end
assert(szY(1) == K, 'X and Y must have the same number of rows.');
assert(XN == 1 || XN == N, 'X must have N=%d columns or 1 column', N);
assert(YN == 1 || YN == N, 'Y must have N=%d columns or 1 column', N);

% Fill last row of C with Y values, and first row with linear slopes
C = repmat(reshape(Y(1:K-1,:), [1, K-1, YN]), [2, 1, N / YN]);
C(1,:,:) = reshape(diff(Y, 1, 1) ./ diff(X, 1, 1), [1, K-1, N]);

% Reshape trailing dimensions of C to be like the input
C = reshape(C, szC);

% Evaluate the piecewise polynomial, if requested
if nargin == 3
    C = ppc_val(X, C, x);
end


