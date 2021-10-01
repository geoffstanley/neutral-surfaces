function C = ppc_deriv(X, C, J) %#codegen
%PPC_DERIV  Differentiate piecewise polynomial.
%
%
% D = ppc_deriv(X, C, J)
% returns the coefficients D for the piecwise polynomial that is the J'th
% derivative of the piecewise polynomial whose coefficients are C.  Both
% piecewise polynomials have knots given by X.
%
%
% --- Input:
% X [    K   x M], knots of the piecewise polynomials
% C [O x K-1 x N], coefficients of the piecewise polynomial
% J [1, 1],        degree of the derivatives 
%
%
% --- Output:
% D [O-J x K-1 x N], coefficients for the derivative
%
%
% --- Notes:
% X(:,n) must be monotonically increasing for all n.
%
% NaN's in X are treated as +Inf (and as such must come at the end of each
% column).
%
% M, the last dimension of X, can be either 1 or N.
%
% Any dimension N can actually be higher-dimensional, so long as it has N
% elements.
%
%
% --- Acknowledgements:
% This code is adapted from PPDIFF by Jonas Lundgren, available as part of
% SPLINEFIT, available on the MATLAB File Exchange.

% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 24/10/2019 - initial release

% Get and check sizes
szC = size(C);
szX = size(X);
O = szC(1);  % Order of the piecewise polynomial
K = szX(1);  % number of knots of the piecewise polynomials
N = prod(szC(3:end)); % Number of parallel problems
XN = prod(szX(2:end));
assert(szC(2) == K-1, 'C must have one fewer columns than X has rows.');
assert(XN == 1 || XN == N,  'X must have %d columns or 1 column', N);
C = reshape(C, [O, K-1, N]); % Just for speed

% Get and check diff order
if nargin == 3 
    assert(isreal(J) && mod(J,1) == 0 && J >= 0, 'J must be a non-negative integer!');
else
    J = 1;
end


if J == 0
    % Do nothing
elseif J < O % J < Order
    % Derivative of order J
    D = [O-J : -1 : 1; ones(J-1, O-J)];
    D = cumsum(D, 1);
    D = prod(D, 1);
    C = C(1:O-J, :);
    for o = 1 : O-J
        C(o,:) = D(o) * C(o,:);
    end
else % J >= Order, so annihilates the piecewise polynomial
    C = zeros(1, K-1, N);
end

C = reshape(C, [max(O-J,1), K-1, szC(3:end)]);