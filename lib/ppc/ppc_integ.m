function C = ppc_integ(X, C, A, B) %#codegen
%PPC_INTEG  Integrate piecewise polynomial
%
%
% I = ppc_integ(X, C, A)
% returns the indefinite integral starting from A of the piecewise
% polynomial with knots X and coefficients C.  The default for A is the
% leftmost knot in X.
%
% D = ppc_integ(X, C, A, B)
% returns the definite integral from A to B of the piecewise
% polynomial with knots X and coefficients C.
%
%
% --- Input:
% X [    K   x M], knots of the piecewise polynomials
% C [O x K-1 x N], coefficients of the piecewise polynomial
% A [    1   x M], lower limit of integral
% B [    L   x M], upper limit of integral
%
%
%
% --- Output:
% I [O+1 x K-1 x N], coefficients for the indefinite integral
%   --or--
% D [L x N], the definite integral
%
%
% --- Notes:
% X(:,n) must be monotonically increasing for all n.
%
% NaN's in X are treated as +Inf (and as such must come at the end of each
% column).
%
% M, the last dimension of X and A and B, can be either 1 or N.
%
% Any dimension N can actually be higher-dimensional, so long as it has N
% elements.
%
% Even if L == 1, A and B do need leading singleton dimensions.
%
% If L == 1, this dimension is squeeze()'d out of D.
%
%
% --- Acknowledgements:
% This code is adapted from PPINT by Jonas Lundgren, available as part of
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
if nargin >= 3
    szA = size(A);
    AM = prod(szA(2:end));
    assert(szA(1) == 1, 'A must have one row');
    assert(AM == 1 || AM == N,  'A must have %d columns or 1 column', N);
end
if nargin == 4
    szB = size(B);
    L = szB(1);      % number of Levels
    BM = prod(szB(2:end));
    assert(BM == 1 || BM == N,  'B must have %d columns or 1 column', N);
end
C = reshape(C, [O, K-1, N]); % Just for speed


% Interval lengths
h = reshape(diff(X, 1, 1), [1, K-1, XN]);

% Integrate each piecewise polynomial
C(1,:) = C(1,:) / O;
y = C(1,:,:) .* h;
for o = 2:O
    C(o,:) = C(o,:) / (O-o+1);
    y = (y + C(o,:,:)) .* h;
end

% Cumulatively sum the integrals pieces, then set the offset
y = cumsum(y,2);
C(O+1,:,:) = cat(2, zeros(1,1,N), y(1,1:K-2,:));

if nargin == 3 && (~isscalar(A) || A ~= X(1))
    % Indefinite integral from A
    C(O+1,:,:) = C(O+1,:,:) - reshape(ppc_val(X,C,A), [1, 1, N]);
    szI = [O+1, K-1, szC(3:end)];
    C = reshape(C, szI);
elseif nargin == 4
    % Definite integral from A to B
    if L == 1
        C = ppc_val(X,C,B) - ppc_val(X,C,A);
    else
        % Add leading singleton dimension as needed
        C = ppc_val(X,C,B) - ppc_val(X,C,A);
    end
    szD = [L, szC(3:end), 1];    % Add a trailing 1, to ensure 2 dimensions at least
    C = squeeze(reshape(C, szD));
end

