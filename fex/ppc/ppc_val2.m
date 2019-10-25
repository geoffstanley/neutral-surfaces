function [y,z] = ppc_val2(X, C, D, x) %#codegen
%PPVALQ2  Piecewise Polynomial Evaluation, quick and twice
%
%
% [y,z] = ppc_val2(X, C, D, x)
% evaluates the piecewise polynomials whose coefficients are C and D and
% whose knots are X, at data sites x.
%
%
% --- Input:
% X [    K   x M], knots of the piecewise polynomials
% C [O x K-1 x N], coefficients of the first piecewise polynomial
% D [O x K-1 x N], coefficients of the second piecewise polynomial
% x [    L   x M], evaluation sites
%
%
% --- Output:
% y [L x N], the first piecewise polynomial evaluated at x
% z [L x N], the second piecewise polynomial evaluated at x
%
%
% --- Notes:
% X(:,n) must be monotonically increasing for all n.
%
% NaN's in X are treated as +Inf (and as such must come at the end of each
% column).
%
% M, the last dimension of X and x, can be either 1 or N.
%
% Any dimension N can actually be higher-dimensional, so long as it has N
% elements.
%
% Even if L == 1, x does need a leading singleton dimension.
%
% If L == 1, this dimension is squeeze()'d out of y and z.
%
%
% --- Acknowledgements:
% This code is adapted from MATLAB's ppval.m

% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 24/10/2019 - initial release

szC = size(C);
szX = size(X);
szx = size(x);
O = szC(1);  % Order of the piecewise polynomial
K = szX(1);  % number of knots of the piecewise polynomials
L = szx(1);  % number of levels to interpolate
N = prod(szC(3:end));
XM = prod(szX(2:end));
xM = prod(szx(2:end));
szy = [L, szC(3:end), 1]; % Add a trailing 1, to ensure 2 dimensions at least

assert(szC(2) == K-1,       'C must have one fewer columns than X has rows.');
assert(xM == 1 || xM == N,  'x must have %d columns or 1 column', N);
assert(XM == 1 || XM == N,  'X must have %d columns or 1 column', N);
assert(all(szC == size(D)), 'C and D must have the same size');

x = reshape(x, L, []);
y = nan([L, N]);
z = nan([L, N]);

% Evaluate each piecewise polynomial
x01 = double(xM > 1);                   % used for linear indexing
MX = K * double(XM > 1);                % used for linear indexing
MC = O * (K - 1) * double(N > 1);       % used for linear indexing
nX = 0 ;                                % used for linear indexing
for n = 0:N-1
    for l = 1:L
        
        xln = x(l + n*L*x01); % x(l,n) or x(l,1) as appropriate
        
        if ~(isnan(xln) || isnan(X(nX + 1)) || xln < X(nX + 1) || xln > X(nX + K))
            
            % Leftmost binary search to find i such that:
            % i = 1                      if xln <= X(1), or
            % i = M                      if X(M) < xln
            % X(i-1,n) < xln <= X(i,n)   otherwise
            % We use leftmost so that NaN's at the end of X are treated as though they are Inf.
            i = 1; % Result will be >= 1 always
            k = K; % Result will be <= M always
            while (i < k)
                j = floor((i + k) / 2);
                if X(nX + j) < xln % [X(j,n) or X(j,1) as appropriate]  < xln
                    i = j + 1;
                else
                    k = j;
                end
            end
            
            ln = l + n * L;
            if i == 1 
                % Note: X(nX + 1) == xln   is guaranteed
                y(ln) = C(O + (i-1)*O + n*MC);      % y(l,n) = C(O, i, n);
                z(ln) = D(O + (i-1)*O + n*MC);
            else
                % Evaluate this piece of the polynomial (see ppval.m)
                t = xln - X(nX + i-1); % Switch to local coordinates
                
                % Overload variable i, to speed with indexing
                % subtract 1 from i so that 1 <= i <= M-1, and X(i) <= xln < X(i+1),
                % subtract another 1 for indexing. 
                i = (i-2)*O + n*MC; 
                y(ln) = C(1 + i);                   % y(l,n) = C(1,i+1,n+1);
                z(ln) = D(1 + i);
                for o = 2:O
                    y(ln) = t * y(ln) + C(o + i);   % y(l,n) = s * y(l,n) + C(o,i+1,n+1);
                    z(ln) = t * z(ln) + D(o + i);   
                end
            end

        end
    end % for l
    
    nX = nX + MX;
end % for n

% Reshape output to be like input.  Also remove leading dimensions if L ==
% 1, but leave row vectors as row vectors.
y = squeeze(reshape(y, szy));
z = squeeze(reshape(z, szy));
