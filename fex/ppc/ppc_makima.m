function C = ppc_makima(X, Y, x) %#codegen
%PPC_MAKIMA  Modified Akima interpolant
%
%
% C = ppc_makima(X,Y)
% builds the coefficients C of a modified Akima piecewise cubic hermite
% interpolating polynomial  that interpolates each column of Y as a
% function of each column of X.
%
% y = ppc_makima(X,Y,x)
% evaluates the modified Akima function formed from C = ppc_makima(X,Y) at
% each column of x. That is, y(l,n) interpolates Y(:,n), as a function of
% X(:,n), at x(l,n). No extrapolation is done: when x(l,n) is out of bounds
% of X(:,n), y(l,n) is NaN.
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
%
% 
% --- Acknowledgements:
% This code is adapted from MATLAB's makima.m and pwch.m

% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 24/10/2019 - initial release

% O = 4;  % Order of the polynomial pieces

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

% Calculations on all X and Y
h = reshape(diff(X, 1, 1), K-1, XN); % distance between data sites
delta = diff(Y, 1, 1) ./ h;          % slope of linear interpolant between data sites
BotK = sum(isfinite(Y));             % Index of bottom valid data point in each column

% Loop over columns, each of which is an interpolation problem
s = zeros(K,1);              % slope of interpolant at data sites
del = zeros(K+3,1);          % working variable for delta
OMm1 = 4 * (K-1);            % used for linear indexing
Kh = (K-1) * double(XN > 1); % used for linear indexing
nh = 0 ;                     % used for linear indexing
for n = 1:N
    k = BotK(n);
    if k > 2
        
        % Calculate Modified Akima slopes (see makima.m)
        % coding note: Unnecessary values for del and w are calculated
        % using data past k, extending to M; this keeps array sizes steady
        % and speeds execution.
        del(3:K+1) = delta(:,n);       
        del(2) = 2*del(3)   - del(4);
        del(1) = 2*del(2) - del(3);
        del(k+2) = 2*del(k+1) - del(k);
        del(k+3) = 2*del(k+2) - del(k+1);
        w = abs(del(2:K+3) - del(1:K+2)) + abs(del(2:K+3) + del(1:K+2))/2;
        for i = 1 : k
            w12 = w(i) + w(i+2);
            if w12 == 0
                s(i) = 0;
            else
                s(i) = (w(i+2) * del(i+1) + w(i) * del(i+2)) / w12;
            end
        end
        
        % Build piecewise cubic Hermite polynomial (see pwch.m)
        nn = (n-1)*OMm1;
        for i = 1 : k-1
            dzzdx = (del(i+2) - s(i)) / h(nh+i);
            dzdxdx = (s(i+1) - del(i+2)) / h(nh+i);
            C(nn + i*4 - 3) = (dzdxdx - dzzdx) / h(nh+i); % C(1, i, n)
            C(nn + i*4 - 2) = 2*dzzdx - dzdxdx;           % C(2, i, n)
            C(nn + i*4 - 1) = s(i);                       % C(3, i, n)
            C(nn + i*4    ) = Y(i + (n-1)*K);             % C(4, i, n),  Y(i,n)
        end

    elseif k == 2
        % Special case n = 2, use linear interpolation.
        nn = (n-1)*OMm1;
        C(nn + 1) = 0;              % C(1, 1, n)
        C(nn + 2) = 0;              % C(2, 1, n)
        C(nn + 3) = del(3);         % C(3, 1, n)
        C(nn + 4) = Y(1 + (n-1)*K); % C(4, 1, n),  Y(i,n)
    end
    
    nh = nh + Kh;
end % for n

% Reshape trailing dimensions of C to be like the input
C = reshape(C, szC);

% Evaluate the piecewise polynomial, if requested
if nargin == 3
    C = ppc_val(X, C, x);
end


