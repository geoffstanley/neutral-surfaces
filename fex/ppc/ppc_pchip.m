function C = ppc_pchip(X, Y, x) %#codegen
%PPC_PCHIP  Piecewise Cubic Hermite Interpolant
%
%
% C = ppc_pchip(X,Y)
% builds the coefficients C of a PCHIP that interpolates each column of Y
% as a function of each column of X.
%
% y = ppc_pchip(X,Y,x)
% evaluates the PCHIP formed from C = ppc_pchip(X,Y) at each column of x. That
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
% This code is adapted from MATLAB's pchip.m and pwch.m

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
h = reshape(diff(X, 1, 1), K-1, XN);            % distance between data sites
delta = diff(Y, 1, 1) ./ h;                     % slope of linear interpolant between data sites
BotK = min(sum(isfinite(X)), sum(isfinite(Y))); % Index of bottom valid data point in each column

% Loop over columns, each of which is an interpolation problem
d = zeros(K,1);              % slope of interpolant at data sites
OKm1 = 4 * (K-1);            % used for linear indexing
Kh = (K-1) * double(XN > 1); % used for linear indexing
KY =  K    * double(YN > 1); % used for linear indexing
nh = 0 ;                     % used for linear indexing
nY = 0 ;                     % used for linear indexing
for n = 1:N
    k = BotK(n);
    if k > 2
        
        % Calculate PCHIP slopes (see pchip.m) 
        %  Slopes at end points: 
        %   Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.
        %  Slopes at interior points:
        %   d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
        %   d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.
        
        del = delta(:,n);
        
        d(1) = ((2*h(nh+1)+h(nh+2))*del(1) - h(nh+1)*del(2))/(h(nh+1)+h(nh+2));
        if sign(d(1)) ~= sign(del(1))
            d(1) = 0;
        elseif (sign(del(1)) ~= sign(del(2))) && (abs(d(1)) > abs(3*del(1)))
            d(1) = 3*del(1);
        end
        
        for i = 2 : k-1
            if sign(del(i-1)) * sign(del(i)) > 0
                w1 = h(nh+i-1) + 2*h(nh+i);
                w2 = 2*h(nh+i-1) + h(nh+i);
                d(i) = (w1+w2) / (w1 / del(i-1) + w2 / del(i));
            else
                d(i) = 0;
            end
        end
        
        d(k) = ((2*h(nh+k-1)+h(nh+k-2))*del(k-1) - h(nh+k-1)*del(k-2))/(h(nh+k-1)+h(nh+k-2));
        if sign(d(k)) ~= sign(del(k-1))
            d(k) = 0;
        elseif (sign(del(k-1)) ~= sign(del(k-2))) && (abs(d(k)) > abs(3*del(k-1)))
            d(k) = 3*del(k-1);
        end
        
        
        % Build piecewise cubic Hermite polynomial (see pwch.m)
        nn = (n-1)*OKm1;
        for i = 1 : k-1
            dzzdx = (del(i) - d(i)) / h(nh+i);
            dzdxdx = (d(i+1) - del(i)) / h(nh+i);
            C(nn + i*4 - 3) = (dzdxdx - dzzdx) / h(nh+i); % C(1, i, n)
            C(nn + i*4 - 2) = 2*dzzdx - dzdxdx;           % C(2, i, n)
            C(nn + i*4 - 1) = d(i);                       % C(3, i, n)
            C(nn + i*4    ) = Y(i + nY);                  % C(4, i, n),  Y(i,n)
        end
        

    elseif k == 2
        % Special case n = 2, use linear interpolation.
        nn = (n-1)*OKm1;
        C(nn + 1) = 0;              % C(1, 1, n)
        C(nn + 2) = 0;              % C(2, 1, n)
        C(nn + 3) = del(1);         % C(3, 1, n)
        C(nn + 4) = Y(1 + nY);      % C(4, 1, n),  Y(i,n)
    end
    
    nh = nh + Kh;
    nY = nY + KY;
end % for n

% Reshape trailing dimensions of C to be like the input
C = reshape(C, szC);

% Evaluate the piecewise polynomial, if requested
if nargin == 3
    C = ppc_val(X, C, x);
end


