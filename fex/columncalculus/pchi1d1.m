function [V,D] = pchi1d1(U,X,Y) %#codegen
%PCHI1D1  Piecewise Cubic Hermite Interpolation and Derivative along all columns.
%
%
% [V,D] = pchi1d1(U,X,Y)
% finds the shape-preserving piecewise cubic interpolants P_n(x) such that
% P_n(X(m,n)) = Y(m,n) for each m = 1 : M, then evaluates them to return
% V(l,n) = P_n(U(l,n)) and D(l,n) = P_n'(U(l,n)). When U(l,n) is outside
% the valid range of X(:,n), V(l,n) and D(l,n) are NaN. Only data in X and
% Y before a NaN is used: e.g. if Y(m+1,n) is NaN, only X(1:m,n) and
% Y(1:m,n) is used.
%
%
% --- Input:
% U [L x N], the interpolant evaluation sites
% X [M x N], the data sites for the function
% Y [M x N], the values of the function at the data sites
%
%
% --- Output:
% V [L x N], the interpolant values
% D [L x N], the derivative
%
%
% --- Notes:
% X(:,n) must be monotonically increasing for all n. 
%
% NaN's in X are treated as +Inf (and as such must come at the end of each
% row).
%
% Any input can be higher-dimensional, so long as the product of its
% dimensions excluding the first dimension is N.
%
% Any input can have a singleton second dimension (N = 1), in which case
% that single value is used for each interpolation problem.
%
% Even if L == 1, U does need a leading singleton dimension.
%
% If L == 1, this dimension is squeeze()'d out of V and D.
%
%
% --- Code generation:
% PCHI1D1 can be compiled into a MEX executible using codegen, as follows:
% >> codegen('pchi1d1', '-args', {zeros(L,N), zeros(M,N), zeros(M,N)}, '-o', 'pchi1d1_mex');
% To call the mex function rather than this, call pchi1d1_mex().
% A more flexible example, using an upper bound for certain dimensions of
% input, is as follows:
% >> type_U = coder.typeof(0, [L N P], [true  true true]);
% >> type_X = coder.typeof(0, [M N P], [false true true]);
% >> codegen('pchi1d1', '-args', {type_U, type_X, type_X}, '-o', 'pchi1d1_mex');
%
%
% --- Higher derivatives and integrals
% Commented-out code is included to compute higher derivatives and the
% first integral (from the lower limit of the domain of each interpolant to
% U).
%
%
% --- Acknowledgements:
% Modified from original code by Cleve Moler called pchiptx.
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 13/12/2018 - initial release

szU = size(U);
szX = size(X);
szY = size(Y);
XN = prod(szX(2:end));
YN = prod(szY(2:end));

% If U given as [N], add leading singleton dimension.
% This has been removed, as a necessity for codegen.
% Un = prod(szU);
% if (Un == XN || Un == YN) && szU(1) ~= 1
%     U = reshape(U, [1, szU]);
%     szU = size(U);
% end
UN = prod(szU(2:end));

M = szX(1);
L = szU(1);
[N,ii] = max([UN, XN, YN]);

assert(szY(1) == M, 'X and Y must have the same number of rows.');
assert(UN == 1 || UN == N, 'U must have N=%d columns or 1 column', N);
assert(XN == 1 || XN == N, 'X must have N=%d columns or 1 column', N);
assert(YN == 1 || YN == N, 'Y must have N=%d columns or 1 column', N);

szV = {szU(2:end), szX(2:end), szY(2:end)};
szV = [L, szV{ii}];
V = nan(szV);
D = nan(szV);

% PCHIP variables for each interpolant:
c = ones(1,M); % Actually c and b are size 1 by M-1, but c(M) and b(M) are
b = ones(1,M); % left as 1 always. This helps indexing when u contains +inf
d = ones(1,M);
h = ones(1,M-1);
delta = ones(1,M-1);

% Variables to index the 1'st or n'th column:
U_01 = double(UN > 1);
X_01 = double(XN > 1);
Y_01 = double(YN > 1);
jU = 1 ;
jX = 1 ;
jY = 1 ;

% Loop over columns, each of which is a pchip problem
for n = 1:N
    
    % Select data for this pchip problem
    u = U(:,jU);
    x = X(:,jX);
    y = Y(:,jY);
    
    % Handle NaN's and calculate first differences. After this loop, both x
    % and y are non-nan from m1 to m2. m1 is is the first place where x and
    % y are non-nan, and m2 is the last place after m1 where x and y are
    % non-nan.
    if isnan(x(1)) || isnan(y(1))
        m1 = 2;
    else
        m1 = 1;
    end
    m2 = M;
    for m = 2:M
        if isnan(x(m)) || isnan(y(m))
            if m > m1
                m2 = m - 1;
                break
            else
                m1 = m + 1;
            end
        end
        h(m-1) = x(m) - x(m-1);
        delta(m-1) = (y(m) - y(m-1)) / h(m-1);
    end
    
    % pchip's require at least 2 data points.
    if m2 - m1 >= 1
        
        % First derivatives:
        if m2 - m1 == 1
            %  Special case m=2, use linear interpolation.
            d(m1:m2) = delta(m1);
        else
            % Slopes for shape-preserving Hermite cubic
            % pchipslopes(h,delta) computes d(k) = P'(x(k)).
            
            %  Slopes at interior points
            %  delta = diff(y)./diff(x).
            %  d(k) = 0 if delta(k-1) and delta(k) have opposites
            %         signs or either is zero.
            %  d(k) = weighted harmonic mean of delta(k-1) and
            %         delta(k) if they have the same sign.
            for m = m1+1 : m2-1
                if sign(delta(m-1)) * sign(delta(m)) > 0
                    w1 = 2*h(m)+h(m-1);
                    w2 = h(m)+2*h(m-1);
                    d(m) = (w1+w2) ./ (w1./delta(m-1) + w2./delta(m));
                else
                    d(m) = 0;
                end
            end
            
            %  Slopes at endpoints
            d(m1) = pchipend(h(m1),h(m1+1),delta(m1),delta(m1+1));
            d(m2) = pchipend(h(m2-1),h(m2-2),delta(m2-1),delta(m2-2));
            
        end
        
        %  Piecewise polynomial coefficients:
        c(m1:m2-1) = (3*delta(m1:m2-1) - 2*d(m1:m2-1) - d(m1+1:m2)) ./ h(m1:m2-1);
        b(m1:m2-1) = (d(m1:m2-1) - 2*delta(m1:m2-1) + d(m1+1:m2)) ./ (h(m1:m2-1).*h(m1:m2-1));
        
        
        
        %  Find subinterval indices j so that x(j) <= u < x(j+1)
        [~,m] = histc(u, x(m1:m2));
        m = m + m1 - 1;
        
        for l = 1:L
            ml = m(l);
            if ml >= m1
                s = u(l) - x(ml);
                
                % Evaluate interpolant
                V(l,n) = y(ml) + s*(d(ml) + s*(c(ml) + s*b(ml)));
                
                % Evaluate derivative
                D(l,n) = d(ml) + s*(2*c(ml) + 3*s*b(ml));
                
                % Higher derivatives can be added as follows:
                % D2(l,n) = 2*c(ml) + 6*s*b(ml); % the second derivative
                % D3(l,n) = 6*b(ml); % the third derivative
            end
        end
        
        % Integrals can also be calculated, as follows for the first integral:
        %{
        % First (cumulatively) integrate each part of the piecewise function
        J = zeros(1,M-1);
        for k = m1:max(m)-1
            % J(k) is the integral of interpolant from x(k-1) to x(k).
            J(k+1) = J(k) + h(k)*(y(k) + h(k)*(d(k)/2 + h(k)*(c(k)/3 + h(k)*b(k)/4)));
        end
        
        % Then refine to the precise location, u
        for l = 1:L
            ml = m(l);
            if ml >= m1
                s = u(l) - x(ml);
                I(l,n) = J(ml) + s*(y(ml) + s*(d(ml)/2 + s*(c(ml)/3 + s*b(ml)/4)));
            end
        end
        %}
        
    end
    
    
    % Augment indices
    jU = jU + U_01;
    jX = jX + X_01;
    jY = jY + Y_01;
end

% Remove leading dimensions if L == 1, but leave row vectors as row vectors
V = squeeze(V);
D = squeeze(D);


% -------------------------------------------------------

function d = pchipend(h1,h2,del1,del2)
%  Noncentered, shape-preserving, three-point formula.
d = ((2*h1+h2)*del1 - h1*del2)/(h1+h2);
if sign(d) ~= sign(del1)
    d = 0;
elseif (sign(del1)~=sign(del2)) && (abs(d)>abs(3*del1))
    d = 3*del1;
end