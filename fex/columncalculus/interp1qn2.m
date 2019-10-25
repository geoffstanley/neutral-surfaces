function [y,z] = interp1qn2(x,X,Y,Z) %#codegen
%INTERP1QN2  Quick linear interpolation along all columns, twice.
%
%
% [y,z] = interp1qn2(x,X,Y,Z) 
% linearly interpolates each column of Y and each column of Z, as functions
% of each column of X, at each column of x. That is, y(l,n) interpolates
% Y(:,n), as a function of X(:,n), at x(l,n), and similarly for z. No
% extrapolation is done: when x(l,n) is out of bounds of X(:,n), y(l,n) and
% z(l,n) are NaN.
%
%
% --- Input:
% x [L x N], the interpolant evaluation sites
% X [K x N], the data sites for the functions
% Y [K x N], the values of the first function at the data sites
% Z [K x N], the values of the second function at the data sites
%
%
% --- Output:
% y [L x N], the interpolant values for the first function
% z [L x N], the interpolant values for the second function
%
%
% --- Notes:
% X(:,n) must be monotonically increasing for all n.
% 
% NaN's in X are treated as +Inf (and as such must come at the end of each
% row).
%
% Any input can have a singleton second dimension (N = 1), in which case
% that single value is used for each interpolation problem.
%
% Any dimension N can actually be higher-dimensional, so long as it has N
% elements.
%
% Even if L == 1, x does need a leading singleton dimension. 
%
% If L == 1, this dimension is squeeze()'d out of y and z. 

% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.1
% History   : 13/12/2018 - initial release
%           : 25/10/2019 - minor speed and documentation updates

szX = size(X);
szY = size(Y);
szZ = size(Z);
assert(all(szY == szZ), 'Y and Z must have identical dimensions.');
XN = prod(szX(2:end));
YN = prod(szY(2:end));

% If x given as [N], add leading singleton dimension.
% This has been removed, as a necessity for codegen.
szx = size(x);
% xn = numel(x);
% if szx(1) ~= 1 && xn == max(XN, YN)
%     x = reshape(x, [1, szx]);
%     szx = size(x);
% end
xN = prod(szx(2:end));
K = szX(1);
L = szx(1);

[N,ii] = max([xN, XN, YN]);
assert(szY(1) == K, 'X and Y must have the same number of rows.');
assert(szZ(1) == K, 'X and Z must have the same number of rows.');
assert(xN == 1 || xN == N, 'x must have N=%d columns or 1 column', N);
assert(XN == 1 || XN == N, 'X must have N=%d columns or 1 column', N);
assert(YN == 1 || YN == N, 'Y must have N=%d columns or 1 column', N);

szy = {szx(2:end), szX(2:end), szY(2:end)};
szy = [L, szy{ii}];
y = nan(szy);
z = nan(szy);
if K == 1
    % Can't really interpolate with one data point.
    % Exit now to prevent out-of-bounds indexing later.
    return
end

% Variables to index the 1'st or n'th column:
Y01 = double(YN > 1);
x01 = double(xN > 1);
KX = K * double(XN > 1);
nX = 0 ;

% Loop over columns, each of which is an interpolation problem
for n = 0:N-1
    for l = 1:L
        xln = x(n*L*x01 + l); % x(l,n) or x(l,1) as appropriate
        
        if ~(isnan(xln) || xln < X(nX + 1) || xln > X(nX + K))
            
            % Binary search, to find i for which X(i-1,n) < xln <= X(i,n)
            i = 2; % Result will be >= 2 always
            k = K; % Result will be <= M always
            while (i < k)
                j = floor((i + k) / 2);
                if X(nX + j) < xln % [X(j,n) or X(j,1) as appropriate]  < xln
                    i = j + 1;
                else
                    k = j;
                end
            end
            
            % X(nX + i-1)  is  X(i-1,n) or X(i-1,1) as appropriate
            % X(nX + i  )  is  X(i  ,n) or X(i  ,1) as appropriate
            d = (xln - X(nX + i-1)) ./ (X(nX + i) - X(nX + i-1));
            idx = n*K*Y01 + i;
            y(n*L + l) = Y(idx-1) * (1-d) + Y(idx) * d;
            z(n*L + l) = Z(idx-1) * (1-d) + Z(idx) * d;
            % The above line is, essentially,
            %V(l,n)    = Y(i-1,n) * (1-d) + Y(i,n) * d;
            
            
        end
    end
    nX = nX + KX;
end

% Remove leading dimensions if L == 1, but leave row vectors as row vectors
y = squeeze(y);
z = squeeze(z);