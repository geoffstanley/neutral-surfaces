function [y,z] = pchipqn2(x,X,Y,Z) %#codegen
%PCHIPQN2  Piecewise Cubic Hermite Interpolation along all columns, twice.
%
%
% [y,z] = pchipqn2(x,X,Y,Z)
% uses a PCHIP to interpolate each column of Y and each column of Z, as
% functions of each column of X, at each column of x. That is, y(l,m)
% interpolates Y(:,m), as a function of X(:,m), at U(l,m), and similarly
% for z. No extrapolation is done: when x(l,n) is out of bounds of X(:,n),
% y(l,n) and z(l,n) are NaN.
%
%
% --- Input:
% x [L x N], the interpolant evaluation sites
% X [M x N], the data sites for the functions
% Y [M x N], the values of the first function at the data sites
% Z [M x N], the values of the second function at the data sites
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
% column).
%
% Any input can be higher-dimensional, so long as the product of its
% dimensions excluding the first dimension is N.
%
% Any input can have a singleton second dimension (N = 1), in which case
% that single value is used for each interpolation problem.
%
% Even if L == 1, x does need a leading singleton dimension.
%
% If L == 1, this dimension is squeeze()'d out of y and z.
%
%
% --- Code generation:
% PCHIPQN2 can be compiled into a MEX executible using codegen, as follows:
% >> codegen('pchipqn2', '-args', {zeros(L,N), zeros(M,N), zeros(M,N), zeros(M,N)}, '-o', 'pchipqn2_mex');
% To call the mex function rather than this, call pchipqn2_mex().
% A more flexible example, using an upper bound for certain dimensions of
% input, is as follows:
% >> type_x = coder.typeof(0, [L P Q], [false true true]);
% >> type_X = coder.typeof(0, [M P Q], [false true true]);
% >> codegen('pchipqn2', '-args', {type_x, type_X, type_X, type_X}, '-o', 'pchipqn2_mex');
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.0
% History   : 05/09/2019 - initial release

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
M = szX(1);
L = szx(1);

[N,ii] = max([xN, XN, YN]);
assert(szY(1) == M, 'X and Y must have the same number of rows.');
assert(xN == 1 || xN == N, 'x must have N=%d columns or 1 column', N);
assert(XN == 1 || XN == N, 'X must have N=%d columns or 1 column', N);
assert(YN == 1 || YN == N, 'Y must have N=%d columns or 1 column', N);

szy = {szx(2:end), szX(2:end), szY(2:end)};
szy = [L, szy{ii}];
y = nan(szy);
z = nan(szy);

% Variables to index the 1'st or n'th column:
x01 = double(xN > 1);
MX = M * double(XN > 1);
MY = M * double(YN > 1);
nX = 0 ;
nY = 0 ;

% Pre-assign sizes for PCHIP variables.
h = zeros(3,1);
deltaY = zeros(3,1);
deltaZ = zeros(3,1);
dY = zeros(2,1);
dZ = zeros(2,1);

% Loop over columns, each of which is an interpolation problem
for n = 0:N-1
    if isnan(Y(nY + 2)) || isnan(X(nX + 2))
        nX = nX + MX;
        nY = nY + MY;
        continue
    end
    for l = 1:L
        xln = x(n*L*x01 + l); % x(l,n) or x(l,1) as appropriate
        
        if isnan(xln) || xln < X(nX + 1) || xln > X(nX + M)
            
            y(n*L + l) = nan;
            z(n*L + l) = nan;
            
        else
            
            % Binary search, to find i for which X(i-1,n) < xln <= X(i,n)
            i = 2; % Result will be >= 2 always
            k = M; % Result will be <= M always
            while (i < k)
                j = floor((i + k) / 2);
                if X(nX + j) < xln % [X(j,n) or X(j,1) as appropriate]  < xln
                    i = j + 1;
                else
                    k = j;
                end
            end
            
            % Check whether xln is adjacent to the start or end of this X
            at_start = (i == 2);
            at_end = (i == M) || isnan(X(nX + i+1)) || isnan(Y(nY + i+1));
            
            if at_start && at_end
                % . X(1) < x < X(2) .   Revert to Linear Interpolation
                d = (xln - X(nX + i-1)) ./ (X(nX + i) - X(nX + i-1));
                y(n*L + l) = Y(nY+i-1) * (1-d) + Y(nY+i) * d;
                z(n*L + l) = Z(nY+i-1) * (1-d) + Z(nY+i) * d;
                
            else
                if at_start
                    % . X(1) < x <= X(2) < X(3) ...
                    h(2) = X(nX+i  ) - X(nX+i-1);
                    h(3) = X(nX+i+1) - X(nX+i  );
                    deltaY(2) = (Y(nY+i  ) - Y(nY+i-1)) / h(2);
                    deltaY(3) = (Y(nY+i+1) - Y(nY+i  )) / h(3);
                    deltaZ(2) = (Z(nY+i  ) - Z(nY+i-1)) / h(2);
                    deltaZ(3) = (Z(nY+i+1) - Z(nY+i  )) / h(3);
                    
                    %  Noncentered, shape-preserving, three-point formula:
                    dY(1) = ((2*h(2) + h(3)) * deltaY(2) - h(2) * deltaY(3)) / (h(2)+h(3));
                    if sign(dY(1)) ~= sign(deltaY(2))
                        dY(1) = 0;
                    elseif (sign(deltaY(2)) ~= sign(deltaY(3))) && (abs(dY(1)) > abs(3*deltaY(2)))
                        dY(1) = 3 * deltaY(2);
                    end
                    dZ(1) = ((2*h(2) + h(3)) * deltaZ(2) - h(2) * deltaZ(3)) / (h(2)+h(3));
                    if sign(dZ(1)) ~= sign(deltaZ(2))
                        dZ(1) = 0;
                    elseif (sign(deltaZ(2)) ~= sign(deltaZ(3))) && (abs(dZ(1)) > abs(3*deltaZ(2)))
                        dZ(1) = 3 * deltaZ(2);
                    end
                    
                    % Standard PCHIP formula
                    if sign(deltaY(2)) * sign(deltaY(3)) > 0
                        w1 = 2*h(3) +   h(2);
                        w2 =   h(3) + 2*h(2);
                        dY(2) = (w1+w2) / (w1/deltaY(2) + w2/deltaY(3));
                    else
                        dY(2) = 0;
                    end
                    if sign(deltaZ(2)) * sign(deltaZ(3)) > 0
                        w1 = 2*h(3) +   h(2);
                        w2 =   h(3) + 2*h(2);
                        dZ(2) = (w1+w2) / (w1/deltaZ(2) + w2/deltaZ(3));
                    else
                        dZ(2) = 0;
                    end
                    
                elseif at_end
                    % ... X(i-2) < X(i-1) < x <= X(i) .
                    h(1) = X(nX+i-1) - X(nX+i-2);
                    h(2) = X(nX+i  ) - X(nX+i-1);
                    deltaY(1) = (Y(nY+i-1) - Y(nY+i-2)) / h(1);
                    deltaY(2) = (Y(nY+i  ) - Y(nY+i-1)) / h(2);
                    deltaZ(1) = (Z(nY+i-1) - Z(nY+i-2)) / h(1);
                    deltaZ(2) = (Z(nY+i  ) - Z(nY+i-1)) / h(2);
                    
                    % Standard PCHIP formula
                    if sign(deltaY(1)) * sign(deltaY(2)) > 0
                        w1 = 2*h(2) +   h(1);
                        w2 =   h(2) + 2*h(1);
                        dY(1) = (w1+w2) / (w1/deltaY(1) + w2/deltaY(2));
                    else
                        dY(1) = 0;
                    end
                    if sign(deltaZ(1)) * sign(deltaZ(2)) > 0
                        w1 = 2*h(2) +   h(1);
                        w2 =   h(2) + 2*h(1);
                        dZ(1) = (w1+w2) / (w1/deltaZ(1) + w2/deltaZ(2));
                    else
                        dZ(1) = 0;
                    end
                    
                    %  Noncentered, shape-preserving, three-point formula:
                    dY(2) = ((h(1) + 2*h(2)) * deltaY(2) - h(2) * deltaY(1)) / (h(1)+h(2));
                    if sign(dY(2)) ~= sign(deltaY(2))
                        dY(2) = 0;
                    elseif (sign(deltaY(2)) ~= sign(deltaY(1))) && (abs(dY(2)) > abs(3*deltaY(2)))
                        dY(2) = 3 * deltaY(2);
                    end
                    dZ(2) = ((h(1) + 2*h(2)) * deltaZ(2) - h(2) * deltaZ(1)) / (h(1)+h(2));
                    if sign(dZ(2)) ~= sign(deltaZ(2))
                        dZ(2) = 0;
                    elseif (sign(deltaZ(2)) ~= sign(deltaZ(1))) && (abs(dZ(2)) > abs(3*deltaZ(2)))
                        dZ(2) = 3 * deltaZ(2);
                    end
                    
                else
                    % ... X(i-2) < X(i-1) < x <= X(i) < X(i+1) ...
                    h(1) = X(nX+i-1) - X(nX+i-2); % Way faster to do this
                    h(2) = X(nX+i  ) - X(nX+i-1); % than  
                    h(3) = X(nX+i+1) - X(nX+i  ); % diff(X(nx+i-2:nX+i+1))
                    deltaY(1) = (Y(nY+i-1) - Y(nY+i-2)) / h(1);
                    deltaY(2) = (Y(nY+i  ) - Y(nY+i-1)) / h(2);
                    deltaY(3) = (Y(nY+i+1) - Y(nY+i  )) / h(3);
                    deltaZ(1) = (Z(nY+i-1) - Z(nY+i-2)) / h(1);
                    deltaZ(2) = (Z(nY+i  ) - Z(nY+i-1)) / h(2);
                    deltaZ(3) = (Z(nY+i+1) - Z(nY+i  )) / h(3);
                    
                    % Standard PCHIP formula
                    for j = 1:2
                        if sign(deltaY(j)) * sign(deltaY(j+1)) > 0
                            w1 = 2*h(j+1) +   h(j);
                            w2 =   h(j+1) + 2*h(j);
                            dY(j) = (w1+w2) / (w1/deltaY(j) + w2/deltaY(j+1));
                        else
                            dY(j) = 0;
                        end
                        if sign(deltaZ(j)) * sign(deltaZ(j+1)) > 0
                            w1 = 2*h(j+1) +   h(j);
                            w2 =   h(j+1) + 2*h(j);
                            dZ(j) = (w1+w2) / (w1/deltaZ(j) + w2/deltaZ(j+1));
                        else
                            dZ(j) = 0;
                        end
                    end
                end
                
                % Polynomial coefficients for this piece
                cY = (3*deltaY(2) - 2*dY(1) - dY(2)) / h(2);
                bY = (dY(1) - 2*deltaY(2) + dY(2)) / h(2)^2;
                
                cZ = (3*deltaZ(2) - 2*dZ(1) - dZ(2)) / h(2);
                bZ = (dZ(1) - 2*deltaZ(2) + dZ(2)) / h(2)^2;
                
                % Evaluate cubic interpolant
                s = xln - X(nX + i-1);
                y(n*L + l) = Y(nY + i-1) + s * (dY(1) + s * (cY + s * bY));
                z(n*L + l) = Z(nY + i-1) + s * (dZ(1) + s * (cZ + s * bZ));
            end
        end
    end % for l
    nX = nX + MX;
    nY = nY + MY;
end

% Remove leading dimensions if L == 1, but leave row vectors as row vectors
y = squeeze(y);
z = squeeze(z);
