function I = int1trapqn(X1,X2,X,Y) %#codegen
%INT1TRAPQN  Quick trapezoidal integration along all columns.
%
%
% I = int1trapqn(X1,X2,X,Y) 
% trapezoidally integrates each column of Y, as a function of each column
% of X, from each column of X1 to each column of X2. That is, I(n) is the
% integral from X1(n) to X2(n) of the (simplest) piecewise linear function
% with values Y(:,n) at X(:,n).
%
%
% --- Input:
% X1 [1 x N], the lower integration limit
% X2 [1 x N], the upper integration limit
% X [M x N], the data sites for the function
% Y [M x N], the values of the function at the data sites
%
% 
% --- Output:
% I [1 x N], the integral
%
% 
% --- Notes:
% X(:,n) must be non-decreasing for all n.
%
% Any input can be higher-dimensional, so long as the product of its
% dimensions excluding the first dimension is N.
%
% Any input can have a singleton second dimension (N = 1), in which case
% that single value is used for each interpolation problem.
%
% X1 and X2 need not have their leading singleton dimensions. 
%
% If I has more than 2 dimensions, it is squeeze()'d.
%
% NaN's are returned conservatively as follows (in 1 dimensional examples):
% If x1 < x(lo) where lo is the first non-NaN value of x, then v is NaN.
% If x2 > x(hi) where hi is the last  non-NaN value of x, then v is NaN.
% If y(n) is NaN and x1 <= x(n) <= x2, then v is NaN.
% If x1 == x2,
%   v is 0 if x(n-1) < x1 < x(n) and y(n-1) and y(n) are both non-NaN,
%   v is 0 if x1 == x(n) and y(n) is non-NaN,
%   v is NaN otherwise.
%
%
% --- Code generation:
% INT1TRAPQN can be compiled into a MEX executible using codegen, as
% follows:
% >> codegen('int1trapqn', '-args', {zeros(1,N), zeros(1,N), zeros(M,N), zeros(M,N)}, '-o', 'int1trapqn_mex');
% To call the mex function rather than this, call int1trapqn_mex().
% A more flexible example, using an upper bound for certain dimensions of
% input, is as follows:
% >> type_XY = coder.typeof(0, [M N P], [false true true]);
% >> type_12 = coder.typeof(0, [1 N P], [false true true]);
% >> codegen('int1trapqn', '-args', {type_12, type_12, type_XY, type_XY}, '-o', 'int1trapqn_mex')
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.1
% History   : 13/12/2018 - initial release
%           : 25/10/2019 - minor speed and documentation updates

szX1 = size(X1);
szX2 = size(X2);
szX = size(X);
szY = size(Y);
XN = prod(szX(2:end));
YN = prod(szY(2:end));
M = szX(1);
X1n = numel(X1);
X2n = numel(X2);
[N,ii] = max([X1n, X2n, XN, YN]);

assert(M == szY(1), 'X and Y must have the same number of rows.');
assert(X1n == 1 || X1n == N, 'X1 must have 1 element or N=%d element.', N);
assert(X2n == 1 || X2n == N, 'X2 must have 1 element or N=%d element.', N);
assert(XN == 1  || XN == N , 'X must have 1 column or N=%d columns.', N);
assert(YN == 1  || YN == N , 'Y must have 1 column or N=%d columns.', N);

szV = {szX1, szX2, szX(2:end), szY(2:end)};
szV = [1, szV{ii}];
I = zeros(szV);

% Variables to index the 1'st or n'th column:
X1_01 = double(X1n > 1);
X2_01 = double(X2n > 1);
X_01 = double(XN > 1); % double is faster than using uint32
Y_01 = double(YN > 1);
jX1 = 1 ;
jX2 = 1 ;
jX = 1 ;
jY = 1 ;

% Loop over columns, each of which is an integral problem
for j = 1:N
    
    x1 = X1(jX1);
    x2 = X2(jX2);
    x = X(:,jX);
    y = Y(:,jY);
    
    % Integral is NaN if either limit is NaN, or if either limit
    % is (definitely) outside the domain of integration.
    if isnan(x2) || isnan(x1) || x1 < x(1) || x2 > x(M)
        I(j) = nan;
        
    elseif x1 == x2
        % Special case where x1 == x2.
        % Return 0 if x(m-1) < x1 < x(m) and y(m-1) and y(m) are both non-NaN
        % Return 0 if x(m) == x1 and y(m) is non-NaN
        % Return NaN otherwise
        m = find(x1 <= x, 1, 'first');
        if isempty(m)
            I(j) = nan;
        else
            m = m(1); % To keep codegen happy. m was a scalar anyway.
            if x1 == x(m) && ~isnan(y(m))
                I(j) = 0;
            elseif ~isnan(x(m-1)) && ~isnan(y(m-1)) && ~isnan(y(m))
                I(j) = 0;
            else
                I(j) = NaN;
            end
        end
        
    else
        
        % Check whether integration limits are reversed.
        if x2 < x1
            tmp = x2;
            x2 = x1;
            x1 = tmp;
            flip = true;
        else
            flip = false;
        end
        
        % Search for m1, smallest index such that x1 < x(m1), and
        % pre-emptively undo (subtract) the integral from x(m-1) to x1.
        m1 = 0; % Just a special flag
        for m = 2 : M
            if x1 < x(m)
                y1 = interp1(x(m), x(m-1), y(m), y(m-1), x1);
                I(j) = -(x1 - x(m-1)) * (y(m-1) + y1) / 2;
                m1 = m;
                break
            end
        end
        
        if m1 == 0
            % x1 >= the biggest (non-NaN) x value.
            % Result can only be non-NaN if x1 == biggest (non-NaN) x value,
            % and x2 == x1, as covered above
            I(j) = nan;
            
        else
            
            % Integrate from x(m-1) to x(m), for all m from m1 to the first
            % index m where x2 <= x(m), then undo (subtract) the
            % integration from x2 to x(m)
            for m = m1 : M
                I(j) = I(j) + (x(m) - x(m-1)) * (y(m) + y(m-1)) / 2;
                if x2 <= x(m)
                    y2 = interp1(x(m-1), x(m), y(m-1), y(m), x2);
                    I(j) = I(j) - (x(m) - x2) * (y(m) + y2) / 2;
                    break
                end
            end
            
            % Negate the final answer if integration limits were reversed
            if flip
                I(j) = -I(j);
            end
            
        end
    end
    
    % Augment indices
    jX1 = jX1 + X1_01;
    jX2 = jX2 + X2_01;
    jX = jX + X_01;
    jY = jY + Y_01;
end

% Remove leading dimensions, but leave row vectors as row vectors:
I = squeeze(I);

end


function yi = interp1(x1, x2, y1, y2, xi)
% Quick linear interpolation, no input checks
u = (xi-x1) ./ (x2-x1);
yi = y1 .* (1-u) + y2 .* u;
end
