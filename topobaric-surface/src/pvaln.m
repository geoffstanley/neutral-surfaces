function y = pvaln(f,x)
%PVALN  Evaluate many polynomials, without checks for domain.
%
%
% y = pvaln(f,x)
% evaluates the polynomial(s) f at location site(s) x. f is a matrix with
% N+2 rows and L columns, with column l representing a polynomial of order
% N, valid on domain f(1,l) to f(2,l). x is a matrix of M rows and L
% columns. Then
% y(m,l) = f(3  ,l) * (x(m,l) - f(1,l))^(N-1)
%        + ...
%        + f(N+1,l) * (x(m,l) - f(1,l))
%        + f(N+2,l)
% for 1 <= m <= M and 1 <= l <= L.
% Either f or x can have a singleton second dimension, in which case that
% polynomial or value is used in conjunction with all values or
% polynomials.
% Though the domain of f is specified, this is ignored for execution speed,
% and the above formula applied regardless of the value of x(m).
%
%
% --- Input:
% f [N+2, L]: L polynomials of order N (may have L = 1)
% x [M  , L]: evaluation sites (may have L = 1)
%
%
% --- Output:
% y [M  , L]: y(m,l) is f(:,l) evaluated at x(m,l).

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% Np2-2 is the polynomial's order. (e.g. Np2==4 for linear)
Np2 = size(f,1);

% Change to local coordinates:
x = x - f(1,:);

% Evaluate polynomial:
xp = 1;
y = f(end,:);
for k = 1:Np2-3
    xp = xp .* x; % xp == x.^k.  xp is same size as x
    y = y + f(Np2-k,:) .* xp;
end