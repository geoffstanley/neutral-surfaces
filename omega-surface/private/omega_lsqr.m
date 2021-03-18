function [x,flag,relres] = omega_lsqr(A,b,tol)
%OMEGA_LSQR  A specialized version of MATLAB's lsqr.
% 
% 
% omega_lsqr is a stripped down version of MATLAB's lsqr, mostly to demand
% A is a matrix rather than a function, and demand exactly three inputs. It
% supports MATLAB's code generation, but further work is needed to make the
% compiled MEX function actually faster than this native MATLAB version;
% this will require getting BLAS or LAPACK to handle the matrix
% multiplication in the compiled MEX function, whereas native MATLAB uses
% these packages naturally.
%
% This function is for use inside omega_surface only; do not use this! Use
% MATLAB's lsqr instead.
%
% For further documentation, see MATLAB's lsqr.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


%LSQR   LSQR Method.
%   X = LSQR(A,B) attempts to solve the system of linear equations A*X=B
%   for X if A is consistent, otherwise it attempts to solve the least
%   squares solution X that minimizes norm(B-A*X). The N-by-P coefficient
%   matrix A need not be square but the right hand side column vector B
%   must have length N.
%
%   X = LSQR(AFUN,B) accepts a function handle AFUN instead of the matrix A.
%   AFUN(X,'notransp') accepts a vector input X and returns the
%   matrix-vector product A*X while AFUN(X,'transp') returns A'*X. In all
%   of the following syntaxes, you can replace A by AFUN.
%
%   X = LSQR(A,B,TOL) specifies the tolerance of the method. If TOL is []
%   then LSQR uses the default, 1e-6.
%
%   X = LSQR(A,B,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then LSQR uses the default, min([N,P,20]).
%
%   X = LSQR(A,B,TOL,MAXIT,M) and LSQR(A,B,TOL,MAXIT,M1,M2) use P-by-P
%   preconditioner M or M = M1*M2 and effectively solve the system
%   A*inv(M)*Y = B for Y, where Y = M*X. If M is [] then a preconditioner
%   is not applied. M may be a function handle MFUN such that
%   MFUN(X,'notransp') returns M\X and MFUN(X,'transp') returns M'\X.
%
%   X = LSQR(A,B,TOL,MAXIT,M1,M2,X0) specifies the P-by-1 initial guess. If
%   X0 is [] then LSQR uses the default, an all zero vector.
%
%   [X,FLAG] = LSQR(A,B,...) also returns a convergence FLAG:
%    0 LSQR converged to the desired tolerance TOL within MAXIT iterations.
%    1 LSQR iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 LSQR stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during LSQR became too
%      small or too large to continue computing.
%
%   [X,FLAG,RELRES] = LSQR(A,B,...) also returns estimates of the relative
%   residual NORM(B-A*X)/NORM(B). If RELRES <= TOL, then X is a
%   consistent solution to A*X=B. If FLAG is 0 but RELRES > TOL, then X is
%   the least squares solution which minimizes norm(B-A*X).
%
%   [X,FLAG,RELRES,ITER] = LSQR(A,B,...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = LSQR(A,B,...) also returns a vector of
%   estimates of the residual norm at each iteration including NORM(B-A*X0).
%
%   [X,FLAG,RELRES,ITER,RESVEC,LSVEC] = LSQR(A,B,...) also returns a vector
%   of least squares estimates at each iteration:
%   NORM((A*inv(M))'*(B-A*X))/NORM(A*inv(M),'fro'). Note the estimate of
%   NORM(A*inv(M),'fro') changes, and hopefully improves, at each iteration.
%
%   Example:
%      n = 100; on = ones(n,1); A = spdiags([-2*on 4*on -on],-1:1,n,n);
%      b = sum(A,2); tol = 1e-8; maxit = 15;
%      M1 = spdiags([on/(-2) on],-1:0,n,n);
%      M2 = spdiags([4*on -on],0:1,n,n);
%      x = lsqr(A,b,tol,maxit,M1,M2);
%   Or, use this matrix-vector product function
%      %-----------------------------------%
%      function y = afun(x,n,transp_flag)
%      if strcmp(transp_flag,'transp')
%         y = 4 * x;
%         y(1:n-1) = y(1:n-1) - 2 * x(2:n);
%         y(2:n) = y(2:n) - x(1:n-1);
%      elseif strcmp(transp_flag,'notransp')
%         y = 4 * x;
%         y(2:n) = y(2:n) - 2 * x(1:n-1);
%         y(1:n-1) = y(1:n-1) - x(2:n);
%      end
%      %-----------------------------------%
%   as input to LSQR:
%      x1 = lsqr(@(x,tflag)afun(x,n,tflag),b,tol,maxit,M1,M2);
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, MINRES, PCG, QMR,
%   SYMMLQ, TFQMR, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2017 The MathWorks, Inc.

n = size(A, 2);

if tol <= eps
    warning(message('MATLAB:lsqr:tooSmallTolerance'));
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:lsqr:tooBigTolerance'));
    tol = 1-eps;
end

maxit = n; % maximum number of iterations


% Set up for the method
n2b = norm(b);                     % Norm of rhs vector, b
flag = 1;
tolb = tol * n2b;                  % Relative tolerance
u = b;

beta = n2b;

% Norm of residual r=b-A*x is estimated well by prod_i abs(sin_i)
normr = beta;
if beta ~= 0
    u = u / beta;
end
c = 1;
s = 0;
phibar = beta;
v = A' * u;
 

% Initialize x
x = zeros(n,1);

alpha = norm(v);
if alpha ~= 0
    v = v / alpha;
end
d = zeros(n,1);

% norm((A*inv(M))'*r) = alpha_i * abs(sin_i * phi_i)
normar = alpha * beta;

% Check for exact solution
if (normar == 0)               % if alpha_1 == 0 | beta_1 == 0
    flag = 0;                  % a valid solution has been obtained
    relres = normr/n2b;
    return
end

% Poorly estimate norm(A*inv(M),'fro') by norm(B_{ii+1,ii},'fro')
% which is in turn estimated very well by
% sqrt(sum_i (alpha_i^2 + beta_{ii+1}^2))
norma = 0;
% norm(inv(A*inv(M)),'fro') = norm(D,'fro')
% which is poorly estimated by sqrt(sum_i norm(d_i)^2)
sumnormd2 = 0;
stag = 0;                      % stagnation of the method
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)
for ii = 1 : maxit
    
    z = v;
    u = A * z - alpha * u;
    beta = norm(u);
    u = u / beta;
    norma = norm([norma alpha beta]);
    thet = - s * alpha;
    rhot = c * alpha;
    rho = sqrt(rhot^2 + beta^2);
    c = rhot / rho;
    s = - beta / rho;
    phi = c * phibar;
    if (phi == 0)              % stagnation of the method
        stag = 1;
    end
    phibar = s * phibar;
    d = (z - thet * d) / rho;
    sumnormd2 = sumnormd2 + (norm(d))^2;
    
    % Check for stagnation of the method
    if abs(phi)*norm(d) < eps*norm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
     % check for convergence in min{|b-A*x|} || in A*x=b
    if normar/(norma*normr) <= tol || normr <= tolb
        flag = 0;
        break
    end
    
    if stag >= maxstagsteps
        flag = 3;
        break
    end
    
    x = x + phi * d;
    normr = abs(s) * normr;
    vt = A' * u;
    
    v = vt - beta * v;
    alpha = norm(v);
    v = v / alpha;
    normar = alpha * abs( s * phi);
    
end                            % for ii = 1 : maxit


%      check for convergence in min{|b-A*x|} || in A*x=b
if flag == 1 && (normar/(norma*normr) <= tol || normr <= tolb) 
    flag = 0;
end

relres = normr/n2b;
