function [rho, rho_s, rho_t, rho_p] = densjmd95_d0ds1dt1dp1(s,t,p,getRHO)
% Compute in-situ density r, and its first partial derivatives
% Based on the densjmd95 function, documented below.
% Geoff Stanley, geoff.stanley@physics.ox.ac.uk.

%
% DENSJMD95    Density of sea water
%=========================================================================
%
% USAGE:  dens = densjmd95(S,Theta,P)
%
% DESCRIPTION:
%    Density of Sea Water using Jackett and McDougall 1995 (JAOT 12)
%    polynomial (modified UNESCO polynomial).
%
% INPUT:  (all must have same dimensions)
%   S     = salinity    [psu      (PSS-78)]
%   Theta = potential temperature [degree C (IPTS-68)]
%   P     = pressure    [dbar]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   dens = density  [kg/m^3]
%
% AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu)
%
% check value
% S     = 35.5 PSU
% Theta = 3 degC
% P     = 3000 dbar
% rho   = 1041.83267 kg/m^3


% Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388

% created by mlosch on 2002-08-09
% $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2 2007/02/17 23:49:43 jmc Exp $
% $Name:  $

narginchk(2,4);
if nargin < 3 || isempty(p)
    % Assume s and t are paired in a cell array.
    p = t;
    t = s{2};
    s = s{1};
end

% default is to get in-situ density, not specific volume.
getRHO = (nargin < 4 || isempty(getRHO));


DP1 = (nargout >= 4);

% convert pressure to bar
p = .1*p;

% coefficients nonlinear equation of state in pressure coordinates for
% 1. density of fresh water at p = 0
eosJMDCFw = zeros(6,1);
eosJMDCFw(1) =  999.842594;
eosJMDCFw(2) =    6.793952e-02;
eosJMDCFw(3) = -  9.095290e-03;
eosJMDCFw(4) =    1.001685e-04;
eosJMDCFw(5) = -  1.120083e-06;
eosJMDCFw(6) =    6.536332e-09;
% 2. density of sea water at p = 0
eosJMDCSw = zeros(9,1);
eosJMDCSw(1) =    8.244930e-01;
eosJMDCSw(2) = -  4.089900e-03;
eosJMDCSw(3) =    7.643800e-05 ;
eosJMDCSw(4) = -  8.246700e-07;
eosJMDCSw(5) =    5.387500e-09;
eosJMDCSw(6) = -  5.724660e-03;
eosJMDCSw(7) =    1.022700e-04;
eosJMDCSw(8) = -  1.654600e-06;
eosJMDCSw(9) =    4.831400e-04;

% coefficients in pressure coordinates for
% 3. secant bulk modulus K of fresh water at p = 0
eosJMDCKFw = zeros(5,1);
eosJMDCKFw(1) =   1.965933e+04;
eosJMDCKFw(2) =   1.444304e+02;
eosJMDCKFw(3) = - 1.706103e+00;
eosJMDCKFw(4) =   9.648704e-03;
eosJMDCKFw(5) = - 4.190253e-05;
% 4. secant bulk modulus K of sea water at p = 0
eosJMDCKSw = zeros(7,1);
eosJMDCKSw(1) =   5.284855e+01;
eosJMDCKSw(2) = - 3.101089e-01;
eosJMDCKSw(3) =   6.283263e-03;
eosJMDCKSw(4) = - 5.084188e-05;
eosJMDCKSw(5) =   3.886640e-01;
eosJMDCKSw(6) =   9.085835e-03;
eosJMDCKSw(7) = - 4.619924e-04;
% 5. secant bulk modulus K of sea water at p
eosJMDCKP = zeros(14,1);
eosJMDCKP( 1) =   3.186519e+00;
eosJMDCKP( 2) =   2.212276e-02;
eosJMDCKP( 3) = - 2.984642e-04;
eosJMDCKP( 4) =   1.956415e-06;
eosJMDCKP( 5) =   6.704388e-03;
eosJMDCKP( 6) = - 1.847318e-04;
eosJMDCKP( 7) =   2.059331e-07;
eosJMDCKP( 8) =   1.480266e-04;
eosJMDCKP( 9) =   2.102898e-04;
eosJMDCKP(10) = - 1.202016e-05;
eosJMDCKP(11) =   1.394680e-07;
eosJMDCKP(12) = - 2.040237e-06;
eosJMDCKP(13) =   6.128773e-08;
eosJMDCKP(14) =   6.207323e-10;

t2 = t.*t;
t3 = t2.*t;
t4 = t3.*t;

s1o2 = sqrt(s);


r = eosJMDCFw(1) + eosJMDCFw(2)*t + eosJMDCFw(3)*t2 + eosJMDCFw(4)*t3 + eosJMDCFw(5)*t4 + eosJMDCFw(6)*t4.*t ;

tmp = eosJMDCSw(1) + eosJMDCSw(2)*t + eosJMDCSw(3)*t2 + eosJMDCSw(4)*t3 + eosJMDCSw(5)*t4;
r = r + s .* tmp;
r_s = tmp;

tmp = s1o2 .* (eosJMDCSw(6) + eosJMDCSw(7)*t + eosJMDCSw(8)*t2);
r = r + s .* tmp;
r_s = r_s + 3/2 * tmp;

tmp = s .* eosJMDCSw(9);
r = r + s .* tmp;
r_s = r_s + 2*tmp;


b = eosJMDCKFw(1) + eosJMDCKFw(2)*t + eosJMDCKFw(3)*t2 + eosJMDCKFw(4)*t3 + eosJMDCKFw(5)*t4 ;

b_s = eosJMDCKSw(1) + eosJMDCKSw(2)*t + eosJMDCKSw(3)*t2 + eosJMDCKSw(4)*t3;
b = b + s .* b_s;

tmp = s1o2 .* (eosJMDCKSw(5) + eosJMDCKSw(6)*t + eosJMDCKSw(7)*t2) ;
b = b + s .* tmp;
b_s = b_s + 3/2 .* tmp;

if DP1
    b_p = eosJMDCKP(1) + eosJMDCKP(2)*t + eosJMDCKP(3)*t2 + eosJMDCKP(4)*t3 ;
    b = b + p .* b_p;
else
    b = b + p .* (eosJMDCKP(1) + eosJMDCKP(2)*t + eosJMDCKP(3)*t2 + eosJMDCKP(4)*t3);
end

tmp = eosJMDCKP(5) + eosJMDCKP(6)*t + eosJMDCKP(7)*t2 ;
b_s = b_s + p .* (tmp + 3/2 * eosJMDCKP(8) * s1o2);

if DP1
    tmp = s .* (tmp + eosJMDCKP(8) * s1o2);
    b_p = b_p + tmp;
    b = b + p .* tmp;
else
    b = b + p .* s .* (tmp + eosJMDCKP(8) * s1o2);
end

tmp = p .* (eosJMDCKP(12) + eosJMDCKP(13)*t + eosJMDCKP(14)*t2) ;
b_s = b_s + p .* tmp;
if DP1
    tmp = s .* tmp + p .* (eosJMDCKP(9) + eosJMDCKP(10)*t + eosJMDCKP(11)*t2);
    b = b + p .* tmp;
    b_p = b_p + 2 * tmp;
else
    b = b + p .* (s .* tmp + p .* (eosJMDCKP(9) + eosJMDCKP(10)*t + eosJMDCKP(11)*t2));
end
clear tmp

r_t = ...
    eosJMDCFw(2) + 2*eosJMDCFw(3)*t + 3*eosJMDCFw(4)*t2 + 4*eosJMDCFw(5)*t3 + 5*eosJMDCFw(6)*t4 ...
    + s .* ( ...
        eosJMDCSw(2) + 2*eosJMDCSw(3)*t + 3*eosJMDCSw(4)*t2 + 4*eosJMDCSw(5)*t3 ...
        + s1o2 .* ( eosJMDCSw(7) + 2*eosJMDCSw(8)*t) ...
    );
clear t4

b_t = ...
    eosJMDCKFw(2) + 2*eosJMDCKFw(3)*t + 3*eosJMDCKFw(4)*t2 + 4*eosJMDCKFw(5)*t3 ...
    + s .* ( ...
        eosJMDCKSw(2) + 2*eosJMDCKSw(3)*t + 3*eosJMDCKSw(4)*t2 ...
        + s1o2 .* ( eosJMDCKSw(6) + 2*eosJMDCKSw(7)*t ) ...
    ) ...
    + p .* ( ...
        eosJMDCKP(2) + 2*eosJMDCKP(3)*t + 3*eosJMDCKP(4)*t2 ...
        + p .* ( eosJMDCKP(10) + 2*eosJMDCKP(11)*t ) ...
        + s .* ( ...
            eosJMDCKP(6) + 2*eosJMDCKP(7)*t ...
            + p .* ( eosJMDCKP(13) + 2*eosJMDCKP(14)*t ) ...
        ) ...
    );
clear t3 t2 s1o2


if getRHO
    
    iBminusP = 1 ./ (b - p);
    rho = r .* b .* iBminusP;
    rho_s = ((r_s .* b + r .* b_s) - rho .*  b_s     ) .* iBminusP;
    rho_t = ((r_t .* b + r .* b_t) - rho .*  b_t     ) .* iBminusP;
    if DP1
        rho_p = ((           r .* b_p) - rho .* (b_p - 1)) .* iBminusP * 1e-1;
    end
    
else % get Specific Volume

    iRB = 1 ./ (r .* b);
    rho = (b - p) .* iRB;
    rho_s = ( b_s     - rho .* (r_s .* b + r .* b_s) ) .* iRB;
    rho_t = ( b_t     - rho .* (r_t .* b + r .* b_t) ) .* iRB;
    if DP1
        rho_p = ( b_p - 1 - rho .* (         + r .* b_p) ) .* iRB * 1e-1;
    end
    
end


