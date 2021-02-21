function rho = densjmd95(S,T,P) %#codegen
%DENSJMD95  Fast JMD95 in-situ density.
%
%
% rho = densjmd95(S,T,P)                                       [ kg / m^3 ]
% computes the JMD95 in-situ density given practical salinity S, potential
% temperature T, and pressure P.
%
%
% This function is derived from densjmd95.m, documented below. Input checks
% and expansion of variables have been removed; instead automatic expansion
% (requiring MATLAB 2016b or later) is used. The calculation has also been
% streamlined by pre-allocating arrays for coeffcients, and modifying the
% coefficients such that pressure need not be internally converted from
% dbar to bar. This function is compatible with MATLAB's codegen.
%
% The input's and output's units are as in densjmd95, documented below.
%
% The input's sizes are more general: they may be arrays of any dimension
% and any size, so long as their sizes match each other excepting that any
% input can have any dimension as a singleton. For example, s and t can be
% 3D arrays of size [nz, nx, ny], while p can be a vector of size [nz,1].
%
%
% Author(s)       : Geoff Stanley
% Email           : g.stanley@unsw.edu.au
% Email           : geoffstanley@gmail.com
% Version   : 1.0
% Distributed with: Topobaric Surface
%
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


% Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388

% created by mlosch on 2002-08-09
% $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2 2007/02/17 23:49:43 jmc Exp $
% $Name:  $
%

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
eosJMDCSw(3) =    7.643800e-05;
eosJMDCSw(4) = -  8.246700e-07;
eosJMDCSw(5) =    5.387500e-09;
eosJMDCSw(6) = -  5.724660e-03;
eosJMDCSw(7) =    1.022700e-04;
eosJMDCSw(8) = -  1.654600e-06;
eosJMDCSw(9) =    4.831400e-04;

% coefficients in pressure coordinates for
% 3. secant bulk modulus K of fresh water at p = 0
eosJMDCKFw = zeros(5,1);
eosJMDCKFw(1) =   1.965933e+05; % == original * 10
eosJMDCKFw(2) =   1.444304e+03; % == original * 10
eosJMDCKFw(3) = - 1.706103e+01; % == original * 10
eosJMDCKFw(4) =   9.648704e-02; % == original * 10
eosJMDCKFw(5) = - 4.190253e-04; % == original * 10
% 4. secant bulk modulus K of sea water at p = 0
eosJMDCKSw = zeros(7,1);
eosJMDCKSw(1) =   5.284855e+02; % == original * 10
eosJMDCKSw(2) = - 3.101089e+00; % == original * 10
eosJMDCKSw(3) =   6.283263e-02; % == original * 10
eosJMDCKSw(4) = - 5.084188e-04; % == original * 10
eosJMDCKSw(5) =   3.886640e+00; % == original * 10
eosJMDCKSw(6) =   9.085835e-02; % == original * 10
eosJMDCKSw(7) = - 4.619924e-03; % == original * 10
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
eosJMDCKP( 9) =   2.102898e-05; % == original / 10
eosJMDCKP(10) = - 1.202016e-06; % == original / 10
eosJMDCKP(11) =   1.394680e-08; % == original / 10
eosJMDCKP(12) = - 2.040237e-07; % == original / 10
eosJMDCKP(13) =   6.128773e-09; % == original / 10
eosJMDCKP(14) =   6.207323e-11; % == original / 10

s1o2 = sqrt(S);

K =              eosJMDCKFw(1) + T.*(eosJMDCKFw(2) + T.*(eosJMDCKFw(3) + T.*(eosJMDCKFw(4) + T.*eosJMDCKFw(5)))) ...
    +   S .*    (eosJMDCKSw(1) + T.*(eosJMDCKSw(2) + T.*(eosJMDCKSw(3) + T.* eosJMDCKSw(4))) ...
    +   s1o2 .* (eosJMDCKSw(5) + T.*(eosJMDCKSw(6) + T.* eosJMDCKSw(7)) )) ...
    + P .* (     eosJMDCKP(1 ) + T.*(eosJMDCKP(2 ) + T.*(eosJMDCKP(3) + T.*eosJMDCKP(4))) ...
    +   S .*    (eosJMDCKP(5 ) + T.*(eosJMDCKP(6 ) + T.* eosJMDCKP(7)) + s1o2.*eosJMDCKP(8)) ...
    + P .* (     eosJMDCKP(9 ) + T.*(eosJMDCKP(10) + T.* eosJMDCKP(11))  ...
    +   S .*    (eosJMDCKP(12) + T.*(eosJMDCKP(13) + T.* eosJMDCKP(14))) ));


rho = ( ...
    +           eosJMDCFw(1) + T.*(eosJMDCFw(2) + T.*(eosJMDCFw(3) + T.*(eosJMDCFw(4) + T.*(eosJMDCFw(5) + T*eosJMDCFw(6))))) ...
    + S .* ( ...
    +           eosJMDCSw(1) + T.*(eosJMDCSw(2) + T.*(eosJMDCSw(3) + T.*(eosJMDCSw(4) + T*eosJMDCSw(5)))) ...
    + s1o2 .* ( eosJMDCSw(6) + T.*(eosJMDCSw(7) + T*eosJMDCSw(8))) ...
    + S    .*   eosJMDCSw(9) ...
    )) ./ (1 - P ./ K); 

end