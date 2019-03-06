function nu = specvoljmd95_tb(s,t,p) %#ok<FNDEF> %#codegen
%SPECVOLJMD95_TB  Fast JMD95 specific volume.
%
%
% nu = specvoljmd95_tb(s,t,p)                                  [ m^3 / kg ]
% computes the JMD95 specific volume given practical salinity s, potential
% temperature t, and pressure p.
%
%
% specvoljmd95_tb is derived from densjmd95.m, documented below. Input
% checks and expansion of variables have been removed; instead automatic
% expansion (requiring MATLAB 2016b or later) is used. The calculation has
% also been streamlined by pre-allocating arrays for coeffcients, and
% modifying the coefficients such that pressure need not be internally
% converted from dbar to bar. This function is compatible with MATLAB's
% codegen.
%
% The inputs' units are as in densjmd95, documented below.
%
% The inputs' sizes are more general: they may be arrays of any dimension
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

s1o2 = sqrt(s);

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

b1 = eosJMDCKP(1) + t.*(eosJMDCKP(2) + t.*(eosJMDCKP(3) + t*eosJMDCKP(4))) ...
    + s .* ( eosJMDCKP(5) + t.*(eosJMDCKP(6) + t*eosJMDCKP(7)) + eosJMDCKP(8)*s1o2 ) ;
    
b2 = p .* ( eosJMDCKP(9) + t.*(eosJMDCKP(10) + t*eosJMDCKP(11))  ...
    + s .*  (eosJMDCKP(12) + t.*(eosJMDCKP(13) + t*eosJMDCKP(14))) ) ;
 
% b == bulkmod
b =   eosJMDCKFw(1) + t.*(eosJMDCKFw(2) + t.*(eosJMDCKFw(3) + t.*(eosJMDCKFw(4) + t*eosJMDCKFw(5)))) ... % secant bulk modulus of fresh water at the surface
    + s .* (   eosJMDCKSw(1) + t.*(eosJMDCKSw(2) + t.*(eosJMDCKSw(3) + t*eosJMDCKSw(4))) ...
        + s1o2 .* (eosJMDCKSw(5) + t.*(eosJMDCKSw(6) + t*eosJMDCKSw(7)) ) ...
    ) ... % secant bulk modulus of sea water at the surface
    + p .* ( b1 + b2 ); % secant bulk modulus of sea water at pressure p

rho =   eosJMDCFw(1) + t.*(eosJMDCFw(2) + t.*(eosJMDCFw(3) + t.*(eosJMDCFw(4) + t.*(eosJMDCFw(5) + t*eosJMDCFw(6))))) ... % density of freshwater at the surface
    + s .* ( eosJMDCSw(1) + t.*(eosJMDCSw(2) + t.*(eosJMDCSw(3) + t.*(eosJMDCSw(4) + t*eosJMDCSw(5)))) ...
        + s1o2 .* ( eosJMDCSw(6) + t.*(eosJMDCSw(7) + t*eosJMDCSw(8)) ) ...
        + s .* eosJMDCSw(9) ...
    ); % density of sea water at the surface

nu = (1 - p ./ b) ./ rho;

end