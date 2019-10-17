function rho_z = eoscg_densmjwf2003_d2001_dz(s, t, z) %#ok<FNDEF> %#codegen
%EOSCG_DENSMJWF2003_D2001_DZ  The partial derivative, with respect to
% depth, of the McDougall, Jackett, Wright, and Feistel (2003) equation of
% state, with the Dukowicz (2001) Boussinesq correction.
%
%
% rho_z = eoscg_densmjwf2003_d2001_dz(s,t,z)                   [ kg / m^4 ]
% computes the partial derivative, with respect to depth, of the D01
% Boussinesq version of the MJWF03 in-situ density given practical salinity
% s, potential temperature t, and depth z.
%
%
% --- Input:
% s: Practical Salinity [PSU]
% t: Potential Temperature with reference pressure of 0 [ITS-90]
% z: Depth [m, positive and increasing downwards, like pressure]
%
%
% --- Output: 
% rho_z: partial derivative of in-situ density with respect to depth [kg / m^4]
%
%
% --- Ackowledgements:
% This MATLAB code has been adapted from the FORTRAN function state_mod.F90
% in the Parallel Ocean Program (POP) model.
%
%
% --- References
% Dukowicz, J.K., 2001. Reduction of Density and Pressure Gradient Errors
% in Ocean Simulations. Journal of Physical Oceanography 31, 1915–1921.
% https://doi.org/10.1175/1520-0485(2001)031<1915:RODAPG>2.0.CO;2
%
% McDougall, t.J., Jackett, D.R., Wright, D.G., Feistel, R., 2003. Accurate
% and Computationally Efficient Algorithms for Potential Temperature and
% Density of Seawater. Journal of Atmospheric and Oceanic Technology 20,
% 730–741. https://doi.org/10.1175/1520-0426(2003)20<730:AACEAF>2.0.CO;2
%
%
% Author(s)       : Geoff Stanley
% Email           : g.stanley@unsw.edu.au
% Email           : geoffstanley@gmail.com
% Version         : 1.0



% The FORTRAN code state_mod.F90 reports the following test value: 
% rho = 1.033213242 for s = 35.0 PSU, theta = 20.0, pressz = 200.0
% However, this MATLAB function has a slightly different test value:
% rho = 1.033213387 for s = 35.0 PSU, theta = 20.0, pressz = 200.0
% Maybe a different result of sqrt() ...?

% The following FORTRAN code from state_mod.F90  is missing from this MATLAB function.
% What are the values for tmin, tmax, smin, and smax?
% case (state_range_enforce)
% TQ = min(TEMPK,tmax(kk))
% TQ = max(TQ,tmin(kk))
% SQ = min(SALTK,smax(kk))
% SQ = max(SQ,smin(kk))

% Get pressure from the Dukowicz (2001) reference pressure as a function of depth
%p = 0.059808 * (exp(-0.025 * z) - 1) + (0.100766 + 2.28405e-7 * z) .* z;  % [p] = bar
p = 0.59808 * (exp(-0.025 * z) - 1) + (1.00766 + 2.28405e-6 * z) .* z;  % [p] = dbar

% Dukowicz (2001) scaling factor
%r = 1.02819 - 2.93161e-4 * exp(-0.05 * p) + 4.4004e-5 * p; % when [p] = bar
r = 1.02819 - 2.93161e-4 * exp(-0.005 * p) + 4.4004e-6 * p; % when [p] = dbar

%dpdz = (-0.025 * 0.059808) * exp(-0.025 * z) + (0.100766 + 2 * 2.28405e-7 * z) ; % when [p] = bar, [z] = m
dpdz = (-0.025 * 0.59808) * exp(-0.025 * z) + (1.00766 + 2 * 2.28405e-6 * z) ; % when [p] = dbar, [z] = m

%drdp = (0.05 * 2.93161e-4) * exp(-0.05 * p) + 4.4004e-5; % when [p] = bar, [z] = m
drdp = (0.005 * 2.93161e-4) * exp(-0.005 * p) + 4.4004e-6; % when [p] = dbar, [z] = m

% these constants will be used to construct the numerator
mwjfnp0s0t0 =   9.99843699e+2 ;
mwjfnp0s0t1 =   7.35212840e+0 ;
mwjfnp0s0t2 =  -5.45928211e-2 ;
mwjfnp0s0t3 =   3.98476704e-4 ;
mwjfnp0s1t0 =   2.96938239e+0 ;
mwjfnp0s1t1 =  -7.23268813e-3 ;
mwjfnp0s2t0 =   2.12382341e-3 ;
mwjfnp1s0t0 =   1.04004591e-2 ;
mwjfnp1s0t2 =   1.03970529e-7 ;
mwjfnp1s1t0 =   5.18761880e-6 ;
mwjfnp2s0t0 =  -3.24041825e-8 ;
mwjfnp2s0t2 =  -1.23869360e-11;

% these constants will be used to construct the denominator
mwjfdp0s0t0 =   1.0e+0        ;
mwjfdp0s0t1 =   7.28606739e-3 ;
mwjfdp0s0t2 =  -4.60835542e-5 ;
mwjfdp0s0t3 =   3.68390573e-7 ;
mwjfdp0s0t4 =   1.80809186e-10;
mwjfdp0s1t0 =   2.14691708e-3 ;
mwjfdp0s1t1 =  -9.27062484e-6 ;
mwjfdp0s1t3 =  -1.78343643e-10;
mwjfdp0sqt0 =   4.76534122e-6 ;
mwjfdp0sqt2 =   1.63410736e-9 ;
mwjfdp1s0t0 =   5.30848875e-6 ;
mwjfdp2s0t3 =  -3.03175128e-16;
mwjfdp3s0t1 =  -1.27934137e-17;


bad = isnan(s);
s = max(min(s, 1000), 0);
s(bad) = nan;
t = max(min(t, 1000), -1000);

S1o2 = sqrt(s);
T2 = t .* t;

%***
%*** first calculate numerator of MWJF density [P_1(s,t,p)]
%***

mwjfnums0t0 = mwjfnp0s0t0 + p .* (mwjfnp1s0t0 + p * mwjfnp2s0t0);
mwjfnums0t1 = mwjfnp0s0t1 ;
mwjfnums0t2 = mwjfnp0s0t2 + p .* (mwjfnp1s0t2 + p * mwjfnp2s0t2);
mwjfnums0t3 = mwjfnp0s0t3;
mwjfnums1t0 = mwjfnp0s1t0 + p * mwjfnp1s1t0;
mwjfnums1t1 = mwjfnp0s1t1;
mwjfnums2t0 = mwjfnp0s2t0;

rho = mwjfnums0t0 + t .* (mwjfnums0t1 + t .* (mwjfnums0t2 + ...
    mwjfnums0t3 * t)) + s .* (mwjfnums1t0 +              ...
    mwjfnums1t1 * t + mwjfnums2t0 * s);

%***
%*** now calculate denominator of MWJF density [P_2(s,t,p)]
%***

mwjfdens0t0 = mwjfdp0s0t0 + p * mwjfdp1s0t0;
mwjfdens0t1 = mwjfdp0s0t1 + p.^3 * mwjfdp3s0t1;
mwjfdens0t2 = mwjfdp0s0t2;
mwjfdens0t3 = mwjfdp0s0t3 + p.*p * mwjfdp2s0t3;
mwjfdens0t4 = mwjfdp0s0t4;
mwjfdens1t0 = mwjfdp0s1t0;
mwjfdens1t1 = mwjfdp0s1t1;
mwjfdens1t3 = mwjfdp0s1t3;
mwjfdensqt0 = mwjfdp0sqt0;
mwjfdensqt2 = mwjfdp0sqt2;

DENOMK =  mwjfdens0t0 ...
    + t .* (mwjfdens0t1 + t .* (mwjfdens0t2 + t .* (mwjfdens0t3 + mwjfdens0t4 * t))) ...
    + s .* (mwjfdens1t0 + t .* (mwjfdens1t1 + T2 * mwjfdens1t3) ...
    + S1o2 .* (mwjfdensqt0 + T2 * mwjfdensqt2)) ...
    ;

rho = rho ./ DENOMK; % the in-situ density, without Dukowicz correction. 

DWORK1DP = (mwjfnp1s0t0 + 2 * mwjfnp2s0t0 * p) ...
    + T2 .* (mwjfnp1s0t2 + 2 * mwjfnp2s0t2 * p) + ...
    + s .* mwjfnp1s1t0;

DWORK2DP = mwjfdp1s0t0 + ...
    t .* p .* (3*mwjfdp3s0t1 * p + (2*mwjfdp2s0t3) .* T2);
rho_z = (DWORK1DP - rho .* DWORK2DP) ./ DENOMK; % This is d(rho)/dp before the Dukowicz correction

rho_z = (rho_z - rho ./ r .* drdp) .* dpdz ./ r; 



