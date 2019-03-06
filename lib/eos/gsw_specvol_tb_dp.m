function v_P = gsw_specvol_tb_dp(SA,CT,p) %#ok<FNDEF> %#codegen
%GSW_SPECVOL_TB_DP  Fast derivative of the GSW specific volume w.r.t. pressure.
%
%
% v_P = gsw_specvol_tb_dp(SA,CT,p)               [ m^3 / (dbar kg) ]
% computes the partial derivative of the GSW specific volume with respect
% to pressure in dbar, given Absolute Salinity SA, Conservative Temperature
% CT, and pressure p.
%
%
% gsw_specvol_tb_dp is derived from gsw_specvol_first_derivatives.m,
% documented below. Input checks and expansion of variables have been
% removed; instead automatic expansion (requiring MATLAB 2016b or later) is
% used. This function is compatible with MATLAB's codegen.
%
% The inputs' units are as in gsw_specvol_first_derivatives, as below.
%
% The inputs' sizes are more general: they may be arrays of any dimension
% and any size, so long as their sizes match each other excepting that any
% input can have any dimension as a singleton. For example, SA and CT can
% be 3D arrays of size [nz, nx, ny], while p can be a vector of size
% [nz,1].
%
%
% Author(s)       : Geoff Stanley
% Email           : g.stanley@unsw.edu.au
% Email           : geoffstanley@gmail.com
% Version   : 1.0
% Distributed with: Topobaric Surface
%
%
% gsw_specvol_first_derivatives                     first order derivatives
%                                     of specific volume (75-term equation)
% =========================================================================
%
% USAGE:
%  [v_SA, v_CT, v_P] = gsw_specvol_first_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three first-order derivatives of specific
%  volume (v),
%   (1) v_SA, first-order derivative with respect to Absolute Salinity 
%       at constant CT & p.
%   (2) v_CT, first-order derivative with respect to CT at 
%       constant SA & p. 
%   (3) v_P, first-order derivative with respect to P at constant SA 
%       and CT. 
%
%  Note that this function uses the using the computationally-efficient
%  75-term expression for specific volume (Roquet et al., 2015).  There is 
%  an alternative to calling this function, namely 
%  gsw_specvol_first_derivatives_CT_exact(SA,CT,p) which uses the full 
%  Gibbs function (IOC et al., 2010).   
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described 
%  in McDougall et al. (2010).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  v_SA  =  The first derivative of specific volume with respect to 
%           Absolute Salinity at constant CT & p.     [ (m^3/kg)(g/kg)^-1 ]
%  v_CT  =  The first derivative of specific volume with respect to 
%           CT at constant SA and p.                         [ m^3/(K kg) ]
%  v_P   =  The first derivative of specific volume with respect to 
%           P at constant SA and CT.                        [ m^3/(Pa kg) ]
%
% AUTHOR:   
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall and P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
offset = 5.971840214030754e-1;                      % offset = deltaS*sfac.

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;
z = p.*1e-4;

c000 = -6.0799143809e-5; 
c001 =  1.99712338438e-5; 
c002 = -3.3928084311e-6; 
c003 =  4.2124612320e-7; 
c004 = -6.3236306430e-8; 
c005 =  1.1768102358e-8; 
c010 =  1.8505765429e-5; 
c011 = -2.3472773462e-6; 
c012 = -1.09581019659e-6; 
c013 =  1.25816399608e-6; 
c020 = -1.1716606853e-5; 
c021 =  4.2610057480e-6; 
c022 =  8.6087715477e-7; 
c030 =  7.9279656173e-6; 
c031 = -9.2265080074e-7; 
c040 = -3.4102187482e-6; 
c041 = -1.26705833028e-7; 
c050 =  5.0736766814e-7; 
c100 =  2.4262468747e-5; 
c101 = -1.16968865968e-6; 
c102 =  1.08930565545e-6; 
c103 = -4.4588501692e-7; 
c110 = -9.5677088156e-6; 
c111 = -1.11398309114e-5; 
c112 = -8.1887088711e-7; 
c120 = -2.3678308361e-7; 
c121 =  7.8274774160e-7; 
c130 = -3.4558773655e-6; 
c131 =  1.55237776184e-8; 
c140 =  1.2956717783e-6; 
c200 = -3.4792460974e-5; 
c201 = -9.6244503194e-6; 
c202 =  5.0238911340e-8; 
c210 =  1.1100834765e-5; 
c211 =  1.09241497668e-5; 
c220 =  2.9283346295e-6; 
c221 = -1.31462208134e-6; 
c230 =  3.1655306078e-7; 
c300 =  3.7470777305e-5; 
c301 =  9.8526213996e-6; 
c310 = -9.8447117844e-6; 
c311 = -2.7088371254e-6; 
c320 = -4.8826139200e-7; 
c400 = -1.7322218612e-5; 
c401 = -3.5623949454e-6; 
c410 =  2.5909225260e-6; 
c500 =  3.0927427253e-6; 

v_p_part = c000 + xs.*(c100 + xs.*(c200 + xs.*(c300 + xs.*(c400 + c500.*xs)))) ...
    + ys.*(c010 + xs.*(c110 + xs.*(c210 + xs.*(c310 + c410.*xs))) + ys.*(c020 ...
    + xs.*(c120 + xs.*(c220 + c320.*xs)) + ys.*(c030 + xs.*(c130 + c230*xs) ...
    + ys.*(c040 + c140.*xs + c050.*ys)))) + z.*(c001 + xs.*(c101 + xs.*(c201 ...
    + xs.*(c301 + c401.*xs))) + ys.*(c011 + xs.*(c111 + xs.*(c211 + c311.*xs)) ...
    + ys.*(c021 + xs.*(c121 + c221.*xs) + ys.*(c031 + c131.*xs + c041.*ys))) ...
    + z.*( c002 + xs.*(c102 + c202.*xs) + ys.*(c012 + c112.*xs + c022.*ys) ...
    + z.*(c003 + c103.*xs + c013.*ys + z.*(c004 + c005.*z))));

v_P = 1e-4.*v_p_part; % Original used 1e-8 for per Pa

end