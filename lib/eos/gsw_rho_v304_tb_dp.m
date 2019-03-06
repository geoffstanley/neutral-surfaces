function drho_dp = gsw_rho_tb_dp(SA,CT,p)  %#ok<FNDEF> %#codegen
%GSW_RHO_TB_DP  Fast derivative of the GSW in-situ density w.r.t. pressure.
%
%
% drho_dp = gsw_rho_tb_dp(SA,CT,p)                      [ kg / (dbar m^3) ]
% computes the partial derivative of the GSW in-situ density with respect
% to pressure, given Absolute Salinity SA, Conservative Temperature CT, and
% pressure p.
%
%
% gsw_rho_tb_dp is derived from gsw_rho_first_derivatives.m, documented
% below. Input checks and expansion of variables have been removed; instead
% automatic expansion (requiring MATLAB 2016b or later) is used. This
% function is compatible with MATLAB's codegen.
%
% The inputs' units are as in gsw_rho_first_derivatives, documented below.
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
% gsw_rho_first_derivatives                SA, CT and p partial derivatives
%                                             of density (48-term equation)
%==========================================================================
% 
% USAGE:  
% [drho_dSA, drho_dCT, drho_dP] = gsw_rho_first_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the three (3) partial derivatives of in situ density with 
%  respect to Absolute Salinity, Conservative Temperature and pressure.  
%  Note that the pressure derivative is done with respect to pressure in 
%  Pa, not dbar.  This function uses the computationally-efficient 48-term 
%  expression for density in terms of SA, CT and p.  
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2010).  The GSW library function 
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
%  drho_dSA  =  partial derivatives of density             [ kg^2/(g m^3) ]
%                 with respect to Absolute Salinity
%  drho_dCT  =  partial derivatives of density               [ kg/(K m^3) ]
%                 with respect to Conservative Temperature
%  drho_dP   =  partial derivatives of density              [ kg/(Pa m^3) ]
%                 with respect to pressure in Pa
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.04 (10th December, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual. 
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

v01 =  9.998420897506056e+2;
v02 =  2.839940833161907;
v03 = -3.147759265588511e-2;
v04 =  1.181805545074306e-3;
v05 = -6.698001071123802;
v06 = -2.986498947203215e-2;
v07 =  2.327859407479162e-4;
v08 = -3.988822378968490e-2;
v09 =  5.095422573880500e-4;
v10 = -1.426984671633621e-5;
v11 =  1.645039373682922e-7;
v12 = -2.233269627352527e-2;
v13 = -3.436090079851880e-4;
v14 =  3.726050720345733e-6;
v15 = -1.806789763745328e-4;
v16 =  6.876837219536232e-7;
v17 = -3.087032500374211e-7;
v18 = -1.988366587925593e-8;
v19 = -1.061519070296458e-11;
v20 =  1.550932729220080e-10;
v21 =  1.0;
v22 =  2.775927747785646e-3;
v23 = -2.349607444135925e-5;
v24 =  1.119513357486743e-6;
v25 =  6.743689325042773e-10;
v26 = -7.521448093615448e-3;
v27 = -2.764306979894411e-5;
v28 =  1.262937315098546e-7;
v29 =  9.527875081696435e-10;
v30 = -1.811147201949891e-11;
v31 = -3.303308871386421e-5;
v32 =  3.801564588876298e-7;
v33 = -7.672876869259043e-9;
v34 = -4.634182341116144e-11;
v35 =  2.681097235569143e-12;
v36 =  5.419326551148740e-6;
v37 = -2.742185394906099e-5;
v38 = -3.212746477974189e-7;
v39 =  3.191413910561627e-9;
v40 = -1.931012931541776e-12;
v41 = -1.105097577149576e-7;
v42 =  6.211426728363857e-10;
v43 = -1.119011592875110e-10;
v44 = -1.941660213148725e-11;
v45 = -1.864826425365600e-14;
v46 =  1.119522344879478e-14;
v47 = -1.200507748551599e-15;
v48 =  6.057902487546866e-17;

c01 = -2.233269627352527e-2;
c02 = -3.436090079851880e-4;
c03 =  3.726050720345733e-6;
c04 = -1.806789763745328e-4;
c05 =  6.876837219536232e-7;
c06 = -6.174065000748422e-7;
c07 = -3.976733175851186e-8;
c08 = -2.123038140592916e-11;
c09 =  3.101865458440160e-10;
c10 = -2.742185394906099e-5;
c11 = -3.212746477974189e-7;
c12 =  3.191413910561627e-9;
c13 = -1.931012931541776e-12;
c14 = -1.105097577149576e-7;
c15 =  6.211426728363857e-10;
c16 = -2.238023185750219e-10;
c17 = -3.883320426297450e-11;
c18 = -3.729652850731201e-14;
c19 =  2.239044689758956e-14;
c20 = -3.601523245654798e-15;
c21 =  1.817370746264060e-16;

sqrtSA = sqrt(SA);

v_hat_denominator = v01 + CT.*(v02 + CT.*(v03 + v04*CT))  ...
             + SA.*(v05 + CT.*(v06 + v07*CT) ...
         + sqrtSA.*(v08 + CT.*(v09 + CT.*(v10 + v11*CT)))) ...
              + p.*(v12 + CT.*(v13 + v14*CT) + SA.*(v15 + v16*CT) ...
              + p.*(v17 + CT.*(v18 + v19*CT) + v20*SA));
          
v_hat_numerator = v21 + CT.*(v22 + CT.*(v23 + CT.*(v24 + v25*CT))) ...
           + SA.*(v26 + CT.*(v27 + CT.*(v28 + CT.*(v29 + v30*CT))) + v36*SA ...
       + sqrtSA.*(v31 + CT.*(v32 + CT.*(v33 + CT.*(v34 + v35*CT)))))  ...
            + p.*(v37 + CT.*(v38 + CT.*(v39 + v40*CT))  ...
           + SA.*(v41 + v42*CT) ...
            + p.*(v43 + CT.*(v44 + v45*CT + v46*SA) ...
            + p.*(v47 + v48*CT)));

dvhatden_dp =  c01 + CT.*(c02 + c03*CT) ...
    + SA.*(c04 + c05*CT) ...
    + p.*(c06 + CT.*(c07 + c08*CT) + c09*SA);

dvhatnum_dp = c10 + CT.*(c11 + CT.*(c12 + c13*CT)) ...
    + SA.*(c14 + c15*CT) ...
    + p.*(c16 + CT.*(c17 + c18*CT + c19*SA) ...
    + p.*(c20 + c21*CT));

rec_num = 1./v_hat_numerator;
       
rho = rec_num.*v_hat_denominator;

drho_dp = (dvhatden_dp - dvhatnum_dp.*rho).*rec_num;

end