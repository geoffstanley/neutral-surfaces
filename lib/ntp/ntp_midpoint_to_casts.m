function dp = ntp_midpoint_to_casts(Sppc_A, Tppc_A, P_A, k_A, Sppc_B, Tppc_B, P_B, k_B, p_A, p_B, tolp) %#codegen
%NTP_MIDPOINT_TO_CASTS  Find the Neutral Tangent Plane from a midpoint
%                       between a pair of water columns
%
%
% dp = ntp_midpoint_to_casts(Sppc_A, Tppc_A, P_A, k_A, Sppc_B, Tppc_B, P_B, k_B, p_A, p_B, tolp)
% finds the P difference, dp, that satisfy, with accuracy tolp,
%     eos(S_A(p_avg + dp / 2), T_A(p_avg + dp / 2), p_avg)
%   = eos(S_B(p_avg - dp / 2), T_B(p_avg - dp / 2), p_avg)
% where
%   S_A, T_A are the S, T profiles on cast A,
%   S_B, T_B are the S, T profiles on cast B,
%   p_avg = (p_A + p_B) / 2  is the midpoint z.
% Here, S_A and T_A have been written as functions of p, but more
% specifically they are piecewise polynomial interpolants with coefficients
% given by Sppc_A and Tppc_A with knots at P_A(1:k_A).  That is,
%   [S_A(z), T_A(z)] = ppc_val2(P_A, Sppc_A, Tppc_A, z).
% Similarly for S_B and T_B.
% The equation of state is given by eos.m on the MATLAB path.
%
%
% --- Input:
% Sppc_A, Sppc_B [O, K-1]: coefficients for piecewise polynomial for
%     practical / Absolute Salinity in terms of P, on casts A and B.  Here,
%     O is the polynomial order and K is the number of data points
%     (including NaN's) on the casts.
% Tppc_A, Tppc_B [O, K-1]: as above, but for potential / Conservative
%     Temperature.
% P_A, P_B [K, 1]: pressure or depth data points on casts A and B.
% k_A, k_B [1, 1]: number of valid (non-NaN) salinity / temperature data 
%   points on casts A and B.  Specifically, if S_A were the salinity data
%   on cast A, then S_A(1:k_A) must be non-NaN.  This translates into
%   requiring Sppc_A(:,1:k_A-1) be non-NaN.  Similarly for the others.
% tolp [1, 1]: tolerance for solving the level of neutral buoyancy (same
%     units as P).
% p_A, p_B [1, 1]: initial guess for the pressure or depth that the NTP
%     intersects casts A and B.  For example, the pressure or depth of a
%     certain potential density surface.  Though this is used as an initial
%     guess for the non-linear root finding solver, the main purpose of p_A
%     and p_B is to define the average pressure or depth, p_avg.
%
% Note: physical units for Sppc_*, Tppc_*, P_*, p_*  are determined by eos.m.
%
% Note: P_A and P_B must increase monotonically.
%
%
% --- Output:
% dp: pressure or depth difference between where NTP intersects casts A and
%     B.  Specifically, the NTP intersects cast A at p_avg + dp / 2 and
%     cast B at p_avg - dp / 2.
%
%
% --- See Also:
% ntp_bottle_to_cast
% ntp_slope
% ntp_slope_error
% ppc_linterp, ppc_pchip


% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com

if k_A > 1 && k_B > 1
  % Casts have at least two data points.  A solution is possible.
  
  p_avg = (p_A + p_B) / 2;
  dp_guess = p_A - p_B;
  
  % Upper and lower bounds of dp must satisfy
  % (a) P_A(1) <= p_avg + dp/2 <= P_A(k_A)
  % (b) P_B(1) <= p_avg - dp/2 <= P_B(k_B)
  % ... therefore ...
  %    max( P_A(1) - p_avg, p_avg - P_B(k_B) )
  % <= dp/2
  % <= min( P_A(k_A) - p_avg, p_avg - P_B(1) )
  lb = 2 * max( P_A(1) - p_avg, p_avg - P_B(k_B) );
  ub = 2 * min( P_A(k_A) - p_avg, p_avg - P_B(1) );
  
  % Search for a sign change, expanding outward from an initial guess
  [lb, ub] = fzero_guess_to_bounds(@diff_fun, dp_guess, lb, ub, ...
    p_avg, Sppc_A, Tppc_A, P_A, Sppc_B, Tppc_B, P_B);
  
  if ~isnan(lb)
    % A sign change was discovered, so a root exists in the interval.
    % Solve the nonlinear root-finding problem using Brent's method
    dp = fzero_brent(@diff_fun, lb, ub, tolp, ...
      p_avg, Sppc_A, Tppc_A, P_A, Sppc_B, Tppc_B, P_B);
  else
    dp = nan;
  end
  
else
  dp = nan;
end

end

function out = diff_fun(dp, p_avg, Sppc_A, Tppc_A, P_A, Sppc_B, Tppc_B, P_B)
% Evaluate difference between (a) eos at location on cast A with
% properties (S_A, T_A, P_A) where the pressure or depth is p_avg + dp/2,
% and (b) eos at location on cast B with properties (S_B, T_B, P_B) where
% the pressure or depth is p_avg - dp/2. Here, eos is always evaluated at
% the mid-pressure or mid-depth, p_avg.
[s_A, t_A] = ppc_val2(P_A, Sppc_A, Tppc_A, p_avg + 0.5 * dp);
[s_B, t_B] = ppc_val2(P_B, Sppc_B, Tppc_B, p_avg - 0.5 * dp);
out = eos(s_A, t_A, p_avg) - eos(s_B, t_B, p_avg);
end
