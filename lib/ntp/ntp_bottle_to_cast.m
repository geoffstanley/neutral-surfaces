function [p, s, t] = ntp_bottle_to_cast(Sppc, Tppc, P, k, sB, tB, pB, tolp) %#codegen
%NTP_BOTTLE_TO_CAST  Find the Neutral Tangent Plane from a bottle to a cast.
%
% [p, s, t] = ntp_bottle_to_cast(Sppc, Tppc, P, k, sB, tB, pB, tolp)
% finds a point on the cast with properties (s, t, p) that is neutrally 
% related to the bottle of (sB, tB, pB), meaning that
%   eos(s, t, p_avg) = eos(sB, tB, p_avg)
% is approximately solved (there is an exact solution within tolp of p),
% where eos is the equation of state given by eos.m in MATLAB's path,
% and   p_avg = (pB + p) / 2 is the average of the bottle's pressure
%                            and the pressure on the cast,
% and   [s,t] are the salinity and temperature on the cast at pressure p
%             given by ppc_val2(P, Sppc, Tppc, p).
% The cast's hydrographic properties are given by piecewise polynomial
% interpolants for salinity and temperature as functions of pressure,
% given by coefficient arrays Sppc and Tppc with knots at P(1:k). If no
% such (s, t, p) is found, each of (s, t, p) is NaN.
%
% For a non-Boussinesq ocean, P, pB, and p are pressure.
% For a Boussinesq ocean, P, pB, and p are depth.
%
%
% --- Input:
% Sppc [O, K-1]: coefficients for piecewise polynomial for practical / 
%     Absolute Salinity in terms of Z on the cast.  Here, O is the
%     polynomial order and K is the number of data points (including NaN's)
%     on the cast.
% Tppc [O, K-1]: as above but for potential / Conservative Temperature.
% P [K, 1]: pressure or depth data points on the cast.
% k [1, 1]: number of valid (non-NaN) data points on the cast.
%          Specifically, Sppc(:,1:k-1) and Tppc(:,1:k-1) must be non-NaN. 
% sB [1 , 1]: practical / Absolute salinity of current bottle
% tB [1 , 1]: potential / Conservative temperature of current bottle
% pB [1 , 1]: pressure or depth of current bottle
% tolp [1, 1]: tolerance for solving the level of neutral buoyancy (same
%              units as P and pB)
%
% Note: physical units for Sppc, Tppc, P, sB, tB, pB, p, s, t are
% determined by eos.m.
%
% Note: P must increase monotonically along its first dimension.
%
%
% --- Output:
% p [1, 1]: pressure or depth in the cast 
%           that is neutrally related to the bottle
% s [1, 1]: practical / Absolute salinity in the cast
%           that is neutrally related to the bottle
% t [1, 1]: potential / Conservative temperature in the cast
%           that is neutrally related to the bottle
%
%
% --- See Also:
% ntp_midpoint_to_casts
% ppc_linterp, ppc_pchip

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


if k > 1

    % Search for a sign-change, expanding outward from an initial guess
    [lb, ub] = fzero_guess_to_bounds(@myfcn, pB, P(1), P(k), ...
      Sppc, Tppc, P, sB, tB, pB);
    
    if ~isnan(lb)
      % A sign change was discovered, so a root exists in the interval.
      % Solve the nonlinear root-finding problem using Brent's method
      p = fzero_brent(@myfcn, lb, ub, tolp, ...
        Sppc, Tppc, P, sB, tB, pB);
      
      % Interpolate S and T onto the updated surface
      [s,t] = ppc_val2(P, Sppc, Tppc, p);

    else
      p = nan;
      s = nan;
      t = nan;
    end

else
    p = nan;
    s = nan;
    t = nan;
end
end

function out = myfcn(p, Sppc, Tppc, P, sB, tB, pB)
% Evaluate difference between (a) eos at location on the cast (S, T, P)
% where the pressure or depth is p, and (b) eos of the bottle (sB, tB, pB);
% here, eos is always evaluated at the average pressure or depth, (p +
% pB)/2.
[s,t] = ppc_val2(P, Sppc, Tppc, p);
p_avg = (pB + p) / 2;
out = eos(sB, tB, p_avg) - eos(s, t, p_avg);
end