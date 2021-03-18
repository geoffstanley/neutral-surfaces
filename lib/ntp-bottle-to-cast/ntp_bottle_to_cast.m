function [p, s, t, success] = ntp_bottle_to_cast(Sppc, Tppc, P, k, sB, tB, pB, tolp) %#codegen
%NTP_BOTTLE_TO_CAST  Find a bottle's level of neutral buoyancy in a water
%                    column, using the Neutral Tangent Plane relationship.
%
% [p, s, t] = ntp_bottle_to_cast(Sppc, Tppc, P, k, sB, tB, pB, tolp)
% finds (s, t, p), with precision in p of tolp, that is at the level of
% neutral buoyancy for a fluid bottle of (sB, tB, pB) in a water column of
% with piecewise polynomial interpolants for S and T given by Sppc and Tppc 
% with knots at P(1:k).  Specifically, s and t are given by
%   [s,t] = ppc_val(P, Sppc, Tppc, p)
% and p satisfies
%      eos(s, t, p') = eos(sB, tB, p')
% where eos is the equation of state given by eos.m in MATLAB's path,
% and   p' is in the range [p_avg - tolp/2, p_avg + tolp/2],
% and   p_avg = (pB + p) / 2 is the average of the fluid bottle's original
%                          and final pressure or depth.
%
% [p, s, t, success] = ntp_bottle_to_cast(...)
% returns a flag value success that is true if a valid solution was found,
% false otherwise.
%
% For a non-Boussinesq ocean, P, pB, and p are pressure.
% For a Boussinesq ocean, P, pB, and p are depth.
%
%
% --- Input:
% Sppc [O, K-1]: coefficients for piecewise polynomial for practical 
%                   / Absolute Salinity in terms of P
% Tppc [O, K-1]: coefficients for piecewise polynomial for potential 
%                   / Conservative Temperature in terms of P
% P [K, 1]: pressure or depth in water column
% k [1, 1]: number of valid (non-NaN) data points in the water column.
%          Specifically, Sppc(end,1:k) and Tppc(end,1:k) must all be valid.
% sB [1 , 1]: practical / Absolute salinity of current bottle
% tB [1 , 1]: potential / Conservative temperature of current bottle
% pB [1 , 1]: pressure or depth of current bottle
% tolp [1, 1]: tolerance for solving the level of neutral buoyancy (same
%             units as P and pB)
%
% Note: physical units for Sppc, Tppc, P, sB, tB, pB, p, s, t  are
% determined by eos.m.
%
% Note: P must increase monotonically along its first dimension.
%
%
% --- Output:
% p [1, 1]: pressure or depth in water column at level of neutral buoyancy
% s [1, 1]: practical / Absolute salinity in water column at level of neutral buoyancy
% t [1, 1]: potential / Conservative temperature in water column at level of neutral buoyancy
% success [1,1]: true if a valid solution was found, false otherwise.

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

      success = true;
    else
      p = nan;
      s = nan;
      t = nan;
      success = false;
    end

    
else
    p = nan;
    s = nan;
    t = nan;
    success = false;
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