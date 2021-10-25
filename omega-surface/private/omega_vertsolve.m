function [p, s, t] = omega_vertsolve(Sppc, Tppc, P, BotK, s, t, p, tolp, phi) %#codegen
%OMEGA_VERTSOLVE  Root finding of a new surface with a specified a density
%                 difference from the current surface
%
%
% [s, t, p] = omega_vertsolve(Sppc, Tppc, P, BotK, s0, t0, p0, tolp, phi)
% determines pressures p that satisfy
%   |eos(S_n(p(n)), T_n(p(n)), p0(n)) + phi(n) - eos(s0(n), t0(n), p0(n))| < tolp
% where S_n(p') and T_n(p') are interpolants whose coefficients are given
% by Sppc(:,:,n) and Tppc(:,:,n) in water column n. These interpolants 
% determine salinities s and temperatures t on the surface, namely
%   s(n) = S_n(p(n)) and t(n) = T_n(p(n)).
% The function eos.m determines either the in-situ density or the specific
% volume.
%
%
% --- Input
% Sppc [O, K-1, N]: coefficients for piecewise polynomial for practical
%                   / Absolute Salinity in terms of P
% Tppc [O, K-1, N]: coefficients for piecewise polynomial for potential
%                   / Conservative Temperature in terms of P
% P [K, N]: knots for the pressure or depth of the casts
% BotK [1, N]: number of valid data points on each cast
% s [1, N]: practical / Absolute salinity on the initial surface
% t [1, N]: potential / Conservative temperature on the initial surface
% p [1, N]: pressure or depth on the initial surface
% tolp [1,1]: precision of solution in pressure or depth
% phi [1,N]: the desired in-situ density or specific volume change of the
%            surface
%
% Note: O is the order of the piecewise polynomials
%       K is the maximum number of knots in these piecewise polynomials,
%           i.e. the maximum number of bottles in any cast
%       N is the number of water columns (possibly including land).
%
% Note: P must increase along its first dimension.
%
% Note: P can have size [K, 1], in which case it is used for each cast.
%
% Note: variables can actually be higher dimensional, e.g. N = [ni, nj],
%       and p can be any dimensional matrix, so long as it has N elements
%       in total.
%
% Note: BotK should be given by
%           BotK = squeeze(sum(isfinite(S), 1));
%
%
% --- Output
% p [same as input p]: pressure or depth of the updated surface
% s [same as input p]: practical / Absolute salinity of the updated surface
% t [same as input p]: potential / Conservative temperature of the updated surface
%
%
% --- Units
% The units of s, t, p, , T, P, tolp, and phi are determined by the
% function eos.m.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com



% Inputs s0, t0, and p0 are named s, t, p so operations are done in-place.

N = numel(p);
Pmat = ~isvector(P);

% Loop over each water column
for n = 1:N
  d = phi(n);
  k = BotK(n);
  if ~isnan(d) && k > 1
    
    % Select this water column
    Sppcn = Sppc(:,1:k-1,n);
    Tppcn = Tppc(:,1:k-1,n);
    if Pmat
      Pn = P(1:k,n);
    else
      Pn = P((1:k).'); % .' is for codegen, so P and (1:k).' both column vectors
    end
    
    rn = eos(s(n), t(n), p(n));
    
    % Search for a sign-change, expanding outward from an initial guess
    %[lb, ub] = fzero_guess_to_bounds(@diff_avgx, p(n), Pn(1), Pn(k), ...
    %  Sppcn, Tppcn, Pn, s(n), t(n), p(n), d);  % DEV:  testing reference p = average of current and point-to-be
    [lb, ub] = fzero_guess_to_bounds(@diff_p0, p(n), Pn(1), Pn(k), ...
      Sppcn, Tppcn, Pn, rn, p(n), d);
    
    if ~isnan(lb)
      % A sign change was discovered, so a root exists in the interval.
      % Solve the nonlinear root-finding problem using Brent's method
      %p(n) = fzero_brent(@diff_avgx, lb, ub, tolp, ...
      %  Sppcn, Tppcn, Pn, s(n), t(n), p(n), d);  % DEV:  testing reference p = average of current and point-to-be
      p(n) = fzero_brent(@diff_p0, lb, ub, tolp, ...
        Sppcn, Tppcn, Pn, rn, p(n), d);
      
      % Interpolate S and T onto the updated surface
      [s(n),t(n)] = ppc_val2(Pn, Sppcn, Tppcn, p(n));
    else
      p(n) = nan;
      s(n) = nan;
      t(n) = nan;
    end
    
  else
    % phi is nan, or only one grid cell so cannot interpolate.
    % This will ensure s,t,p all have the same nan structure
    p(n) = nan;
    s(n) = nan;
    t(n) = nan;
  end
  
end

end


function out = diff_p0(p, Sppc, Tppc, P, r0, p0, d)
% Evaluate difference between (a) eos at location on the cast where the
% pressure or depth is p, and (b) eos at location on the cast where the
% pressure or depth is p0 (where the surface currently is) plus the density
% perturbation d.  Part (b) is precomputed as r0.  Here, eos always
% evaluated at the pressure or depth of the original position, p0; this is
% to calculate locally referenced potential density with reference pressure
% p0.

% Interpolate S and T to the current pressure or depth
[s,t] = ppc_val2(P, Sppc, Tppc, p);

% Calculate the potential density or potential specific volume difference
%out =  eos(s, t, p0) - d - r0 ;
out =  eos(s, t, p0) - r0 - d ;

end

%{
function out = diff_avgx(p, Sppc, Tppc, P, s0, t0, p0, d)
% Evaluate difference between (a) eos at location on the cast where the
% pressure or depth is p, and (b) eos at location on the cast where the
% pressure or depth is p0 (where the surface currently is) plus the density
% perturbation d.  Here, eos is always evaluated at the average pressure or
% depth, (p + p0)/2; this is to calculate locally referenced potential
% density with reference pressure (p + p0)/2.

% Interpolate S and T to the current pressure or depth
[s,t] = ppc_val2(P, Sppc, Tppc, p);

% Average the current pressure or depth and the original pressure or depth
p_avg = (p + p0) / 2;

% Calculate the potential density or potential specific volume difference
out =  eos(s, t, p_avg) - d - eos(s0, t0, p_avg) ;
%out =  eos(s, t, 1500) - d - eos(s0, t0, 1500) ; % DEV: testing omega software to find potential density surface

end
%}