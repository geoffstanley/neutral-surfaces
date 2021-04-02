function [p,s,t] = pot_dens_surf_vertsolve(Sppc, Tppc, P, BotK, p, pref, val, tolp) %#codegen
%POT_DENS_SURF_VERTSOLVE  Helper function for pot_dens_surf, solving 
%                         non-linear root finding problem in each water column
%
%
% Note: to ensure the generated MEX function gives the same output as
% running this in native MATLAB, the MEX function must be generated with K
% specified as an integer class, not double: we use uint16.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


N = numel(p);
Pmat = ~isvector(P);


s = nan(size(p));
t = nan(size(p));

% Loop over each cast
for n = 1:N
    k = BotK(n);
    if k > 1
        
        % Select this water column
        Sppcn = Sppc(:,1:k-1,n);
        Tppcn = Tppc(:,1:k-1,n);
        if Pmat
          Pn = P(1:k,n);
        else
          Pn = P((1:k).'); % .' is for codegen, so X and (1:k).' both column vectors
        end
        
        % Initial guess could be nan.  In this case, try initial guess at mid-depth.
        if isnan(p(n))
            p(n) = (P(1) + P(k)) * 0.5;
        end
        
        % Search for a sign-change, expanding outward from an initial guess 
        [lb, ub] = fzero_guess_to_bounds(@myfcn, p(n), Pn(1), Pn(k), ...
          Sppcn, Tppcn, Pn, pref, val);
        
        if ~isnan(lb)
          % A sign change was discovered, so a root exists in the interval.
          % Solve the nonlinear root-finding problem using Brent's method
          p(n) = fzero_brent(@myfcn, lb, ub, tolp, ...
            Sppcn, Tppcn, Pn, pref, val);
          
          % Interpolate S and T onto the updated surface
          [s(n),t(n)] = ppc_val2(Pn, Sppcn, Tppcn, p(n));
        else
          p(n) = nan;
          s(n) = nan;
          t(n) = nan;
        end
        
    end
end

end

function out = myfcn(p, Sppc, Tppc, P, pref, val)
[s,t] = ppc_val2(P, Sppc, Tppc, p);
out = eos(s, t, pref) - val;
end