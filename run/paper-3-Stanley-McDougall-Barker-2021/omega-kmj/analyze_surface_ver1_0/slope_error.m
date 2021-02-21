function [ss,sx,sy,curl_s,ee,ex,ey,curl_e,fdd] = slope_error(p,g,n2,sns,ctns,pns,e1t,e2t,keyword,wrap,Z2P)  % GJS added Z2P

%           Calculate slope errors
%
% Usage:    [ss,sx,sy,curl_s,ee,ex,ey,curl_e,ver] = slope_error(p,g,n2,sns,ctns,pns,e1t,e2t,keyword,wrap)
%
%           Calculate slope errors, density gradient errors, their curl and
%           the fictitious diapycnal diffusivity
%
% Input:    p           pressure
%           g           gravitational accelleration
%           sns         salinity on neutral surface
%           ctns        conservative temperature on neutral surface
%           pns         pressure on neutral surface
%           n2_ns       buoyancy frequency
%           e1t         zonal scale-factor at tracer points
%           e2t         meridional scale-factor at tracer points
%           lats        latitude
%           keyword     'op' for values on points
%                       'bp' for values in between points
%           wrap        'none'
%                       'long'
%
% Output:   ss         slope error
%           sx         x-component of slope error
%           sy         y-component of slope error
%           curl_s     curl of slope error
%           ee         epsilon
%           ex         x-component of epsilon
%           ey         y-component of epsilon
%           curl_e     curl of epsilon
%           fdd        fictitious diapycnal diffusivity
%
%
% Calls:    grad_surf.m, ab_from_ct.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           buoyancy frequency          s^-1
%           gravitational acceleration  m/s^2
%           scale factors               m
%           fictitious diapycnal diff.  m^2/s
%           
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

%if ~(nargin == 10) % GJS removed
%    error('slope_error.m: requires 10 input arguments') % GJS removed
%end % GJS removed
if nargin < 11; Z2P = 1; end  % GJS

%% initialize

[gi,dummy,dummy] = size(sns); %#ok
[zi,yi,xi] = size(p);

%% make g correct size

if (length(size(g)) == 2)
    g_nan = nan(size(sns));
    g_nan(1,:,:) = g;
    g = g_nan;
end

%% calculate slope errors and epsilon

switch keyword

    case 'op'

        [gradx_ct,grady_ct] = grad_surf(ctns,e1t,e2t,'op',wrap);
        [gradx_s,grady_s] = grad_surf(sns,e1t,e2t,'op',wrap);
        [alpha,beta] = ab_from_ct(sns,ctns,pns,Z2P);  % GJS added Z2P

        % calculate density gradient errors (epsilon)

        ex = ((beta .* gradx_s) - (alpha .* gradx_ct));
        ey = ((beta .* grady_s) - (alpha .* grady_ct));
        ee = ex + ey; 

        % calculate slope errors

        p_mid = (p(2:zi,:,:) + p(1:zi-1,:,:)) ./ 2;
        n2_ns = var_on_surf(pns,p_mid,n2);
        n2_ns = change(n2_ns,'==',0,nan);
        fac = g ./ n2_ns;

        sx = fac .* ex;
        sy = fac .* ey;
        ss = sx + sy; 

        fdd = 1000 * (ss .* ss);

    case 'bp'

        %disp('ss, curl_s, ee, curl_e and fictitious diapycnal diffusivity are not calculated due to sx, sy, ex, ey being at different locations');

        [gradx_ct,grady_ct] = grad_surf(ctns,e1t,e2t,'bp',wrap);
        [gradx_s,grady_s] = grad_surf(sns,e1t,e2t,'bp',wrap);
        [alpha_tmp,beta_tmp] = ab_from_ct(sns,ctns,pns,Z2P);  % GJS added Z2P

        alpha_x(1:gi,1:yi,1:xi-1) = 0.5 * (alpha_tmp(1:gi,1:yi,1:xi-1) + alpha_tmp(1:gi,1:yi,2:xi));
        alpha_y(1:gi,1:yi-1,1:xi) = 0.5 * (alpha_tmp(1:gi,1:yi-1,1:xi) + alpha_tmp(1:gi,2:yi,1:xi));
        beta_x(1:gi,1:yi,1:xi-1) = 0.5 * (beta_tmp(1:gi,1:yi,1:xi-1) + beta_tmp(1:gi,1:yi,2:xi));
        beta_y(1:gi,1:yi-1,1:xi) = 0.5 * (beta_tmp(1:gi,1:yi-1,1:xi) + beta_tmp(1:gi,2:yi,1:xi));

        p_mid = (p(2:zi,:,:) + p(1:zi-1,:,:)) ./ 2;
        n2_ns = var_on_surf(pns,p_mid,n2);
        n2_ns = change(n2_ns,'==',0,nan);

        fac_tmp = g ./ n2_ns;

        fac_x(1:gi,1:yi,1:xi-1) = 0.5 * (fac_tmp(1:gi,1:yi,1:xi-1) + fac_tmp(1:gi,1:yi,2:xi));
        fac_y(1:gi,1:yi-1,1:xi) = 0.5 * (fac_tmp(1:gi,1:yi-1,1:xi) + fac_tmp(1:gi,2:yi,1:xi));

        switch wrap

            case 'none'

                alpha_x(1:gi,1:yi,xi) = nan;
                beta_x(1:gi,1:yi,xi) = nan;
                fac_x(1:gi,1:yi,xi) = nan;
                alpha_y(1:gi,yi,1:xi) = nan;
                beta_y(1:gi,yi,1:xi) = nan;
                fac_y(1:gi,yi,1:xi) = nan;

            case 'long'

                alpha_x(1:gi,1:yi,xi) = 0.5 * (alpha_tmp(1:gi,1:yi,xi) + alpha_tmp(1:gi,1:yi,1));
                beta_x(1:gi,1:yi,xi) = 0.5 * (beta_tmp(1:gi,1:yi,xi) + beta_tmp(1:gi,1:yi,1));
                fac_x(1:gi,1:yi,xi) = 0.5 * (fac_tmp(1:gi,1:yi,xi) + fac_tmp(1:gi,1:yi,1));
                alpha_y(1:gi,yi,1:xi) = nan;
                beta_y(1:gi,yi,1:xi) = nan;
                fac_y(1:gi,yi,1:xi) = nan;

        end

        % calculate density gradient errors (epsilon)

        ex = ((beta_x .* gradx_s) - (alpha_x .* gradx_ct));
        ey = ((beta_y .* grady_s) - (alpha_y .* grady_ct));
        ee = nan; % ex and ey are at different locations

        % calculate slope errors

        sx = fac_x .* ex;
        sy = fac_y .* ey;
        ss = nan; % sx and sy are at different locations

        fdd = nan;

    otherwise
        error('slope_error.m: keyword must be "op" or "bp"')
end

%% calculate curl of slope errors and epsilon

switch keyword

    case 'op'

        % calculate curl of epsilon and curl of slope errors

        [dummy,ex_diff] = grad_surf(ex,e1t,e2t,'op',wrap); %#ok
        [ey_diff,dummy] = grad_surf(ey,e1t,e2t,'op',wrap); %#ok

        [dummy,sx_diff] = grad_surf(sx,e1t,e2t,'op',wrap); %#ok
        [sy_diff,dummy] = grad_surf(sy,e1t,e2t,'op',wrap); %#ok

        curl_e = ey_diff-ex_diff;
        curl_s = sy_diff-ex_diff;

    otherwise

        curl_e = [];
        curl_s = [];
end