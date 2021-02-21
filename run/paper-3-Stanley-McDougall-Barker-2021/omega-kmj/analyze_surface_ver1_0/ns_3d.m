function [sns,ctns,pns,dsns,dctns,dpns] = ns_3d(s,ct,p,rho,rholevels)

%           Find s, ct and p on a density surface
%
% Usage:    [sns,ctns,pns,dsns,dctns,dpns] = ns_3d(s,ct,p,rho,rholevels)
%
%           For a 3-dim field of hydrographic data (s,ct,p) that has been 
%           labelled with a density variable, find the salinities, 
%           conservative temperatures and pressures of density surfaces 
%           specified by rholevels.
%
%
% Input:    s           salinity:
%                       either [depth lats longs] or [time depth lats longs] 
%           ct          conservative temperature: 
%                       either [depth lats longs] or [time depth lats longs] 
%           p           pressure:
%                       either [lats longs] or [depth lats longs] or 
%                       [time depth lats longs]
%           rho         density values:
%                       either [depth lats longs] or [time depth lats longs]
%           rholevels   array of density values
%           lats        latitude
%
% Output:   sns         salinity on the density surfaces       
%           ctns        conservative temperature on the density surfaces  
%           pns         pressure on the density surfaces  
%           dsns        surface salinity errors
%           dctns       surface temperature errors
%           dpns        surface pressure errors
%  
% Calls:    neutral_surfaces.m, neutral_surfaces0.m, t_from_ct.m,
%           ct_from_t.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           rho                         kg/m^3
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 5)
  error('ns_3d.m: requires 5 input arguments')
end 

[dummy,n] = size(rholevels);

%% start calculation of approximate neutral surfaces

time_check = length(size(s));

switch time_check
    
    case 3 % time-independent data

        [zi,yi,xi] = size(s);
        sns = nan(n,yi,xi);
        ctns = nan(n,yi,xi);
        pns = nan(n,yi,xi);
        dsns = nan(n,yi,xi);
        dctns = nan(n,yi,xi);
        dpns = nan(n,yi,xi);

        % find indices of surface saltinity that are finite
        
        inds = find(isfinite(squeeze(s(1,:,:))));
        
        % if p 2-dim then make p 3-dim
        
        if (length(size(p)) == 2)
            p_exp = repmat(p, [1 yi xi]);
            pp = p_exp(:,inds);
        else
            pp = p(:,inds);
        end

        % pack finite hydrography into one section of data
        
        ss = s(:,inds); tt = ct(:,inds);

        gg = rho(:,inds);
        gg2=change(gg,'<',0,NaN);
        
        [s_ns,ct_ns,p_ns,ds_ns,dct_ns,dp_ns] = neutral_surfaces(ss,tt,pp,gg2,rholevels);

        % pack output into 3-d g array

        sns(:,inds) = s_ns;
        ctns(:,inds) = ct_ns;
        pns(:,inds) = p_ns;
        dsns(:,inds) = ds_ns;
        dctns(:,inds) = dct_ns;
        dpns(:,inds) = dp_ns;

    case 4 % time-dependent data

        [ti,zi,yi,xi] = size(s);
        sns = nan(ti,n,yi,xi);
        ctns = nan(ti,n,yi,xi);
        pns = nan(ti,n,yi,xi);
        dsns = nan(ti,n,yi,xi);
        dctns = nan(ti,n,yi,xi);
        dpns = nan(ti,n,yi,xi);

        for l = 1:ti

            % find indices of surface saltinity that are finite
            inds = find(isfinite(squeeze(s(l,1,:,:))));
            
            % pack finite hydrography into one section of data
            ss = squeeze(s(l,:,inds)); tt = squeeze(ct(l,:,inds)); pp = squeeze(p(l,:,inds));

            gg = squeeze(rho(l,:,inds));
            gg2 = change(gg,'<',0,NaN);
            
            [s_ns,ct_ns,p_ns,ds_ns,dct_ns,dp_ns] = neutral_surfaces(ss,tt,pp,gg2,rholevels);

            % pack output into 4d g array

            sns(l,:,inds) = s_ns;
            ctns(l,:,inds) = ct_ns;
            pns(l,:,inds) = p_ns;
            dsns(l,:,inds) = ds_ns;
            dctns(l,:,inds) = dct_ns;
            dpns(l,:,inds) = dp_ns;
        end
        
    otherwise
        
        error('ns_3d.m: wrong size of input arguments')
        
end