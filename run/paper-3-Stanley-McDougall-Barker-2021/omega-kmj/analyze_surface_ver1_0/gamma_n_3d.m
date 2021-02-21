function [gg,dggl,dggh] = gamma_n_3d(s,ct,p,lats,longs)

%           Calculate neutral density (gamma^n)
%
% Usage:    [gg,dggl,dggh] = gamma_n_3d(s,ct,p,lats,longs)
%
%           Label 3-dim s,ct,p field with neutral density (gamma^n) according 
%           to 'A neutral density variable for the world's oceans, Jackett and 
%           McDougall, JPO, 1997'.
%
%
% Input:    s           salinity: 
%                       either [depth lat long] or [time lat long] 
%           ct          conservative temperature: 
%                       either [depth lat long] or [time lat long] 
%           p           pressure:
%                       either [lat long] or [depth lat long] or 
%                       [time depth lat long]
%           lats        latitude
%           longs       longitude
%
% Output:   gg          matrix of gamma_n values       
%           dggl        matrix of gamma_n lower error estimates
%           dggh        matrix of gamma_n upper error estimates
%  
% Calls:    neutral density software
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           gg                          kg/m^3
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 5)
  error('gamma_n_3d.m: requires 5 input arguments')
end  

%% preallocate memory

gg = nan(size(s));
dggl = nan(size(s));
dggh = nan(size(s));

%% start calculation of approximate neutral surfaces

if (length(size(s)) == 3) % time-independent data
    
    [zi,yi,xi] = size(s);
    
    disp('preparing s, ct, p for gamma_n...');
    
    % find indices of surface salinity that are finite for lats in [-80,64]
    
    inds = find(isfinite(squeeze(s(1,:,:)))&-80<=lats&lats<=64); 
    
    % if p 2-dim then make p 3-dim
    
    if (length(size(p)) == 2)
        p_exp = repmat(p, [1 yi xi]);
        pp = squeeze(p_exp(:,inds));
    else
        pp = squeeze(p(:,inds));
    end

    % pack finite hydrography into one section of data
    
    ss = squeeze(s(:,inds)); tt = squeeze(ct(:,inds));
    longss = longs(inds); latss = lats(inds);                         
 
    disp('converting ct to t...');
    tt = t_from_ct(ss,tt,pp); 
    
    % run gamma_n
    
    disp('running gamma_n...');
    [g,dgl,dgh] = gamma_n(ss,tt,pp,longss,latss); 
    
    % pack output into 3-d array
    
    disp('writing gg, dggl, dggh into 3-d arrays...');
    gg(:,inds) = g; 
    dggl(:,inds) = dgl;
    dggh(:,inds) = dgh;

else % time-dependent data

    [ti,zi,yi,xi] = size(s); %#ok
    
    for i = 1:ti
        
        disp(['preparing s, ct, p at time ',int2str(i),' for gamma_n...']);
        
        % find indices of surface saltinity that are finite for lat in [-80,64]
        
        inds = find(isfinite(squeeze(s(i,1,:,:)))&-80<=lats&lats<=64); 
        
        % pack finite hydrography into one section of data
        
        ss = squeeze(s(i,:,inds)); ctt = squeeze(ct(i,:,inds)); pp = squeeze(p(i,:,inds)); 
        longss = longs(inds); latss = lats(inds);      
        
        disp(['converting ct to t for time ',int2str(i),'...']);
        tt = t_from_ct(ss,ctt,pp); 
 
        % run gamma_n
        
        disp(['running gamma_n for time ',num2str(i),'...']);
        [g,dgl,dgh] = gamma_n(ss,tt,pp,longss,latss); 

        % pack output into 4-d array
        
        disp(['writing gg, dggl, dggh for time ',num2str(i),' into 4-d arrays...']);
        gg(i,:,inds) = g; 
        dggl(i,:,inds) = dgl;
        dggh(i,:,inds) = dgh;
    end
end