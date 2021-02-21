function [mld_dpth] = mld(s,ct,p)

%           Calculate mixed-layer depth
%
% Usage:    [mld_dpth] = mld(s,ct,p)
%
%           Calculate the mixed-layer depth as in 'Mixed layer depth over the 
%           global ocean: An examination of profile data and a profile-based 
%           climatology, de Boyer Montégut et al., JGR, 2004'. The criterion
%           selected is a threshold value of density from a surface value 
%           (Δrho = 0.03 kg m−3). 
%
% Input:    s         salinity 
%           ct        conservative temperature 
%           p         pressure 
%
% Output:   mld       mixed-layer depth 
%  
% Calls:    gpoly16ct.m
%
% Units:    salinity                 psu (IPSS-78)
%           conservative temperature degrees C (IPS-90)
%           pressure                 dbar
%           mld                      dbar
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%


%% check input arguments

if ~(nargin == 3)
  error('mld.m: requires 3 input arguments')
end

%% preallocate memory

[zi,yi,xi] = size(s);
mld_dpth = nan(yi,xi);

%% calculate gamma^rf

dens = gpoly16ct(s,ct);

%% calculate mixed layer depth

for j = 1:yi
    for i = 1:xi
        min_dens = dens(1,j,i);
        if isnan(min_dens)
            mld_dpth(j,i) = nan;
        else
            diff = min_dens + 0.3 - dens(:,j,i);
            inds = find(diff(:) > 0);
            mld_dpth(j,i) = p(inds(end),j,i);
        end
    end
end