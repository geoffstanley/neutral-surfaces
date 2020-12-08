function [n2,p_mid] = bfrq_ct(s,ct,p,g)

%           Calculate the buoyancy frequency, N^2 
%
% Usage:    [n2,p_mid] = bfrq_ct(s,ct,p,g)
%
%           Calculate Brunt-Vaisala frequency squared (N^2) at the mid depths
%           from the equation,
%
%               -g      d(pdens)
%           N2 =  ----- x --------
%               pdens     d(z)
%
% Input:    s           salinity 
%           ct          conservative temperature 
%           p           pressure 
%           g           gravitational acceleration [z y x] or [y x]
%           (if no g is given g = 9.8 will be used)
%
% Output:   n2          Brunt-Vaisala frequency squared (N^2)
%           p_mid       mid pressure between p-points
%  
% Calls:    z_from_p.m, rho_from_ct.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin==3 || nargin==4)
    error('bfrq_ct.m: requires 3 or 4 input arguments')
end

%% if g is not given set g=9.8

if nargin == 3
    g = 9.8*ones(size(p));
end

%% get depth from pressure

z = z_from_p(p,g);

%% make g 3-dim.

if (length(size(g)) == 2)
    [k,dummy,dummy] = size(p); %#ok
    g_rep = repmat(g,[1 1 k]);
    g = permute(g_rep,[3 1 2]);
end

%% calculate buoyancy frequency

[zi,dummy,dummy] = size(p); %#ok

p_up = p(1:zi-1,:,:);
p_lo = p(2:zi,:,:);
p_mid = (p_lo + p_up) ./ 2;

rho_up = rho_from_ct(s(1:zi-1,:,:),ct(1:zi-1,:,:),p_mid);
rho_lo = rho_from_ct(s(2:zi,:,:),ct(2:zi,:,:),p_mid);

rho_mid = (rho_up + rho_lo) ./ 2;
rho_diff = rho_up - rho_lo;
g_mid = (g(1:zi-1,:,:) + g(2:zi,:,:)) ./ 2;
z_diff = diff(z);

n2 = -g_mid .* rho_diff ./ (z_diff .* rho_mid);