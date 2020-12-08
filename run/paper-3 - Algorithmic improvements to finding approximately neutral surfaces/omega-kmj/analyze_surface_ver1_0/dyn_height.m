function [dyn_h] = dyn_height(s,ct,p,sref,ctref,pref)

%           Calculate dynamic height
%
% Usage:    [dyn_h] = dyn_height(s,ct,p,sref,ctref,pref)
%
%           Calculate dynamic height at a particular pressure level in a 
%           3-dimensional hydrography
%
% Input:    s                     salinity
%           ct                    conservative temperature
%           p                     pressure
%           sref                  reference salinity 
%           ctref                 reference conservative temperature
%           pref                  reference pressure
%                 
% Output:   dyn_height            dynamic height
%
% Calls:    EOS06 (rho_from_ct.m), var_on_surf.m
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

if ~(nargin == 6) 
    error('dyn_height.m: requires 6 input arguments')
end

%% change units from dbar to Pa

pref_pa = pref .* 10000;
p_pa = p .* 10000;

%% calculate steric anomaly

v_all = 1./rho_from_ct(s,ct,p);
v_ref = 1./rho_from_ct(sref*ones(size(p)),ctref*ones(size(p)),p);

steric_cast = v_all - v_ref;

%% calculate integral of steric anomaly

steric_cast_mid  = 0.5 * (steric_cast(2:end,:,:) + steric_cast(1:end-1,:,:));
top = steric_cast(1,:,:).*p_pa(1,:,:);
delta_steric = (steric_cast_mid .* diff(p_pa,1));
steric_int  = - cumsum([top;delta_steric]);

%% vertically interpolate steric anomaly integral onto pressure surface

p_mid_pa = 0.5 * (p_pa(2:end,:,:) + p_pa(1:end-1,:,:));
dyn_h = var_on_surf(pref_pa*ones(size(steric_int)),p_mid_pa,steric_int);

%% replace zeros with NaNs

dyn_h = change(dyn_h,'==',0,nan);