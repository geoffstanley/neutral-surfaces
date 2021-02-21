function [streamfunc] = geo_streamfunc(s,ct,p,sns,ctns,pns,sref,ctref,pref)

%           Calculate the geostrophic streamfunction on approximately
%           neutral surfaces or potential density surfaces
%
% Usage:    [streamfunc] = geo_streamfunc(s,ct,p,sns,ctns,pns,sref,ctref,pref)
%
%           Calculate the geostrophic streamfunction on approximately
%           neutral surfaces according to 'An approximate geostrophic 
%           streamfunction for use on density surfaces, McDougall, 
%           Trevor J. and A. Klocker, Journal of Physical Oceanography, 
%           in prep.'
%
% Input:    s                     absolute salinity
%           ct                    conservative temperature
%           p                     pressure
%           sns                   absolute salinity on density surface
%           ctns                  conservative temperature on density surface
%           pns                   pressure on density surface
%           sref                  reference absolute salinity
%           ctref                 reference conservative temperature
%           pref                  reference pressure
%
%                 
% Output:   streamfunc            geostrophic streamfunction 
%
% Calls:    EOS06 (rho_from_ct.m, ab_from_ct.m, eosall_from_ct.m),
%           var_on_surf.m
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

if ~(nargin == 10) && ~(nargin == 11)
    error('geo_stream.m: requires 10 or 11 input arguments')
end

[zi,dummy,dummy] = size(p); %#ok

%% set reference values

sref = sref * ones(size(sns));
ctref = ctref * ones(size(sns));
pref = pref * ones(size(sns));

%% change units from dbar to Pa

pns_pa = pns .* 10000; 
p_pa = p .* 10000;
pref_pa = pref .* 10000;

%% calculate term 1 

% calculate steric anomaly on surface

v_tilde = 1./rho_from_ct(sns,ctns,pns);
v_tilde_2 = 1./rho_from_ct(sref,ctref,pns);
steric_tilde = v_tilde - v_tilde_2;

% pressure difference of pressure on surface and reference pressure

p_diff = pns_pa-pref_pa;

% calculate term 1 

term1 = (1/2) .* p_diff .* steric_tilde;

%% calculate term 2 

ct_diff = ctns - ctref;

% calculate thermobaric parameter

[alpha,beta,aonb,T_b] = ab_from_ct(sns,ctns,pns); 

% change from dbar into Pa

T_b = T_b ./ 10000; 

% calculate specific volume anomaly on surface

v_surf = 1./rho_from_ct(sns,ctns,pns);

% calculate term 2 

term2 = (1/12) .* (T_b .*v_surf) .* ct_diff .* p_diff .*  p_diff;

%% calculate term 3 

% calculate steric anomaly

v_all = 1./rho_from_ct(s,ct,p);

sref_rep = repmat(sref,[zi 1 1]);
ctref_rep = repmat(ctref,[zi 1 1]);

v_ref = 1./rho_from_ct(sref_rep,ctref_rep,p);
steric_cast = v_all - v_ref;

% calculate integral of steric anomaly

steric_cast_mid  = 0.5 * (steric_cast(2:end,:,:) + steric_cast(1:end-1,:,:));
top = steric_cast(1,:,:).*p_pa(1,:,:);
delta_steric = (steric_cast_mid .* diff(p_pa,1));
steric_int  = cumsum([top;delta_steric]);

% vertically interpolate steric anomaly integral onto pressure surface

p_mid_pa = 0.5 * (p_pa(2:end,:,:) + p_pa(1:end-1,:,:));

% calculate term 3

term3 = var_on_surf(pns_pa,p_mid_pa,steric_int);

%% calculate streamfunction

streamfunc = term1 - term2 - term3;