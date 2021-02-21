function [streamfunc] = bernoulli_streamfunc(s,ct,p,sns,ctns,pns)

%           Calculate the Bernoulli streamfunction
%
% Usage:    [streamfunc] = bernoulli_streamfunc(s,ct,p,sns,ctns,pns)
%
%           Calculate the Bernoulli streamfunction according to 'Approximate 
%           geostrophic streamfunctions in potential density and approximately 
%           neutral density surfaces, McDougall, Trevor J. and A. Klocker, 
%           Journal of Physical Oceanography, in prep.'
%
% Input:    s                 salinity
%           ct                conservative temperature
%           p                 pressure
%           sns               salinity on density surface
%           ctns              conservative temperature on density surface
%           pns               pressure on density surface
%                 
% Output:   streamfunc        Bernoulli streamfunction 
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
    error('bernoulli_streamfunc.m: requires 6 input arguments')
end

%% get size of data set

[zi,dummy,dummy] = size(p); %#ok

%% change units from dbar to Pa

pns_pa = pns .* 10000; 
p_pa = p .* 10000;

%% calculate steric anomaly on casts at tracer grid points

sns_rep = repmat(sns,[zi 1 1]);
ctns_rep = repmat(ctns,[zi 1 1]);

v_hat = 1./rho_from_ct(s,ct,p);
v_hat_2 = 1./rho_from_ct(sns_rep,ctns_rep,p);

steric_hat = v_hat - v_hat_2;

%% calculate pressure integrals of steric anomaly

steric_hat_mid  = 0.5 * (steric_hat(2:end,:,:) + steric_hat(1:end-1,:,:));

top = steric_hat(1,:,:).*p_pa(1,:,:);
delta_steric = (steric_hat_mid .* diff(p_pa,1));

% calculate cumulative pressure integral of steric anomaly 

steric_int  = - cumsum([top;delta_steric]);

% interpolate cumulative pressure integral of steric anomaly onto density 
% surface 

p_mid_pa = 0.5 * (p_pa(2:end,:,:) + p_pa(1:end-1,:,:));
streamfunc = var_on_surf(pns_pa,p_mid_pa,steric_int);

streamfunc = change(streamfunc,'==',0,nan);