function [streamfunc] = montgomery_streamfunc(s,ct,p,sns,ctns,pns,sref,ctref,pref,method)

%           Calculate the Montgomery/Zhang & Hogg  streamfunction
%
% Usage:    [mont_streamfunc] = montgomery_streamfunc(s,ct,p,sns,ctns,pns,...
%               sref,ctref,pref,method)
%
%           Calculate the Montgomery streamfunction (the acceleration potential)
%           in a steric anomaly surface according to 'A suggested method
%           for representing gradient flow in isentropic surfaces, Montgomery, 
%           Bull. Amer. Meteor. Soc., 18, 210-212.' or according to 
%           'Circulation and water mass balance in the Brazil Basin, Zhang, H-M 
%           and N. G. Hogg, J. Marine Research, 50, 385-420.'
%
% Input:    s                     salinity
%           ct                    conservative temperature
%           p                     pressure
%           sns                   salinity on density surface
%           ctns                  conservative temperature on density surface
%           pns                   pressure on density surface
%           sref                  reference salinity
%           ctref                 reference temperature
%           pref                  reference pressure
%           method                'mont' streamfunction according to
%                                  Montgomery
%                                 'zhang_hogg' streamfunction
%                                 according to Zhang and Hogg
%                 
% Output:   streamfunc            Montgomery/Zhang & Hogg streamfunction 
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

if ~(nargin == 10)
    error('montgomery_streamfunc.m: requires 10 input arguments')
end

%% get size of data set

[zi,dummy,dummy] = size(p); %#ok

%% make reference values into a matrix the same size as surface values

s_ref = sref * ones(size(sns));
ct_ref = ctref * ones(size(sns));
p_ref = pref * ones(size(sns));

%% change units from dbar to Pa

pns_pa = pns .* 10000; 
p_pa = p .* 10000;
p_ref_pa = p_ref .* 10000;

%% calculate steric anomaly on surface 

v = 1./rho_from_ct(sns,ctns,pns);
v_2 = 1./rho_from_ct(s_ref,ct_ref,pns);

steric_surf = v - v_2;

%% calculate steric anomaly at tracer grid points above surface

v_all = 1./rho_from_ct(s,ct,p);

s_ref = repmat(s_ref,[zi 1 1]);
ct_ref = repmat(ct_ref,[zi 1 1]);

v_cast = v_all;
v_2_cast = 1./rho_from_ct(s_ref,ct_ref,p);
steric_cast = v_cast - v_2_cast;

%% calculate pressure integrals of steric anomaly

steric_cast_mid  = 0.5 * (steric_cast(2:end,:,:) + steric_cast(1:end-1,:,:));

top = steric_cast(1,:,:).*p_pa(1,:,:);
delta_steric = (steric_cast_mid .* diff(p_pa,1));

% calculate cumulative pressure integral of steric anomaly 

steric_int  = cumsum([top;delta_steric]);

% interpolate cumulative pressure integral of steric anomaly onto density
% surface 

p_mid_pa = 0.5 * (p_pa(2:end,:,:) + p_pa(1:end-1,:,:));
steric_int_surf = var_on_surf(pns_pa,p_mid_pa,steric_int);

steric_int_surf = change(steric_int_surf,'==',0,nan);

%% calculate streamfunction

switch method
    case 'mont'
        streamfunc = pns_pa .* steric_surf - steric_int_surf;
    case 'zhang_hogg'
        streamfunc = (pns_pa - p_ref_pa) .* steric_surf - steric_int_surf;
end