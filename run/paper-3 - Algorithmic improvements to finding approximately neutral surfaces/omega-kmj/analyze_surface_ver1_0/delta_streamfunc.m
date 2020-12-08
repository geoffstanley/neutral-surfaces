function [delta_stream_ew,delta_stream_ns] = delta_streamfunc(s,ct,p,sns,ctns,pns,wrap)

%           Calculate local difference in streamfunctions
%
% Usage:    [delta_stream_ew,delta_stream_ns] =
%               delta_streamfunc(s,ct,p,sns,ctns,pns,density)
%
%           Calculate local differences in locally referenced streamfunctions 
%           according to 'An approximate geostrophic streamfunction for use 
%           on density surfaces, McDougall, Trevor J. and A. Klocker, 
%           Journal of Physical Oceanography, in prep.'
%
% Input:    s                     salinity
%           ct                    conservative temperature
%           p                     pressure
%           sns                   salinity on density surface 
%           ctns                  conservative temperature on density
%                                 surface
%           pns                   pressure on density surface
%           wrap                  'none'
%                                 'long'
%                 
% Output:   delta_stream_ew       east-west difference of streamfunction
%           delta_stream_ns       north-south difference of streamfunction
%
% calls:    EOS06 (rho_from_ct.m, ab_from_ct.m, eosall_from_ct.m), var_on_surf.m
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

if ~(nargin == 7)
    error('delta_streamfunc.m: requires 7 input arguments')
end

%% get size of hydrography

[zi,yi,xi] = size(p);

%% change units from dbar to Pa

pns_pa = pns .* 10000; 
p_pa = p .* 10000;

%% preallocate memory

delta_p_ew = nan(size(pns));
delta_ct_ew = nan(size(pns));
delta_p_ns = nan(size(pns));
delta_ct_ns = nan(size(pns));
v_tilde_e = nan(size(pns));
v_2_tilde_e = nan(size(pns));
v_hat_e = nan(size(pns));
v_2_hat_e = nan(size(pns));
v_tilde_n = nan(size(pns));
v_2_tilde_n = nan(size(pns));
v_hat_n = nan(size(pns));
v_2_hat_n = nan(size(pns));
v_cast_e  = nan(size(pns));
v_2_cast_e = nan(size(pns));
v_cast_n  = nan(size(pns));
v_2_cast_n = nan(size(pns));

delta_e_steric = nan(zi-1,yi,xi);
steric_int_e_surf = nan(1,yi,xi);
delta_n_steric = nan(zi-1,yi,xi);
steric_int_n_surf = nan(1,yi,xi);

%% calculate term 1 

% calculate delta_p and delta_ct

% for east-west links
    
delta_p_ew(1,1:yi,1:xi-1) = pns_pa(1,1:yi,2:xi) - pns_pa(1,1:yi,1:xi-1); 
delta_ct_ew(1,1:yi,1:xi-1) = ctns(1,1:yi,2:xi) - ctns(1,1:yi,1:xi-1);

switch wrap
    
    case 'long'
        
        delta_p_ew(1,1:yi,xi) = pns_pa(1,1:yi,1) - pns_pa(1,1:yi,xi); 
        delta_ct_ew(1,1:yi,xi) = ctns(1,1:yi,1) - ctns(1,1:yi,xi); 
        
end

% for north-south links

delta_p_ns(1,1:yi-1,1:xi) = pns_pa(1,2:yi,1:xi) - pns_pa(1,1:yi-1,1:xi);
delta_ct_ns(1,1:yi-1,1:xi) = ctns(1,2:yi,1:xi) - ctns(1,1:yi-1,1:xi);

% calculate locally referenced specific volume anomaly

v_all = 1./rho_from_ct(s,ct,p);
v_surf = 1./rho_from_ct(sns,ctns,pns);

% calculate steric anomaly on density surface of eastern point

v_tilde_e(1,1:yi,1:xi-1) = 1./rho_from_ct(sns(1,1:yi,2:xi),ctns(1,1:yi,2:xi),pns(1,1:yi,2:xi));
v_2_tilde_e(1,1:yi,1:xi-1) = 1./rho_from_ct(sns(1,1:yi,1:xi-1),ctns(1,1:yi,1:xi-1),pns(1,1:yi,2:xi));

v_hat_e(1,1:yi,1:xi-1) = 1./rho_from_ct(sns(1,1:yi,1:xi-1),ctns(1,1:yi,1:xi-1),pns(1,1:yi,1:xi-1));
v_2_hat_e(1,1:yi,1:xi-1) = 1./rho_from_ct(sns(1,1:yi,2:xi),ctns(1,1:yi,2:xi),pns(1,1:yi,1:xi-1));

switch wrap
    
    case 'long'
        
        v_tilde_e(1,1:yi,xi) = 1./rho_from_ct(sns(1,1:yi,1),ctns(1,1:yi,1),pns(1,1:yi,1));
        v_2_tilde_e(1,1:yi,xi) = 1./rho_from_ct(sns(1,1:yi,xi),ctns(1,1:yi,xi),pns(1,1:yi,1));
        
        v_hat_e(1,1:yi,xi) = 1./rho_from_ct(sns(1,1:yi,xi),ctns(1,1:yi,xi),pns(1,1:yi,xi));
        v_2_hat_e(1,1:yi,xi) = 1./rho_from_ct(sns(1,1:yi,1),ctns(1,1:yi,1),pns(1,1:yi,xi));
        
end
        
steric_tilde_e = v_tilde_e - v_2_tilde_e;
steric_hat_e = v_hat_e - v_2_hat_e;

% calculate steric anomaly on density surface of northern point

v_tilde_n(1,1:yi-1,1:xi) = 1./rho_from_ct(sns(1,2:yi,1:xi),ctns(1,2:yi,1:xi),pns(1,2:yi,1:xi));
v_2_tilde_n(1,1:yi-1,1:xi) = 1./rho_from_ct(sns(1,1:yi-1,1:xi),ctns(1,1:yi-1,1:xi),pns(1,2:yi,1:xi));

steric_tilde_n = v_tilde_n - v_2_tilde_n;

v_hat_n(1,1:yi-1,1:xi) = 1./rho_from_ct(sns(1,1:yi-1,1:xi),ctns(1,1:yi-1,1:xi),pns(1,1:yi-1,1:xi));
v_2_hat_n(1,1:yi-1,1:xi) = 1./rho_from_ct(sns(1,2:yi,1:xi),ctns(1,2:yi,1:xi),pns(1,1:yi-1,1:xi));

steric_hat_n = v_hat_n - v_2_hat_n;

% calculate term1

term1_e = (1/4) .* delta_p_ew .* (steric_tilde_e + steric_hat_e);
term1_n = (1/4) .* delta_p_ns .* (steric_tilde_n + steric_hat_n);

%% calculate term 2 

% calculate thermobaric parameter

[alpha,beta,aonb,T_b] = ab_from_ct(sns,ctns,pns);
T_b_pa = T_b ./ 10000;

% calculate term 2

term2_e = (1/12) .* (T_b_pa .* v_surf) .* delta_ct_ew .* delta_p_ew .* delta_p_ew;
term2_n = (1/12) .* (T_b_pa .* v_surf) .* delta_ct_ns .* delta_p_ns .* delta_p_ns;

%% calculate term 3 

% calculate steric anomaly on casts above density surface

sref = repmat(sns,[zi 1 1]);
ctref = repmat(ctns,[zi 1 1]);

v_cast = v_all;
v_2_cast = 1./rho_from_ct(sref,ctref,p);
steric_cast = v_cast - v_2_cast;

% for eastern points

v_cast_e(1:zi,1:yi,1:xi-1) = v_all(1:zi,1:yi,2:xi);
v_2_cast_e(1:zi,1:yi,1:xi-1) = 1./rho_from_ct(sref(1:zi,1:yi,1:xi-1),ctref(1:zi,1:yi,1:xi-1),p(1:zi,1:yi,2:xi));

switch wrap
    
    case 'long'
        
        v_cast_e(1:zi,1:yi,xi) = v_all(1:zi,1:yi,1);
        v_2_cast_e(1:zi,1:yi,xi) = 1./rho_from_ct(sref(1:zi,1:yi,xi),ctref(1:zi,1:yi,xi),p(1:zi,1:yi,1));

end
        
steric_cast_e = v_cast_e - v_2_cast_e;
        
% for northern points

v_cast_n(1:zi,1:yi-1,1:xi) = v_all(1:zi,2:yi,1:xi);
v_2_cast_n(1:zi,1:yi-1,1:xi) = 1./rho_from_ct(sref(1:zi,1:yi-1,1:xi),ctref(1:zi,1:yi-1,1:xi),p(1:zi,2:yi,1:xi));
steric_cast_n(1:zi,1:yi-1,1:xi) = v_cast_n(1:zi,1:yi-1,1:xi) - v_2_cast_n(1:zi,1:yi-1,1:xi);

% calculate pressure integrals of steric anomaly

% for central point (western or southern point depending on link)

steric_cast_mid  = 0.5 * (steric_cast(2:end,:,:) + steric_cast(1:end-1,:,:));
top = steric_cast(1,:,:).*p_pa(1,:,:);

delta_steric = (steric_cast_mid .* diff(p_pa,1));
steric_int  = cumsum([top;delta_steric]);

p_mid_pa = 0.5 * (p_pa(2:end,:,:) + p_pa(1:end-1,:,:));
pp = [p_pa(1,:,:);p_mid_pa];
steric_int_surf = var_on_surf(pns_pa(1,:,:),pp,steric_int);

steric_int_surf = change(steric_int_surf,'==',0,nan);

% for eastern point
 
top = nan(size(pns));

steric_cast_e_mid  = 0.5 * (steric_cast_e(2:end,:,:) + steric_cast_e(1:end-1,:,:));
top(1,1:yi,1:xi-1) = steric_cast_e(1,1:yi,1:xi-1).*p_pa(1,1:yi,2:xi);
delta_e_steric(:,1:yi,1:xi-1) = (steric_cast_e_mid(:,1:yi,1:xi-1) .* diff(p_pa(:,1:yi,2:xi),1));

switch wrap
    
    case 'long'
        
        top(1,1:yi,xi) = steric_cast_e(1,1:yi,xi).*p_pa(1,1:yi,1);
        delta_e_steric(:,1:yi,xi) = (steric_cast_e_mid(:,1:yi,xi) .* diff(p_pa(:,1:yi,1),1));
        
end


steric_e_int  = cumsum([top;delta_e_steric]);

steric_int_e_surf(1,1:yi,1:xi-1) = var_on_surf(pns_pa(1,1:yi,2:xi),pp(:,1:yi,2:xi),steric_e_int);

switch wrap
    
    case 'long'
        
        steric_int_e_surf(1,1:yi,xi) = var_on_surf(pns_pa(1,1:yi,1),pp(:,1:yi,1),steric_e_int(:,1:yi,xi));
        
end
        

term3_e = change(steric_int_e_surf,'==',0,nan);

% for northern point

top = nan(size(pns));

steric_cast_n_mid  = 0.5 * (steric_cast_n(2:end,:,:) + steric_cast_n(1:end-1,:,:));
top(1,1:yi-1,1:xi) = steric_cast_n(1,1:yi-1,1:xi).*p_pa(1,2:yi,1:xi);

delta_n_steric(:,1:yi-1,1:xi) = (steric_cast_n_mid(:,1:yi-1,1:xi) .* diff(p_pa(:,2:yi,1:xi),1));
steric_n_int  = cumsum([top;delta_n_steric]);

steric_int_n_surf(1,1:yi-1,1:xi) = var_on_surf(pns_pa(1,2:yi,1:xi),pp(:,2:yi,1:xi),steric_n_int);

term3_n = change(steric_int_n_surf,'==',0,nan);
 
%% calculate streamfuction difference 

delta_stream_ew = term1_e - term2_e - term3_e + steric_int_surf;
delta_stream_ns = term1_n - term2_n - term3_n + steric_int_surf;
        