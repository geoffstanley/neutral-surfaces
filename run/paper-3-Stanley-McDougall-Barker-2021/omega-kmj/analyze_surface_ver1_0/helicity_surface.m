function [hel] = helicity_surface(sns,ctns,pns,n2_ns,g,e1t,e2t,wrap)

%           Calculate neutral helicity on a density surface
%
% Usage:    [hel] = helicity_surface(sns,ctns,pns,n2_ns,g,e1t,e2t,wrap)
%
%           Calculate neutral helicity on a density surface according to 'The
%           Thinness of the ocean in S-Theta-p Space and the Implications 
%           for Mean Diapycnal Advection, Trevor J McDougall and David R 
%           Jackett,  Journal of Physical Oceanography, 2007'. 
%
%
% Input:    sns         salinity on density surface
%           ctns        conservative temperature on density surface
%           pns         pressure on density surface
%           n2_ns       buoyancy frequency on density surface
%           g           gravitational acceleration
%           e1t         zonal scale-factor at tracer points 
%           e2t         meridional scale-factor at tracer points
%           wrap        'none'
%                       'long'
%
% Output:   hel         neutral helicity
%  
% Calls:    grad_surf.m, ab_from_ct.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           buoyancy frequency          s^-1
%           gravitational acceleration  m/s^2
%           scale factors               m
%           neutral helicity            m^-3                      
%           
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 8)
  error('helicity_surface.m: requires 8 input arguments')
end 

%% calculate thermobaric parameter

[alpha,beta,aonb,T_b] = ab_from_ct(sns,ctns,pns);

%% calculate gradients of ct and p on density surfaces

[slope_ctx,slope_cty] = grad_surf(ctns,e1t,e2t,'op',wrap);
[slope_px,slope_py] = grad_surf(pns,e1t,e2t,'op',wrap);

%% change size of g

g_tmp = nan(size(n2_ns));
g_tmp(1,:,:) = g;

%% calculate neutral helicity on density surface

hel(1,:,:) = (n2_ns ./ g_tmp) .* T_b .* (slope_px .* slope_cty - slope_py .* slope_ctx);