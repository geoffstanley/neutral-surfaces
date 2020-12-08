function [hel] = helicity_pressure(s,ct,p,e1t,e2t,wrap)

%           Calculate neutral helicity on pressure levels
%
% Usage:    [hel] = helicity_pressure(s,ct,p,e1t,e2t,wrap)
%
%           Calculate neutral helicity on a pressure level according to 'The
%           Thinness of the ocean in S-Theta-p Space and the Implications 
%           for Mean Diapycnal Advection, Trevor J McDougall and David R 
%           Jackett,  Journal of Physical Oceanography, 2007'. 
%
%
% Input:    s           salinity
%           ct          conservative temperature
%           p           pressure
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
%           scale factors               m
%           hel                         m^-3
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 6)
  error('helicity_pressure.m: requires 6 input arguments')
end 

%% calculate thermobaric parameter

[alpha,beta,aonb,T_b] = ab_from_ct(s,ct,p);

%% calculate gradient of s and ct on pressure levels

[slope_sx,slope_sy] = grad_surf(s,e1t,e2t,'op',wrap);
[slope_ctx,slope_cty] = grad_surf(ct,e1t,e2t,'op',wrap);

%% calculate neutral helicity

hel = -1 * beta .* T_b .* (slope_sx .* slope_cty - slope_sy .* slope_ctx);