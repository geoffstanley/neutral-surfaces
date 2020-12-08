function [geo_vel_x,geo_vel_y] = geo_vel(streamfunc,s,ct,p,sref,ctref,lats,e1t,e2t,p_dynh,wrap)

%           Calculate geostrophic velocities
%
% Usage:    [geo_vel_x,geo_vel_y] = geo_vel(streamfunc,s,ct,p,lats,e1t,e2t,pref,wrap)
%
%           Calculate geostrophic velocities from a geostrophic streamfunction
%           reltaive to the dynamic height at a specified pressure.
%
% Input:    streamfunc        geostrophic streamfunction
%           s                 salinity
%           ct                conservative temperature
%           p                 pressure
%           sref              reference salinity for steric anomaly
%                             calculation
%           ctref             reference conservative temperature for steric 
%                             anomaly calculation
%           lats              latitude
%           e1t               zonal scale factor
%           e2t               meridional scale factor
%           p_dynh            pressure of reference level
%           wrap              'none'
%                             'long' 
%                 
% Output:   geo_vel_x         geostrophic velocity in zonal direction
%           geo_vel_y         geostrophic velocity in meridional direction
% 
% Calls:    grad_surf.m, dyn_height.m, coriolis.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           scale factors               m
%           geo_vel                     m/s
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 11)
    error('geo_vel.m: requires 11 input arguments')
end

[zi,yi,xi] = size(streamfunc);

%% calculate gradients of streamfunction

[grad_streamfunc_x,grad_streamfunc_y] = grad_surf(streamfunc,e1t,e2t,'bp',wrap);

%% calculate steric height and its gradients  

[dynh] = dyn_height(s,ct,p,sref,ctref,p_dynh);

[grad_streamfunc_ref_x,grad_streamfunc_ref_y] = grad_surf(dynh,e1t,e2t,'bp',wrap);

%% calculating relative gradients

grad_streamfunc_rel_x = grad_streamfunc_x - grad_streamfunc_ref_x;
grad_streamfunc_rel_y = grad_streamfunc_y - grad_streamfunc_ref_y;

%% calculate coriolis parameter

f_tmp = coriolis(lats);
f = nan(1,yi,xi);
f(1,1:yi,1:xi) = f_tmp;

%% calculate geostrophic velocities

geo_vel_x = -grad_streamfunc_rel_y ./ f; 
geo_vel_y = grad_streamfunc_rel_x ./ f; 

%% disregard geostrophic velocities +/- 2 deg from the equator

for j = 1:yi
    for i = 1:xi
        if ((lats(j,i) <= 2) && (lats(j,i) >= -2))
            geo_vel_x(1,j,i) = nan;
            geo_vel_y(1,j,i) = nan;
        end
    end
end