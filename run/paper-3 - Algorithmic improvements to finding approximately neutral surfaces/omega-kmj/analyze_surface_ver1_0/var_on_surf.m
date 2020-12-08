function [var1_ns,var2_ns] = var_on_surf(pns,p,var1,var2)

%           Vertically interpolate a variable onto a density surface
%
% Usage:    [var1_ns,var2_ns] = var_on_surf(pns,p,var1,var2)
%
%           Vertically interpolate one or two variable(s) onto a density surface
%
% Input:    var1          variable to be interpolated onto
%                         surface
%           var2          variable 2 to be interpolated onto
%                         surface (optional)
%           pns           pressure on density surface      
%           p             pressure 
%
% Output:   var1_ns       variable 1 interpolated onto density surface
%           var2_ns       variable 2 interpolated onto density surface
%                 
% Calls:           
%
% Units:    pressure    dbar
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 3) && ~(nargin == 4)
    error('var_surf.m: requires 3 or 4 input arguments')
end

if (nargin == 3) && (nargout == 2)
    error('var_surf.m: number of input and output variables is not the same')
elseif (nargin == 4) && (nargout == 1)
    error('var_surf.m: number of input and output variables is not the same')
end

%% initialize and preallocate memory

[zi,yi,xi] = size(p);

var1_ns = nan(1,yi,xi);
var2_ns = nan(1,yi,xi);

%% make pns 3-dimensional

pns_squeeze = squeeze(pns(1,:,:));
pns_rep = repmat(pns_squeeze,[1 1 zi]);
pns_rep = permute(pns_rep,[3 1 2]);

%% find zero-crossing of density surface

if xi ~= 1
    p_diff = p - pns_rep;
else
    p_diff = squeeze(p) - squeeze(pns_rep);
end

%% interpolation of values onto zero-crossing

if (nargin == 3) % if interpolating one variable
    
    for j = 1:yi
        for i = 1:xi
            if ~isnan(p_diff(1,j,i))
                inds = find(p_diff(:,j,i) >= 0);
                if (numel(inds) ~= 0)
                    inds_l = inds(1);
                    inds_u = inds(1) - 1;
                    if (inds_l ~= 1)
                        r = (pns(1,j,i)-p(inds_u,j,i))/(p(inds_l,j,i)-p(inds_u,j,i));
                        var1_ns(1,j,i) = var1(inds_u,j,i) + r*(var1(inds_l,j,i)-var1(inds_u,j,i));
                    else
                        var1_ns(1,j,i) = var1(1,j,i);
                    end
                end
            end
        end
    end
    
elseif (nargin == 4) % if interpolating a second variable
    
    for j = 1:yi
        for i = 1:xi
            if ~isnan(p_diff(1,j,i))
                inds = find(p_diff(:,j,i) >= 0);
                if (numel(inds) ~= 0)
                    inds_l = inds(1);
                    inds_u = inds(1) - 1;
                    if (inds_l ~= 1)
                        r = (pns(1,j,i)-p(inds_u,j,i))/(p(inds_l,j,i)-p(inds_u,j,i));
                        var1_ns(1,j,i) = var1(inds_u,j,i) + r*(var1(inds_l,j,i)-var1(inds_u,j,i));
                        var2_ns(1,j,i) = var2(inds_u,j,i) + r*(var2(inds_l,j,i)-var2(inds_u,j,i));
                    else
                        var1_ns(1,j,i) = var1(1,j,i);
                        var2_ns(1,j,i) = var2(1,j,i);
                    end
                end
            end
        end
    end
    
end