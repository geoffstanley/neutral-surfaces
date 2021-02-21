function [e2t,e1t] = scale_fac(lats,longs,wrap)

%           Find distances between gridpoints of given latitude/longitude 
%
% Usage:    [e2t,e1t] = scale_fac(lats,longs,wrap)
%
%           Find distances between gridpoints of given latitude/longitude 
%           in a 2-dim dataset 
%
% Input:    lats        latitude
%           longs       longitude
%           wrap        'none'
%                       'long'
%
% Output:   e2t         distance between gridpoints in N-S direction             
%           e1t         distance between gridpoints in E-W direction 
%
% Units:    e2t         m
%           e1t         m
%
% Calls:    sw_dist.m
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 3)
    error('scale_fac.m: requires 3 input arguments')
end

%% initialize and preallocate memory

[yi,xi] = size(lats);
e1t = nan(yi,xi);
e2t = nan(yi,xi);

%% calculate distances

switch wrap
    case 'none'
        for j = 1:yi
            for i = 1:xi-1
                e1t(j,i) = sw_dist([lats(j,i) lats(j,i)],[longs(j,i) longs(j,i+1)],'km');
            end
        end
        for j = 1:yi-1
            for i = 1:xi
                e2t(j,i) = sw_dist([lats(j,i) lats(j+1,i)],[longs(j,i) longs(j,i)],'km');
            end
        end
    case 'long'
        for j = 1:yi
            for i = 1:xi-1
                e1t(j,i) = sw_dist([lats(j,i) lats(j,i)],[longs(j,i) longs(j,i+1)],'km');
            end
        end
        for j = 1:yi
            e1t(j,xi) = sw_dist([lats(j,i) lats(j,i)],[longs(j,xi) longs(j,1)],'km');
        end
        for j = 1:yi-1
            for i = 1:xi
                e2t(j,i) = sw_dist([lats(j,i) lats(j+1,i)],[longs(j,i) longs(j,i)],'km');
            end
        end
end

%% convert distances from kilometers to meters

e1t = 1000 * e1t;
e2t = 1000 * e2t;