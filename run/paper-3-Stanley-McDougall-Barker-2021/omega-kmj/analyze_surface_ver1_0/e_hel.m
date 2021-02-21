function [e_hel,e_hel_x,e_hel_y] = e_hel(slope_x,slope_y,uns,vns,grid,wrap)

%           Calculate the diapycnal velocity due to neutral helicity
%
% Usage:    [e_hel,e_hel_x,e_hel_y] = e_hel(slope_x,slope_y,uns,vns,grid,wrap) 
%
%           Calculate the diapycnal velocity (v*s) due to neutral helicity 
%           (v is lateral velocity, s is the slope-difference between a 
%           density surface and a neutral tangent plane). 
%
% Input:    slope_x           zonal component of slope error
%           slope_y           meridional component of slope error
%           uns               zonal component of lateral velocity
%           vns               meridional component of lateral velocity
%           grid              'bgrid' for Arakawa B-grid
%                             'cgrid' for Arakawa C-grid
%           wrap              'none'
%                             'long'
%
% Output:   e_hel             diapycnal velocity due to neutral helicity
%           e_hel_x           x-component of diapycnal velocity due
%                             to neutral helicity
%           e_hel_y           y-component of diapycnal velocity due
%                             to neutral helicity
%  
% Calls:    var_on_surf.m
%
% Units:    lateral velocities          m/s
%           e_hel                       m/s
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 6)
  error('e_hel.m: requires 6 input arguments')
end

%% initialize and preallocate memory

[gi,yi,xi] = size(uns);

e_hel = nan(gi,yi,xi);
e_hel_x = nan(gi,yi,xi);
e_hel_y = nan(gi,yi,xi);
uns_tracer = nan(gi,yi,xi);
vns_tracer = nan(gi,yi,xi);

%% calculate u,v on tracer-points;

switch grid

    case 'bgrid'

        switch wrap

            case {'none'}

                uns_tracer(1:gi,2:yi,2:xi) = 0.25 * (uns(1:gi,1:yi-1,1:xi-1) + uns(1:gi,1:yi-1,2:xi-1) + uns(1:gi,2:yi-1,1:xi-1) + uns(1:gi,2:yi-1,2:xi-1));
                vns_tracer(1:gi,2:yi,2:xi) = 0.25 * (vns(1:gi,1:yi-1,1:xi-1) + vns(1:gi,1:yi-1,2:xi-1) + vns(1:gi,2:yi-1,1:xi-1) + vns(1:gi,2:yi-1,2:xi-1));

            case {'long'}

                uns_tracer(1:gi,2:yi,1) = 0.25 * (uns(1:gi,1:yi-1,xi) + uns(1:gi,2:yi,xi) + uns(1:gi,1:yi-1,1) + uns(1:gi,2:yi,1));
                uns_tracer(1:gi,2:yi,2:xi) = 0.25 * (uns(1:gi,1:yi-1,1:xi-1) + uns(1:gi,1:yi-1,2:xi) + uns(1:gi,2:yi,1:xi-1) + uns(1:gi,2:yi,2:xi));
                vns_tracer(1:gi,2:yi,1) = 0.25 * (vns(1:gi,1:yi-1,xi) + vns(1:gi,2:yi,xi) + vns(1:gi,1:yi-1,1) + vns(1:gi,2:yi,1));
                vns_tracer(1:gi,2:yi,2:xi) = 0.25 * (vns(1:gi,1:yi-1,1:xi-1) + vns(1:gi,1:yi-1,2:xi) + vns(1:gi,2:yi,1:xi-1) + vns(1:gi,2:yi,2:xi));

        end

    case 'cgrid'

        switch wrap

            case 'none'

                uns_tracer(1:gi,1:yi,2:xi) = 0.5 * (uns(1:gi,1:yi,1:xi-1) + uns(1:gi,1:yi,2:xi));
                vns_tracer(1:gi,2:yi,1:xi) = 0.5 * (vns(1:gi,1:yi-1,1:xi) + vns(1:gi,2:yi,1:xi));

            case {'long'}

                uns_tracer(1:gi,1:yi,2:xi) = 0.5 * (uns(1:gi,1:yi,1:xi-1) + uns(1:gi,1:yi,2:xi));
                uns_tracer(1:gi,1:yi,1) = 0.5 * (uns(1:gi,1:yi,xi) + uns(1:gi,1:yi,1));
                vns_tracer(1:gi,2:yi,1:xi) = 0.5 * (vns(1:gi,1:yi-1,1:xi) + vns(1:gi,2:yi,1:xi));

        end

end

%% calculate diapycnal velocity due to neutral helicity

e_hel(1:gi,1:yi,1:xi) = (slope_x(1:gi,1:yi,1:xi) .* uns_tracer(1:gi,1:yi,1:xi)) + (slope_y(1:gi,1:yi,1:xi) .* vns_tracer(1:gi,1:yi,1:xi));
e_hel_x(1:gi,1:yi,1:xi) = (slope_x(1:gi,1:yi,1:xi) .* uns_tracer(1:gi,1:yi,1:xi));
e_hel_y(1:gi,1:yi,1:xi) = (slope_y(1:gi,1:yi,1:xi) .* vns_tracer(1:gi,1:yi,1:xi));