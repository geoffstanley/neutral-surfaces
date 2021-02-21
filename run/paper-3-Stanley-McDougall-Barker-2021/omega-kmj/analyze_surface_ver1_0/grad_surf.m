function [grad_x,grad_y] = grad_surf(var_surf,e1t,e2t,keyword,wrap)

%           Find the gradient of a variable in x and y-direction.
%
% Usage:    [grad_x,grad_y] = grad_surf(var_surf,e1t,e2t,keyword,wrap)
%
%
% Input:    var                   variable 
%           e1t                   zonal scale-factor at var_surf points 
%           e2t                   meridional scale-factor at var_surf 
%                                 points
%           keyword               'op' for slopes on var_surf points
%                                 'bp' for slopes between var_surf
%                                 points
%           wrap                  'none'
%                                 'long'
%
% Output:   grad_x                gradient in zonal direction 
%           grad_y                gradient in meridional direction 
%
% Units:    scale factors         m         
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 5)
    error('grad_surf.m: requires 5 input arguments')
end

%% calculate gradients

test = length(size(var_surf));

if  (test == 3) % time-independent data

    [gi,yi,xi] = size(var_surf);
    
    grad_x = nan(gi,yi,xi);
    grad_y = nan(gi,yi,xi);

    e1t = repmat(e1t,[1 1 gi]);
    e1t = permute(e1t,[3 1 2]);

    e2t = repmat(e2t,[1 1 gi]);
    e2t = permute(e2t,[3 1 2]);

    % calculate slope in longitude direction

    grad_x_tmp(1:gi,1:yi,1:xi-1) = (var_surf(1:gi,1:yi,2:xi) - var_surf(1:gi,1:yi,1:xi-1)) ./ e1t(1:gi,1:yi,1:xi-1);

    switch wrap
        
        case {'none'}
            
            grad_x_tmp(1:gi,1:yi,xi) = nan;
            
        case {'long'}
            
            grad_x_tmp(1:gi,1:yi,xi) = (var_surf(1:gi,1:yi,1) - var_surf(1:gi,1:yi,xi)) ./ e1t(1:gi,1:yi,1);
            
    end

    % calculate slope in latitude direction

     grad_y_tmp(1:gi,1:yi-1,1:xi) = (var_surf(1:gi,2:yi,1:xi) - var_surf(1:gi,1:yi-1,1:xi)) ./ e2t(1:gi,1:yi-1,1:xi);
     grad_y_tmp(1:gi,yi,1:xi) = nan;
 
     switch keyword

         case 'op' % averages slopes onto var_surf points

             switch wrap

                 case 'none'

                     grad_x(1:gi,1:yi,2:xi-1) = 0.5 * (grad_x_tmp(1:gi,1:yi,1:xi-2) + grad_x_tmp(1:gi,1:yi,2:xi-1));
                     grad_y(1:gi,2:yi-1,1:xi) = 0.5 * (grad_y_tmp(1:gi,1:yi-2,1:xi) + grad_y_tmp(1:gi,2:yi-1,1:xi));
                     
                 case 'long'

                     grad_x(1:gi,1:yi,1) = 0.5 * (grad_x_tmp(1:gi,1:yi,xi) + grad_x_tmp(1:gi,1:yi,1));
                     grad_x(1:gi,1:yi,2:xi) = 0.5 * (grad_x_tmp(1:gi,1:yi,1:xi-1) + grad_x_tmp(1:gi,1:yi,2:xi));
                     grad_y(1:gi,2:yi-1,1:xi) = 0.5 * (grad_y_tmp(1:gi,1:yi-2,1:xi) + grad_y_tmp(1:gi,2:yi-1,1:xi));

             end

         otherwise

             grad_x = grad_x_tmp;
             grad_y = grad_y_tmp;

     end
    
else % time-dependent data

    [ti,gi,yi,xi] = size(var_surf);

    grad_x = nan(1:ti,1:gi,1:yi,1:xi);
    grad_y = nan(1:ti,1:gi,1:yi,1:xi);
    
    e1t = repmat(e1t,[1 1 gi ti]);
    e1t = permute(e1t,[4 3 1 2]);

    e2t = repmat(e2t,[1 1 gi ti]);
    e2t = permute(e2t,[4 3 1 2]);

    % calculate slope in longitude direction

    grad_x_tmp(1:ti,1:gi,1:yi,1:xi-1) = (var_surf(1:ti,1:gi,1:yi,2:xi) - var_surf(1:ti,1:gi,1:yi,1:xi-1)) ./ e1t(1:ti,1:gi,1:yi,1:xi-1);

    switch wrap
        
        case {'long'}
            
            grad_x_tmp(1:ti,1:gi,1:yi,xi) = (var_surf(1:ti,1:gi,1:yi,1) - var_surf(1:ti,1:gi,1:yi,xi)) ./ e1t(1:ti,1:gi,1:yi,1);
            
    end

    % calculate slope in latitude direction

    grad_y_tmp(1:ti,1:gi,1:yi-1,1:xi) = (var_surf(1:ti,1:gi,2:yi,1:xi) - var_surf(1:ti,1:gi,1:yi-1,1:xi)) ./ e2t(1:ti,1:gi,1:yi-1,1:xi);
    grad_y_tmp(1:ti,1:gi,yi,1:xi) = nan;

    switch keyword
        
        case 'op' % averages slopes onto var_surf points
            
            switch wrap
                
                case 'none'
                    
                    grad_x(1:ti,1:gi,1:yi,2:xi-1) = 0.5 * (grad_x_tmp(1:ti,1:gi,1:yi,1:xi-2) + grad_x_tmp(1:ti,1:gi,1:yi,2:xi-1));
                    grad_y(1:ti,1:gi,2:yi-1,1:xi) = 0.5 * (grad_y_tmp(1:ti,1:gi,1:yi-2,1:xi) + grad_y_tmp(1:ti,1:gi,2:yi-1,1:xi));
                    
                case 'long'
                    
                    grad_x(1:ti,1:gi,1:yi,1) = 0.5 * (grad_x_tmp(1:ti,1:gi,1:yi,xi) + grad_x_tmp(1:ti,1:gi,1:yi,1));
                    grad_x(1:ti,1:gi,1:yi,2:xi) = 0.5 * (grad_x_tmp(1:ti,1:gi,1:yi,1:xi-1) + grad_x_tmp(1:ti,1:gi,1:yi,2:xi));
                    grad_y(1:ti,1:gi,2:yi-1,1:xi) = 0.5 * (grad_y_tmp(1:ti,1:gi,1:yi-2,1:xi) + grad_y_tmp(1:ti,1:gi,2:yi-1,1:xi));
                    
            end
            
        otherwise
            
            grad_x = grad_x_tmp;
            grad_y = grad_y_tmp;
            
    end
    
end