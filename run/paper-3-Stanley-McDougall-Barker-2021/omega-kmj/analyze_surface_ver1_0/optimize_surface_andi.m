function [sns_i,ctns_i,pns_i] = optimize_surface(s,ct,p,g,n2,sns,ctns,pns,e1t,e2t,nit,choice,wrap)

%           Optimize density surfaces to minimise the fictitious diapycnal diffusivity
%
% Usage:    [sns_i,ctns_i,pns_i] =
%           Optimize_surface(s,ct,p,g,n2,sns,ctns,pns,e1t,e2t,nit,choice,wrap)
%
%           Optimize density surface using an iterative method to minimise
%           fictitious diapycnal diffusivity according 'A new method of
%           forming approximately neutral surfaces, Klocker et al., Ocean
%           Science, 5, 155-172, 2009.'
%
%
% Input:    s           salinity
%           ct          conservative temperature
%           p           pressure
%           g           gravitational acceleration
%           n2          buoyancy frequency
%           sns         salinity on initial density surface
%           ctns        conservative temperature on initial density
%                       surface
%           pns         pressure on initial density surface
%           e1t         grid length in meters in x-direction
%           e2t         grid length in meters in y-direction
%           nit         number of iterations
%           choice      choice between
%                       's' slope error
%                       'epsilon' density gradient error
%           wrap        'none'
%                       'long'
%
% Output:   sns_i       salinity on optimized surface
%           ctns_i      conservative temperature on optimized surface
%           pns_i       pressure on optimized surface
%
% Calls:    cut_off.m, mld.m, slope_error.m, var_on_surf.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           gravitational acceleration  m/s^2
%           buoyancy frequency          s^-1
%           scale factors               m
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% get size of 3-dim hydrography

[zi,yi,xi] = size(s);

%% calculate buoyancy frequency on density surface

p_mid = (p(2:zi,:,:) + p(1:zi-1,:,:)) ./ 2;
n2_ns = var_on_surf(pns,p_mid,n2);

%% calculate slope errors/density gradient errors on initial density surface

[ss,sx,sy,curl_s,ee,ex,ey,curl_e,ver] = slope_error(p,g,n2,sns,ctns,pns,e1t,e2t,'bp',wrap); %#ok

%% choice between minimizing slope errors or density gradient errors

switch choice
    case 's'
        xx = sx;
        yy = sy;
    case 'epsilon'
        xx = ex;
        yy = ey;
end

%% disregard data above mixed layer depth

cut_off_choice(1,1:yi,1:xi) = mld(s,ct,p);

[pns] = cut_off(pns,pns,cut_off_choice);
[n2_ns] = cut_off(n2_ns,pns,cut_off_choice);

%% prepare data

slope_square = nan(nit,1);

pns_squeeze = squeeze(pns);
n2_ns_squeeze = squeeze(n2_ns);
xx_squeeze = squeeze(xx);
yy_squeeze = squeeze(yy);

%% iterations of inversion

for it = 1:nit

    disp(['iteration nr.',int2str(it)]); % print number of iteration

    if (it > 1)

        % disregard data above mixed layer depth

        [pns_i] = cut_off(pns_i,pns_i,cut_off_choice);
        [n2_ns_i] = cut_off(n2_ns_i,pns_i,cut_off_choice);

        % prepare data for next iteration

        pns_squeeze = squeeze(pns_i);
        n2_ns_squeeze = squeeze(n2_ns_i);

        switch choice
            case 'epsilon'
                [ex_i] = cut_off(ex_i,pns_i,cut_off_choice);
                [ey_i] = cut_off(ey_i,pns_i,cut_off_choice);
                xx_squeeze = squeeze(ex_i);
                yy_squeeze = squeeze(ey_i);
            case 's'
                [sx_i] = cut_off(sx_i,pns_i,cut_off_choice);
                [sy_i] = cut_off(sy_i,pns_i,cut_off_choice);
                xx_squeeze = squeeze(sx_i);
                yy_squeeze = squeeze(sy_i);
        end
    end

    [yi,xi] = size(pns_squeeze); % find dimensions in lats and longs

    % preallocate memory

    depth_change = nan(yi,xi);
    depth_change_e = nan(yi,xi);

    b = zeros(3*xi*yi,1,1);
    s1 = zeros(3*xi*yi,1);
    s2 = zeros(3*xi*yi,1);
    s3 = zeros(3*xi*yi,1);

    % find grid points which are not continents

    inds = find(~isnan(n2_ns_squeeze(:)));

    % build matrix where the ocean gridpoints have their indices and
    % continet gridpoints nans

    ng = nan(yi,xi);
    ng(inds) = inds;

    % select regions

    % build neighbour matrix -> build a matrix which is like a look-up
    % table to see which gridpoints communicate with each other

    neighbour = nan(5,length(inds));

    switch wrap

        case 'none'

            for i = 1:length(inds)

                [jj,ii] = ind2sub([yi,xi],inds(i));

                if (jj == 1) && (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = nan;
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (jj == yi) && (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = nan;
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == 1) && (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = nan;
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (jj == yi) && (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = nan;
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = nan;
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == yi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = nan;
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                else
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                end
            end

        case 'long'

            for i = 1:length(inds)

                [jj,ii] = ind2sub([yi,xi],inds(i));

                if (jj == 1) && (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,xi);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (jj == yi) && (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == 1) && (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (jj == yi) && (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,xi);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,xi);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == yi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                else
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                end
            end
    end

    nregion = 0;
    region_matrix = nan(yi,xi);

    % find independent regions -> a least-squares problem is solved for
    % each of these regions

   while (length(find(~isnan(region_matrix))) ~= length(neighbour))

        pos = 0;
        region = [];

        nregion = nregion + 1;

        % find starting point of region

        for i = 1:length(inds)

            if ~isnan(neighbour(1,i)) && isnan(region_matrix(neighbour(1,i)))
                pos = pos+1;
                region(pos) = neighbour(1,i);
                region_matrix(neighbour(1,i)) = nregion;
                neighbour(1,i) = nan;

                if ~isnan(neighbour(2,i))
                    if isnan(region_matrix(neighbour(2,i)))
                        pos = pos+1;
                        region(pos) = neighbour(2,i);
                        region_matrix(neighbour(2,i)) = nregion;
                        neighbour(2,i) = nan;
                    end
                end

                if ~isnan(neighbour(3,i))
                    if isnan(region_matrix(neighbour(3,i)))
                        pos = pos+1;
                        region(pos) = neighbour(3,i);
                        region_matrix(neighbour(3,i)) = nregion;
                        neighbour(3,i) = nan;
                    end
                end

                if ~isnan(neighbour(4,i))
                    if isnan(region_matrix(neighbour(4,i)))
                        pos = pos+1;
                        region(pos) = neighbour(4,i);
                        region_matrix(neighbour(4,i)) = nregion;
                        neighbour(4,i) = nan;
                    end
                end

                if ~isnan(neighbour(5,i))
                    if isnan(region_matrix(neighbour(5,i)))
                        pos = pos+1;
                        region(pos) = neighbour(5,i);
                        region_matrix(neighbour(5,i)) = nregion;
                        neighbour(5,i) = nan;
                    end
                end

                break

            end
        end

        region_old = [];

        while(length(region) ~= length(region_old))

            region_old = region;

            % find all other points of region

            for i = 1:length(inds)

                if ~isnan(neighbour(1,i))
                    
                    if (region_matrix(neighbour(1,i)) == nregion) %#ok

                        if ~isnan(neighbour(1,i))
                            if isnan(region_matrix(neighbour(1,i)))
                                pos = pos+1;
                                region(pos) = neighbour(1,i);
                                region_matrix(neighbour(1,i)) = nregion;
                                neighbour(1,i) = nan;
                            end
                        end

                        if ~isnan(neighbour(2,i))
                            if isnan(region_matrix(neighbour(2,i)))
                                pos = pos+1;
                                region(pos) = neighbour(2,i);
                                region_matrix(neighbour(2,i)) = nregion;
                                neighbour(2,i) = nan;
                            end
                        end

                        if ~isnan(neighbour(3,i))
                            if isnan(region_matrix(neighbour(3,i)))
                                pos = pos+1;
                                region(pos) = neighbour(3,i);
                                region_matrix(neighbour(3,i)) = nregion;
                                neighbour(3,i) = nan;
                            end
                        end

                        if ~isnan(neighbour(4,i))
                            if isnan(region_matrix(neighbour(4,i)))
                                pos = pos+1;
                                region(pos) = neighbour(4,i);
                                region_matrix(neighbour(4,i)) = nregion;
                                neighbour(4,i) = nan;
                            end
                        end

                        if ~isnan(neighbour(5,i))
                            if isnan(region_matrix(neighbour(5,i)))
                                pos = pos+1;
                                region(pos) = neighbour(5,i);
                                region_matrix(neighbour(5,i)) = nregion;
                                neighbour(5,i) = nan;
                            end
                        end
                    end
                end
            end
        end

        % only use points of the region improved in this loop

        xx_inds = xx_squeeze(region);
        yy_inds = yy_squeeze(region);
        e1t_inds = e1t(region);
        e2t_inds = e2t(region);

        % build a matrix where all the points of the region are labelled
        % with 1,2,....length(region), other regions and continents are
        % filled with nans

        ng_region = nan(yi,xi);
        ng_region(region) = 1:length(region);

        %% set up east-west equations for weighted inversion

        neq = 0; % set equation number to 0
        ieq = 0;

        for i = 1:length(region)

            [jj,ii] = ind2sub([yi,xi],region(i));

            switch wrap

                case {'none'}

                    if (ii+1 <= xi) && (~isempty(find(region == ng(jj,ii)))) && (~isempty(find(region == ng(jj,ii+1)))) && (~isnan(xx_inds(ng_region(jj,ii)))) %#ok

                        neq = neq + 1;
                        ieq = ieq + 1;
                        s1(ieq) = neq;
                        s2(ieq) = ng_region(jj,ii);
                        s3(ieq) = -1;
                        ieq = ieq+1;
                        s1(ieq) = neq;
                        s2(ieq) = ng_region(jj,ii+1);
                        s3(ieq) = 1;
                        b(neq,1) = xx_inds(ng_region(jj,ii)) * e1t_inds(ng_region(jj,ii));

                    end

                case {'long'}

                    if (ii+1 <= xi)

                        if (~isempty(find(region == ng(jj,ii)))) && (~isempty(find(region == ng(jj,ii+1)))) && (~isnan(xx_inds(ng_region(jj,ii)))) %#ok

                            neq = neq + 1;
                            ieq = ieq + 1;
                            s1(ieq) = neq;
                            s2(ieq) = ng_region(jj,ii);
                            s3(ieq) = -1;
                            ieq = ieq+1;
                            s1(ieq) = neq;
                            s2(ieq) = ng_region(jj,ii+1);
                            s3(ieq) = 1;
                            b(neq,1) = xx_inds(ng_region(jj,ii)) * e1t_inds(ng_region(jj,ii));

                        end

                    elseif (ii+1 > xi)

                        if (~isempty(find(region == ng(jj,ii))) && ~isempty(find(region == ng(jj,1)))) && (~isnan(xx_inds(ng_region(jj,ii))))%#ok

                            neq = neq + 1;
                            ieq = ieq + 1;
                            s1(ieq) = neq;
                            s2(ieq) = ng_region(jj,ii);
                            s3(ieq) = -1;
                            ieq = ieq+1;
                            s1(ieq) = neq;
                            s2(ieq) = ng_region(jj,1);
                            s3(ieq) = 1;
                            b(neq,1) = xx_inds(ng_region(jj,ii)) * e1t_inds(ng_region(jj,ii));

                        end

                    end
            end
        end

        %% set up north-south equations for weighted inversion

        for i = 1:length(region)

            [jj,ii] = ind2sub([yi,xi],region(i));

            if (jj+1 <= yi) && (~isempty(find(region == ng(jj,ii)))) && ~isempty(find(region == ng(jj+1,ii))) && (~isnan(yy_inds(ng_region(jj,ii))))%#ok

                neq = neq + 1;
                ieq = ieq + 1;
                s1(ieq) = neq;
                s2(ieq) = ng_region(jj,ii);
                s3(ieq) = -1;
                ieq = ieq + 1;
                s1(ieq) = neq;
                s2(ieq) = ng_region(jj+1,ii);
                s3(ieq) = 1;
                b(neq,1) = yy_inds(ng_region(jj,ii)) * e2t_inds(ng_region(jj,ii));
            end
        end

        % make the average of all density changes zero -> this should keep
        % the surface from drifting away from the initial condition

        neq = neq + 1;

        for i = 1:length(region)
            [jj,ii] = ind2sub([yi,xi],region(i));
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ng_region(jj,ii);
            s3(ieq) = 1;
        end

        b(neq,1) = 0;

        % cut A and b to appropriate size

        s1 = s1(1:ieq);
        s2 = s2(1:ieq);
        s3 = s3(1:ieq);
        b = b(1:neq,1);

        % make matrix sparse and invert

        A = sparse(s1,s2,s3);
        b = sparse(b);

        disp(['solving for region ',int2str(nregion)]);
        
        x = lsqr(A,b,1e-6,10000);
           
        x = full(x)'; %#ok

        % put density changes calculated by the least-squares solver into
        % their appropriate position in the matrix

        switch choice
            case 's'
                depth_change(region) = - x;
            case 'epsilon'
                depth_change_e(region) = x;
        end

        eval(['region_',int2str(nregion),' = region;']);

    end

    disp(['calculating ',int2str(nregion),' regions']);

    % calculate density change into depth change

    switch choice
        case 'epsilon'
            if (it == 1)
                depth_change = (-9.81 * (depth_change_e)) ./ (3e-6 + squeeze(n2_ns));
            else
                depth_change = (-9.81 * (depth_change_e)) ./ (3e-6 + squeeze(n2_ns_i));
            end
    end

    if (it > 1)
        pns_old = pns_i;
    else
        pns_old = pns;
    end

    % calculate new depth of approximate neutral surface

    % damp damps the depth change - if the algorithm goes unstable decrease damp
    
    damp = 0.2;
    
    pns_i = nan(1,yi,xi);
    pns_i(1,:,:) = squeeze(pns_old) + damp .* depth_change;

    % calculate ct, s, slope errors and buoyancy frequency on new
    % approximate neutral surface

    [ctns_i,sns_i] = var_on_surf(pns_i,p,ct,s);
    n2_ns_i = var_on_surf(pns_i,p_mid,n2);

    [ss_i,sx_i,sy_i,curl_s_i,ee_i,ex_i,ey_i,curl_e_i] = ...
        slope_error(p,g,n2,sns_i,ctns_i,pns_i,e1t,e2t,'bp',wrap); %#ok

    % calculate epsilon^2 to estimate quality of approximate neutral
    % surface

    switch choice
        case 'epsilon'
            square = ex_i .* ex_i + ey_i .* ey_i;
            slope_square(it,1) = nansum(square(:));
        case 's'
            square = sx_i .* sx_i + sy_i .* sy_i;
            slope_square(it,1) = nansum(square(:));
    end

    % plot depth change and evolution of Veronis error

    if (mod(it,10) == 0) % plot every 10th iteration
        figure('Position',[10, 10, 1000, 1000])

        subplot(2,1,1)
        fpcolor(depth_change)
        xlabel('longitude','fontsize',30,'fontweight','bold')
        ylabel('latitude','fontsize',30,'fontweight','bold')
        clear depth_change
        colorbar
        hold on

        subplot(2,1,2)
        semilogy(slope_square)
        xlabel('iterations','fontsize',30,'fontweight','bold')

        switch choice
            case 'epsilon'
                ylabel('\epsilon ^2','fontsize',30,'fontweight','bold')
            case 's'
                ylabel('s^2','fontsize',30,'fontweight','bold')
        end

        hold on
        grid on
    end
end
