function [opt_streamfunc] = optimize_streamfunc(s,ct,p,sns,ctns,pns,wrap)

%           Optimize streamfunction on density surfaces 
%
% Usage:    [stream_func] = optimize_streamfunc(s,ct,p,sns,ctns,pns,wrap,density,pref)
%
%           Optimize streamfunction on density surfaces according to 
%           'An approximate geostrophic streamfunction for use 
%           on density surfaces, McDougall, Trevor J. and A. Klocker, 
%           Journal of Physical Oceanography, in prep.'
%
% Input:    s           salinity
%           ct          conservative temperature
%           p           pressure
%           sns         salinity on initial density surface
%           ctns        conservative temperature on initial density
%                       surface
%           pns         pressure on initial density surface
%           wrap        'none'
%                       'long'
%
% Output:   streamfunc          streamfunction on density surface
%           streamfunc_x        x-component of streamfunction
%           streamfunc_y        y-component of streamfunction
%
% Calls:    delta_stream_func.m, dyn_height.m, grad_surf.m
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
    error('optimize_stream_func.m: requires 7 input arguments')
end

%% get size of 3-dim hydrography

[zi,yi,xi] = size(s);

%% calculate differrence in neutral surface geostrophic streamfunctions

[delta_stream_ew,delta_stream_ns] = delta_streamfunc(s,ct,p,sns,ctns,pns,wrap);

%% prepare data

pns_squeeze = squeeze(pns);
delta_stream_ew_squeeze = squeeze(delta_stream_ew);
delta_stream_ns_squeeze = squeeze(delta_stream_ns);

%% preallocate memory

opt_streamfunc = nan(1,yi,xi);
streamfunc_tmp = nan(yi,xi);

b = zeros(3*xi*yi,1,1);
s1 = zeros(3*xi*yi,1);
s2 = zeros(3*xi*yi,1);
s3 = zeros(3*xi*yi,1);

%% find connected regions in the ocean

% find grid points which are not continents

inds = find(~isnan(pns_squeeze(:)));

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

%% find connected regions and solve least-squares problem for each 

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

    delta_stream_ew_inds = delta_stream_ew_squeeze(region);
    delta_stream_ns_inds = delta_stream_ns_squeeze(region);

    % build a matrix where all the points of the region are labelled
    % with 1,2,....length(region), other regions and continents are
    % filled with nans

    ng_region = nan(yi,xi);
    ng_region(region) = 1:length(region);

    % set up east-west equations for weighted inversion

    neq = 0; % set equation number to 0
    ieq = 0;

    for i = 1:length(region)

        [jj,ii] = ind2sub([yi,xi],region(i));

        switch wrap

            case {'none'}

                if (ii+1 <= xi) && (~isempty(find(region == ng(jj,ii)))) && (~isempty(find(region == ng(jj,ii+1)))) && (~isnan(delta_stream_ew_inds(ng_region(jj,ii)))) %#ok

                    neq = neq + 1;
                    ieq = ieq + 1;
                    s1(ieq) = neq;
                    s2(ieq) = ng_region(jj,ii);
                    s3(ieq) = -1;
                    ieq = ieq+1;
                    s1(ieq) = neq;
                    s2(ieq) = ng_region(jj,ii+1);
                    s3(ieq) = 1;
                    b(neq,1) = delta_stream_ew_inds(ng_region(jj,ii));

                end

            case {'long'}

                if (ii+1 <= xi)

                    if (~isempty(find(region == ng(jj,ii)))) && (~isempty(find(region == ng(jj,ii+1)))) && (~isnan(delta_stream_ew_inds(ng_region(jj,ii)))) %#ok

                        neq = neq + 1;
                        ieq = ieq + 1;
                        s1(ieq) = neq;
                        s2(ieq) = ng_region(jj,ii);
                        s3(ieq) = -1;
                        ieq = ieq+1;
                        s1(ieq) = neq;
                        s2(ieq) = ng_region(jj,ii+1);
                        s3(ieq) = 1;
                        b(neq,1) = delta_stream_ew_inds(ng_region(jj,ii));

                    end

                elseif (ii+1 > xi)

                    if (~isempty(find(region == ng(jj,ii))) && ~isempty(find(region == ng(jj,1)))) && (~isnan(delta_stream_ew_inds(ng_region(jj,ii)))) %#ok

                        neq = neq + 1;
                        ieq = ieq + 1;
                        s1(ieq) = neq;
                        s2(ieq) = ng_region(jj,ii);
                        s3(ieq) = -1;
                        ieq = ieq+1;
                        s1(ieq) = neq;
                        s2(ieq) = ng_region(jj,1);
                        s3(ieq) = 1;
                        b(neq,1) = delta_stream_ew_inds(ng_region(jj,ii));

                    end

                end
        end
    end

    % set up north-south equations for weighted inversion

    for i = 1:length(region)

        [jj,ii] = ind2sub([yi,xi],region(i));

        if (jj+1 <= yi) && (~isempty(find(region == ng(jj,ii)))) && ~isempty(find(region == ng(jj+1,ii))) && (~isnan(delta_stream_ns_inds(ng_region(jj,ii)))) %#ok

            neq = neq + 1;
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ng_region(jj,ii);
            s3(ieq) = -1;
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ng_region(jj+1,ii);
            s3(ieq) = 1;
            b(neq,1) = delta_stream_ns_inds(ng_region(jj,ii));
        end
    end

    % constrain the average of the streamfunction to be zero

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

    % put streamfunction field calculated witht the least-squares solver into
    % their appropriate position in the matrix

    streamfunc_tmp(region) = x;

    eval(['region_',int2str(nregion),' = region;']);

end

disp(['calculating ',int2str(nregion),' regions']);

%% calculate streamfunction and its gradients

opt_streamfunc(1,1:yi,1:xi) = streamfunc_tmp;