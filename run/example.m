% Examples of topobaric surface, modified topobaric surface, and several
% geostrophic stream functions (gsf) including the topobaric gsf and the
% orthobaric Montgomery potential.

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2019  Geoff Stanley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com
% Version   : 2.1.0
%
% Modified by : --
% Date        : --
% Changes     : --

%% --- BEGIN SETUP --------------------------------------------------------
warning('off', 'MATLAB:nargchk:deprecated')
set(0, 'defaultfigurecolor', [1 1 1]); % white figure background
V = filesep(); % /  or  \  depending on OS.

PATH_LOCAL = [fileparts(mfilename('fullpath')) V]; % Get path to this file.
%PATH_LOCAL = '~/work/projects-gfd/neutral-surfaces/run/'; % Manually set path to this file.

% Get path to one directory up, containing all of Topobaric Surface
PATH_PROJECT = PATH_LOCAL(1 : find(PATH_LOCAL(1:end-1) == V, 1, 'last'));

% Add Neutral Surfaces to MATLAB's path
run([PATH_PROJECT 'ns_add_to_path.m']);

% Make a directory for figures
PATH_FIGS = [PATH_LOCAL 'figs' V];
if ~exist(PATH_FIGS, 'dir')
    mkdir(PATH_FIGS)
end

% Read path to ECCO2 data from PATH_ECCO2.txt
file_id = fopen([PATH_PROJECT 'lib' V 'dat' V 'PATH_ECCO2.txt'], 'rt');
PATH_ECCO2 = textscan(file_id, '%s');
PATH_ECCO2 = PATH_ECCO2{1}{1};
fclose(file_id);

%fileID = 1; % For standard output to the screen
fileID = fopen([PATH_LOCAL 'run ' datestr(now, 'yyyy mm dd hh MM ss') '.txt'], 'wt'); % For output to a file

db2Pa = 1e4; % dbar to Pa conversion
Pa2db = 1e-4; % Pa to dbar conversion

%% Load ECCO
TIMESTEP = '20021223';
g = load_ECCO2(PATH_ECCO2, 'grid');
ni = g.nx;
nj = g.ny;
nk = g.nz;
[S, T, ATMP, ETAN, SAP] = load_ECCO2(PATH_ECCO2, 'casts', TIMESTEP);
grav = g.grav;
rho_c = g.rho_c;
Z = -g.RC(:); % We are going to be Z > 0 people!
cori_V = g.cori_YG;

% Re-order data so water-columns are contiguous data:
S = permute(S, [3 1 2]); % [nz,nx,ny]. depth  by  long  by  lat
T = permute(T, [3 1 2]);

land = isnan(squeeze(S(1,:,:))); % true where water column is land, false where ocean

% C-grid functions for 2D data (longitude = i = rows, latitude = j = columns)
jm1 = @(F) subsasgn(circshift(F, [0 +1]), struct('type','()', 'subs', {{':',   1}}), nan); %#ok<SUBSASGN>

%% Set alias functions
% Choose the Boussinesq densjmd95 and set grav and rho_c in eos.m and eos_x.m
eoscg_set_bsq_param([PATH_PROJECT 'lib' V 'eos' V 'eoscg_densjmd95_bsq.m'   ] , [PATH_PROJECT 'eos.m'  ], grav, rho_c);
eoscg_set_bsq_param([PATH_PROJECT 'lib' V 'eos' V 'eoscg_densjmd95_bsq_dz.m'] , [PATH_PROJECT 'eos_x.m'], grav, rho_c);
eoscg_set_bsq_param([PATH_PROJECT 'lib' V 'eos' V 'eoscg_densjmd95_bsq_s_t.m'], [PATH_PROJECT 'eos_s_t.m'], grav, rho_c);

% Choose vertical interpolation method
interpfn = @ppc_linterp;

%% Pre-compute in-situ density and interpolants
R = eos(S, T, Z);

SppZ = interpfn(Z, S);
TppZ = interpfn(Z, T);


%% Select reference cast and depth
i0 = round(ni / 2); % reference cast zonal grid index
j0 = round(nj / 2); % reference cast meridional grid index
ztar = 2000; % Target depth for surfaces at (i0,j0), [m]

%% Prepare options, OPTS
OPTS.WRAP = [true, false];   % Periodic in longitude, not latitude
OPTS.VERBOSE = 1;            % Display modest info during computation
OPTS.X_TOL = 1e-4;           % error tolerance when updating the surface
OPTS.X_EXPN = 500;           % expansion of domain to search for solutions in each water column
OPTS.ITER_MAX = 3;           % Not too many iterations, for quick gratification
OPTS.ITER_START_WETTING = 1; % Start wetting right away
OPTS.REF_IJ = [i0, j0];      % Specify reference cast
OPTS.SppX = SppZ;            % Use pre-computed interpolant
OPTS.TppX = TppZ;            % Use pre-computed interpolant


%% Compute potential density surface
zref = ztar; % Reference depth for potential density [m]
z_PDENS = pot_dens_surf(S, T, Z, zref, [i0, j0, ztar], OPTS);

% The usual approach is to interpolate a 3D potential density field. 
% Check the differences. 
SIGMA = sort(eos(S, T, zref), 1); % Make potential density monotonic.
val_PDENS = interp1(Z, SIGMA(:,i0,j0), ztar);
z = ppc_linterp(SIGMA, Z, val_PDENS);
fig_map(g.XCvec, g.YCvec, z_PDENS - z, land); colorbar
title('depth difference between more and less accurate potential density surface')
clear SIGMA z



%% Compute specific volume anomaly surface
[s0,t0] = ppc_val2(Z, SppZ(:,:,i0,j0), TppZ(:,:,i0,j0), ztar);
z_DELTA = delta_surf(S, T, Z, s0, t0, [i0, j0, ztar], OPTS);

% The usual approach is to interpolate a 3D specific volume anomaly field. 
% Check the differences. 
DELTA = sort(R - eos(s0, t0, Z), 1); % Make monotonic.
val_DELTA = interp1(Z, DELTA(:,i0,j0), ztar);
z = ppc_linterp(DELTA, Z, val_DELTA);
fig_map(g.XCvec, g.YCvec, z_DELTA - z, land); colorbar
title('depth difference between more and less accurate in-situ density anomaly surface')
clear DELTA z

%% --- Compute the Reeb graph of the depth of a surface, without updating the surface

% Select the above potential density surface
z = z_PDENS;  

% Find linear index to reference cast
I0 = sub2ind([ni nj], i0, j0);

% Find the connected component containing the reference cast I0
[~,~,~,wet] = bfs_conncomp(isfinite(z), grid_adjacency([ni,nj], OPTS.WRAP), I0);

% Remove other regions from the surface
z(~wet) = nan;

% Calculate the Reeb graph of z
[z, RG] = calc_reeb_graph(z, OPTS);

%% Map each region of the Reeb graph:
% First, find the arc that each pixel in the one connected region (each
% vertex in the  simplical mesh) belongs to:
arc_segment = RG.arc_segment; % dereference for speed
arc_ = zeros(1, length(RG.n_casts));
for e = 1 : RG.nArcs
    arc_(arc_segment{e}) = e;
end

% Now build a map that sends pixel (i,j) to its arc in the Reeb graph:
arc = nan(ni, nj);
arc(RG.wet) = arc_;

% Map the regions associated with each arc:
fig_map(g.XCvec, g.YCvec, arc, land);
colormap(colorcube(RG.nArcs));
title('Map of all regions associated to arcs of the Reeb graph')

%% Show the Reeb graph -- in a small region
% (otherwise the graphics part of this will be very slow and very confusing)

% Longitude limits and latitude limits for the zoom region
x1 = 195; x2 = 223; y1 = 5; y2 = 25; % Just east of Hawaii

% Grid index limits for the zoom region
i1 = interp1(g.XCvec,1:ni,x1,'next');
i2 = interp1(g.XCvec,1:ni,x2,'prev');
j1 = interp1(g.YCvec,1:nj,y1,'next');
j2 = interp1(g.YCvec,1:nj,y2,'prev');

% ij2v maps from (i,j) coordinates to vertices in the simplical mesh
ij2v = reshape(cumsum(RG.wet(:)), ni, nj);
ij2v(~RG.wet) = 0;

% v2I maps from vertices in the simplical mesh to linear index for 2D data
% that you'd normally index by (i,j)
v2I = find(ij2v);

% Convert the linear indices I into 2D indices i and j.
[i,j] = ind2sub([ni, nj], v2I(RG.node_v));

figure('Position', [0 0 800 800]);
for panel = 1:2
    if panel == 1
        % Map depth of surface
        ax = axes('Position', [.05 .05 .9 .4]);
        imagesc(ax, g.XCvec(i1:i2), g.YCvec(j1:j2), -z(i1:i2,j1:j2).');
        colormap(ax,flipud(parula));
        title('Depth [m] of the surface')
    else
        % Map the regions associated to each arc
        ax = axes('Position', [.05 .55 .9 .4]);
        arczoom = arc(i1:i2,j1:j2);
        bad = isnan(arczoom);
        arczoom(bad) = -1;
        [~,~,arczoom] = unique(arczoom); %re-index them, from 1 to however many there are in the box
        arczoom = reshape(arczoom, i2-i1+1, j2-j1+1);
        arczoom(bad) = -99; % just so they appear different in imagesc
        imagesc(ax, g.XCvec(i1:i2), g.YCvec(j1:j2), arczoom.');
        colormap(ax,jet);
        title('Associated regions')
    end
    colorbar; ax.YDir = 'normal';
    hold(ax, 'on');
    
    
    % Plot nodes that are inside the zoomed region
    x = g.XCvec(i); x = x(:);
    y = g.YCvec(j); y = y(:);
    nodesub = x > x1 & x < x2 & y > y1 & y < y2;
    scatter(x(nodesub), y(nodesub), 24, 'ko');
    
    % Plot arcs in black, only those incident to at least one node inside the zoomed region
    arcsub = unique([vertcat(RG.node_next{nodesub}); vertcat(RG.node_prev{nodesub})]);
    xx = [x(RG.arc_from(arcsub)).'; x(RG.arc_to(arcsub)).'];
    yy = [y(RG.arc_from(arcsub)).'; y(RG.arc_to(arcsub)).'];
    plot(xx, yy, '-', 'Color', [1 1 1]*.5, 'LineWidth', .5);
    
    % Plot arcs in grey, only those between two nodes that are both inside the zoomed region
    xx = [x(RG.arc_from).'; x(RG.arc_to).'];
    yy = [y(RG.arc_from).'; y(RG.arc_to).'];
    arcsub = all(xx > x1,1) & all(xx < x2,1) & all(yy > y1,1) & all(yy < y2,1);
    plot(xx(:,arcsub), yy(:,arcsub), '-', 'Color', [1 1 1]*0);
    
    % Reset zoom
    ax.XLim = [x1 x2];
    ax.YLim = [y1 y2];
end


%% --- Compute an "orthobaric" surface
% initialized from the potential density surface
OPTS.REEB = false;              % Make single-valued functions
z_ORTHO = topobaric_surface(S, T, Z, z_PDENS, OPTS);
fig_map(g.XCvec, g.YCvec, z_ORTHO, land); colorbar;
title('Depth of an "orthobaric" surface')

% We could get more outputs from topobaric_surface, including the
% single-valued function for in-situ density anomaly as a function of
% depth. But let's carry on to topobaric surfaces...

%% --- Compute a topobaric surface
% initialized from the potential density surface
OPTS.REEB = true;     % Make multi-valued functions. (Default)
OPTS.GEOSTRF = false; % We won't (yet) compute a geostrophic stream function, so let it be ill-defined
[z_TOPOB, ~, ~, RG, s0, t0, dfn] = topobaric_surface(S, T, Z, z_PDENS, OPTS);
fig_map(g.XCvec, g.YCvec, z_TOPOB, land); colorbar;
title('Depth of a topobaric surface')

%% Scatter plot depth vs in-situ density anomaly on the surface
% Get S and T on the surface
lead1 = @(x) reshape(x, [1 size(x)]);
[s,t] = ppc_val2(Z, SppZ, TppZ, lead1(z_TOPOB));

% Calculate in-situ density anomaly on the surface
d = eos(s, t, z_TOPOB) - eos(s0, t0, z_TOPOB); % delta

% Select just those water columns in the one connected region
z_ = z_TOPOB(RG.wet);
d_ = d(RG.wet);

% Determine which arc each water column in the region belongs to:
arc_segment = RG.arc_segment; % dereference for speed
arc_ = zeros(1, RG.n_casts);
for e = 1 : RG.nArcs
    arc_(arc_segment{e}) = e;
end

% Scatter plot d vs z, coloured according to which arc each data point is in
figure;
scatter(z_, d_, 8, arc_);
colormap(colorcube(RG.nArcs));
xlabel('Depth [m]')
ylabel('in-situ density anomaly [kg m^{-3}]')


%% How well does the density function match the density?
% Map the difference between the in-situ density anomaly on the surface,
% and the multivalued function for the in-situ density anomaly evaluated at
% the depth of the surface.  They are almost identical. That they are not
% exactly identical is because the root-finding procedure in each water
% column that matches these quantities stops at a certain tolerance.

% Evaluate the multivalued delta function at each water column in the
% connected region
dfn_at_z_ = nan(RG.n_casts, 1);
for e = 1 : RG.nArcs
    dfn_at_z_(arc_segment{e}) = pvallin(dfn(:,e), z_(arc_segment{e}));
end

% ... and transfer that onto a map of physical space
dfn_at_z = nan(ni,nj);
dfn_at_z(RG.wet) = dfn_at_z_;

% Map the difference in deltas (!)
fig_map(g.XCvec, g.YCvec, d - dfn_at_z, land);
caxis([-1 1]*1e-6); colorbar;
title('Difference between $\rho$ and $\hat{\rho}(p)$ on surface', 'Interpreter', 'latex');

%% --- Compute an omega surface
% initialized from the potential density surface
z_OMEGA = omega_surface(S, T, Z, z_PDENS, OPTS);
fig_map(g.XCvec, g.YCvec, z_TOPOB, land); colorbar;
title('Depth of omega surface')

%% Compare the neutral error between the potential density surface, topobaric surface, and omega surface
figure('Position', [0 0 600 800]);

npanel = 3; % number of vertical panels
mv = .05; % vertical margin
hv = (1-(npanel+1)*mv) / npanel;

z = z_PDENS;
[s,t] = ppc_val2(Z, SppZ, TppZ, lead1(z));
eps_x = ntp_epsilon_r_x(s, t, z, g.DXCvec, g.DYCsc, false, g.WRAP);
ax = axes('Position', [.1 3*mv+2*hv .87 hv]);
fig_map(ax, g.XCvec, g.YCvec, log10(abs(eps_x)), land, struct('CLIM', [-12 -5]));
colorbar
title('Zonal neutral error on potential density surface [log_{10}]')

z = z_TOPOB;
[s,t] = ppc_val2(Z, SppZ, TppZ, lead1(z));
eps_x = ntp_epsilon_r_x(s, t, z, g.DXCvec, g.DYCsc, false, g.WRAP);
ax = axes('Position', [.1 2*mv+hv .87 hv]);
fig_map(ax, g.XCvec, g.YCvec, log10(abs(eps_x)), land, struct('CLIM', [-12 -5]));
colorbar
title('Zonal neutral error on topobaric surface [log_{10}]')

z = z_OMEGA;
[s,t] = ppc_val2(Z, SppZ, TppZ, lead1(z));
eps_x = ntp_epsilon_r_x(s, t, z, g.DXCvec, g.DYCsc, false, g.WRAP);
ax = axes('Position', [.1 mv .87 hv]);
fig_map(ax, g.XCvec, g.YCvec, log10(abs(eps_x)), land, struct('CLIM', [-12 -5]));
colorbar
title('Zonal neutral error on omega surface [log_{10}]')


%% --- Geostrophic Stream Function comparison
% on a potential density surface
% mapping errors in predicting the zonal geostrophic velocity
Y = hsap3(Z, ATMP, ETAN, R, grav, rho_c); % Pre-compute hydrostatic acceleration potential

z = z_PDENS;

[s,t] = ppc_val2(Z, SppZ, TppZ, lead1(z));
[y,r] = hsap2(s, t, z, Z, R, Y, grav, rho_c);

% "True" zonal geostrophic velocities, in Boussinesq mode, computed from
% gradients sloping in the surface. (-) before z undoes that we've been z>0
% people! g.cori_YG is the Coriolis parameter on V grid (for non zonal
% grids, should average this)
ugs = -((y - jm1(y)) + (grav / (2*rho_c)) * (-z + jm1(z)) .* (r + jm1(r))) ./ (g.DYCsc * cori_V);

% Reference S and T values for various geostrophic stream functions:
s0 = s(i0,j0);
t0 = t(i0,j0);
z0 = z(i0,j0);

hf = figure('Position', [0 0 800 800]);
OPTS_FIGS = struct('CLIM', [-8 0]);
ax = subplot(3,2,1);
ZH = zhanghogg92(s, t, z, Z, R, Y, s0, t0, z0, grav, rho_c);
uzh = -(ZH - jm1(ZH)) ./ (g.DYCsc * cori_V); % Zonal geostrophic velocity by the Zhang and Hogg gsf
fig_map(ax, g.XCvec, g.YCvec, log10(abs(uzh - ugs)), land, OPTS_FIGS); colorbar;
title(ax, 'log_{10} |u_g error| - Zhang and Hogg (1992)');

ax = subplot(3,2,2);
MK = mcdougallklocker10(s, t, z, Z, R, Y, s0, t0, z0, grav, rho_c);
umk = -(MK - jm1(MK)) ./ (g.DYCsc * cori_V);
fig_map(ax, g.XCvec, g.YCvec, log10(abs(umk - ugs)), land, OPTS_FIGS); colorbar;
title(ax, 'log_{10} |u_g error| - McDougall and Klocker (2010)');

ax = subplot(3,2,3);
CU = cunningham00(s, t, z, Z, R, Y, z0, grav, rho_c);
ucu = -(CU - jm1(CU)) ./ (g.DYCsc * cori_V);
fig_map(ax, g.XCvec, g.YCvec, log10(abs(ucu - ugs)), land, OPTS_FIGS); colorbar;
title(ax, 'log_{10} |u_g error| - Cunningham (2000)');

ax = subplot(3,2,4);
OM = orthobaric_montgomery(s, t, z, Z, R, Y, [], [], struct('WRAP', g.WRAP), grav, rho_c);
uom = -(OM - jm1(OM)) ./ (g.DYCsc * cori_V);
fig_map(ax, g.XCvec, g.YCvec, log10(abs(uom - ugs)), land, OPTS_FIGS); colorbar;
title(ax, 'log_{10} |u_g error| - Orthobaric Montgomery');

% Compare with orthobaric geostrophic stream function:
ax = subplot(3,2,5);
OPTS.REEB = false; % Turn off the Reeb graph => a single-valued function is fit
OPTS.SPLINE_BREAKS = [0 200 6000]; % [m]
OPTS.SPLINE_ORDER = 4; % cubic splines
OB = topobaric_geostrf(s, t, z, Z, R, Y, s0, t0, OPTS, grav, rho_c);
uob = -(OB - jm1(OB)) ./ (g.DYCsc * cori_V);
fig_map(ax, g.XCvec, g.YCvec, log10(abs(uob - ugs)), land, OPTS_FIGS); colorbar;
title(ax, 'log_{10} |u_g error| - Orthobaric gsf');

ax = subplot(3,2,6);
OPTS.REEB = true;    % Use multivalued functions for empirical fits
OPTS.GEOSTRF = true; % Force well-defined geostrophic stream function
TB = topobaric_geostrf(s, t, z, Z, R, Y, s0, t0, OPTS, grav, rho_c);
utb = -(TB - jm1(TB)) ./ (g.DYCsc * cori_V);
fig_map(ax, g.XCvec, g.YCvec, log10(abs(utb - ugs)), land, OPTS_FIGS); colorbar;
title(ax, 'log_{10} |u_g error| - Topobaric gsf');


%% --- Compute a modified topobaric surface and its exact gsf
OPTS.REEB = true;
OPTS.ITER_MAX = 6;   % A good number of iterations for convergence, so the gsf will be very accurate
OPTS.GEOSTRF = true; % Add second integral constraint so that dfn will integrate to give a well-defined gsf
[z_MTOPO, ~, ~, RG, s0, t0, dfn] = topobaric_surface(S, T, Z, z_PDENS, OPTS);
z = z_MTOPO;

% --- Compute zonal geostrophic velocity on the modified topobaric surface
[s,t] = ppc_val2(Z, SppZ, TppZ, lead1(z));
[y,r] = hsap2(s, t, z, Z, R, Y, grav, rho_c);
ugs = -((y - jm1(y)) + (grav / (2*rho_c)) * (-z + jm1(z)) .* (r + jm1(r))) ./ (g.DYCsc * cori_V);

% --- Compute topobaric geostrophic stream function on a modified topobaric surface
% This geostrophic stream function is exact on this surface.
OPTS.REEB = true; % Use multivalued functions for empirical fits
OPTS.RG    = RG;  % taken from Modified Topobaric Surface
OPTS.dfn   = dfn; % taken from Modified Topobaric Surface
TB = topobaric_geostrf(s, t, z, Z, R, Y, s0, t0, OPTS, grav, rho_c);
OPTS = rmfield(OPTS, {'RG', 'dfn'});

utb = -(TB - jm1(TB)) ./ (g.DYCsc * cori_V); % Zonal geostrophic velocity by the topobaric gsf
fig_map(g.XCvec, g.YCvec, log10(abs(utb - ugs)), land, OPTS_FIGS); colorbar;
title('log_{10} geostrophic velocity error, using topobaric gsf on modified topobaric surface');
