% Make most of the figures and results, using ECCO2 data, in
% Stanley, G.J., 2019a. Neutral surface topology. Ocean Modelling 138,
% 88–106. https://doi.org/10.1016/j.ocemod.2019.01.008
% and
% Stanley, G.J., 2019b. The exact geostrophic streamfunction for neutral
% surfaces. Ocean Modelling 138, 107–121.
% https://doi.org/10.1016/j.ocemod.2019.04.002
%
% Note: To reproduce exactly these papers' results, an earlier version of
% this code must be used, namely v1.0 (commit cae0c75) of Topobaric-Surface
% available at
% https://github.com/geoffstanley/Topobaric-Surface/releases/tag/v1.0

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
% Version   : 2.0.0
%
% Modified by : --
% Date        : --
% Changes     : --

%% --- BEGIN SETUP --------------------------------------------------------
%#ok<*UNRCH>
warning('off', 'MATLAB:nargchk:deprecated')
set(0, 'defaultfigurecolor', [1 1 1]); % white figure background
V = filesep(); % /  or  \  depending on OS.

OMEGA_FRESH = true;  % Use this to compute omega surfaces from scratch
%OMEGA_FRESH = false; % Use this to load pre-computed omega surfaces

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

%% Load ECCO grid
TIMESTEP = '20021223';
g = load_ECCO2(PATH_ECCO2, 'grid');
ni = g.nx;
nj = g.ny;
nk = g.nz;
g.RAC = repmat(g.RACvec, ni, 1);
g.RAS = repmat(g.RASvec, ni, 1);
g.RAW = repmat(g.RAWvec, ni, 1);
g.YC = repmat(g.YCvec, ni, 1);
g.YG = repmat(g.YGvec, ni, 1);
rho_c = g.rho_c;
grav = g.grav;
Z = -g.RC(:); % We are going to be Z > 0 people!
Z2P = Pa2db * rho_c * grav ; % Note > 0
lead1 = @(x) reshape(x, [1 size(x)]);

neigh = grid_adjacency([ni,nj], g.WRAP);

%% Set alias functions
% Choose the Boussinesq densjmd95 and set grav and rho_c in eos.m and eos_x.m
eos_set_bsq_param([PATH_PROJECT 'lib' V 'eos' V 'eoscg_densjmd95_bsq.m'   ], [PATH_PROJECT 'lib' V 'alias' V 'eos.m'  ], grav, rho_c);
eos_set_bsq_param([PATH_PROJECT 'lib' V 'eos' V 'eoscg_densjmd95_bsq_dz.m'], [PATH_PROJECT 'lib' V 'alias' V 'eos_x.m'], grav, rho_c);

% Choose vertical interpolation method
copyfile([PATH_PROJECT 'fex' V 'columncalculus' V 'interp1qn2.m'], [PATH_PROJECT 'lib' V 'alias' V 'interp_firstdim_twovars.m']); % linear interpolation
%copyfile([PATH_PROJECT 'fex' V 'columncalculus' V 'pchipqn2.m'], [PATH_PROJECT 'lib' V 'alias' V 'interp_firstdim_twovars.m']);  % PCHIP interpolation
clear interp_firstdim_twovars % Make sure new file gets used

%% C-grid functions for 2D data (longitude = i = rows, latitude = j = columns)
im1 = @(F) circshift(F, [+1 0]);
ip1 = @(F) circshift(F, [-1 0]);
jm1 = @(F) subsasgn(circshift(F, [0 +1]), struct('type','()', 'subs', {{':',   1}}), nan); %#ok<SUBSASGN>
jp1 = @(F) subsasgn(circshift(F, [0 -1]), struct('type','()', 'subs', {{':',nj}}), nan); %#ok<SUBSASGN>

avg_i_ip1 = @(X) (X + ip1(X)) / 2;   % Average (i,j) and (i+1,j).
avg_j_jp1 = @(X) (X + jp1(X)) / 2;   % Average (i,j) and (i,j+1).
avg_i_im1 = @(X) (X + im1(X)) / 2;   % Average (i,j) and (i-1,j).
avg_j_jm1 = @(X) (X + jm1(X)) / 2;   % Average (i,j) and (i,j-1).
avg_T_onto_U = avg_i_im1;
avg_T_onto_V = avg_j_jm1;
dTdx_on_U = @(T) (T - im1(T)) ./ g.DXCvec; % Derivative between (i,j) and (i-1,j)
dTdy_on_V = @(T) (T - jm1(T)) / g.DYCsc ; % Derivative between (i,j) and (i,j-1)

y_to_j = @(y) floor((y - g.YGvec(1)) * g.resy) + 1;
x_to_i = @(x) mod(floor((x - g.XGvec(1)) * g.resx), ni) + 1;

%% Load ECCO data
[S, T, ATMP, ETAN] = load_ECCO2(PATH_ECCO2, 'casts', TIMESTEP);
GAMMA = load_ECCO2(PATH_ECCO2, 'gamma', TIMESTEP);

% Re-order data so water-columns are contiguous data:
S = permute(S, [3 1 2]); % [nz,nx,ny]. depth  by  long  by  lat
T = permute(T, [3 1 2]);
GAMMA = permute(GAMMA, [3 1 2]);

% Compute in-situ density
R = eos(S, T, Z);

% Pre-compute mixed layer
%MLZ = mixed_layer(S, T, Z); % [m].  Use this to remove the mixed layer
MLZ = [];  % Use this to leave the mixed layer in

%% Selecting surfaces and depths
list_surf = {'\sigma', '\delta', '\gamma^n', '\sigma_\nu', '\omega', '\tau', '\tau'''};
LIST_SURF = {'SIGMA', 'DELTA', 'GAMMA', 'ORTHO', 'OMEGA', 'TOPOB', 'MTOPO'};
iSURF = containers.Map(LIST_SURF, 1:length(list_surf));
nS = length(list_surf);

list_zref = [1000 2000];
nZ = length(list_zref);

z_i0j0 = nan(nZ,1); % depth of omega surface at (i0,j0)

% Will save extra output from topobaric_surface() for tau' surface
mtopo_data = struct('RG', cell(nZ,1), 's0', cell(nZ,1), 't0', cell(nZ,1), 'dfn', cell(nZ,1));
zs = nan(ni,nj,nS,nZ); % depth for each surface

%% Setup options for (regular / modified) topobaric and "orthobaric" surfaces

% Start in the Pacific Ocean
x0 = 180; y0 = 0;
i0 = x_to_i(x0);
j0 = y_to_j(y0);

% Tolerance for bisection to find surfaces
tolx = 1e-4;

OPTS = struct();
OPTS.WRAP = g.WRAP;    % Periodic in longitude, not in latitude
OPTS.REF_IJ = [i0 j0]; % Pacific Ocean
OPTS.GEOSTRF = false;  % Don't add extra criteria that make geostrophic streamfunction well-defined
OPTS.SIMPLIFY_ARC_REMAIN = Inf; % No simplification
OPTS.FILL_IJ = [];     % No filling
OPTS.FILL_PIX = 0;     % No filling.
%OPTS.MLX = MLZ;        % Mixed Layer Depth
OPTS.MLX = [];         % Leave the Mixed Layer in.
OPTS.X_TOL = tolx;     % error tolerance when updating the surface
OPTS.X_EXPN = 500;     % expansion of domain to search for solutions in each water column
OPTS.ITER_L2_CHANGE = 1e-3; % Stop iterations when the L2 change of pressure on the surface is less than this
OPTS.ITER_MAX = 6;     % The maximum number of iterations
OPTS.ITER_START_WETTING = 1; % Start wetting on first iteration
OPTS.VERBOSE = 1;      % Some output while executing topobaric_surface()
OPTS.FILE_ID = fileID; % Write output to this file
OPTS.FIGS_SHOW = false; % Don't show figures (for omega surface)

%% Figure setup:
land = squeeze(isnan(S(1,:,:)));

OPTS_FIGS = struct();
OPTS_FIGS.FONTSIZE = 10 + 2*ismac(); % 10 normally, 12 for mac...
OPTS_FIGS.XSH = 20;  % Shift so 20 E is at edges of figure.
OPTS_FIGS.LATTICK = -60 : 30 : 90;
OPTS_FIGS.LONLIM = [0 360];
OPTS_FIGS.LATLIM = [-80 +72]; % This will change
OPTS_FIGS.GRIDCOL = []; % no grid
OPTS_FIGS.LANDCOL = [1 1 1]*0; % Black
OPTS_FIGS.NANCOL = [1 1 1]* .75; % Grey
OPTS_FIGS.VIS = 'on'; % show figures on the screen
%OPTS_FIGS.VIS = 'off'; % hide figures from the screen

%% --- DONE SETUP -------------------------------------------------------------

%% Calculate surfaces
fprintf(fileID, '--- Calculate surfaces ---\n');
for iZ = 1 : nZ
    
    zref = list_zref(iZ);
    
    %% Omega surface
    
    if OMEGA_FRESH
        
        % Compute a potential density surface to initialize omega_surface
        [z_sigma, val] = pot_dens_surf(S, T, Z, zref, [i0, j0, zref], tolx, OPTS);
        fprintf(fileID, '(%d,%d,%.4fm): SIGMA = %.4f, z_ref = %.4f\n', ...
            i0, j0, zref, val, zref);

        % Calculate omega surface
        mytic = tic;
        zs(:,:,iSURF('OMEGA'),iZ) = omega_surface(S, T, Z, z_sigma, OPTS);
        fprintf(fileID, '(%d,%d,%.4fm): OMEGA done in time %.2f\n', ...
            i0, j0, z_i0j0(iZ), toc(mytic));
        
        clear z_sigma
    else
        zs(:,:,iSURF('OMEGA'),iZ) = load_ECCO2(PATH_ECCO2, 'omega', TIMESTEP, zref); % Load pre-computed omega surface
    end
    
    % Select reference depth. All surfaces shall intersect this point:
    z_i0j0(iZ) = zs(i0,j0,iSURF('OMEGA'),iZ);
    
    %% Potential density Surface
    [zs(:,:,iSURF('SIGMA'),iZ), val] = pot_dens_surf(S, T, Z, zref, [i0, j0, z_i0j0(iZ)], tolx, OPTS);
    fprintf(fileID, '(%d,%d,%.4fm): SIGMA = %.4f, z_ref = %.4f\n', ...
        i0, j0, z_i0j0(iZ), val, zref);
    
    % Figure showing depth of isopycnal
    %{
    hf = figure;
    ax = axes;
    fig_map(ax, g.XCvec, g.YCvec, -zs(:,:,iSURF('SIGMA'),iZ), land, setfield(setfield(OPTS_FIGS, 'CLIM', []), 'CLIM', []));
    caxis([-zref-500 0]);
    colorbar
    fn = sprintf('z_SIGMA%d_(%d,%d,%.4fm)', list_zref(iZ), x0, y0, z_i0j0(iZ));
    export_fig(hf, [PATH_FIGS fn], '-jpg', '-m2');
    %}
    
    % Figure showing depth of omega surface - depth of isopycnal
    %{
    hf = figure;
    ax = axes;
    fig_map(ax, g.XCvec, g.YCvec, -zs(:,:,iSURF('SIGMA'),iZ) - -zs(:,:,iSURF('OMEGA'),iZ), land, setfield(setfield(OPTS_FIGS, 'CLIM', []), 'CLIM', []));
    caxis([-1 1] * 200);
    colorbar
    % colormap(ax, bluewhitered)
    fn = sprintf('z_SIGMA_minus_OMEGA_%d_(%d,%d,%.4fm)', list_zref(iZ), x0, y0, z_i0j0(iZ));
    export_fig(hf, [PATH_FIGS fn], '-jpg', '-m2');
    %}
    
    
    %% Density anomaly surface
    
    % Get reference S and T as Southern Ocean average:
    [s,t] = interp_firstdim_twovars(lead1(zs(:,:,iSURF('OMEGA'),iZ)), Z, S, T);
    jj = y_to_j(-55) : y_to_j(-50);
    S_delta = nanmean(reshape(s(:,jj), [], 1));
    T_delta = nanmean(reshape(t(:,jj), [], 1));
    clear s t
    
    [zs(:,:,iSURF('DELTA'),iZ),val] = delta_surf(S, T, Z, S_delta, T_delta, [i0, j0, z_i0j0(iZ)], tolx, OPTS);
    fprintf(fileID, '(%d,%d,%.4fm): DELTA = %.4e, S_delta = %.4f, T_delta = %.4f\n', ...
        i0, j0, z_i0j0(iZ), val, S_delta, T_delta);
    
    
    %% Neutral Density surface
    val = interp1qn(z_i0j0(iZ), Z, squeeze(GAMMA(:,i0,j0)));
    zs(:,:,iSURF('GAMMA'),iZ) = interp1qn(val, GAMMA, Z);
    
    fprintf(fileID, '(%d,%d,%.4fm): GAMMA = %.4f\n', ...
        i0, j0, z_i0j0(iZ), val);
    if iZ == nZ
        clear GAMMA
    end
    
    %% Calculate "orthobaric" surface
    OPTS.REEB = false;
    OPTS.GEOSTRF = false;
    if zref == 1000
        OPTS.SPLINE_BREAKS = [0 200 6000];
        OPTS.SPLINE_ORDER = 4;
    elseif zref == 2000
        OPTS.SPLINE_BREAKS = [0 200 1500 1800 6000];
        OPTS.SPLINE_ORDER = 4;
    end
    
    mytic = tic;
    zs(:,:,iSURF('ORTHO'),iZ) = topobaric_surface(S, T, Z, zs(:,:,iSURF('SIGMA'),iZ), OPTS);
    fprintf(fileID, '(%d,%d,%.4fm): ORTHOBARIC done in time %.2f\n', ...
        i0, j0, z_i0j0(iZ), toc(mytic));
    
    %% Calculate topobaric surface
    % TBS does depend on initial surface, but using potential density or
    % orthobaric surfaces or gamma^n surfaces, is essentially irrelevant.
    % Starting from Omega surface does produce a small improvement (less
    % than 10%).
    
    OPTS.REEB = true;
    OPTS.GEOSTRF = false;
    
    mytic = tic;
    zs(:,:,iSURF('TOPOB'),iZ) = topobaric_surface(S, T, Z, zs(:,:,iSURF('SIGMA'),iZ), OPTS);
    fprintf(fileID, '(%d,%d,%.4fm): TOPOBARIC done in time %.2f\n', ...
        i0, j0, z_i0j0(iZ), toc(mytic));
    
    % depth difference between topobaric and omega surfaces:
    %{
    if iZ == 2
        hf = figure;
        ax = axes;
        fig_map(ax, g.XCvec, g.YCvec, -zs(:,:,iSURF('TOPOB'),iZ) - -zs(:,:,iSURF('OMEGA'),iZ), land, setfield(setfield(OPTS_FIGS, 'CLIM', []), 'CLIM', []));
        caxis([-1 1] * 60);
        colorbar
        % colormap(ax, bluewhitered)
        fn = sprintf('z_TOPOB_minus_OMEGA_%d_(%d,%d,%.4fm)', list_zref(iZ), x0, y0, z_i0j0(iZ));
        export_fig(hf, [PATH_FIGS fn], '-jpg', '-m2');
        close(hf)
    end
    %}
    %% Calculate modified topobaric surface
    OPTS.REEB = true;
    OPTS.GEOSTRF = true;
    
    mytic = tic;
    [zs(:,:,iSURF('MTOPO'),iZ), ~, ~, mtopo_data(iZ).RG, mtopo_data(iZ).s0, mtopo_data(iZ).t0, mtopo_data(iZ).dfn] ...
        = topobaric_surface(S, T, Z, zs(:,:,iSURF('SIGMA'),iZ), OPTS);
    fprintf(fileID, '(%d,%d,%.4fm): MODIFIED TOPOBARIC done in time %.2f\n', ...
        i0, j0, z_i0j0(iZ), toc(mytic));
    
end % iZ

%% Figure: Map associated regions for sigma_2 surface
iZ = 2;
z = zs(:,:,iSURF('SIGMA'),iZ); 
OPTS_FIGS.LATLIM = [-80, 72]; % The z(180,0) = 2000 surfaces exclude the Arctic.

% Find the connected component containing the reference cast
I0 = sub2ind([ni,nj], i0, j0); % Linear index to reference cast (i0,j0)
[~,qts,~,wet] = bfs_conncomp(isfinite(z), neigh, I0);
z(~wet) = nan;
n_casts = qts(end)-1; % == sum(wet(:))

% Calculate the Reeb graph
[~, RG] = calc_reeb_graph(z, OPTS);

arc_ = zeros(n_casts,1);
for e = 1:RG.nArcs
    arc_(RG.arc_segment{e}) = e;
end
if any(arc_ == 0)
    fprintf(fileID, '*** WARNING *** %d vertices lost.\n', sum(arc_ == 0));
end
arc = nan(ni,nj);
arc(wet) = arc_;

hf = figure; ax = axes;
fig_map(ax, g.XCvec, g.YCvec, arc, land, OPTS_FIGS);
cm = distinguishable_colors(512, OPTS_FIGS.NANCOL);
rp = randperm(size(cm,1));
cm = cm(rp,:);
colormap(ax,cm);
colorbar(ax, 'off');
ax.CLim = [0 RG.nArcs];

clear z RG arc arc_ cm rp

%% Save figure
fn = sprintf('arc(%d,%d,%.4fm)', i0, j0, z_i0j0(iZ));
export_fig(hf, [PATH_FIGS fn], '-jpg', '-m3');
close(hf);


%% --- BEGIN NEUTRAL CALCULATIONS ---------------------------------------------
%% Make a common mask for neutral calculations
I0 = sub2ind([ni,nj], i0, j0);
good_mask = false(ni,nj,1,nZ);
for iZ = 1:nZ
    
    if ~isempty(MLZ)
        % Will only sample from data where all surfaces are valid & below mixed layer
        good_mask(:,:,1,iZ) = all(zs(:,:,:,iZ) > MLZ, 3);
    else
        good_mask(:,:,1,iZ) = all(isfinite(zs(:,:,:,iZ)), 3);
    end
    
    % Select only one connected ocean region
    [~,~,~,good_mask(:,:,1,iZ)] = bfs_conncomp(good_mask(:,:,1,iZ), neigh, I0);
    
end

%% --- Calculate neutral errors -----------------------------------------------
mytic = tic;

eps = struct();
fdd = struct();
eps.mapx = nan(ni,nj,nS,nZ); % neutral error, zonal
eps.mapy = nan(ni,nj,nS,nZ); % neutral error, meridional
fdd.map  = nan(ni,nj,nS,nZ); % fictitious diapycnal diffusivity

K_iso = 1000; % Isopycnal diffusivity (useful ball park estimate)

for iZ = 1:nZ
    for iS = 1:nS
        surfname = LIST_SURF{iS};
        
        % Select surface and apply common mask
        z = zs(:,:,iS,iZ);
        z(~good_mask(:,:,1,iZ)) = nan;
        
        if strcmp(surfname, 'GAMMA') || (strcmp(surfname, 'OMEGA') && ~OMEGA_FRESH)
            % OMEGA-surface was pre-computed using linear interpolation
            % GAMMA-surface was linearly interpolated in this script
            interpfn = @interp1qn2;
        else
            interpfn = @interp_firstdim_twovars;
        end
        [s, t] = interpfn(lead1(z), Z, S, T);
        
        % Get neutrality errors on the U,V grid
        [eps.mapx(:,:,iS,iZ), eps.mapy(:,:,iS,iZ)] = ntp_epsilon_r_x(s,t,z,g.DXCvec,g.DYCsc,false,g.WRAP);
        
        % Get slope errors on the Tracer grid
        [~, ~, sx, sy] = ntp_epsilon_r_x(s,t,z,g.DXCvec,g.DYCsc,true,g.WRAP,[],S,T,Z);
        
        % Fictitious Diapycnal Diffusivity
        fdd.map(:,:,iS,iZ) = K_iso * (sx.^2 + sy.^2);
    end
end
fprintf(fileID, 'Neutral errors calculated in time %.2f\n', toc(mytic));
clear sx sy

%% Calculate histograms and norms of neutral errors and FDD

eps.L1 = nan(nS,nZ);
eps.L2 = nan(nS,nZ);
fdd.L1 = nan(nS,nZ);
fdd.L2 = nan(nS,nZ);

eps.edges = linspace(-14.5, -5, 100);
eps.bins = (eps.edges(1:end-1) + eps.edges(2:end)) / 2;
eps.hc = nan(length(eps.bins),nS,nZ);
eps.XTick = -14 : -5;
eps.table_symb = '||\epsilon||';
eps.XLabel = '$\log_{10}(|\epsilon|)$';

fdd.edges = linspace(-15.5, -2, 100);
fdd.bins = (fdd.edges(1:end-1) + fdd.edges(2:end)) / 2;
fdd.hc = nan(length(fdd.bins),nS,nZ);
fdd.XTick = -15 : -2;
fdd.table_symb = '||D^f||';
fdd.XLabel = '$\log_{10}(D^f)$';

for iZ = 1:nZ
    
    for iS = 1:nS
        % Calculate error norms and histogram counts for epsilon errors
        data = [reshape(eps.mapx(:,:,iS,iZ), [], 1); reshape(eps.mapy(:,:,iS,iZ), [], 1)];
        if iS == 1
            good = isfinite(data);
            AREA = [g.RAW(:); g.RAS(:)];
            AREA = AREA(good);
        end
        data = data(good);
        eps.L1(iS,iZ) = sum(     AREA .* abs(data)) / sum(AREA);
        eps.L2(iS,iZ) = sqrt(sum(AREA .*   data.^2) / sum(AREA));
        
        eps.hc(:,iS,iZ) = histcounts(log10(abs(data)), eps.edges);
        eps.hc(:,iS,iZ) = eps.hc(:,iS,iZ) / sum(eps.hc(:,iS,iZ)); % Normalize
        
    end
    
    for iS = 1:nS
        % Calculate error norms and histogram counts for FDD errors
        data = fdd.map(:,:,iS,iZ);
        if iS == 1
            good = isfinite(data);
            AREA = g.RAC(good);
        end
        data = data(good);
        
        fdd.L1(iS,iZ) = sum(     AREA .* abs(data)) / sum(AREA);
        fdd.L2(iS,iZ) = sqrt(sum(AREA .*   data.^2) / sum(AREA));
        
        fdd.hc(:,iS,iZ) = histcounts(log10(abs(data)), fdd.edges);
        fdd.hc(:,iS,iZ) = fdd.hc(:,iS,iZ) / sum(fdd.hc(:,iS,iZ)); % Normalize
        
    end
end
clear AREA data good

%% Report: norms of errors
fprintf(fileID, '--- Norms of errors ---\n');
fprintf(fileID, '|eps|_1, |eps|_2, |fdd|_1, |fdd|_2,\n');
for iZ = 1:nZ
    fprintf(fileID, 'Surfaces through (%d,%d,%.4fm)\n', i0, j0, z_i0j0(iZ));
    for iS = 1:nS
        fprintf(fileID, '%.4e; %.4e; %.4e; %.4e;  %s\n', ...
            eps.L1(iS,iZ), ...
            eps.L2(iS,iZ), ...
            fdd.L1(iS,iZ), ...
            fdd.L2(iS,iZ), ...
            list_surf{iS} );
    end
end

%% Report: norms of errors, as ratio
fprintf(fileID, '--- Norms of errors, as ratio ---\n');
iS0 = iSURF('TOPOB');
fprintf(fileID, '|eps|_1, |eps|_2, |fdd|_1, |fdd|_2,\n  all as ratios from given surface to %s\n', list_surf{iS0});
for iZ = 1:nZ
    fprintf(fileID, 'Surfaces through (%d,%d,%.4fm)\n', i0, j0, z_i0j0(iZ));
    for iS = 1:nS
        fprintf(fileID, '%.4e; %.4e; %.4e; %.4e;  %s\n', ...
            eps.L1(iS,iZ) / eps.L1(iS0,iZ), ...
            eps.L2(iS,iZ) / eps.L2(iS0,iZ), ...
            fdd.L1(iS,iZ) / fdd.L1(iS0,iZ), ...
            fdd.L2(iS,iZ) / fdd.L2(iS0,iZ), ...
            list_surf{iS} );
    end
end


%% Report: Fraction of area with FDD > 1e-5
fprintf(fileID, '--- Fraction of area with FDD > 1e-5 ---\n');
for iZ = 1:nZ
    fprintf(fileID, 'Column 1: # grid points with FDD > 1e-5 m2 s-1\n');
    fprintf(fileID, 'Column 2: areal fraction of grid points with FDD > 1e-5 m2 s-1\n');
    fprintf(fileID, 'Column 3: surface through (%d,%d,%.4fm)\n', i0, j0, z_i0j0(iZ));
    for iS = 1:nS
        bigfdd = fdd.map(:,:,iS,2);
        good = isfinite(bigfdd);
        bigfdd = bigfdd(good) > 1e-5;
        fprintf(fileID, '%06d ; %.6e ;  %s\n', sum(bigfdd), sum(g.RAC(good) .* bigfdd) / sum(g.RAC(good)), list_surf{iS});
    end
end
clear good bigfdd

%% --- Figure: histogram of neutral errors
SURFsubset = {'SIGMA', 'DELTA', 'GAMMA', 'ORTHO', 'OMEGA', 'TOPOB'};
iSURFsubset = cellfun(@(c) iSURF(c), SURFsubset);
nSs = length(SURFsubset);
lw = 2; % line width
fs = 14 + double(ismac) * 4; % fontsize

subaxopts = {'marginleft',.03,'marginright',.01,'marginbottom',.07,'margintop',.01, ...
    'spacinghorizontal', .03, 'spacingvertical', .01};

hf = figure('Position', [0 0 600*2 400*nZ], 'Visible', OPTS_FIGS.VIS);

cm = [ ... % Colormap for lines
    0.0       0.0       1.0
    1.0       0.0       0.0
    0.85      0.6       0.0
    0.5       0.2       0.7
    0.0       0.5       0.0
    0.0       0.0       0.0 ];
cm = cm([1 6 2 4 3 5],:);

panels = {eps, fdd};
for iZ = 1:nZ
    zref = list_zref(iZ);
    
    for iP = 1 : 2
        pan = panels{iP};
        
        ax = subaxis(nZ,2,iP,iZ,subaxopts{:});
        hold(ax,'on'); box(ax,'on'); grid(ax,'on')
        ax.XTick = pan.XTick;
        ax.XLim = pan.edges([1 end]);
        ax.FontSize = fs;
        
        x = ax.Position(1) + .005;
        w = ax.Position(3) * .325;
        h = ax.Position(4) * .4;
        y = ax.Position(2) + ax.Position(4) - .005 - h;
        ax2 = axes('Position', [x y w h]);
        ax2.XTick = []; ax2.YTick = [];
        box(ax2,'on');
        taby = .92;
        text(ax2, .35, taby, ['$' pan.table_symb '_1$'], 'units','n', 'HorizontalAlignment', 'center', 'FontSize', fs, 'interpreter', 'latex');
        text(ax2, .8 , taby, ['$' pan.table_symb '_2$'], 'units','n', 'HorizontalAlignment', 'center', 'FontSize', fs, 'interpreter', 'latex');
        taby = .85;
        for iS = iSURFsubset
            
            taby = taby - .85 / nSs;
            if iS == iSURF('SIGMA')
                txt = sprintf('_{%g}', zref/1000);
            else
                txt = '';
            end
            text(ax2, .03 ,taby, ['$' list_surf{iS} txt '$']  , 'units','n','Color',cm(iS,:), 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'left' , 'FontSize', fs, 'Interpreter', 'latex');
            text(ax2, .53, taby, sprintf('%.2e',pan.L1(iS,iZ)), 'units','n','Color',cm(iS,:), 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'right', 'FontSize', fs, 'Interpreter', 'latex');
            text(ax2, .97, taby, sprintf('%.2e',pan.L2(iS,iZ)), 'units','n','Color',cm(iS,:), 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'right', 'FontSize', fs, 'Interpreter', 'latex');
            
            plot(ax, pan.bins, pan.hc(:,iS,iZ), '-', 'Color', cm(iS,:), 'LineWidth', lw);
            
        end
        ax.YLim = [0 max(max(pan.hc(:,iSURFsubset,iZ)))];
        ax.YTick = [];
        ax.XMinorTick = 'on';
        ax.XMinorGrid = 'off';
        if iZ < nZ
            ax.XTickLabel = {};
        else
            ax.XLabel.Interpreter = 'latex';
            ax.XLabel.String = pan.XLabel;
        end
        
        if iP == 1
            ax.YLabel.String = 'Frequency';
            ax.YAxisLocation = 'right';
            txt = sprintf('$z$($180^\\circ$E, $0^\\circ$N) = $-%.2f$m', z_i0j0(iZ));
            text(ax, -.03, .2, txt, 'Units','n', 'FontSize', fs, 'Interpreter','latex', 'Rotation', 90);
            text(ax, -.03+.015, .2, '~', 'Units', 'norm', 'FontSize', fs-3, 'Rotation', 90);     % under z
        end
        
    end
end

%% Save figure
fn = sprintf('%shist_Eps,D_%sm_%dsurf', PATH_FIGS, strrep(num2str(list_zref), ' ', '_'), nSs);
export_fig(hf, fn, '-pdf');
close(hf)

%% --- Figure: maps of neutral errors
SURFsubset = {'SIGMA', 'DELTA', 'GAMMA', 'ORTHO', 'OMEGA', 'TOPOB'};
nSs = length(SURFsubset);
fs1 = 15;
fs2 = 12;
OPTS_FIGS.CLIM = [-13 -2];
OPTS_FIGS.CM = parula(256);
cbticks = unique([OPTS_FIGS.CLIM(1), ceil(OPTS_FIGS.CLIM(1)) : floor(OPTS_FIGS.CLIM(2)), OPTS_FIGS.CLIM(2)]);
cblabel = arrayfun(@(x) sprintf('10^{%d}', x), cbticks, 'UniformOutput', false);
for iZ = 1 : nZ
    zref = list_zref(iZ);
    if zref <= 1000
        OPTS_FIGS.LATLIM = [-80, 90]; % The z(180,0) = 1000 surfaces include the Arctic.
    else
        OPTS_FIGS.LATLIM = [-80, 72]; % The z(180,0) = 2000 surfaces exclude the Arctic.
    end
    
    hf = figure('Position', [0 0 800 600], 'Visible', OPTS_FIGS.VIS);
    
    margin_v = .05;
    margin_l = .07;
    margin_r = .12;
    subaxopts = {'MarginLeft',margin_l,'MarginRight', margin_r, 'MarginTop', margin_v, 'MarginBottom', margin_v, ...
        'SpacingVertical', 0.035, 'SpacingHorizontal', .01};
    
    for iS = 1 : nSs
        
        surfname = SURFsubset{iS};
        
        ax = subaxis(ceil((nS-1)/2),2,iS,subaxopts{:});
        fig_map(ax, g.XCvec, g.YCvec, log10(abs(fdd.map(:,:,iSURF(surfname),iZ))), land, OPTS_FIGS);
        colorbar(ax, 'off');
        ax.FontSize = fs2;
        
        if iS <= nSs - 2
            ax.XTickLabel = {};
        end
        if mod(iS,2) == 0
            ax.YTickLabel = {};
        end
        
        if iS == 1
            cb = colorbar('Position', [1-margin_r+.01, margin_v, margin_r/2.5, 1-2*margin_v], ...
                'FontSize', fs2, 'Ticks', cbticks, 'TickLabels', cblabel);
        end
        
        if strcmp(surfname, 'SIGMA')
            txt = sprintf('_{%g}', zref/1000);
        else
            txt = '';
        end
        text(ax, 0, 1.06, ...
            sprintf('$(%s) \\; %s%s$', char('a' + iS-1), list_surf{iSURF(surfname)}, txt), ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', fs1, 'Interpreter', 'latex');
        text(ax, .98, 1.06, ...
            sprintf('$||D^f||_1=$%.2e, $||D^f||_2=$%.2e.', fdd.L1(iSURF(surfname),iZ), fdd.L2(iSURF(surfname),iZ)), ...
            'Units', 'norm', 'HorizontalAlignment', 'right', 'FontSize', fs2, 'Interpreter', 'latex');
        
    end
    
    %% Save figure
    fn = sprintf('FDD_%dsurfs_(%d,%d,%.4fm)', nSs, i0, j0, z_i0j0(iZ));
    export_fig(hf, [PATH_FIGS fn], '-jpg', '-m3');
    close(hf);
    
end

%% --- END NEUTRAL CALCULATIONS -----------------------------------------------
% clear eps fdd

%% --- BEGIN GEOSTROPHIC STREAMFUNCTION --------------------------------------

SURFsubset = {'SIGMA', 'MTOPO', 'OMEGA'}; % only do geostrf analysis on these surfaces
nSs = length(SURFsubset); % number of Surfaces in this subset

EQBAND = 1; % Ignore 1 deg band either side of equator

% For orthobaric Montgomery and topobaric gsf
OPTS_STRF = struct();
OPTS_STRF.WRAP = g.WRAP;

list_geos_utilde = {...
    '\utilde{\mathbf{u}}', ...
    '\utilde{\mathbf{u}_g}^{(ZH)}', ...
    '\utilde{\mathbf{u}_g}^{(Cu)}', ...
    '\utilde{\mathbf{u}_g}^{(MK)}', ...
    '\utilde{\mathbf{u}_g}^{(oM)}', ...
    '\utilde{\mathbf{u}_g}^{(tb)}', ...
    '\utilde{\mathbf{u}_g}^{(\nabla_z p)}', ...
    };
list_geos = {...
    '\mathbf{u}', ...
    '\mathbf{u}_g^{(ZH)}', ...
    '\mathbf{u}_g^{(Cu)}', ...
    '\mathbf{u}_g^{(MK)}', ...
    '\mathbf{u}_g^{(oM)}', ...
    '\mathbf{u}_g^{(tb)}', ...
    '\mathbf{u}_g^{(\nabla_z p)}', ...
    };
list_geosu = {...
    'u', ...
    'u_g^{(ZH)}', ...
    'u_g^{(Cu)}', ...
    'u_g^{(MK)}', ...
    'u_g^{(oM)}', ...
    'u_g^{(tb)}', ...
    'u_g^{(\nabla_z p)}', ...
    };
iSTRF = containers.Map({'AG', 'ZH', 'CU', 'MK', 'OM', 'TB', 'ZS'}, 1:length(list_geos));

cori_V = avg_i_ip1(g.cori_YG); % Coriolis parameter on V grid
cori_U = avg_j_jp1(g.cori_YG); % Coriolis parameter on U grid

nF = length(list_geos); % number of streamFunctions
ex = nan(ni,nj,nF,nSs,nZ); % zonal geostrophic velocity error
ey = nan(ni,nj,nF,nSs,nZ); % meridional geostrophic velocity error

% Load true velocities, then average them onto opposite grid
% (U to V, V to U) as the model does to calculate the Coriolis force,
% then permute them so the data for each water column is contiguous
UVEL = load_ECCO2(PATH_ECCO2, 'U', TIMESTEP);
VVEL = load_ECCO2(PATH_ECCO2, 'V', TIMESTEP);
UVEL = UVEL .* g.DYGsc;
UVEL = UVEL + jm1(UVEL);
UVEL = UVEL + ip1(UVEL);
UVEL = 0.25 * UVEL ./ g.DYCsc;
VVEL = VVEL .* g.DXGvec;
VVEL = VVEL + im1(VVEL);
VVEL = VVEL + jp1(VVEL);
VVEL = 0.25 * VVEL ./ g.DXCvec;
UVEL = permute(UVEL, [3 1 2]); % [nz,nx,ny]. depth  by  long  by  lat
VVEL = permute(VVEL, [3 1 2]);

Y = hsap3(Z, ATMP, ETAN, R, grav, rho_c); % Pre-compute hydrostatic acceleration potential

%% Calculate geostrophic streamfunctions
mytic = tic;
fprintf(fileID, '--- Calculate geostrophic streamfunctions ---\n');
for iZ = 1 : nZ
    for iS = 1 : nSs
        % --- Prepare properties on the surface
        surfname = SURFsubset{iS};
        z = zs(:,:,iSURF(surfname),iZ);
        
        if strcmp(surfname, 'GAMMA') || (strcmp(surfname, 'OMEGA') && ~OMEGA_FRESH)
            % OMEGA-surface was pre-computed using linear interpolation
            % GAMMA-surface was linearly interpolated in this script
            interpfn = @interp1qn2;
        else
            interpfn = @interp_firstdim_twovars;
        end
        [s, t] = interpfn(lead1(z), Z, S, T);
        
        [y,r] = hsap2(s, t, z, Z, R, Y, grav, rho_c);
        
        % Geostrophic velocities from in-surface gradients. Boussinesq. (-) on z undoes that we've been z>0 people!
        ugs = -((y - jm1(y)) + (grav / (2*rho_c)) * (-z + jm1(z)) .* (r + jm1(r))) ./ (g.DYCsc * cori_V);
        vgs =  ((y - im1(y)) + (grav / (2*rho_c)) * (-z + im1(z)) .* (r + im1(r))) ./ (g.DXCvec .* cori_U);
        
        % Zonal velocities:
        z_avg = (z + jm1(z)) / 2; % on V grid
        u = interp1qn(lead1(z_avg), Z, UVEL); % Full velocity, interpolated onto the surface
        y     =     hsap2(S, T,     z_avg , Z, R, Y, grav, rho_c, interpfn) ;
        y_adj = jm1(hsap2(S, T, jp1(z_avg), Z, R, Y, grav, rho_c, interpfn));
        ugz = -(y - y_adj) ./ (g.DYCsc * cori_V); % Geostrophic velocity from z level gradient
        
        % Meridional velocities:
        z_avg = (z + im1(z)) / 2; % on U grid
        v = interp1qn(lead1(z_avg), Z, VVEL); % Full velocity, interpolated onto the surface
        y     =     hsap2(S, T,     z_avg , Z, R, Y, grav, rho_c, interpfn) ;
        y_adj = im1(hsap2(S, T, ip1(z_avg), Z, R, Y, grav, rho_c, interpfn));
        vgz =  (y - y_adj) ./ (g.DXCvec .* cori_U); % Geostrophic velocity from z level gradient
        
        % Geostrophic streamfunctions:
        s0 = s(i0,j0);
        t0 = t(i0,j0);
        z0 = z(i0,j0);
        ZH = zhanghogg92(       s, t, z, Z, R, Y, s0, t0, z0, grav, rho_c);
        MK = mcdougallklocker10(s, t, z, Z, R, Y, s0, t0, z0, grav, rho_c);
        CU = cunningham00(      s, t, z, Z, R, Y,         z0, grav, rho_c);
        
        % Orthobaric Montgomery potential:
        % Ensure orthobaric_montgomery chooses its own reference S and T!
        % Do this by passing [] for s0 and t0.
        OPTS_STRF.SPLINE_PIECES = 12; % 12 pieces in spline
        OPTS_STRF.SPLINE_ORDER = 4; % cubic splines
        OM = orthobaric_montgomery(s, t, z, Z, R, Y, [], [], OPTS_STRF, grav, rho_c);
        
        % Topobaric geostrophic streamfunction:
        OPTS_STRF.REEB = true; % use multivalued functions (topobaric not orthobaric)
        OPTS_STRF.GEOSTRF = true; % ensure topobaric gsf is well-defined
        OPTS_STRF.REF_IJ = OPTS.REF_IJ;
        if strcmp(surfname, 'MTOPO')
            OPTS_STRF.RG = mtopo_data(iZ).RG;
            OPTS_STRF.dfn = mtopo_data(iZ).dfn;
            s0 = mtopo_data(iZ).s0;
            t0 = mtopo_data(iZ).t0;
        end
        TB = topobaric_geostrf(s, t, z, Z, R, Y, s0, t0, OPTS_STRF, grav, rho_c);
        if strcmp(surfname, 'MTOPO')
            OPTS_STRF = rmfield(OPTS_STRF, {'RG', 'dfn'});
        end
        
        %{
        % Map a geostrophic streamfunction, just as a simple example
        if strcmp(surfname, 'OMEGA') && iZ == 2
            hf = figure;
            ax = axes;
            fig_map(ax, g.XCvec, g.YCvec, TB, land, ...
                setfield(setfield(rmfield(OPTS_FIGS, {'CLIM', 'CM'}), 'XSH', 0), 'LATLIM', [-80 72]));
            colorbar;
            cl = caxis(ax);
            contour(ax, g.XCvec, g.YCvec, TB.', [-5.3, -7.3 -9.3 -11.3 -13.3 -15.3], '-k');
            fn = sprintf('TB_on_OMEGA_%d_(%d,%d,%.4fm)', list_zref(iZ), x0, y0, z_i0j0(iZ));
            export_fig(hf, fn, '-jpg', '-m2');
            close(hf);
        end
        %}
        
        % --- Calculate error from "true" geostrophic velocity (ugs, vgs)
        ex(:,:,iSTRF('AG'),iS,iZ) = u  -  ugs;
        ex(:,:,iSTRF('ZH'),iS,iZ) = -dTdy_on_V(ZH) ./ cori_V  -  ugs;
        ex(:,:,iSTRF('CU'),iS,iZ) = -dTdy_on_V(CU) ./ cori_V  -  ugs;
        ex(:,:,iSTRF('MK'),iS,iZ) = -dTdy_on_V(MK) ./ cori_V  -  ugs;
        ex(:,:,iSTRF('OM'),iS,iZ) = -dTdy_on_V(OM) ./ cori_V  -  ugs;
        ex(:,:,iSTRF('TB'),iS,iZ) = -dTdy_on_V(TB) ./ cori_V  -  ugs;
        ex(:,:,iSTRF('ZS'),iS,iZ) = ugz  -  ugs;
        
        ey(:,:,iSTRF('AG'),iS,iZ) = v  -  vgs;
        ey(:,:,iSTRF('ZH'),iS,iZ) = +dTdx_on_U(ZH) ./ cori_U  -  vgs;
        ey(:,:,iSTRF('CU'),iS,iZ) = +dTdx_on_U(CU) ./ cori_U  -  vgs;
        ey(:,:,iSTRF('MK'),iS,iZ) = +dTdx_on_U(MK) ./ cori_U  -  vgs;
        ey(:,:,iSTRF('OM'),iS,iZ) = +dTdx_on_U(OM) ./ cori_U  -  vgs;
        ey(:,:,iSTRF('TB'),iS,iZ) = +dTdx_on_U(TB) ./ cori_U  -  vgs;
        ey(:,:,iSTRF('ZS'),iS,iZ) = vgz  -  vgs;
        
    end %iS
end % iZ
clear ZH CU MK TB OM ugs ugz u vgs vgz v z_avg y y_adj r s t
fprintf(fileID, 'Streamfunctions calculated in time %.2f\n', toc(mytic));

%% Build a common mask all errors, for each surface
% (The orthobaric and topobaric geostrophic streamfunctions have different
% masks to the others, since they internally select a single connected
% component of the surface.)
good_mask_ex = false(ni,nj,1,nSs,nZ);
good_mask_ey = false(ni,nj,1,nSs,nZ);
for iZ = 1:nZ
    for iS = 1:nSs
        goodX = all(isfinite(ex(:,:,:,iS,iZ)), 3) | cori_V == 0;
        goodY = all(isfinite(ey(:,:,:,iS,iZ)), 3);
        [~,~,~,goodX] = bfs_conncomp(goodX, neigh, I0);
        [~,~,~,goodY] = bfs_conncomp(goodY, neigh, I0);
        goodX(abs(g.YG) < EQBAND) = false;
        goodY(abs(g.YC) < EQBAND) = false;
        good_mask_ex(:,:,1,iS,iZ) = goodX;
        good_mask_ey(:,:,1,iS,iZ) = goodY;
    end
end

%% Calculate norms and histograms of geostrophic velocity errors
geos.L1 = nan(nF,nSs,nZ);
geos.L2 = nan(nF,nSs,nZ);
geos.L1x = nan(nF,nSs,nZ);
geos.L2x = nan(nF,nSs,nZ);
geos.edges = linspace(-11, 0.5, 100);
geos.bins = (geos.edges(1:end-1) + geos.edges(2:end)) / 2;
geos.hc = nan(length(geos.bins),nF,nSs,nZ);
geos.wp = nan(5,nF,nSs,nZ);
for iZ = 1:nZ
    for iS = 1:nSs
        
        qq = zeros(5, nF);
        
        AREA = [g.RAS(good_mask_ex(:,:,1,iS,iZ)); g.RAW(good_mask_ey(:,:,1,iS,iZ))];
        for iF = 1 : nF
            x = ex(:,:,iF,iS,iZ);
            y = ey(:,:,iF,iS,iZ);
            err = [x(good_mask_ex(:,:,1,iS,iZ)); y(good_mask_ey(:,:,1,iS,iZ))];
            geos.L1(iF,iS,iZ) = sum(     AREA .* abs(err)) / sum(AREA);
            geos.L2(iF,iS,iZ) = sqrt(sum(AREA .*   err.^2) / sum(AREA));
            
            % Un-weighted histogram
            %geos.hc(:,iF,iS,iZ) = histcounts(log10(abs(err)), geos.edges);
            
            % Area-weighted histogram
            [~,~,bins] = histcounts(log10(abs(err)), [-inf, geos.edges(2:end-1), inf]);
            geos.hc(:,iF,iS,iZ) = accumarray(bins, AREA, [length(geos.bins), 1]);
            
            % Area-weighted percentiles
            geos.wp(:,iF,iS,iZ) = wprctile(log10(abs(err)), [5 25 50 75 95], AREA);
        end
        
        AREA = g.RAS(good_mask_ex(:,:,1,iS,iZ));
        for iF = 1 : nF
            x = ex(:,:,iF,iS,iZ);
            err = x(good_mask_ex(:,:,1,iS,iZ));
            geos.L1x(iF,iS,iZ) = sum(     AREA .* abs(err)) / sum(AREA);
            geos.L2x(iF,iS,iZ) = sqrt(sum(AREA .*   err.^2) / sum(AREA));
        end
    end
end
geos.hc = geos.hc ./ sum(geos.hc,1); % Normalize

%% --- Table (for Latex): L1 and L2 norms of geostrophic velocity error
fprintf(fileID, '--- Table of L1 and L2 norms of geostrophic velocity error ---\n');
fprintf(fileID, '\\tabcolsep=0.1cm\n');
for l = 1:2
    fprintf(fileID, '\\begin{center}\n');
    fprintf(fileID, '\\begin{tabular}{|l|ccc|ccc|}\n');
    fprintf(fileID, '\\hline\n');
    fprintf(fileID, ' (%s) & \\multicolumn{3}{c|}{$\\zn(\\bm{x}_\\rv) = \\SI{-%.2f}{m}$} & \\multicolumn{3}{c|}{$\\zn(\\bm{x}_\\rv) = \\SI{-%.2f}{m}$} \\\\ \n', ...
        char('a' + l - 1), z_i0j0(1), z_i0j0(2));
    fprintf(fileID, '$||\\mathbf{\\epsilon}||_%d$', l);
    for iZ = 1 : nZ
        zref = list_zref(iZ);
        for iS = 1:nSs
            surfname = SURFsubset{iS};
            if strcmp(SURFsubset{iS}, 'SIGMA')
                txt = sprintf('_{%g}', zref/1000);
            else
                txt = '';
            end
            fprintf(fileID, '& $%s%s$ ', list_surf{iSURF(surfname)}, txt);
        end
    end
    fprintf(fileID, '\\\\[5pt]\n');
    fprintf(fileID, '\\hline\n');
    fprintf(fileID, '& & & & & & \\\\[-7pt]\n');
    for iF = 1 : nF
        fprintf(fileID, '$%s$ ', list_geos_utilde{iF});
        for iZ = 1 : nZ
            for iS = 1:nSs
                fprintf(fileID, '& %.1e ', geos.(['L' num2str(l)])(iF,iS,iZ));
            end
        end
        fprintf(fileID, '\\\\[5pt]\n');
    end
    fprintf(fileID, '\\hline\n');
    fprintf(fileID, '\\end{tabular}\n');
    fprintf(fileID, '\\end{center}\n');
end

%% --- Table (for Latex): Fraction of area with geostrophic error > ageostrophic velocity
fprintf(fileID, '--- Fraction of area with geostrophic error > ageostrophic velocity ---\n');
fprintf(fileID, '\\tabcolsep=0.28cm\n');
fprintf(fileID, '\\begin{center}\n');
fprintf(fileID, '\\begin{tabular}{|l|ccc|ccc|}\n');
fprintf(fileID, '\\hline\n');
fprintf(fileID, '& \\multicolumn{3}{c|}{$\\zn(\\bm{x}_\\rv) = \\SI{-%.2f}{m}$} & \\multicolumn{3}{c|}{$\\zn(\\bm{x}_\\rv) = \\SI{-%.2f}{m}$} \\\\ \n', z_i0j0(1), z_i0j0(2));
for iZ = 1 : nZ
    zref = list_zref(iZ);
    for iS = 1:nSs
        surfname = SURFsubset{iS};
        if strcmp(SURFsubset{iS}, 'SIGMA')
            txt = sprintf('_{%g}', zref/1000);
        else
            txt = '';
        end
        fprintf(fileID, '& $%s%s$ ', list_surf{iSURF(surfname)}, txt);
    end
end
fprintf(fileID, '\\\\[5pt]\n');
fprintf(fileID, '\\hline\n');
fprintf(fileID, '& & & & & & \\\\[-7pt]\n');
for iF = cellfun(@(c) iSTRF(c), {'ZH', 'CU', 'MK', 'OM', 'TB', 'ZS'}) % all but AG
    fprintf(fileID, '$%s$ ', list_geos_utilde{iF});
    for iZ = 1 : nZ
        for iS = 1:nSs
            x = abs(ex(:,:,iF,iS,iZ)) > abs(ex(:,:,iSTRF('AG'),iS,iZ));
            y = abs(ey(:,:,iF,iS,iZ)) > abs(ey(:,:,iSTRF('AG'),iS,iZ));
            err = [x(good_mask_ex(:,:,1,iS,iZ)); y(good_mask_ey(:,:,1,iS,iZ))];
            AREA = [g.RAS(good_mask_ex(:,:,1,iS,iZ)); g.RAW(good_mask_ey(:,:,1,iS,iZ))];
            
            err = sum(AREA .* err) / sum(AREA);
            fprintf(fileID, '& %.2f ', err*100);
        end
    end
    fprintf(fileID, '\\\\[5pt]\n');
end
fprintf(fileID, '\\hline\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{center}\n');
clear good bigerr


%% --- Figure: box plots of errors, and norms
bw = .3; % box width
fs = 12; % font size
ml = .1; % margin left
mr = .01; % margin right
mb = .06; % margin bottom
mt = .03; % margin top
sh = .015; % spacing horiz
sv = .04; % spacing vert
hf = figure('Position', [0 0 400*nZ 266*nSs], 'Visible', OPTS_FIGS.VIS);

%yl = [-7.5 -1.5; -7.5 -1.5; -9 -1.5];
%yl = [-8.5 -1; -8.5 -1; -8.5 -1];
yl = [-9 -1; -11 -1; -9 -1];
w = (1-ml-mr-sh*nZ ) / nZ ; % axis width
h = nan(nSs,1); % axis height, variable. Vert grid spacing will be uniform.
for i = 1:nSs
    h(i) = (1-mb-mt-sv*nSs) * diff(yl(i,:)) / sum(diff(yl,[],2));
end

for iZ = 1:nZ
    zref = list_zref(iZ);
    
    for iS = 1:nSs
        surfname = SURFsubset{iS};
        
        l = ml + (sh + w) * (iZ - 1);
        b = mb + (sv    ) * (nSs - iS) + sum(h(iS+1:nSs));
        ax = axes('Position', [l b w h(iS)]);
        
        hold(ax,'on'); box(ax,'on'); %grid(ax,'on')
        ax.YAxis.FontSize = fs;
        ax.XLim = [1-2*bw, nF+2*bw];
        ax.YLim = yl(iS,:);
        ax.YTick = ceil(yl(iS,1)) : floor(yl(iS,2));
        
        % Faint horiz grid lines
        line(ax, repmat(ax.XLim()', 1, length(ax.YTick)), ax.YTick + [0;0], 'Color', [1 1 1]*.8, 'LineWidth', .5);
        
        line(ax, (1:nF) + [-1; 1]*bw, geos.wp(2,:,iS,iZ) + [0;0], 'Color','k'); % Horiz line at .25
        line(ax, (1:nF) + [-1; 1]*bw, geos.wp(4,:,iS,iZ) + [0;0], 'Color','k'); % Horiz line at .75
        line(ax, (1:nF) + [-1;-1]*bw, geos.wp([2 4],:,iS,iZ), 'Color','k'); % Left  vert line .25 to .75
        line(ax, (1:nF) + [ 1; 1]*bw, geos.wp([2 4],:,iS,iZ), 'Color','k'); % Right vert line .25 to .75
        line(ax, (1:nF) + [-1; 1]*bw, geos.wp(3,:,iS,iZ) + [0;0], 'Color','k'); % Horiz line at .5 (median)
        line(ax, (1:nF) + [0; 0], geos.wp(1:2,:,iS,iZ), 'Color', 'k'); % lower whisker: .05 to .25
        line(ax, (1:nF) + [0; 0], geos.wp(4:5,:,iS,iZ), 'Color', 'k'); % upper whisker: .75 to .95
        
        scatter(ax, (1:nF), log10(geos.L1(:,iS,iZ)'), 48, 'o', 'MarkerEdgeColor', [1 1 1]*0); % L1 norm
        scatter(ax, (1:nF), log10(geos.L2(:,iS,iZ)'), 72, 'x', 'MarkerEdgeColor', [0 0 0]  ); % L2 norm
        
        if iZ == 1
            ax.YLabel.Interpreter = 'latex';
            ax.YLabel.String = '$\mathbf{\epsilon} \quad [\mathrm{m} \ \mathrm{s}^{-1}]$';
            ax.YLabel.FontSize = 17;
            ax.YTickLabel = arrayfun(@(x) sprintf('10^{%d}', x), ax.YTick, 'UniformOutput', false);
        else
            ax.YTickLabel = {};
        end
        if iS == iSURF('SIGMA')
            txt = sprintf('$%s_{%g}$', list_surf{iSURF(surfname)}, zref/1000);
        else
            txt = ['$' list_surf{iSURF(surfname)} '$'];
        end
        text(ax, .02, .05, txt, 'Units', 'n', 'FontSize', fs, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
        
        ax.XTickLabel = {};
        if iS == 1
            txt = sprintf('$z$($180^\\circ$E, $0^\\circ$N) = $-%.2f$m', z_i0j0(iZ));
            text(ax, .3, 1.06, txt, 'Units','n', 'FontSize', fs, 'Interpreter','latex');
            text(ax, .3, 1.06-.07, '$\widetilde{\phantom{z}}$', 'Units', 'norm', 'FontSize', fs, 'Interpreter','latex');     % under z
        elseif iS == nSs
            y = ax.YLim(1) - diff(ax.YLim)*.13;
            y_ = ax.YLim(1) - diff(ax.YLim)*.17;
            for iF = cellfun(@(c) iSTRF(c), {'ZH', 'CU', 'MK', 'OM', 'TB', 'ZS'}) % all but AG
                text(iF, y, ['$' list_geos{iF} '$'], 'FontSize', fs, 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                text(iF - .22, y_, '$\widetilde{\phantom{\mathbf{u_g}}}$', 'FontSize', fs', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
            end
            iF = iSTRF('AG');
            text(iF, y, '$\mathbf{u}$', 'FontSize', fs, 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(iF, y_, '$\widetilde{\phantom{\mathbf{u}}}$', 'FontSize', fs, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
        end
    end
end

%% Save figure
fn = sprintf('%sboxplot%d%d_GeosVelErr_EQBAND%d', PATH_FIGS, nSs, nF, EQBAND);
export_fig(hf, fn, '-pdf');
close(hf)

%% --- Figures: Maps of zonal geostrophic velocity error
deg = char(176);
fs = 10;   % font size
ml = .06;  % margin left
mr = .09;  % margin right
mb = .04;  % margin bottom
mt = .03;  % margin top
sh = .01;  % spacing horiz
sv = .035; % spacing vert

subaxopts = {'marginleft',ml,'marginright',mr,'marginbottom',mb,'margintop',mt, ...
    'spacinghorizontal', sh, 'spacingvertical', sv};

OPTS_FIGS.CM = parula(128);
OPTS_FIGS.CLIM = [-8 -1];
OPTS_FIGS.LATLIM = [-80, 90]; % The z(180,0) = 1000 surfaces include the Arctic.

cbticks = unique([OPTS_FIGS.CLIM(1), ceil(OPTS_FIGS.CLIM(1)) : floor(OPTS_FIGS.CLIM(2)), OPTS_FIGS.CLIM(2)]);
cblabel = arrayfun(@(x) sprintf('10^{%d}', x), cbticks, 'UniformOutput', false);

for iS = 1:nSs
    surfname = SURFsubset{iS};
    for iZ = 1:nZ
        zref = list_zref(iZ);
        
        hf = figure('Position', [0 0 800 round(200 * ceil(nF/2))], 'Visible', OPTS_FIGS.VIS);
        
        for iF =  1 : nF
            iFp = iF + (iF == nF);
            err = ex(:,:,iF,iS,iZ);
            err(~good_mask_ex(:,:,1,iS,iZ)) = nan;
            
            ax = subaxis(ceil(nF/2), 2, iFp, subaxopts{:});
            fig_map(ax, g.XCvec, g.YCvec, log10(abs(err)), land, OPTS_FIGS);
            colorbar(ax, 'off');
            
            if mod(iFp,2) == 0 && iFp < nF % left panels
                ax.YTickLabel = {};
            end
            if iFp <= 4 || iFp == 6  % non-bottom panels
                ax.XTickLabel = {};
            end
            
            text(ax, -.01, 1.06, ['$(' char('a' + iF-1) ')$'], ...
                'Units', 'norm', 'HorizontalAlignment', 'left', 'FontSize', fs, 'Interpreter', 'latex');
            text(ax, .05, 1.06, ['$\epsilon^{(x)} = ' list_geosu{iF} ' - u_g.$'], ...
                'Units', 'norm', 'HorizontalAlignment', 'left', 'FontSize', fs, 'Interpreter', 'latex');
            text(ax, .05, .99, '$\phantom{\epsilon^{(x)} = }\, \widetilde{\phantom{u_g.}}$', ...
                'Units', 'norm', 'HorizontalAlignment', 'left', 'FontSize', fs, 'Interpreter', 'latex');
            text(ax, .05, .99, ['$\phantom{\epsilon^{(x)} = ' list_geosu{iF} ' - } \, \widetilde{\phantom{u_g}}$'], ...
                'Units', 'norm', 'HorizontalAlignment', 'left', 'FontSize', fs, 'Interpreter', 'latex');
            text(ax, 1, 1.06, ...
                sprintf('$||\\epsilon^{(x)}||_1=$%.2e, $||\\epsilon^{(x)}||_2=$%.2e.', geos.L1(iF,iS,iZ), geos.L2(iF,iS,iZ)), ...
                'Units', 'norm', 'HorizontalAlignment', 'right', 'FontSize', fs, 'Interpreter', 'latex');
        end
        colorbar(ax, 'Position', [1-mr+.01 mb .0333 1-mb-mt], 'FontSize', fs, ...
            'Ticks', cbticks, 'TickLabels', cblabel);
        
        
        % --- Scatter plot of delta vs z, only in the main region
        % (connected to OPTS.REF_IJ)
        z = zs(:,:,iSURF(surfname),iZ);
        if strcmp(surfname, 'GAMMA') || (strcmp(surfname, 'OMEGA') && ~OMEGA_FRESH)
            % OMEGA-surface was pre-computed using linear interpolation
            % GAMMA-surface was linearly interpolated in this script
            interpfn = @interp1qn2;
        else
            interpfn = @interp_firstdim_twovars;
        end
        [s, t] = interpfn(lead1(z), Z, S, T);
        
        [~,~,~,ROI] = bfs_conncomp(isfinite(z), neigh, I0);
        s_ = s(ROI); clear s
        t_ = t(ROI); clear t
        z_ = z(ROI);
        lat_ = g.YC(ROI);
        
        % Re-order data so first data is not always overlain by later data
        rp = randperm(length(s_));
        s_ = s_(rp);
        t_ = t_(rp);
        z_ = z_(rp);
        lat_ = lat_(rp);
        
        [~,Ideep] = max(z_);
        S_delta = s_(Ideep);
        T_delta = t_(Ideep);
        
        d_ = eos(s_, t_, z_) - eos(S_delta, T_delta, z_);
        
        num_pieces = 12;
        order = 4;
        breaks = linspace(min(d_), max(d_), num_pieces + 1);
        fnOM = splinefit(d_, -z_, breaks, order);
        
        ax = subaxis(ceil(nF/2), 2, nF, subaxopts{:});
        ax.Position = ax.Position + [0.02, .02, -.07, -.04];
        
        scatter(ax, d_, -z_, 4, lat_, 'MarkerFaceColor', 'flat');
        hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
        colormap(ax, jet)
        ax.XLim = [min( d_), max( d_)];
        ax.YLim = [min(-z_), 0];
        ax.XLim = [ax.XLim(1), ax.XLim(2) + diff(ax.XLim)*.28];
        ax.CLim = [-1 1] * max(abs(lat_));
        hcb = colorbar(ax, 'Location', 'east', 'FontSize', fs-1);
        SNstr = 'S N';
        hcb.Ticks = -60 : 30 : 60;
        hcb.TickLabels = arrayfun(@(a) [num2str(abs(a)) deg SNstr(sign(a)+2)], hcb.Ticks, 'UniformOutput', false);
        ax.FontSize = 12;
        
        text(ax, .4, -.2, '$\delta$  [kg m$^{-3}$]', 'Interpreter', 'latex', 'units', 'norm', 'fontsize', fs);
        text(ax, .4, -.3, '$\widetilde{\phantom{\delta}}$', 'units', 'norm', 'interpreter', 'latex', 'fontsize', fs)
        text(ax, -.15, .8, '$z$  [m]', 'units', 'norm', 'interpreter', 'latex', 'fontsize', fs, 'Rotation', 90)
        text(ax, -.11, .8, '$\widetilde{\phantom{z}}$', 'units', 'norm', 'interpreter', 'latex', 'fontsize', fs, 'Rotation', 90)
        
        text(ax, -.125, 1.1, '$(*)$', ...
            'Units', 'norm', 'HorizontalAlignment', 'left', 'FontSize', fs, 'Interpreter', 'latex');
        
        dd = linspace(min(d_), max(d_), 200);
        zz = ppval(fnOM, dd);
        plot(ax, dd, zz, '-', 'LineWidth', 1.5, 'Color', [1 1 1]*.0);
        
        %% Save figure
        fn = sprintf('GeosVelUErr_%d_on_%s_%d_EQBAND%d', nF, list_surf{iSURF(surfname)}(2:end), zref, EQBAND);
        export_fig(hf, [PATH_FIGS fn] , '-jpg', '-m3');
        % Close figure
        close(hf);
        
    end
end
clear rp s_ t_ z_ lat_ S_delta T_delta d_ num_pieces order breaks fnOM dd zz

%% Report: Geostrophic velocity errors when streamfunction is ill-defined
fprintf(fileID, '--- Geostrophic velocity errors when streamfunction is ill-defined ---\n');
iZ = 2;
surfname = 'OMEGA';
z = zs(:,:,iSURF(surfname),iZ);
if strcmp(surfname, 'GAMMA') || (strcmp(surfname, 'OMEGA') && ~OMEGA_FRESH)
    % OMEGA-surface was pre-computed using linear interpolation
    % GAMMA-surface was linearly interpolated in this script
    interpfn = @interp1qn2;
else
    interpfn = @interp_firstdim_twovars;
end
[s, t] = interpfn(lead1(z), Z, S, T);
s0 = s(i0,j0);
t0 = t(i0,j0);

OPTS_STRF.REEB = true; % use multivalued functions (topobaric not orthobaric)
OPTS_STRF.GEOSTRF = false; % Make geostrophic streamfunction is ill-defined
[~,TBdiff] = topobaric_geostrf(s, t, z, Z, R, Y, s0, t0, OPTS_STRF, grav, rho_c);

errX = -dTdy_on_V(TBdiff) ./ cori_V;
errY = +dTdx_on_U(TBdiff) ./ cori_U;

goodX = isfinite(errX) | cori_V == 0;
[~,~,~,goodX] = bfs_conncomp(goodX, neigh, I0);
goodX(abs(g.YG) < EQBAND) = false;

goodY = isfinite(errY);
[~,~,~,goodY] = bfs_conncomp(goodY, neigh, I0);
goodY(abs(g.YC) < EQBAND) = false;

jumpX = errX(goodX & errX ~= 0);
jumpY = errY(goodY & errY ~= 0);
nJump = length(jumpX) + length(jumpY);
nGood = sum(goodX(:)) + sum(goodY(:));
L1jump = (sum(abs(jumpX)) + sum(abs(jumpY))) / nJump;
L2jump = sqrt( (sum(jumpX.^2) + sum(jumpY.^2)) / nJump);
L1jumpX = sum(abs(jumpX)) / length(jumpX);
L2jumpX = sqrt( sum(jumpX.^2) / length(jumpX));

fprintf(fileID, 'On %s at (%d,%d,%.4fm), the geostrophic velocity from\n', surfname, i0, j0, z_i0j0(iZ));
fprintf(fileID, '  from the ill-defined topobaric geostrophic streamfunction\n');
fprintf(fileID, '  has %d discontinuous data points in the good mask which has %d points (ratio: %.4f) \n', nJump, nGood, nJump / nGood);
fprintf(fileID, '  and of these points, the L1 error is %.4e, and the L2 error is %.4e\n', L1jump, L2jump);

%figure; scatter(g.YG(goodX & errX ~= 0), abs(jumpX), 12, 'k', '.')

%% Compute "orthobaric" geostrophic streamfunction
% which is like the topobaric geostrophic streamfunction but with a single
% region, and a single-valued functional relationship fit between density
% and pressure / depth. Compare it with the orthobaric Montgomery
% potential.
fprintf(fileID, '--- Orthobaric geostrophic streamfunction vs orthobaric Montgomery potential ---\n');
surfname = 'OMEGA';
if strcmp(surfname, 'GAMMA') || (strcmp(surfname, 'OMEGA') && ~OMEGA_FRESH)
    % OMEGA-surface was pre-computed using linear interpolation
    % GAMMA-surface was linearly interpolated in this script
    interpfn = @interp1qn2;
else
    interpfn = @interp_firstdim_twovars;
end
iS = find(cellfun(@(c) strcmp(c, 'OMEGA'), SURFsubset));
assert(~isempty(iS), 'Failed to select OMEGA surface.');
for iZ = 1 : nZ
    
    % --- Prepare properties on the surface
    z = zs(:,:,iSURF(surfname),iZ);
    [s,t] = interpfn(lead1(z), Z, S, T);
    [y,r] = hsap2(s, t, z, Z, R, Y, grav, rho_c);
    
    % Geostrophic velocities from in-surface gradients. Boussinesq. (-) on z undoes that we've been z>0 people!
    ugs = -((y - jm1(y)) + (grav / (2*rho_c)) * (-z + jm1(z)) .* (r + jm1(r))) ./ (g.DYCsc * cori_V);
    vgs =  ((y - im1(y)) + (grav / (2*rho_c)) * (-z + im1(z)) .* (r + im1(r))) ./ (g.DXCvec .* cori_U);
    
    % calculate "Orthobaric" geostrophic streamfunction
    OPTS_STRF.REEB = false;
    zref = list_zref(iZ);
    if zref == 1000
        OPTS_STRF.SPLINE_BREAKS = [0 200 6000];
        OPTS_STRF.SPLINE_ORDER = 4;
    elseif zref == 2000
        OPTS_STRF.SPLINE_BREAKS = [0 200 1500 1800 6000];
        OPTS_STRF.SPLINE_ORDER = 4;
    end
    %OB = topobaric_geostrf(ATMP, ETAN, S, T, Z, z, OPTS_STRF, R);
    s0 = s(i0,j0);
    t0 = t(i0,j0);
    OB = topobaric_geostrf(s, t, z, Z, R, Y, s0, t0, OPTS_STRF, grav, rho_c);
    
    % Compare OB and OM:
    exOB = -dTdy_on_V(OB) ./ cori_V  -  ugs;
    eyOB = +dTdx_on_U(OB) ./ cori_U  -  vgs;
    exOM = ex(:,:,iSTRF('OM'),iS,iZ);
    eyOM = ey(:,:,iSTRF('OM'),iS,iZ);
    
    goodX = good_mask_ex(:,:,1,iS,iZ);
    goodY = good_mask_ey(:,:,1,iS,iZ);
    AREA = [g.RAS(good_mask_ex(:,:,1,iS,iZ)); g.RAW(good_mask_ey(:,:,1,iS,iZ))];
    L2OM = sqrt( ...
        ( sum(g.RAS(goodX) .* exOM(goodX).^2) + sum(g.RAW(goodY) .* eyOM(goodY).^2) ) ...
        / ...
        ( sum(g.RAS(goodX)) + sum(g.RAW(goodY)) ) );
    L2OB = sqrt( ...
        ( sum(g.RAS(goodX) .* exOB(goodX).^2) + sum(g.RAW(goodY) .* eyOB(goodY).^2) ) ...
        / ...
        ( sum(g.RAS(goodX)) + sum(g.RAW(goodY)) ) );
    
    fprintf(fileID, 'On %s through (%d,%d,%.2f), \n', surfname, i0, j0, z_i0j0(iZ));
    fprintf(fileID, 'for OM, ||eps||_2 = %.4e\n', L2OM);
    fprintf(fileID, 'for OB, ||eps||_2 = %.4e\n', L2OB);
end

%% --- Figure: Histogram of geostrophic errors
lw = 2; % line width
fs = 10; % fontsize
subaxopts = {'marginleft',.01,'marginright',.01,'marginbottom',.06,'margintop',.03, ...
    'spacinghorizontal', .03, 'spacingvertical', .01};

cm = [ ... % RGB color for lines
    0.5       0.5       0.5
    0.0       0.0       1.0
    0.85      0.6       0.0
    1.0       0.0       0.0
    1.0       0.0       1.0
    0.5       0.2       0.7
    0.0       0.5       0.0
    0.0       0.0       0.0
    ];

hf = figure('Position', [0 0 400*nZ 266*nSs], 'Visible', OPTS_FIGS.VIS);
for iZ = 1:nZ
    for iS = 1:nSs
        surfname = SURFsubset{iS};
        
        ax = subaxis(nSs,nZ,iZ,iS,subaxopts{:});
        hold(ax,'on'); box(ax,'on'); grid(ax,'on')
        ax.XLim = geos.edges([1 end]);
        ax.XTick = -10:0;
        ax.FontSize = fs;
        
        if iZ == 1
            x = ax.Position(1) - .005;
        else
            x = ax.Position(1) + .005;
        end
        w = ax.Position(3) * .37;
        h = ax.Position(4) * .7;
        y = ax.Position(2) + ax.Position(4) - .005 - h;
        ax2 = axes('Position', [x y w h]);
        ax2.XTick = []; ax2.YTick = [];
        box(ax2,'on');
        taby = .92;
        
        txt = sprintf('(%s,%s)', repelem('I', 1, iS), repelem('i', 1, iZ));
        text(ax2, .05 , taby, txt, 'units','n', 'HorizontalAlignment', 'left', 'FontSize', fs, 'interpreter', 'latex');
        text(ax2, .45, taby,'$||\mathbf{\epsilon}||_1$' , 'units','n', 'HorizontalAlignment', 'center', 'FontSize', fs, 'interpreter', 'latex');
        text(ax2, .8 , taby,'$||\mathbf{\epsilon}||_2$' , 'units','n', 'HorizontalAlignment', 'center', 'FontSize', fs, 'interpreter', 'latex');
        taby = .1;
        for iF = [nF nF-1,  1:nF-2]
            
            taby = .85 * (nF-iF)/(nF) + .02;
            text(ax2, .03, taby, ['$' list_geos{iF} '$'], 'units', 'n','Color', cm(iF,:), 'HorizontalAlignment', 'left' , 'FontSize', fs, 'Interpreter', 'latex', 'VerticalAlignment', 'bottom');
            text(ax2, .035, taby-.07, '$\widetilde{\phantom{\mathbf{u}}}$', 'units', 'n','Color', cm(iF,:), 'HorizontalAlignment', 'left' , 'FontSize', fs, 'Interpreter', 'latex', 'VerticalAlignment', 'bottom');
            text(ax2, .625, taby, sprintf('%.2e',geos.L1(iF,iS,iZ)), 'units','n','Color',cm(iF,:), 'HorizontalAlignment', 'right', 'FontSize', fs, 'Interpreter', 'latex', 'VerticalAlignment', 'bottom');
            text(ax2, .97 , taby, sprintf('%.2e',geos.L2(iF,iS,iZ)), 'units','n','Color',cm(iF,:), 'HorizontalAlignment', 'right', 'FontSize', fs, 'Interpreter', 'latex', 'VerticalAlignment', 'bottom');
            
            plot(ax, geos.bins, geos.hc(:,iF,iS,iZ), '-', 'Color', cm(iF,:), 'LineWidth', lw);
            
        end
        ax.YLim = [0 max(max(geos.hc(:,:,iS,iZ)))];
        ax.YTick = [];
        if iZ == 2
            ax.YLabel.String = 'Frequency';
        end
        if iS == 1
            txt = sprintf('$z$($180^\\circ$E, $0^\\circ$N) = $-%.2f$m', z_i0j0(iZ));
            text(ax, .3, 1.04, txt, 'Units','n', 'FontSize', fs, 'Interpreter','latex');
            text(ax, .3, 1.04-.045, '$\widetilde{\phantom{z}}$', 'Units', 'norm', 'FontSize', fs, 'Interpreter','latex');     % under z
        elseif iS == nSs
            ax.XLabel.Interpreter = 'latex';
            ax.XLabel.String = '$\log_{10}(|\epsilon|)$';
        end
        if iS < nSs
            ax.XTickLabel = {};
        end
        if iZ == 1
            if iS == iSURF('SIGMA')
                txt = sprintf('$%s_{%g}$', list_surf{iSURF(surfname)}, list_zref(1)/1000);
                text(ax, 1.02, .95, txt, 'Units', 'n', 'FontSize', fs+2, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
                txt = 'or';
                text(ax, 1.02, .88, txt, 'Units', 'n', 'FontSize', fs+2, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
                txt = sprintf('$%s_{%g}$', list_surf{iSURF(surfname)}, list_zref(2)/1000);
                text(ax, 1.02, .81, txt, 'Units', 'n', 'FontSize', fs+2, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
                
            else
                txt = ['$' list_surf{iSURF(surfname)} '$'];
                text(ax, 1.02, .92, txt, 'Units', 'n', 'FontSize', fs+2, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
            end
        end
        ax.XMinorTick = 'on';
        ax.XMinorGrid = 'off';
        
    end
end
%% Save figure
fn = sprintf('%shist%d%d_GeosVelErr_EQBAND%d', PATH_FIGS, nSs, nF, EQBAND);
export_fig(hf, fn, '-pdf');
close(hf)

%% --- END GEOSTROPHIC STREAMFUNCTION ----------------------------------------

%% --- Quick test of non-Boussinesq form for geostrophic streamfunctions
% Switch equation of state to the jmd95 specific volume
copyfile([PATH_PROJECT 'lib' V 'eos' V 'eoscg_specvoljmd95.m'   ], [PATH_PROJECT 'lib' V 'alias' V 'eos.m'  ]);
copyfile([PATH_PROJECT 'lib' V 'eos' V 'eoscg_specvoljmd95_dp.m'], [PATH_PROJECT 'lib' V 'alias' V 'eos_x.m']);
clear eos eos_x % Make sure the copied file gets used

%% Get the hydrostatic pressure (and pretend it's the actual pressure)
P = Y * (rho_c * Pa2db);

%% Get the specific volume
A = eos(S, T, P);

% Pre-compute hydrostatic acceleration potential
Y = hsap3(P, ATMP, ETAN, A, grav);

% Get a potential density surface to work on
pref = 2000; % [dbar]
SIGMA = sort(-eos(S, T, pref), 1); % (-) so that increasing with depth. sort to ensure monotonic
isoval = interp1qn(pref, P(:,i0,j0), SIGMA(:,i0,j0));
p = interp1qn(isoval, SIGMA, P);
clear SIGMA

% Properties on the surface
[s,t] = interp1qn2(lead1(p), P, S, T);
[y,a] = hsap2(s, t, p, P, A, Y);

% Geostrophic velocities from in-surface gradients. Non-Boussinesq.
ugs = -( (db2Pa/2) * (a + jm1(a)) .* (p - jm1(p)) + (y - jm1(y)) ) ./ (g.DYCsc * cori_V);
vgs =  ( (db2Pa/2) * (a + im1(a)) .* (p - im1(p)) + (y - im1(y)) ) ./ (g.DXCvec .* cori_U);

% Zonal velocities:
p_avg = (p + jm1(p)) / 2; % on V grid
u = interp1qn(lead1(p_avg), P, UVEL); % Full velocity, interpolated onto the surface
y     =     hsap2(S, T,     p_avg , P, A, Y) ;
y_adj = jm1(hsap2(S, T, jp1(p_avg), P, A, Y));
ugz = -(y - y_adj) ./ (g.DYCsc * cori_V); % Geostrophic velocity from z level gradient

% Meridional velocities:
p_avg = (p + im1(p)) / 2; % on U grid
v = interp1qn(lead1(p_avg), P, VVEL); % Full velocity, interpolated onto the surface
y     =     hsap2(S, T,     p_avg , P, A, Y) ;
y_adj = im1(hsap2(S, T, ip1(p_avg), P, A, Y));
vgz =  (y - y_adj) ./ (g.DXCvec .* cori_U); % Geostrophic velocity from z level gradient

% Geostrophic streamfunctions:
s0 = s(i0,j0);
t0 = t(i0,j0);
p0 = p(i0,j0);
ZH = zhanghogg92(       s, t, p, P, A, Y, s0, t0, p0);
MK = mcdougallklocker10(s, t, p, P, A, Y, s0, t0, p0);
CU = cunningham00(      s, t, p, P, A, Y,         p0);

% Orthobaric Montgomery potential:
% Ensure orthobaric_montgomery chooses its own reference S and T!
% Do this by passing [] for s0 and t0.
OM = orthobaric_montgomery(s, t, p, P, A, Y, [], [], OPTS_STRF);

% Topobaric geostrophic streamfunction:
OPTS_STRF.REEB = true;
OPTS_STRF.GEOSTRF = true;
TB = topobaric_geostrf(s, t, p, P, A, Y, s0, t0, OPTS_STRF);

% --- Calculate error from "true" geostrophic velocity (ugs, vgs)
ex = nan(ni,nj,nF);
ex(:,:,iSTRF('AG')) = u  -  ugs;
ex(:,:,iSTRF('ZH')) = -dTdy_on_V(ZH) ./ cori_V  -  ugs;
ex(:,:,iSTRF('CU')) = -dTdy_on_V(CU) ./ cori_V  -  ugs;
ex(:,:,iSTRF('MK')) = -dTdy_on_V(MK) ./ cori_V  -  ugs;
ex(:,:,iSTRF('OM')) = -dTdy_on_V(OM) ./ cori_V  -  ugs;
ex(:,:,iSTRF('TB')) = -dTdy_on_V(TB) ./ cori_V  -  ugs;
ex(:,:,iSTRF('ZS')) = ugz  -  ugs;

ey = nan(ni,nj,nF);
ey(:,:,iSTRF('AG')) = v  -  vgs;
ey(:,:,iSTRF('ZH')) = +dTdx_on_U(ZH) ./ cori_U  -  vgs;
ey(:,:,iSTRF('CU')) = +dTdx_on_U(CU) ./ cori_U  -  vgs;
ey(:,:,iSTRF('MK')) = +dTdx_on_U(MK) ./ cori_U  -  vgs;
ey(:,:,iSTRF('OM')) = +dTdx_on_U(OM) ./ cori_U  -  vgs;
ey(:,:,iSTRF('TB')) = +dTdx_on_U(TB) ./ cori_U  -  vgs;
ey(:,:,iSTRF('ZS')) = vgz  -  vgs;

goodX = all(isfinite(ex), 3) | cori_V == 0;
goodY = all(isfinite(ey), 3);
[~,~,~,goodX] = bfs_conncomp(goodX, neigh, I0);
[~,~,~,goodY] = bfs_conncomp(goodY, neigh, I0);
goodX(abs(g.YG) < EQBAND) = false;
goodY(abs(g.YC) < EQBAND) = false;

%% Figure: visually check that errors are similar to the Boussinesq case.
OPTS_FIGS.CM = parula(128);
OPTS_FIGS.CLIM = [-8 -1];
OPTS_FIGS.LATLIM = [-80, 90]; % The z(180,0) = 1000 surfaces include the Arctic.
cbticks = unique([OPTS_FIGS.CLIM(1), ceil(OPTS_FIGS.CLIM(1)) : floor(OPTS_FIGS.CLIM(2)), OPTS_FIGS.CLIM(2)]);
cblabel = arrayfun(@(x) sprintf('10^{%d}', x), cbticks, 'UniformOutput', false);
fs = 10;   % font size
ml = .06;  % margin left
mr = .09;  % margin right
mb = .04;  % margin bottom
mt = .03;  % margin top
sh = .01;  % spacing horiz
sv = .035; % spacing vert
subaxopts = {'marginleft',ml,'marginright',mr,'marginbottom',mb,'margintop',mt, ...
    'spacinghorizontal', sh, 'spacingvertical', sv};
hf = figure('Position', [0 0 800 round(200 * ceil(nF/2))], 'Visible', OPTS_FIGS.VIS);
for iF =  1 : nF
    iFp = iF + (iF == nF);
    err = ex(:,:,iF);
    err(~goodX) = nan;
    ax = subaxis(ceil(nF/2), 2, iFp, subaxopts{:});
    fig_map(ax, g.XCvec, g.YCvec, log10(abs(err)), land, OPTS_FIGS);
    colorbar(ax,'off');
    if mod(iFp,2) == 0 && iFp < nF % left panels
        ax.YTickLabel = {};
    end
    if iFp <= 4 || iFp == 6  % non-bottom panels
        ax.XTickLabel = {};
    end
end
colorbar(ax, 'Position', [1-mr+.01 mb .0333 1-mb-mt], 'FontSize', fs, ...
    'Ticks', cbticks, 'TickLabels', cblabel);

%% Save figure
fn = sprintf('GeosVelUErr_%d_on_SIGMA_%d_NONBSQ', nF, pref);
export_fig(hf, [PATH_FIGS fn] , '-jpg', '-m3');
close(hf)