% Make figure showing pitch of neutral trajectory around an island,
% plus that for the same shape neutral trajectory but shifted in longitude
% Repeat for several islands, for several starting depths.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


%% --- BEGIN SETUP --------------------------------------------------------
warning('off', 'MATLAB:nargchk:deprecated')
set(0, 'defaultfigurecolor', [1 1 1]); % white figure background
V = filesep(); % /  or  \  depending on OS.

% Add Neutral Surfaces to MATLAB's path
PATH_NS = '~/work/projects-gfd/neutral-surfaces/'; % adjust as needed.
run([PATH_NS 'ns_add_to_path.m']);

% Make a directory for figures and other output
PATH_OUT = '~/work/dphil/projects/topobaric/output/'; % adjust as needed.
PATH_OUT = [PATH_OUT datestr(now, 'dd-mm-yyyy') '/'];
if ~exist(PATH_OUT, 'dir'); mkdir(PATH_OUT); end

% Folder containing the functions eos.m, eos_p.m, and eos_s_t.m
PATH_EOS = '~/work/MATLAB/eos/eos/';
addpath(PATH_EOS);  % Doing this last to ensure at top of MATLAB's path, above eos fcns in other dirs

% Read path to OCCA data from PATH_OCCA.txt
file_id = fopen([PATH_NS 'lib' V 'dat' V 'PATH_ECCO2.txt']);
PATH_ECCO2 = textscan(file_id, '%s');
PATH_ECCO2 = PATH_ECCO2{1}{1};
fclose(file_id);

%fileID = 1; % For standard output to the screen
fileID = fopen([PATH_OUT 'run ' datestr(now, 'yyyy mm dd hh MM ss') '.txt'], 'wt'); % For output to a file

db2Pa = 1e4; % dbar to Pa conversion
Pa2db = 1e-4; % Pa to dbar conversion

%% Load ECCO
TIMESTEP = '20021223';
g = load_ECCO2(PATH_ECCO2, 'grid');
ni = g.nx;
nj = g.ny;
nk = g.nz;
g.XC = repmat(g.XCvec, 1, nj);
g.YC = repmat(g.YCvec, ni, 1);
Z = -g.RC; % We are going to be Z > 0 people!
WRAP = g.WRAP;
[S, T, ATMP, ETAN, SAP] = load_ECCO2(PATH_ECCO2, 'casts', TIMESTEP);

% Re-order data so water-columns are contiguous data:
S = permute(S, [3 1 2]); % [nz,nx,ny]. depth  by  long  by  lat
T = permute(T, [3 1 2]);

tol = 1e-6; % Tolerance on depth [m] when solving neutral trajectories
BotK = squeeze(sum(isfinite(S),1)); % vertical index for deepest ocean grid point
BotK(BotK == 0) = 1; % Just for finite indexing

y_to_j = @(y) floor((y - g.YGvec(1)) * g.resy) + 1;
x_to_i = @(x) mod(floor((x - g.XGvec(1)) * g.resx), ni) + 1;

grav = g.grav;
rho_c = g.rho_c;

%% Set alias functions.  << ENSURE THIS GETS DONE >>
% Choose the Boussinesq densjmd95 and set grav and rho_c in eos.m and eos_p.m
if false % only need to do this once!  Doing every time makes codegen run every time
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq.m'   ] , [PATH_EOS 'eos.m'  ], grav, rho_c); %#ok<UNRCH>
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq_dz.m'] , [PATH_EOS 'eos_p.m'], grav, rho_c);
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq_s_t.m'], [PATH_EOS 'eos_s_t.m'], grav, rho_c);
end
clear eos eos_p eos_s_t % Make sure the copied files get used


%% Build interpolants for S and T in terms of Z
interpfn = @ppc_linterp; % linear interpolation

Sppc = interpfn(Z, S);
Tppc = interpfn(Z, T);

%% Run codegen on ntp_bottle_to_cast
ntp_bottle_to_cast_codegen(nk);

%% Select Islands and depths to study
island_data = { ...
    'madagascar', 46, -20; ...
    'sandwich', 323.1, -54.12; ...
    'kerguelen', 69, -49; ...
    'new zealand', 169, -45; ...
    'australia and new guinea', 135, -25; ...
    'hawaii', 204.6, 19.62; ...
    'fiji', 178, -17.75; ...
    'iceland', 340, 65; ...
    'antarctica', 180, -89; ...
    'galapagos', 268.9, -.625; ...
    'flemish cap', 315.4, 47.12; ...
    'japan', 140.4, 37; ...
    'svalbard', 15, 79; ...
    'corsica', 9, 42; ...
    };
island_x = containers.Map(island_data(:,1), island_data(:,2));
island_y = containers.Map(island_data(:,1), island_data(:,3));

island_names = {'Flemish Cap', 'Australia and New Guinea', 'Fiji', 'Madagascar', 'Kerguelen'};
island_names = repelem(island_names, 1, 2); % do each twice
z0s = [400 500, 400 700, 500 650, 500 1000, 750 2000]; % starting depths

%% Run neutral trajectories around each island
H = length(island_names);
island_pitch = nan(1,H);
island_pos = nan(2,H);
other_pitch = nan(ni,H);
neigh = grid_adjacency([ni,nj], 4, WRAP);
LMAX = 2000;
for h = 1 : H
    
    z0 = z0s(h);
    k0 = find(Z > z0, 1, 'first');
    
    % Don't use ground = g.Depth < z0, which can differ wildly from
    % Z(BotK)... because ECCO2 lat-lon data is actually interpolations of
    % cube-sphere data
    ground = isnan(squeeze(S(k0,:,:)));
    
    % Determine shape of islands/holes
    [~,~,~,~,LM] = bfs_conncomp(ground, neigh, []);
    
    % Determine shape of the chosen island/hole
    island_name = lower(island_names{h});
    i = x_to_i(island_x(island_name));
    j = y_to_j(island_y(island_name));
    island = (LM == LM(i,j));
    
    % Find most easterly part aligned with central latitude of chosen island/hole
    j0 = y_to_j(mean(g.YC(island)));
    i0 = find(island(:,j0),1,'last');
    i0 = mod(i0,ni) + 1; % Take 1 step east
    
    % Prepare to track a neutral trajectory around the chosen island/hole
    z = z0;
    ij = nan(LMAX,2);
    ij(1,:) = [i0, j0];
    
    i1 = ij(1,1); j1 = ij(1,2);
    k = BotK(i1,j1);
    [s,t] = ppc_val2(Z(1:k), Sppc(:,1:k-1,i1,j1), Tppc(:,1:k-1,i1,j1), z);
    
    dir = [0 +1; +1 0; 0 -1; -1 0]; % N E S W
    d = 3; % first step will try to move west
    l = 2;
    
    % Begin neutral trajectory, one step at a time
    while true
        
        % Test each of four directions until a step can be made
        for turn = 1:4
            d = mod(d,4) + 1;
            i2 = i1 + dir(d,1);
            j2 = j1 + dir(d,2);
            k = BotK(i2,j2);
            if k < 2
                % Only 1 grid point in this water column. Don't even bother
                continue
            end
            
            % Test a neutral trajectory from current location to neighbouring location
            [z2, s2, t2, success] = ntp_bottle_to_cast_mex(Sppc(:,:,i2,j2), Tppc(:,:,i2,j2), Z, k, s, t, z, tol);
            
            if success
                % A successful step! We did not ground. This direction worked -- we grounded.
                d = d-2;
                z = z2;
                s = s2;
                t = t2;
                break
            end
        end
        if isnan(z2)
            error('Stuck at (%d,%d) at step %d. Current depth: %.2f\n', ...
                i2, j2, l, z);
        end
        
        ij(l,:) = [i2, j2];
        i1 = i2;
        j1 = j2;
        
        if i2 == i0 && j2 == j0 || l > LMAX
            % What if it has spiralled in or out? More advanced code needed for that.
            if l > LMAX
                z = nan;
            end
            break
        end
        l = l + 1;
    end
    STEPS = l;
    ij = ij(1:l,:);
    
    % Record the island's pitch and starting position:
    island_pitch(h) = -(z - z0);  % (-) undoes that z > 0
    island_pos(1,h) = g.XCvec(ij(1,1));
    island_pos(2,h) = g.YCvec(ij(1,2));
    
    fprintf('%s: Change = %.4fcm. s = %03d. z0=%d.\n', ...
        island_name, island_pitch(h) * 100, l, z0)
    
    % --- Now repeat, shifted in longitude, with trajectory replicating the
    % shape of the one around the island/hole
    j1 = ij(1,2);
    %for i1 = find(botK(:,j1) >= 2).'
    for i1 = 1 : ni
        k = BotK(i1,j1);
        % Only check helix if starting column has >= 2 grid points, and
        % it's not the helix around the island already checked.
        if k >= 2 && i1 ~= ij(1,1)
            
            z = z0;
            [s,t] = ppc_val2(Z(1:k), Sppc(:,1:k-1,i1,j1), Tppc(:,1:k-1,i1,j1), z);
            
            for l = 2 : STEPS
                
                i2 = mod(ij(l,1) - i0 + i1 - 1, ni) + 1;
                j2 = ij(l,2);
                k = BotK(i2,j2);
                
                [z, s, t, success] = ntp_bottle_to_cast_mex(Sppc(:,:,i2,j2), Tppc(:,:,i2,j2), Z, k, s, t, z, tol);
                
                if ~success
                    break
                end
            end
            other_pitch(i1,h) = -(z - z0); % (-) undoes that z > 0
        end
    end
    
end

%% Make figure
deg = char(176);
OPTS_FIGS.FONTSIZE = 10 + 2*(ismac);
OPTS_FIGS.LONLAB = {['0' deg], ['60' deg 'E'], ['120' deg 'E'], ['180' deg], ['120' deg 'W'], ['60' deg 'W'], ['0' deg '']};

width = 360;
subaxopts = {'Margin', .05, 'MarginTop', .02, 'MarginLeft', 0.1, 'MarginRight', .13 'SpacingVertical', .03};

hf = figure('Position', [0 0 width round(width*sqrt(2)*1.6)]);
for h = 1:H
    island_name = lower(island_names{h});
    
    ax = subaxis(H,1,h,subaxopts{:});
    hold(ax, 'on'); box(ax, 'off');
    plot(ax, [0 360], [0 0], ':', 'Color', [1 1 1]*0); % Zero line
    
    scatter(ax, g.XCvec, other_pitch(:,h), 12, '.', 'k');
    plot(ax, island_pos(1,h), island_pitch(h), 'ok');
    
    horiz_lines = [min(other_pitch(:,h)), island_pitch(h), max(other_pitch(:,h))];
    ax.YTick = [min(horiz_lines), max(horiz_lines)];
    ax.YTickLabel = arrayfun(@(c) sprintf('%.2f', c), ax.YTick, 'UniformOutput', false);
    
    ax.YLim = ax.YTick([1 end]);
    ax.XLim = [0 360];
    ax.XTick = 0 : 60 : 360;
    if h < H
        ax.XTickLabel = {};
    else
        
        ax.XTickLabel = OPTS_FIGS.LONLAB;
        ax.XTickLabelRotation = 33;
        
        for i = 0 : 6
            x = i * 60;
            text(ax, x, ax.YTick(1), OPTS_FIGS.LONLAB{i+1}, 'FontSize', OPTS_FIGS.FONTSIZE, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        end
        
    end
    ax.XColor = 'none';
    ax.YColor = 'none';
    
    myspace = {'', ' '};
    for i = 1 : 3
        y = horiz_lines(i);
        plot(ax, [0 360], [1 1]*y, '-k');
        if h == 6 % Fudge some text placement for Fiji at -650m
            if i == 2 % the island's pitch
                dy = .03;
            elseif i == 3 % the maximum pitch of island duplicates
                dy = -.03;
            end
        else
            dy = 0;
        end
        text(ax, 365, y+dy, sprintf('%s%.2f', myspace{(y>0)+1}, y), 'FontSize', OPTS_FIGS.FONTSIZE, 'HorizontalAlignment', 'left');
    end
    
    plot(ax, (0:60:360) .* [1;1], ax.YLim, '-', 'Color', [1 1 1] * .8);
    
    txt = sprintf('(%s)  %s at %dm', char('a' + h - 1), island_names{h}, -z0s(h)); % (-)z0s(h) undoes that z0 > 0
    text(.02,1, txt, 'FontSize', OPTS_FIGS.FONTSIZE, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Units', 'n');
    
    assert(sum(other_pitch(:,1) == island_pitch(1)) == 0, 'Found other longitudes that match island''s pitch')
    
end

% ghost axis
axFIG = axes('Position', [0 0 1 1], 'visible','off');
text(axFIG, 0.04, .5, 'Neutral Trajectory Pitch [m]', 'FontSize', OPTS_FIGS.FONTSIZE, 'rotation', 90, 'units', 'n', 'HorizontalAlignment', 'center')

%% Save figure
fn = 'Island_Pitch';
export_fig(hf, [PATH_OUT fn], '-pdf');
close(hf);
return

%% Examine contours of pressure on approx neutral surface around each island
% This requires some human exploration!
H = length(island_names); %#ok<UNRCH>
for h = 1:H
    %%
    % -- vv copied from above
    z0 = z0s(h);
    k0 = find(Z > z0, 1, 'first');
    ground = isnan(squeeze(S(k0,:,:)));
    [~,~,~,~,LM] = bfs_conncomp(ground, neigh, []);
    island_name = lower(island_names{h});
    x = island_x(island_name);
    y = island_y(island_name);
    i = x_to_i(x);
    j = y_to_j(y);
    island = (LM == LM(i,j));
    j0 = y_to_j(mean(g.YC(island)));
    i0 = find(island(:,j0),1,'last');
    i0 = mod(i0,nx) + 1; % Take 1 step east
    % -- ^^ copied from above
    
    % Use locally referenced potential density as approx neutral surface
    SIGMA = sort(eos(S, T, z0), 1); % Make monotonic
    sigval = interp1(Z, SIGMA(:,i0,j0), z0);
    z = squeeze(interp1qn(sigval, SIGMA, Z));
    clear SIGMA
    
    % Plot depth of surface, and some contours. You'll want to zoom in yourself.
    figure;
    ax = axes;
    imagesc(ax, g.XCvec, g.YCvec, z'); ax.YDir = 'n';
    colorbar; hold on
    plot(ax, x, y, 'xr'); % Pinpoint the island
    contour(ax, g.XCvec, g.YCvec, z', z0 - 200 : 50 : z0 + 200, 'k');
    title(sprintf('(%s)  %s at %dm', char('a' + h - 1), island_names{h}, -z0s(h))); % (-)z0s(h) undoes that z0 > 0
    
end

