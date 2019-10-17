% Produces figures on the illustrative omega surface using OCCA data:
%  (a) map of the pressure
%  (b) scatter plot of pressure vs. density anomaly
%  (c) simplified Reeb graph of the pressure
%  (d) associated regions of the Reeb graph of the pressure
% Also produces a table showing the location, pressure value, and incident
%   arcs for each node.
%
% The figure produced is similar to Fig 3 in Stanley (2019a) but has
% cosmetic differences.  That figure can be reproduced exactly using the
% original Topobaric-Surface package, at
% https://github.com/geoffstanley/Topobaric-Surface/
%
% Stanley, G.J., 2019a. Neutral surface topology. Ocean Modelling 138,
% 88–106. https://doi.org/10.1016/j.ocemod.2019.01.008

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

% Read path to OCCA data from PATH_OCCA.txt
file_id = fopen([PATH_PROJECT 'lib' V 'dat' V 'PATH_OCCA.txt']);
PATH_OCCA = textscan(file_id, '%s');
PATH_OCCA = PATH_OCCA{1}{1};
fclose(file_id);

%fileID = 1; % For standard output to the screen
fileID = fopen([PATH_LOCAL 'run ' datestr(now, 'yyyy mm dd hh MM ss') '.txt'], 'wt'); % For output to a file

db2Pa = 1e4; % dbar to Pa conversion
Pa2db = 1e-4; % Pa to dbar conversion

%% Load OCCA grid & data
[g, S, T, P, ETAN, ATMP] = load_OCCA(PATH_OCCA);
[ni,nj,nk] = size(T);

% Re-order data so water-columns are contiguous data:
S = permute(S, [3 1 2]); % [nk,ni,nj]. depth  by  long  by  lat
T = permute(T, [3 1 2]);
P = permute(P, [3 1 2]);
ETAN = permute(ETAN, [3 1 2]);

land = squeeze(isnan(S(1,:,:)));

X = g.XCvec;
Y = g.YCvec;
XG = g.XGvec;
YG = g.YGvec;
WRAP = [true false];

%% Set alias functions
% Ensure the equation of state is for the densjmd95 in-situ density:
copyfile([PATH_PROJECT 'lib' V 'eos' V 'eoscg_densjmd95.m'   ], [PATH_PROJECT 'lib' V 'alias' V 'eos.m'  ]);
copyfile([PATH_PROJECT 'lib' V 'eos' V 'eoscg_densjmd95_dp.m'], [PATH_PROJECT 'lib' V 'alias' V 'eos_x.m']);
clear eos eos_x % Make sure the copied file gets used

% Choose vertical interpolation method
copyfile([PATH_PROJECT 'fex' V 'columncalculus' V 'interp1qn2.m'], [PATH_PROJECT 'lib' V 'alias' V 'interp_firstdim_twovars.m']); % linear interpolation
%copyfile([PATH_PROJECT 'fex' V 'columncalculus' V 'pchipqn2.m'], [PATH_PROJECT 'lib' V 'alias' V 'interp_firstdim_twovars.m']);  % PCHIP interpolation
clear interp_firstdim_twovars % Make sure new file gets used

%% Selecting surfaces and depths
xref = 180; yref = 0;
[~,iref] = min(abs(g.XCvec - xref));
[~,jref] = min(abs(g.YCvec - yref));
I0 = sub2ind([ni,nj], iref, jref);

% Load the illustrative omega surface:
pref = 1000;
LOAD = load(sprintf('%somega.0406annclim.from_SIGMA%04d_through_(180,0,%04d).nonBoussinesq.mat', [PATH_OCCA V 'omega_v1.1gjs' V], pref, pref), 'pns_i');
p = squeeze(LOAD.pns_i).'; % transpose so that data is long by lat
clear LOAD

surfname = sprintf('OMEGA%04d_(%03d,%02d,%.0f)', pref, xref, yref, p(iref,jref));

%% Compute the Simplified Reeb Graph
OPTS = struct();
OPTS.WRAP = WRAP;

% Fill in everything except land connected to these locations:
FILL_LATLON = [ ...
    103    -69    ; ... % Antarctica
    130    -25    ; ... % Australia
    20      0     ; ... % Africa
    260     40    ; ... % North America
    320     75    ; ... % Greenland
    170    -45    ; ... % New Zealand South Island
    176    -39    ; ... % New Zealand North Island
    46     -20    ; ... % Madagascar
    140    -5     ; ... % Papua New Guinea
    113     0     ; ... % Borneo
    101     0     ; ... % Sumatra
    138     36    ; ... % Japan
    341     65    ; ... % Iceland
    303     48.5  ; ... % Newfoundland
    358     52    ; ... % Great Britain
    352     53    ; ... % Ireland
    290     66.38 ; ... % Baffin Island
    278     78    ; ... % Ellesmere Island
    250     70    ; ... % Victoria Island
    ];
% Convert above longitudes & latitudes to i & j grid indices:
lon_to_i = @(x) interp1(X, 1:ni, x, 'nearest');
lat_to_j = @(y) interp1(Y, 1:nj, y, 'nearest');
OPTS.FILL_IJ = FILL_LATLON;
for i = 1:size(OPTS.FILL_IJ,1)
    OPTS.FILL_IJ(i,:) = [lon_to_i(FILL_LATLON(i,1)), lat_to_j(FILL_LATLON(i,2))];
end

% Choose simplification level
OPTS.SIMPLIFY_ARC_REMAIN = 43;
OPTS.SIMPLIFY_WEIGHT_PERSIST = 0.5;


% Find the connected component containing the reference cast I0
[~,~,~,wet] = bfs_conncomp(isfinite(p), grid_adjacency([ni,nj], OPTS.WRAP), I0);
p(~wet) = nan;

% Calculate the Reeb graph
[p, RG] = calc_reeb_graph(p, OPTS);

% Interpolate water properties on the surface
lead1 = @(x) reshape(x, [1 size(x)]);
[s, t] = interp1qn2(lead1(p), P, S, T);

% Get in-situ density on the surface
r = eos(s, t, p);

%% 1 Big Figure: In situ density vs. Pressure, pressure map, regions, and Reeb graph.
fontsize = 10 + 2*(ismac);

ALPHABET = 'ABCDEFGHJKMNPQRSTUVWXYZ'; % Skip I L O
alphabet = 'abcdefghjkmnpqrstuvwxyz'; % Skip i l o
node_label = repmat('0', 1, RG.nNodes);
L1 = 0; l1 = 0;
for n = 1 : RG.nNodes
    if RG.node_type(n) == 2
        L1 = L1 + 1;
        node_label(n) = ALPHABET(L1);
    else
        l1 = l1 + 1;
        node_label(n) = alphabet(l1);
    end
end

I = 1:RG.nArcs;


cm = distinguishable_colors(RG.nArcs, {[1 1 1]});
cm = cm(I,:);

I_cm = I;

% jig enables moving the node labels for better visuals
% Columns 1,2 for Associated Region map [lon shift, lat shift]
% Columns 3,4 for Reeb Graph [pressure shift, lon shift]
% Column 5 for the P vs R scatter plot [density shift]
jig = zeros(5, RG.nNodes);


% nodelonjig enables actually shifting the Reeb Graph's nodes for better visuals
nodelonjig = zeros(1, RG.nNodes);


% Just the casts indexed by RG.arc_segment. Gets rid of disconnected
% components on the surface, and ground.
p_ = p(wet);

inds_ocean = find(wet); % inverse to ij2v. a map from 2D space onto the 1D data vectors (p -> p_)
v2ij = @(v) ind2sub([ni nj], inds_ocean(v));
%v2i = @(v) mod1(inds_ocean(v), ni);  % function from a 1D data index to the i (longitude) in 2D space
%v2j = @(v) ceil(inds_ocean(v) / ni); % function from a 1D data index to the j (latitude) in 2D space

arc_ = zeros(length(p_),1);
for e = 1:RG.nArcs
    arc_(RG.arc_segment{e}) = e;
end
if any(arc_ == 0)
    fprintf('%d vertices lost.\n', sum(arc_ == 0));
end
arc = nan(ni,nj);
arc(wet) = arc_;

xsh = 100; % put 100 deg E on the edges
ish = interp1(XG, 0:ni-1, mod(xsh,360), 'nearest');
xsh = XG(ish+1);

OPTS_FIGS.XSH = xsh;
OPTS_FIGS.LATLIM = [-77 77];
OPTS_FIGS.LONLIM = [0 360]; % prior to shifting
OPTS_FIGS.LATTICK = (-60:20:60);


rgn = zeros(5, RG.nNodes); % Each column is a node: [x; y; p; x; r].
for n = 1 : RG.nNodes
    % If a min or max, recompute the min or max in that arc segment.
    % node_fn can't be used for it because the leaf pruning process could
    % have changed it to some other leaf.
    if RG.node_type(n) == 1 % a min
        e = RG.node_next{n};
        seg = vertcat(RG.arc_segment{e}); % vertcat because possible to have min node but 2 up arcs
        
        if node_label(n) == 'h'
            % Manually move this node! It has a data point in an annoying
            % spot, very close to other nodes. Instead of the minimum p in
            % the region, take the p between the 1% and the 99%
            [p_sort,I] = sort(p_(seg));
            J = interp1(linspace(0,100,length(p_sort)), 1:length(p_sort), 1, 'nearest');
            rgn(3,n) = p_sort(J);
            I = I(J);
            clear p_sort J
        else
            [rgn(3,n),I] = min(p_(seg));
        end
        [i,j] = v2ij(seg(I));
    elseif RG.node_type(n) == 3 % a max
        e = RG.node_prev{n};
        seg = vertcat(RG.arc_segment{e}); % vertcat because possible to have max node but 2 down arcs
        [rgn(3,n),I] = max(p_(seg));
        [i, j] = v2ij(seg(I));
    else % saddle. Can use node_fn here, because leaf pruning never got these.
        [i,j] = v2ij(RG.node_v(n));
        rgn(3,n) = RG.node_fn(n);
    end
    
    rgn(1,n) = X(i);
    rgn(2,n) = Y(j);
    rgn(4,n) = X(i) + nodelonjig(n);
    rgn(5,n) = r(i,j);
    
end

% Shift into desired longitude range which is: [0, 360]+xsh
rgn([1 4],:) = mod(rgn([1 4],:) - xsh, 360) + xsh;

% Define the in-situ density reference profile
sref = s(iref,jref);
tref = t(iref,jref);
rreffn = @(p) eos(sref,tref,p);

% Find dupliate arcs. They'll be made thick.
[~, iA, iC] = unique([RG.arc_from, RG.arc_to], 'rows', 'stable');
iAiC = iA(iC);
nondup = (iAiC - (1:length(iC))') == 0;
dup = find(~nondup);

WEstr = 'W E';
SNstr = 'S N';

PC = [500 1000 1250 1350 ]; % for OMEGA 1000
COL = repmat([1; 0; 1; 0], 1, 3);
LS = {'-', '-', '-', '-'};
LW = [2 1 1 2];
LWb = [3 1 2 2];

plim = [0 max(p_)+15];
pbreak = 750;

% Set up figure and axes
hf1 = figure('Position', [0 0 800 900]);

ml = .07; % margin left
mr = .07; % margin right
mb = .05; % margin bottom
mt = .06; % margin top
sv = .05; % spacing vertical
sh = .01; % spacing horizontal
h = (1-mb-mt-sv); % height of each panel, were things evenly spaced
w = (1-ml-mr-sh)/2; % width of each panel
hb = h * (2/3) * .85; % height for bottom part of split panel
hm = h * (2/3) * (1-.85); % height for middle panel (top part of split panel)
ht = h * (1/3); % height for top panel

% Pressure map: Top left
axP   = axes('Position', [ml     , mb+hb+hm+sv , w, ht]); % Bottom left

% Regions: Top right
axMAP = axes('Position', [ml+w+sh, mb+hb+hm+sv , w, ht]); % Bottom right

% Density vs pressure: Bottom left
axRP(2) = axes('Position', [ml, mb   , w, hb]); % lower
axRP(1) = axes('Position', [ml, mb+hb, w, hm]); % upper

% Reeb graph: Bottom right
axRG(2) = axes('Position', [ml+w+sh, mb   , w, hb]); % lower
axRG(1) = axes('Position', [ml+w+sh, mb+hb, w, hm]); % upper

% little axes for the serrated line
hw = .0025;
axWIGGLE = axes('Position', [ml, mb+hb-hw/2, w*2+sh, hw]);
xx = 0:150;
yy = mod(xx, 2);
plot(axWIGGLE, xx, yy, '-k');
axWIGGLE.XAxis.Visible = 'off';
axWIGGLE.YAxis.Visible = 'off';

hcb = colorbar(axP, 'Location', 'NorthOutside');
hcb.Position = [ml+0.005, mb+hb+hm+sv+ht-.018, .2, .025];


% PANEL: Pressure on the surface, with same mask as the arc map
% Prep colormap
nc = 512;
b0 = 0; b1 = 750; b2 = 1150; b3 = plim(2); % for OMEGA 1000
axP.CLim = [b0 b3];
db = nc/(b3 - b0);
cm2 = flipud(jet(round((b2-b1)*db)));
cm2 = cm2(10:end-1,:);
cm3 = [cm2(end,:); .6 .35 1; 0 0 .25];
cm3 = interp1([0 .5 1], cm3, linspace(0,1,(b3-b2)*db), 'pchip');
cm1 = [1 1 1; cm2(1,:)];
cm1 = interp1([0 1], cm1, linspace(0,1,(b1-b0)*db).^.75, 'pchip');
cmP = vertcat(cm1, cm2, cm3);

OPTS_FIGS.CLIM = [b0 b3];
OPTS_FIGS.CM = cmP;
fig_map(axP, X, Y, p, land, OPTS_FIGS);

% Colorbar for the pressure map
haxcb = axes('position', hcb.Position, 'xlim', hcb.Limits, 'color', 'none', 'visible','off');
hold(haxcb, 'on');
for i = 1:length(PC)
    contour(axP, X + xsh, Y, circshift(p, [-ish 0]).', PC(i)*[1 1] , LS{i}, 'linewidth', LW(i), 'Color', COL(i,:));
    plot(haxcb, [1 1]*PC(i), haxcb.YLim, LS{i}, 'linewidth', LW(i), 'Color', COL(i,:))
end
hcb.Ticks = 0 : 500 : plim(2);

% PANEL: Regions associated to each arc
OPTS_FIGS.CLIM = [0 RG.nArcs];
OPTS_FIGS.CM = cm;
fig_map(axMAP, X, Y, arc - 0.5, land, OPTS_FIGS);
axMAP.YAxisLocation = 'right';


% PANEL: Scatter plot of density vs pressure
yy = p_;
xx = (r(wet) - rreffn(yy)) * 1e3;
ll = arc_;
I = randperm(length(xx));
xx = xx(I); yy = yy(I); ll = ll(I);
for i = 1:2
    scatter(axRP(i), xx, yy, 12, ll, '.');
end

axRP(1).XLim = [min(xx)-10, max(xx)+10];
for i = 1:2
    hold(axRP(i), 'on');
    axRP(i).CLim = [1 RG.nArcs];
    axRP(i).YDir = 'reverse';
    axRP(i).XGrid = 'on';
    axRP(i).Color = [1 1 1]*.85;
    colormap(axRP(i),cm);
    
    hold(axRG(i), 'on');
    axRG(i).XGrid = 'on';
    axRG(i).XLim = [0 360] + xsh;
    axRG(i).YDir = 'reverse';
    axRG(i).XTick = axMAP.XTick;
    axRG(i).YTick = axRP(i).YTick;
    axRG(i).Color = axRP(i).Color;
    axRG(i).XTickLabel = axMAP.XTickLabel;
end


axRP(1).YTick = 0 : 250 : pbreak;
axRP(2).YTick = (pbreak+50) : 50 : plim(2) ;
axRP(2).XLim = [min(xx(yy > pbreak)), max(xx(yy > pbreak))];
axRP(1).YLim = [plim(1) pbreak];
axRP(2).YLim = [pbreak plim(2)];
axRP(1).XAxisLocation = 'top';

for i = 1:2
    axRG(i).YLim = axRP(i).YLim;
    axRG(i).YTick = axRP(i).YTick;
    axRG(i).YAxisLocation = 'right';
end
axRG(1).XAxisLocation = 'top';
axRG(1).XTickLabel = {};

% A way to get "box on" but only for 3 out of 4 sides...
for ax = [axRP axRG]
    axStealth = axes();
    axStealth.Position = ax.Position;
    axStealth.XLim = ax.XLim;
    axStealth.YLim = ax.YLim;
    axStealth.XTick = ax.XTick;
    axStealth.YTick = ax.YTick;
    axStealth.XTickLabel = {};
    axStealth.YTickLabel = {};
    axStealth.XAxisLocation = ax.XAxisLocation;
    if strcmp(ax.YAxisLocation, 'left')
        axStealth.YAxisLocation = 'right';
    else
        axStealth.YAxisLocation = 'left';
    end
    axStealth.YDir = ax.YDir;
    axStealth.Color = 'none';
end

% Horizonal lines indicating which levels are contoured in the map
for l = 1 : length(PC)
    i = (PC(l) > pbreak) + 1;
    plot(axRP(i), axRP(i).XLim, PC(l)*[1 1], LS{l}, 'linewidth', LWb(l), 'Color', COL(l,:));
    plot(axRG(i), axRG(i).XLim, PC(l)*[1 1], LS{l}, 'linewidth', LWb(l), 'Color', COL(l,:));
end

% PANEL: The Reeb Graph
for e = 1:RG.nArcs
    col = cm(e,:);
    lw = 3 + 2 * ismember(e, dup) ;
    n = [RG.arc_from(e), RG.arc_to(e)];
    pp = rgn(3,n); % row vector
    xx = rgn(4,n); % row vector
    [pp,idx] = sort(pp.'); % col vector
    xx = xx(idx).';        % col vector
    per = abs(diff(xx)) > 180; % wrap around (0/360)-xsh;
    
    if all(pp < pbreak)
        p3 = [pp; nan];
        x3 = [xx; nan];
    elseif all(pp >= pbreak)
        p3 = [nan; pp];
        x3 = [nan; xx];
    else
        xbreak = xx(1);
        p3 = [pp(1); pbreak; pp(2)];
        x3 = [xx(1); xbreak; xx(2)];
    end
    
    for i = 1:2
        xx = x3(i:i+1);
        pp = p3(i:i+1);
        if per && abs(diff(xx)) > 180 % longitude wrapping in this panel
            [xx,idx] = sort(xx);
            pp = pp(idx);
            pmid = interp1qn(xsh, [xx(2)-360; xx(1)], pp([2 1]));
            plot(axRG(i), [xx(1)     xsh], [pp(1) pmid], '-', 'Color', col, 'LineWidth', lw);
            plot(axRG(i), [xx(2) 360+xsh], [pp(2) pmid], '-', 'Color', col, 'LineWidth', lw);
        else % no wrap
            plot(axRG(i), xx, pp, '-', 'Color', col, 'LineWidth', lw);
        end
    end
end

% Set fontsize for all axes
for ax = [axRP axRG axP axMAP hcb]
    ax.FontSize = fontsize;
end

% --- Now label nodes:
label_nodes = 1 : RG.nNodes;
RPlabelnodes = 'abcdefhkpvwxy';
L1 = 0; l1 = 0;
for n = label_nodes
    
    node_x = rgn(1,n);
    node_y = rgn(2,n);
    
    if RG.node_type(n) == 2
        L1 = L1 + 1;
        letter = ALPHABET(L1);
        plot(axP, node_x, node_y, 'xk');
    else
        l1 = l1 + 1;
        letter = alphabet(l1);
        if RG.node_type(n) == 1
            plot(axP, node_x, node_y, '^k'); % min
        else
            plot(axP, node_x, node_y, 'vk'); % max
        end
    end
    
    node_x = node_x + jig(1,n);
    node_y = node_y + jig(2,n);
    
    textborder(axMAP, node_x, node_y, letter, 'w', 'k', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontWeight', 'bold', 'FontSize', fontsize);
    
    
    node_p = rgn(3,n);
    node_x = rgn(4,n) + jig(4,n);
    node_ra = (rgn(5,n) - rreffn(node_p)) * 1e3 + jig(5,n);
    i = (node_p > pbreak)+1;
    
    if ismember(letter, RPlabelnodes)
        text(axRP(i),node_ra,node_p,letter, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'Color', [1 1 1]*0, 'FontWeight', 'normal', 'FontSize', fontsize)
    end
    
    node_p = node_p + jig(3,n);
    i = (node_p > pbreak)+1;
    
    text(axRG(i),node_x,node_p,letter, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'Color', [1 1 1]*0, 'FontWeight', 'normal', 'FontSize', fontsize)
    
end

%% --- Now add labels to the panels:

axFIG = axes('Position', [0 0 1 1], 'visible','off'); %'color', 'none',
hold(axFIG, 'on');

text(axFIG, .01, mb+hb+hm+sv+ht+.01, '(a) $p$', 'FontSize', fontsize+4, 'Interpreter', 'latex');
text(axFIG, .043, mb+hb+hm+sv+ht, '~', 'Units', 'norm', 'FontSize', fontsize+2); % under p
text(axFIG, .24, mb+hb+hm+sv+ht+.0175, '[dbar]', 'FontSize', fontsize+2, 'Interpreter', 'latex');

text(axFIG, .5, mb+hb+hm+sv+ht+.01, '(d) Associated regions of each arc', 'FontSize', fontsize+4, 'Interpreter', 'latex');

text(axFIG, .01, mb-.033, '(b)', 'FontSize', fontsize+4, 'Interpreter', 'latex');
text(axFIG, .2, mb-.033, '$\rho - R(S_0, \theta_0, p)$  [g m$^{-3}$]', ...
    'FontSize', fontsize+4, 'Interpreter', 'latex');
text(axFIG, .2    , mb-.033-.01, '~', 'Units', 'norm', 'FontSize', fontsize+2); % under \rho
text(axFIG, .2+.11, mb-.033-.01, '~', 'Units', 'norm', 'FontSize', fontsize+2); % under p
text(axFIG, .02, .4, '$p$ [dbar]', 'Rotation', 90, ...
    'Units', 'norm', 'FontSize', fontsize+4, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
text(axFIG, .03, .4, '~', 'Rotation', 90, 'Units', 'norm', 'FontSize', fontsize+2); % under p

text(axFIG, .5, mb-.033, '(c) Simplified Reeb Graph of $p$', 'FontSize', fontsize+4, 'Interpreter', 'latex');
text(axFIG, .762, mb-.033-.01, '~', 'Units', 'norm', 'FontSize', fontsize+2); % under p
text(axFIG, .975, .4, '$p$ [dbar]', 'Rotation', 90, ...
    'Units', 'norm', 'FontSize', fontsize+4, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
text(axFIG, .985, .4, '~', 'Rotation', 90, 'Units', 'norm', 'FontSize', fontsize+2); % under p



%% Save figure
fn = sprintf('illustrative_%s', surfname);
%saveas(hf1, [figsdir fn '.fig']);
export_fig(hf1, [PATH_FIGS fn], '-jpg', '-m2.5');

%% Table of RG, formatted for Latex
westr = 'w e';
snstr = 's n';
L1 = 0; l1 = 0;
lines = cell(RG.nNodes,1);
for n = 1 : RG.nNodes
    node_x = rgn(1,n);
    node_y = rgn(2,n);
    node_p = rgn(3,n);
    node_x = mod(node_x+180,360)-180; % Make between [-180, 180].
    neigh_up = RG.arc_to(RG.node_next{n});
    neigh_dn = RG.arc_from(RG.node_prev{n});
    line = sprintf('%s & %.0f\\deg%s & %.0f\\deg%s & %.0f & %s%s', ...
        node_label(n), abs(node_x), westr(sign(node_x)+2), abs(node_y), snstr(sign(node_y)+2), node_p, node_label(neigh_dn), node_label(neigh_up));
    if RG.node_type(n) == 2
        L1 = L1 + 1;
        lines{L1,2} = line;
    else
        l1 = l1 + 1;
        lines{l1,1} = line;
    end
end
lines = lines(1:max(l1,L1), :);
for n = 1 : size(lines,1)
    fprintf('%s & %s \\\\\n', lines{n,1}, lines{n,2});
end
