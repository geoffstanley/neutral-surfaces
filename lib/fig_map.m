function hf = fig_map(varargin)
%FIG_MAP  Produce figure maps
%
%
% hf = fig_map(X, Y, data)
% creates a new figure, with handle hf, and using pcolor() to display data
% with horizontal coordinate X and vertical coordinate Y.
%
% hf = fig_map(X, Y, data, mask)
% also overlays a solid colour where mask is true.
%
% hf = fig_map(X, Y, data, mask, OPTS)
% overrides the default options according to OPTS.
%
% hf = fig_map(ax, ...)
% plots into the axes object ax.
%
%
% --- Input:
% ax [matlab.graphics.axis.Axes]: axes to plot into
% X [nx, 1]: horizontal coordinate
% Y [1, ny]: vertical coordinate
% data [nx, ny]: data to display in colour
% mask [nx, ny]: mask to overlay in a solid colour [logical]
% OPTS [struct]: user-defined options, containing some subset of the
% following (case-sensitive) fields:
%   CLIM [1, 2]: limits of the color bar (ax.CLim)
%   CM: colour map, passed into colormap(ax, cm)
%   MASKCOL[1, 3]: RGB colour for mask
%   NANCOL [1, 3]: RGB colour for NaN values in data
%   VIS [string]: 'on' to display the figure, 'off' to keep it invisible
%   XSH [1, 1]: Longitudal shift: puts xsh degrees east at the edges
%   LONLIM [1, 2]: Longitude limits
%   LATLIM [1, 2]: Latitude limits
%   LATTICK [vector]: Latitude tick marks
%   LONTICK [vector]: Longitude tick marks
%   GRIDCOL [] or [1, 3]: No grid, or RGB colour for grid lines across tick marks
%
%
% --- Output:
% hf [1, 1]: handle to the figure

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

ax = axescheck(varargin{:});
arg0 = double(~isempty(ax)); % 1 if no ax passed in, 0 otherwise
X    = varargin{arg0+1};
Y    = varargin{arg0+2};
data = varargin{arg0+3};
if length(varargin) >= arg0+4
    mask = varargin{arg0+4};
else
    mask = [];
end
[nx,ny] = size(data);
assert(isvector(X) && length(X) == nx, 'X must be a vector of length size(data,1) == %d', nx);
assert(isvector(Y) && length(Y) == ny, 'Y must be a vector of length size(data,2) == %d', ny);
X = X(:);
Y = Y(:);

% --- Set default options (where empty, values are determined from the data, later)
clim = []; % Color axis limits
cm = []; % Color map
maskcol = [0 0 0]; % Color where mask. Black
nancol = [1 1 1]* .75; % Color where nan data. Light grey
vis = 'on'; % Figure visibility
xsh = 20; % Degrees longitude to shift. This puts 270 deg E at the edges
lonlim = [0 360]; % Longitude limits
latlim = [-80 90]; % Latitude limits
lattick = -60 : 30 : 90; % Latitude tick marks
lontick = []; % Longitude tick marks
gridcol = []; % No grid. For dark grey, e.g. use [1 1 1]*.25;
deg = char(176);
if ismac
    fontsize = 12;
else
    fontsize = 10;
end

% --- Override default options with user options
if length(varargin) >= arg0+5  &&  isstruct(varargin{arg0+5})
    OPTS = varargin{arg0+5};
    if isfield(OPTS, 'CLIM') && ~isempty(OPTS.CLIM)
        clim = OPTS.CLIM;
    end
    if isfield(OPTS, 'CM') && ~isempty(OPTS.CM)
        cm = OPTS.CM;
    end
    if isfield(OPTS, 'LANDCOL') && ~isempty(OPTS.LANDCOL)
        maskcol = OPTS.LANDCOL;
    end
    if isfield(OPTS, 'NANCOL') && ~isempty(OPTS.NANCOL)
        nancol = OPTS.NANCOL;
    end
    if isfield(OPTS, 'GRIDCOL') && ~isempty(OPTS.GRIDCOL)
        gridcol = OPTS.GRIDCOL;
    end
    if isfield(OPTS, 'VIS') && ~isempty(OPTS.VIS)
        vis = OPTS.VIS;
    end
    if isfield(OPTS, 'XSH') && ~isempty(OPTS.XSH)
        xsh = OPTS.XSH;
    end
    if isfield(OPTS, 'LONTICK') && ~isempty(OPTS.LONTICK)
        lontick = OPTS.LONTICK;
    end
    if isfield(OPTS, 'LATTICK') && ~isempty(OPTS.LATTICK)
        lattick = OPTS.LATTICK;
    end
    if isfield(OPTS, 'LONLIM') && ~isempty(OPTS.LONLIM)
        lonlim = OPTS.LONLIM;
    end
    if isfield(OPTS, 'LATLIM') && ~isempty(OPTS.LATLIM)
        latlim = OPTS.LATLIM;
    end
    if isfield(OPTS, 'FONTSIZE') && ~isempty(OPTS.FONTSIZE)
        fontsize = OPTS.FONTSIZE;
    end
end

% --- Zonal shift
if xsh == 0
    ish = 0; % No shifting
else
    XG = (X + circshift(X, [+1 0])) / 2; % Assuming a regular grid
    XG(1) = 1.5*X(1) - 0.5*X(2); % Assuming a regular grid, extrapolate backwards
    ish = interp1(XG, 0:nx-1, mod(xsh, 360), 'nearest');
    xsh = XG(ish+1);
end
X = X + xsh;
data = circshift(data, [-ish 0]);
if isempty(lontick) && lonlim(1) == 0 && lonlim(2) == 360
    lontick = (0:60:300) + 60*ceil(xsh/60);
end
lonlim = lonlim + xsh;

% --- Set up figure
if isempty(ax)
    hf = figure('Position', [0 0 1000 600], 'Visible', vis);
    ax = axes;
else
    hf = [];
end
hold(ax,'on'); % Do this before calling pcolor, or else it will destroy pre-existing colorbars.
grid(ax, 'off');

% --- Use surf() to display the mask
if ~isempty(mask)
    mask = circshift(mask, [-ish 0]);
    mask = mask.';
    C = nan(nx*ny,3);
    C(mask,1) = maskcol(1);
    C(mask,2) = maskcol(2);
    C(mask,3) = maskcol(3);
    C = reshape(C, [ny nx 3]);
    C = [C, zeros(ny,1,3); zeros(1,nx+1,3)];
    hm = surf(ax, [X; Inf], [Y; Inf], zeros(ny+1,nx+1), C);
    hm.EdgeColor = 'none';
    view(ax,0,90);
end

% --- Plot grid lines inside mask aligned with the tick marks
if ~isempty(gridcol) && ~isempty(lontick)
    plot(ax, [1;1].*lontick, latlim, '-', 'LineWidth', 1, 'Color', gridcol);
end
if ~isempty(gridcol) && ~isempty(lattick)
    plot(ax, lonlim, [1;1]*lattick, '-', 'LineWidth', 1, 'Color', gridcol);
end


% --- Use pcolor() to display the data
X = repmat([X; nan], 1, ny+1);
Y = repmat([Y.', nan], nx+1, 1);
hp = pcolor(ax, X, Y, [data nan(nx,1); nan(1,ny+1)]);
hp.EdgeColor = 'none';
box(ax, 'on');
ax.Layer = 'top';  % bring the axes (e.g. grid, and bounding box) to top.

% --- Handle colors
ax.Color = nancol; % Set nan colour
if isempty(clim)
    clim = [min(data(:)), max(data(:))];
    if clim(1) == 0 && clim(2) == 0
        clim = [-1 1];
    end
elseif isscalar(clim)
    clim = [-1 1] * clim;
end
ax.CLim = clim;
if ~isempty(cm)
    colormap(ax,cm);
end

% --- Label axes
WEstr = 'W E';
SNstr = 'S N';
ax.XLim = lonlim;
ax.YLim = latlim;
if ~isempty(lontick)
    ax.XTick = lontick;
end
if ~isempty(lattick)
    ax.YTick = lattick;
end
xt = mod(ax.XTick, 360);
xt(xt > 180) = xt(xt > 180) - 360;
ax.XTickLabel = arrayfun(@(a) [num2str(abs(a)) deg WEstr(sign(a-180*(a==180))+2)], xt, 'UniformOutput', false);
ax.YTickLabel = arrayfun(@(a) [num2str(abs(a)) deg SNstr(sign(a)+2)], ax.YTick, 'UniformOutput', false);
ax.TickDir = 'out';
ax.FontSize = fontsize;
