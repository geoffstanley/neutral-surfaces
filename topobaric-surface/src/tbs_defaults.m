function OPTS = tbs_defaults(ni, nj)
%TBS_DEFUALTS  Default options for topobaric_surface.
%
%
% OPTS = tbs_defaults(ni,nj)
% returns a struct OPTS containing default options for use in
% topobaric_surface, appropriate for a grid with ni longitude points and nj
% latitude points.

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

% For topobaric surfaces, a multivalued function is fit to the surface
% according to the Reeb graph. Set REEB to true.
% For "orthobaric" surfaces, a single-valued function is fit to the surface
% using a spline. Set REEB to false and set the SPLINE parameters below.
OPTS.REEB = true;

% Spine parameters, only used for "orthobaric" surfaces
OPTS.SPLINE_BREAKS = [0 200 1500 1800 6000]; % [dbar] or [m]
OPTS.SPLINE_ORDER = 4; % Cubic

% Do not remove the Mixed Layer
OPTS.MLX = [];

% Pixel indices for reference water column. Trying to select (180E,0N)
OPTS.REF_IJ = round([ni/2, nj/2]);

% Reference pressure or depth.  If a scalar is provided, the pressure or
% depth on the topobaric surface at the reference water column will be this
% value (with precision of OPTS.TOL).
OPTS.REF_X = []; % [dbar] or [m, positive]

% Reference practical / Absolute salinity and potential / Conservative
% temperature. If provided as scalars, these determine the way delta
% (specific volume anomaly or in-situ density anomaly) is caluclated.
% Otherwise, these values are taken from the S and T values where the
% reference water column intersects the reference pressure.
OPTS.REF_S = [];
OPTS.REF_T = [];

% Simplical Decomposition: Split each rectangle into 2 triangles
% ('diagonal'). Another option is 'cross', which does 4 way averaging to
% add an extra point in the middle of each rectangle of valid data. This
% option is disabled however, as it can produce arcs with 0 or 1 data
% points in a segment, which makes fitting affine linear functions
% under-determined.
OPTS.DECOMP = 'diagonal';

% Pre-processing (before Reeb graph computation) of the simplical mesh:
OPTS.FILL_IJ = []; % No filling
OPTS.FILL_PIX = 0; % No filling

% Use linear interpolation in the vertical dimension.
OPTS.INTERPFN = @ppc_linterp;

% Post-processing of Reeb Graph -- graph simplification parameters:
OPTS.SIMPLIFY_WEIGHT_PERSIST = 0.5; % Equal weighting between area and persistence
OPTS.SIMPLIFY_ARC_REMAIN = Inf;     % No leaf pruning simplification

% Error tolerance when root-finding to update surface, in the same units as
% X [dbar] or [m].
OPTS.X_TOL = 1e-4;

% Solutions to the root-finding problem in each water column are sought in
% the domain of the local branch of the multivalued function expanded
% outwards by this amount, in the same units as X [dbar] or [m].
OPTS.X_EXPN = 500;

% Pre-computed interpolation functions.  None given here.
OPTS.SppX = [];
OPTS.TppX = [];

% Conditions to terminate iterations:
OPTS.ITER_MAX = 6;  % maximum number of iterations

% quit when the L2 change of pressure on surface exceeds this value (set to
% 0 to deactivate), in the same units as X [dbar] or [m].
OPTS.ITER_L2_CHANGE = 1e-3;

% Iteration number at which to begin wetting (use inf for no wetting)
OPTS.ITER_START_WETTING = 1; % Start wetting immediately

% When true, adjust the empirical function to ensure exact geostrophic
% stream function is well-defined
OPTS.GEOSTRF = false;

% Verbosity level
OPTS.VERBOSE = 1;

% File ID to write output
OPTS.FILE_ID = 1; % standard output to MATLAB terminal

end