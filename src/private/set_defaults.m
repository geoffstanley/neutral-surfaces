function OPTS = set_defaults(nx, ny)
%SET_DEFUALTS  Default options for topobaric_surface.
%
%
% OPTS = set_defaults(nx,ny) 
% returns a struct OPTS containing default options for use in
% topobaric_surface, appropriate for a grid with nx longitude points and ny
% latitude points.

% --- Copyright:
% Copyright 2019 Geoff Stanley
%
% This file is part of Topobaric Surface.
% 
% Topobaric Surface is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
% 
% Topobaric Surface is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with Topobaric Surface.  If not, see
% <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au 
% Email     : geoffstanley@gmail.com
% Version   : 1.0
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
OPTS.SPLINE_BREAKS = [0 200 1500 1800 6000];
OPTS.SPLINE_ORDER = 4; % Cubic

% Calculate the Mixed Layer Pressure using default values from
% mixed_layer_pressure()
OPTS.MLP = struct(); 

% Pixel indices for reference water column. Trying to select (180E,0N)
OPTS.REF_IJ = round([nx/2, ny/2]);

% Reference pressure (if OPTS.RHOB is not provided). If a scalar is
% provided, the pressure on the topobaric surface at the reference water
% column will be this value (with precision of OPTS.TOL).
OPTS.REF_P = []; % [dbar]

% Reference depth (if OPTS.RHOB is provided). If a scalar is provided, the
% depth on the topobaric surface at the reference water column will be this
% value (with precision of OPTS.TOL).
OPTS.REF_Z = []; % [m, positive]

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
OPTS.FILL_PIX = 0; % No filling
OPTS.FILL_IJ = []; % No filling

% Post-processing of Reeb Graph -- graph simplification parameters:
OPTS.SIMPLIFY_ARC_REMAIN = Inf;     % No leaf pruning simplification
OPTS.SIMPLIFY_WEIGHT_PERSIST = 0.5; % Equal weighting between area and persistence

% Error tolerance when root-finding to update surface. If OPTS.RHOB is
% empty, provide as [dbar]. If OPTS.RHOB is given, provide as [m]
% (topobaric_surface will internally convert it to dbar).
OPTS.TOL = 1e-4; 

% Damping when updating to the new surface. 
% No damping at 0. Set > 0 for damping. Set < 0 for overrelaxation.
OPTS.DAMP = 0; 

% Conditions to terminate iterations:
OPTS.ITER_MAX = 6;  % maximum number of iterations

% quit when the L2 change of pressure on surface exceeds this value (set to
% 0 to deactivate). If OPTS.RHOB is empty, provide as [dbar]. If OPTS.RHOB
% is given, provide as [m] (topobaric_surface will internally convert it to
% dbar).
OPTS.ITER_L2_CHANGE = 1e-3; 

% When true, adjust the empirical function to ensure exact geostrophic
% stream function is well-defined
OPTS.GEOSTRF = false;

% Verbosity level 
OPTS.VERBOSE = 0;

% File ID to write output
OPTS.FILE_ID = 1; % standard output to MATLAB terminal

end