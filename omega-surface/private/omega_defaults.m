function OPTS = omega_defaults()
%OMEGA_DEFUALTS  Default options for omega_surface
%
%
% OPTS = omega_defaults()
% returns a struct OPTS containing default options for use in
% omega_surface.

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

OPTS = struct();

OPTS.MLX = []; % Do not remove the Mixed Layer

OPTS.FIGS_SHOW = false; % do not show figures

OPTS.INTERPFN = @ppc_linterp; % Use linear interpolation in the vertical dimension.

OPTS.SppX = [];  % Pre-computed interpolation functions.  None given here.
OPTS.TppX = [];  % Pre-computed interpolation functions.  None given here.

OPTS.ITER_MAX = 10; % maximum number of iterations

OPTS.ITER_START_WETTING = 1; % start wetting immediately

OPTS.TOL_DENS = 1e-7; % Density tolerance [kg m^-3]

% Relative tolerance for LSQR. Since the matrix problem is overdetermined,
% the relative residual will, in general, exceed this tolerance bound. As
% such, it is challenging to relate LSQR's relative tolerance to physical
% tolerances on density or pressure.  From several numerical tests,
% the default relative tolerance of 1e-6 seems to work well.
OPTS.TOL_LSQR_REL = 1e-6;

OPTS.VERBOSE = 1; % show a moderate level of information

% The matrix has 1's and 0's everywhere, except this value in the final row
OPTS.FINAL_ROW_VALUES = 1e-2; % chosen empirically from tests on 1x1deg OCCA data

OPTS.INTEGRATING_FACTOR = [];  % No integrating factor

OPTS.FILE_ID = 1; % standard output to MATLAB terminal

OPTS.REF_IJ = []; % No reference water column at which the surface is pinned.  Instead, maintain the mean density. 