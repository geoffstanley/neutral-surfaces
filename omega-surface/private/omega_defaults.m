function OPTS = omega_defaults()
%OMEGA_DEFUALTS  Default options for omega_surface
%
%
% OPTS = omega_defaults()
% returns a struct OPTS containing default options for use in
% omega_surface.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


OPTS = struct();

OPTS.DIST1_iJ = 1;   % Distance [m] in 1st dimension centred at (I-1/2, J)
OPTS.DIST2_Ij = 1;   % Distance [m] in 2nd dimension centred at (I, J-1/2)
OPTS.DIST2_iJ = 1;   % Distance [m] in 2nd dimension centred at (I-1/2, J)
OPTS.DIST1_Ij = 1;   % Distance [m] in 1st dimension centred at (I, J-1/2)

OPTS.ML = []; % Do not remove the Mixed Layer

OPTS.FIGS_SHOW = false; % do not show figures

OPTS.INTERPFN = @ppc_linterp; % Use linear interpolation in the vertical dimension.

OPTS.Sppc = [];  % Pre-computed interpolation coefficients.  None given here.
OPTS.Tppc = [];  % Pre-computed interpolation coefficients.  None given here.

OPTS.ITER_MIN = 1;  % minimum number of iterations
OPTS.ITER_MAX = 10; % maximum number of iterations

OPTS.ITER_START_WETTING = 1; % start wetting immediately
OPTS.ITER_STOP_WETTING  = 5; % stop wetting after this many iterations (to avoid flip-flopping on adding then removing some nuisance casts)

% Exit iterations when the L1 change of the Locally Referenced Potential
% Density perturbation is less than this value [kg m^-3].  Set to 0 to deactivate. 
OPTS.TOL_LRPD_L1 = 1e-7;

% Exit iterations when the L2 change of pressure (or depth) on the surface
% is less than this value. Set to 0 to deactivate. Units are the same as P [dbar or m].
OPTS.TOL_P_CHANGE_L2 = 0;

% Error tolerance when root-finding to update surface, in the same units as
% P [dbar] or [m].
OPTS.TOL_P_UPDATE = 1e-4;

OPTS.VERBOSE = 1; % show a moderate level of information

OPTS.FILE_ID = 1; % standard output to MATLAB terminal

OPTS.POISSON = true; % Whether to solve the square, symmetric Poisson matrix problem, or a rectangular gradient problem