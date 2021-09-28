function [p,s,t,d0,diags] = pot_dens_surf(S, T, P, pref, var, WRAP, OPTS)
%POT_DENS_SURF  Potential density surface by nonlinear solution in each water column.
%
%
% [p,s,t] = pot_dens_surf(S, T, P, pref, d0, WRAP)
% finds the pressure or depth p -- and its salinity s and temperature t --
% of the isosurface d0 of potential density referenced to pref, in the
% ocean with practical / Absolute salinity S and potential / Conservative
% temperature T at data sites where the and pressure or depth is P.  The
% equation of state is given by eos.m in MATLAB's path, which accepts S, T,
% pref as its 3 inputs. For a non-Boussinesq ocean, p, P, and pref should
% be pressure [dbar].  For a Boussinesq ocean, p, P, and pref should be
% depth [m], positive and increasing down.  The domain is periodic in the
% i'th horizontal dimension iff WRAP(i) is true; this is relevant only to
% measuring neutrality errors in diags output, not to constructing the
% surface.
%
% [p,s,t,d0] = pot_dens_surf(S, T, P, pref, [i0, j0, p0], WRAP)
% as above but finds the surface intersecting a reference cast given by
% grid indices (i0,j0) at pressure or depth p0.
%
% ... = pot_dens_surf(..., OPTS)
% sets algorithmic parameters (see "Options" below for further details).
%
%
% --- Input:
% S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% P [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
% pref [1, 1]: reference pressure [dbar] or depth [m]
% var [1, 1] or [1, 3]: isovalue of potential density, or target location
% WRAP [2 element array]: determines which dimensions are treated periodic
%                         [logical].  Set WRAP(i) to true when periodic in 
%                         the i'th lateral dimension(i=1,2).
% OPTS [struct]: options (see below)
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: physical units for S, T, P, pref, p, d0 are determined by eos.m.
%
% Note: P must increase monotonically along its first dimension.
%
%
% --- Output:
% p [ni, nj]: pressure [dbar] or depth [m] of the surface
% s [ni, nj]: practical / Absolute salinity on the surface
% t [ni, nj]: potential / Conservative temperature the surface
% d0 [1, 1]: potential density [kg m^-3] or specific volume [m^3 kg^-1] of
%            the surface
% diags [struct]: diagnostics of the solution and computation time
%
%
% --- Options:
% OPTS is a struct containing a subset of the following fields.
%   FILE_ID [1, 1]: 1 to write any output to MATLAB terminal, or a file
%       identifier as returned by fopen() to write to a file. Default: 1.
%   INTERPFN [function handle]: vertical interpolation function, used to
%       evaluate Sppc and Tppc if those are not provided.  E.g. INTERPFN =
%       @ppc_linterp.
%   Sppc [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose 
%       knots are P that interpolate S as a function of P in each water 
%       column.  E.g. Sppc = ppc_linterp(P, S);
%   Tppc [O, nk-1, ni, nj]: Coefficients for piecewise polynomials whose 
%       knots are P that interpolate T as a function of P in each water 
%       column.  E.g. Tppc = ppc_linterp(P, T);
%   VERBOSE [scalar]: 0 for no output; 1 for summary of each iteration;
%                     2 for detailed information on each iteration.
%   TOL_P_UPDATE [1, 1]: error tolerance, in the same units as P [dbar] or
%       [m], when root-finding to update the surface.
%   DIST1_iJ [ni, nj]: Distance [m] in 1st dimension centred at (I-1/2, J)
%   DIST2_Ij [ni, nj]: Distance [m] in 2nd dimension centred at (I, J-1/2)
%   DIST2_iJ [ni, nj]: Distance [m] in 2nd dimension centred at (I-1/2, J)
%   DIST1_Ij [ni, nj]: Distance [m] in 1st dimension centred at (I, J-1/2)
%
%
% --- Equation of State:
% The MATLAB path* must contain the function eos.m. This must accept 3
% inputs: S, T, and P. eos(S, T, P) returns the specific volume [m^3 kg^-1]
% or the in-situ density [kg m^-3].
% *Note: It is not sufficient to simply have these eos functions in the
% current working directory, because the compiled MEX functions will not be
% able to find them there.  They must be in the MATLAB path.  If they are
% in the current working directory, use `addpath(pwd)` to add the current
% working directory to the top of MATLAB's path.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


[nk, ni, nj] = size(S);
[~,PN] = size(P);
KP = nk * double(PN > 1);

% Set up default options:
DEFS = struct();
DEFS.TOL_P_UPDATE = 1e-4;  % error tolerance in the vertical dimension
DEFS.INTERPFN = @ppc_linterp;  % linear interpolation in the vertical
DEFS.FILE_ID = 1;  % standard output to MATLAB terminal
DEFS.VERBOSE = 1;  % show a moderate level of information
DEFS.DIST1_iJ = 1;   % Distance [m] in 1st dimension centred at (I-1/2, J)
DEFS.DIST2_Ij = 1;   % Distance [m] in 2nd dimension centred at (I, J-1/2)
DEFS.DIST2_iJ = 1;   % Distance [m] in 2nd dimension centred at (I-1/2, J)
DEFS.DIST1_Ij = 1;   % Distance [m] in 1st dimension centred at (I, J-1/2)

% Override any options with user-specified OPTS
if nargin < 6 || isempty(OPTS)
  OPTS = DEFS;
else
  OPTS = catstruct(DEFS, OPTS); 
end

DIAGS = nargout >= 5;

% Run codegen to create MEX function handling the main computation
ni_ = max(ni, 4096); % using variable size code generation and avoiding
nj_ = max(nj, 4096); % recompiling all the time
pot_dens_surf_vertsolve_codegen(nk, ni_, nj_, isvector(P), OPTS);

% Interpolate S and T as piecewise polynomials of P, or use pre-computed interpolants in OPTS.
if isfield(OPTS, 'Sppc') && isfield(OPTS, 'Tppc')
  [~, kIS, iIS, jIS] = size(OPTS.Sppc);
  [~, kIT, iIT, jIT] = size(OPTS.Tppc);
  if nk-1 == kIS && kIS == kIT && ...
      ni  == iIS && iIS == iIT && ...
      nj  == jIS && jIS == jIT
    Sppc = OPTS.Sppc;
    Tppc = OPTS.Tppc;
  else
    Sppc = OPTS.INTERPFN(P, S);
    Tppc = OPTS.INTERPFN(P, T);
  end
else
  Sppc = OPTS.INTERPFN(P, S);
  Tppc = OPTS.INTERPFN(P, T);
end

% Count number of valid bottles per cast
BotK = squeeze(sum(isfinite(S), 1));

% Decide on the isosurface value
if isscalar(var)
    
    d0 = var;
    
else
    % var should be a 3 element vector
    i0 = var(1);
    j0 = var(2);
    p0 = var(3);
    
    % Get linear index to reference cast
    ij0 = sub2ind([ni,nj], i0, j0);
    n = (ij0-1) * KP;
    
    % Select the reference cast
    P0 = P((n+1:n+nk).');
    Sppc0 = Sppc(:,:,ij0);
    Tppc0 = Tppc(:,:,ij0);
    
    % Choose iso-value that will intersect (i0,j0,p0).
    [s0,t0] = ppc_val2(P0, Sppc0, Tppc0, p0);
    d0 = eos(s0, t0, pref);
    
end

% Calculate 3D field for vertical interpolation
if eos(34.5,3,1000) > 1
    sortdir = 'ascend'; % eos is in-situ density, increasing with 3rd argument
else
    sortdir = 'descend'; % eos is specific volume, decreasing with 3rd argument
end
D = sort(eos(S, T, pref), 1, sortdir);

% Get started with the discrete version (and linear interpolation)
p = interp1qn(d0, D, P);

% Start timer after all the setup has been done.
iter_tic = tic();

% Solve non-linear root finding problem in each cast
[p, s, t] = pot_dens_surf_vertsolve_mex(Sppc, Tppc, P, BotK, p, pref, d0, OPTS.TOL_P_UPDATE);

if DIAGS
    diags = struct();
    diags.clocktime = toc(iter_tic);
    [epsL2, epsL1] = eps_norms(s, t, p, true, WRAP, {}, OPTS.DIST1_iJ, OPTS.DIST2_Ij, OPTS.DIST2_iJ, OPTS.DIST1_Ij);
    diags.epsL1 = epsL1;
    diags.epsL2 = epsL2;
end
