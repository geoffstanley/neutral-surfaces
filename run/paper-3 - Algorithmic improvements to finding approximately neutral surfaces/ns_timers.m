% Comparison of the various surfaces in the Neutral-Surfaces toolbox

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2020  Geoff Stanley
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
%
% Modified by : --
% Date        : --
% Changes     : --

%% --- BEGIN SETUP --------------------------------------------------------
%#ok<*UNRCH>

warning('off', 'MATLAB:nargchk:deprecated')
set(0, 'defaultfigurecolor', [1 1 1]); % white figure background
V = filesep(); % /  or  \  depending on OS.

%PATH_LOCAL = [fileparts(mfilename('fullpath')) V]; % Get path to this file.
PATH_LOCAL = '~/work/unsw/projects/omega/'; % Manually set path to this file.
cd(PATH_LOCAL);

PATH_OMEGA_KMJ = '~/work/unsw/projects/omega/omega_ver1_0/';
addpath(genpath(PATH_OMEGA_KMJ));

% Add Neutral Surfaces to MATLAB's path
PATH_NS = '~/work/projects-gfd/neutral-surfaces/';
run([PATH_NS 'ns_add_to_path.m']);

% Folder containing the functions eos.m, eos_x.m, and eos_s_t.m
PATH_EOS = '~/work/MATLAB/eos/eos/';
addpath(PATH_EOS);  % Doing this last to ensure at top of MATLAB's path, above eos fcns in other dirs

% Folder containing .nc data files
PATH_OCCA = '~/work/data/OCCA/';
PATH_ECCO2 = '~/work/data/ECCO2/latlon/';

% Make a directory for figures
PATH_FIGS = [PATH_LOCAL 'figs' V];
PATH_FIGS = [PATH_FIGS V 'shootout3_' datestr(now, 'dd-mm-yyyy') V];
if ~exist(PATH_FIGS, 'dir')
  mkdir(PATH_FIGS)
end

fileID = 1; % For standard output to the screen
%fileID = fopen([PATH_LOCAL 'run ' datestr(now, 'yyyy mm dd hh MM ss') '.txt'], 'wt'); % For output to a file

maxNumCompThreads(1);   % Just use one CPU, for apples-to-apples comparison

REPEAT = true; % for paper
% REPEAT = false; nrep = 1; % for quicker tests

% RES = 2.^(4:7); %MODELS = {};  % for quicker tests
% RES = 2.^(4:8); %MODELS = {'OCCA'};  % for quicker tests
RES = 2.^(4:11); MODELS = {'OCCA', 'ECCO2'};  % Everything. For paper

RES = RES + 1; % Add 1 for the wall in each dimension... so # ocean grid points is a power of 2
nRES = length(RES);
nMODELS = length(MODELS);
nDATA = nRES + nMODELS;
NWC = nan(nDATA,1); % Number of Water Columns

rho_c = 1027.5; % Boussinesq reference density [kg m-3]
grav = 9.81;    % gravitational acceleration [m s-2]

db2Pa = 1e4; % dbar to Pa conversion
Pa2db = 1e-4; % Pa to dbar conversion

lead1 = @(x) reshape(x, [1 size(x)]);

% Choose vertical interpolation method
INTERPFN = @ppc_linterp;

%% Set alias functions
% Choose the Boussinesq densjmd95 and set grav and rho_c in eos.m and eos_x.m
if false  % <-- only need to do this once!  Doing every time makes codegen run every time
  eoscg_set_bsq_param([PATH_NS 'lib' V 'eos' V 'eoscg_densjmd95_bsq.m'   ] , [PATH_EOS 'eos.m'  ], grav, rho_c);
  eoscg_set_bsq_param([PATH_NS 'lib' V 'eos' V 'eoscg_densjmd95_bsq_dz.m'] , [PATH_EOS 'eos_x.m'], grav, rho_c);
  eoscg_set_bsq_param([PATH_NS 'lib' V 'eos' V 'eoscg_densjmd95_bsq_s_t.m'], [PATH_EOS 'eos_s_t.m'], grav, rho_c);
end

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


%% Selecting surfaces and depths

list_surf = {'\sigma', '\delta', '\omega_{KMJ}', '\omega_{\nabla}', '\omega_+', '\tau'};
LIST_SURF = {'SIGMA', 'DELTA', 'KMJ', 'OMEGAGRAD', 'OMEGA', 'TOPOB'};

nS = length(list_surf);
iSURF = containers.Map(LIST_SURF, 1:nS);

% Save depths and diagnostic outputs from various surfaces
zs    = cell(nDATA, nS);
diags = cell(nDATA, nS);
timer = nan( nDATA, nS);

%% Begin looping over data
for iDATA = 1 : nDATA
  %% Load grid and hydrographic data
  fprintf(1, '--- Starting grid # %d\n', iDATA);
  
  if iDATA <= nRES
    % do Synthetic data first
    DATA_SOURCE = 'SYNTHRAND';
    
    ni = RES(iDATA);
    nj = ni;
    nk = 50;
    
    [S, T, Z, g] = synthocean_rand(ni, nj, nk);
    g.DXC = repmat(g.DXCvec, ni, 1);
    g.DYC = repmat(g.DYCsc, ni, nj);
    g.DXG = repmat(g.DXGvec, ni, 1);
    g.DYG = repmat(g.DYGsc, ni, nj);
    
    %S(:,2:end,:) = repmat(S(:,2,:), [1, ni-1, 1]); % DEV: Zero helicity option
    %T(:,2:end,:) = repmat(T(:,2,:), [1, ni-1, 1]); % DEV: Zero helicity option
    
    z0 = 1500;
    
    i0 = ceil(ni/2); % middle i
    j0 = ceil(nj/2); % middle j
    
  else
    
    DATA_SOURCE = MODELS{iDATA - nRES};
    
    if DATA_SOURCE(1) == 'O'
      [g, S, T, ~, ETAN, ATMP] = load_OCCA(PATH_OCCA);
    elseif DATA_SOURCE(1) == 'E' % ECCO2
      TIMESTEP = '20021223';
      g = load_ECCO2(PATH_ECCO2, 'grid');
      [S, T, ATMP, ETAN] = load_ECCO2(PATH_ECCO2, 'casts', TIMESTEP);
    end
    ni = g.nx;
    nj = g.ny;
    nk = g.nz;
    
    RES(iDATA) = sqrt(ni * nj);  % If this were a square, it'd have this side-length... add to RES vector
    
    assert(rho_c == g.rho_c, 'rho_c differs between this model and this script');
    assert(grav == g.grav, 'grav differs between this model and this script');
    Z = -g.RC(:); % We are going to be Z > 0 people!
    
    
    % Re-order data so water-columns are contiguous data:
    S = permute(S, [3 1 2]); % [nz,ny,nx]. depth  by  lat  by  long
    T = permute(T, [3 1 2]);
    if exist('GAMMA', 'var')
      GAMMA = permute(GAMMA, [3 1 2]);
      GAMMA(GAMMA < 0) = nan;  % change bad values (should be -99) to nan
    end
    
    z0 = 1500;
    
    y_to_j = @(y) floor((y - g.YGvec(1)) * g.resy) + 1;
    x_to_i = @(x) mod(floor((x - g.XGvec(1)) * g.resx), ni) + 1;
    
    % Start in the Pacific Ocean
    x0 = 180; y0 = 0;
    i0 = x_to_i(x0);
    j0 = y_to_j(y0);
    
  end
  
  g.RAC = repmat(g.RACvec, ni, 1);
  g.RAS = repmat(g.RASvec, ni, 1);
  g.RAW = repmat(g.RAWvec, ni, 1);
  g.RAZ = repmat(g.RAZvec, ni, 1);
  g.YC = repmat(g.YCvec, ni, 1);
  g.YG = repmat(g.YGvec, ni, 1);
  g.DXC = repmat(g.DXCvec, ni, 1);
  g.DYC = repmat(g.DYCsc, ni, nj);
  g.DXG = repmat(g.DXGvec, ni, 1);
  g.DYG = repmat(g.DYGsc, ni, nj);
  
  nij = ni * nj;
  I0 = sub2ind([ni nj], i0, j0);
  
  BotK = squeeze(sum(isfinite(S), 1));
  
  SppZ = INTERPFN(Z, S);
  TppZ = INTERPFN(Z, T);
  
  Z2P = Pa2db * rho_c * grav ; % Note > 0
  
  Szyx = permute(S, [1 3 2]); % KMJ omega needs dimensions to be [z,y,x].
  Tzyx = permute(T, [1 3 2]);
  Zzyx = repmat(Z(:), [1, nj, ni]);
  
  NWC(iDATA) = sum(flat(isfinite(S(1,:,:))));
  
  
  % --- Setup options neutral surfaces
  OPTS = struct();
  OPTS.WRAP = WRAP;       % omega's Poisson code now requires doubly periodic domain; we have nan walls to prevent wrapping though, as needed. 
  OPTS.TOL_X_CHANGE_L2 = 1e-3; % Stop iterations when the L2 change of pressure on the surface is less than this
  OPTS.ITER_MIN = 5;      % .. but do at least 5 iterations
  OPTS.ITER_MAX = 30;     % The maximum number of iterations
  %OPTS.ITER_START_WETTING = 1; % Start wetting on first iteration
  OPTS.ITER_START_WETTING = inf; % No Wetting
  OPTS.REF_IJ = [i0 j0]; % Pacific Ocean
  OPTS.REEB = true;       % Do topobaric, not orthobaric, surfaces
  OPTS.GEOSTRF = false;  % Don't add extra criteria that make geostrophic streamfunction well-defined
  OPTS.SIMPLIFY_ARC_REMAIN = Inf; % No simplification
  OPTS.FILL_IJ = [];     % No filling
  OPTS.FILL_PIX = 0;     % No filling.
  OPTS.MLX = [];         % Leave the Mixed Layer in.
  OPTS.TOL_X_UPDATE = 1e-4; % error tolerance when updating the surface
  OPTS.X_EXPN = 500;     % expansion of domain to search for solutions in each water column
  OPTS.TOL_LRPD_L1 = 0; % Don't use the LRPD stopping criterion.
  OPTS.VERBOSE = 1;      % Some output while executing functions
  OPTS.FILE_ID = fileID; % Write output to this file
  OPTS.FIGS_SHOW = false; % Don't show figures (for omega surface)
  OPTS.INTERPFN = INTERPFN; % set interpolation function to match this script
  OPTS.SppX = SppZ;      % Use pre-computed piecewise polynomial for S in terms of X=Z
  OPTS.TppX = TppZ;      % Use pre-computed piecewise polynomial for T in terms of X=Z
  
  OPTS.DIST1_iJ = g.DXC;   % Distance [m] in 1st dimension centred at (I-1/2, J)
  OPTS.DIST2_Ij = g.DYC;   % Distance [m] in 2nd dimension centred at (I, J-1/2)
  OPTS.DIST2_iJ = g.DYG;   % Distance [m] in 2nd dimension centred at (I-1/2, J)
  OPTS.DIST1_Ij = g.DXG;   % Distance [m] in 1st dimension centred at (I, J-1/2)
  
%   OPTS.DIST1_iJ = 1;   % Distance [m] in 1st dimension centred at (I-1/2, J)
%   OPTS.DIST2_Ij = 1;   % Distance [m] in 2nd dimension centred at (I, J-1/2)
%   OPTS.DIST2_iJ = 1;   % Distance [m] in 2nd dimension centred at (I-1/2, J)
%   OPTS.DIST1_Ij = 1;   % Distance [m] in 1st dimension centred at (I, J-1/2)
  

  % Call eos() functions to do Just-in-Time Compilation, so that it's not done in one or another function later
  s = squeeze(S(1,:,:));
  t = squeeze(T(1,:,:));
  z = Z(1);
  eos(s, t, z);
  eos_x(s, t, z);
  eos_s_t(s, t, z);
  
  fprintf(1, 'Loaded %s (%d,%d)\n', DATA_SOURCE, ni, nj);
  
  % --- Calculate surfaces
  %% Potential Density Surface
  if 1
  zref = z0;
  
  mytime = 0;
  if REPEAT
    nrep = max(1, ceil(log2(1024^2 ./ NWC(iDATA)))); % so N x N data gets done once where N=1024.  Others get about same CPU time, if they scale as N^1.
  end
  for rep = 1 : nrep
    [z, s, t, isoval, d] = pot_dens_surf(S, T, Z, zref, [i0, j0, z0], OPTS);
    
    zs{iDATA, iSURF('SIGMA')} = z;
    
    mytime = mytime + sum(d.clocktime);
  end
  
  diags{iDATA, iSURF('SIGMA')} = d;
  timer(iDATA, iSURF('SIGMA')) = mytime / nrep;
  fprintf(fileID, '(%d,%d,%.4fm): SIGMA = %.4f, z_ref = %.4f in time %.2f\n', ...
    i0, j0, z0, isoval, zref, mytime/nrep);
  
  z_sigma = z; % another alias
  s_sigma = s;
  t_sigma = t;
  d_sigma = d;
  end
  %% In-situ density anomaly
  mytime = 0;
  if REPEAT
    nrep = max(1, ceil(log2(1024^2 ./ NWC(iDATA)))); % so N x N data gets done once where N=1024.  Others get about same CPU time, if they scale as N^1.
  end
  for rep = 1 : nrep
    
    % For OCCA or ECCO2, consider choosing the reference s and t values from the Southern Ocean, the nexus of the other oceans.
    %s0 = s_sigma(270,25); % <-- for OCCA
    %t0 = t_sigma(270,25); % <-- for OCCA
    %[z, s, t, isoval, s0, t0, d] = delta_surf(S, T, Z, s0, t0, [i0, j0, z0], OPTS);
    % ... But actually, the results for OCCA are marginally better if we take s0, t0 from the reference Pacific cast.

    [z, s, t, isoval, s0, t0, d] = delta_surf(S, T, Z, [], [], [i0, j0, z0], OPTS);
    
    zs{iDATA, iSURF('DELTA')} = z;
    
    mytime = mytime + sum(d.clocktime);
  end
  
  diags{iDATA, iSURF('DELTA')} = d;
  timer(iDATA, iSURF('DELTA')) = mytime / nrep;
  fprintf(fileID, '(%d,%d,%.4fm): DELTA = %.4f, in time %.2f\n', ...
    i0, j0, z0, isoval, mytime/nrep);
  
  z_delta = z; % another alias
  s_delta = s;
  t_delta = t;
  d_delta = d;
  
  
  % label, internally consistent with this dataset:
  l0 = eos(S(1,i0,j0), T(1,i0,j0), 0) - 1000;
  label = int_grad_lrpd(Z(1), l0, S(:,i0,j0), T(:,i0,j0), Z, z_delta(i0,j0));
  
  %% KMJ's Omega surface
  
  if 1 && RES(iDATA) < 256  % Do up to 128^2 and do OCCA.
    if REPEAT
      nrep = max(1, ceil(log2(128^2 ./ NWC(iDATA)))); % so N x N data gets done once where N=1024.  Others get about same CPU time, if they scale as N^1.
    end
    mytime = 0;
    for rep = 1 : nrep
      
      % The following are all [ny, nx]
      lats = repmat(g.YCvec(:), 1, ni);
      lons = repmat(g.XCvec(:)', nj, 1);
      e1t = OPTS.DIST1_iJ.';
      e2t = OPTS.DIST2_Ij.';
      e1g = OPTS.DIST1_Ij.';
      e2g = OPTS.DIST2_iJ.';
      
      S_mid = 0.5*(Szyx(2:nk,:,:) + Szyx(1:nk-1,:,:));
      T_mid = 0.5*(Tzyx(2:nk,:,:) + Tzyx(1:nk-1,:,:));
      Z_mid = 0.5*(Zzyx(2:nk,:,:) + Zzyx(1:nk-1,:,:));
      r = eos(S_mid, T_mid, Z_mid * Z2P);
      [rs, rt] = eos_s_t(S_mid, T_mid, Z_mid * Z2P);
      clear S_mid T_mid Z_mid
      dS = diff(Szyx,[],1);
      dT = diff(Tzyx,[],1);
      dZ = diff(Zzyx,[],1);
      n2 = grav * (rs .* dS + rt .* dT) ./ (dZ .* r);
      clear dS dT dZ r rs rt
      
      nit = 200; % max iterations
      
      [s, t, z, d] = optimize_surface(Szyx, Tzyx, Zzyx, repmat(grav,nj,ni), n2, lead1(s_delta'), lead1(t_delta'), lead1(z_delta'), e1t, e2t, nit, 'epsilon', 'long', Z2P, e1g, e2g);
      s = squeeze(s).';
      t = squeeze(t).';
      z = squeeze(z).';
      z_omega0 = z;
      s_omega0 = s;
      t_omega0 = t;
      
      mytime = mytime + sum(d.clocktime);
    end
    timer(iDATA, iSURF('KMJ')) = mytime / nrep;
    diags{iDATA, iSURF('KMJ')} = d;
    zs{iDATA, iSURF('KMJ')} = z;
    
    fprintf(fileID, '(%d,%d,%.4fm): KMJ done in time %.2f\n', ...
      i0, j0, z0, mytime / nrep);
  end
  
  %% Omega surface - updated but still gradient equations
  if 1 && ~strcmp(DATA_SOURCE, 'ECCO2')
    % Gradient equations stall on LSQR solution for ECCO2, due to the large
    % disparity of non-zero entries in the matrix, which are the square root
    % of a zonal distance divided by a meridional distance, or vice versa.
    % The ECCO2 grid reaches near the poles, so this ranges from
    % 7.8251096e-09  to 21.409498.  It may be possible to obtain a solution
    % by turning down the relative tolerance of the LSQR solver, but then
    % we'll need to assess the quality of that solution, which would require
    % a separate figure.  Thus, we sidestep the issue by not computing the
    % ECCO2 solution for the gradient omega equations.
    
    % NOTE:  Timing (and accuracy) depends on on TOL_LSQR_REL in omega_matsolve_grad.m
    
    OPTS.POISSON = false;
    mytime = 0;
    if REPEAT
      nrep = max(1, ceil(log2(1024^2 ./ NWC(iDATA)))); % so N x N data gets done once where N=1024.  Others get about same CPU time, if they scale as N^1.
    end
    for rep = 1 : nrep
      [z, ~, ~, d] = omega_surface(S, T, Z, z_delta, OPTS);
      mytime = mytime + sum(d.clocktime(1:OPTS.ITER_MIN));
    end
    timer(iDATA, iSURF('OMEGAGRAD')) = mytime / nrep;
    diags{iDATA, iSURF('OMEGAGRAD')} = d;
    zs{iDATA, iSURF('OMEGAGRAD')} = z;
    
    fprintf(fileID, '(%d,%d,%.4fm): OMEGAGRAD done in time %.2f\n', ...
      i0, j0, z0, mytime / nrep);
  end
  
  
  %% Omega surface -- Poisson version
  if 1
    
    OPTS.POISSON = true;
    mytime = 0;
    if REPEAT
      nrep = max(1, ceil(log2(1024^2 ./ NWC(iDATA)))); % so N x N data gets done once where N=1024.  Others get about same CPU time, if they scale as N^1.
    end
    for rep = 1 : nrep
      [z, ~, ~, d] = omega_surface(S, T, Z, z_delta, OPTS);
      mytime = mytime + sum(d.clocktime(1:OPTS.ITER_MIN));
    end
    timer(iDATA, iSURF('OMEGA')) = mytime / nrep;
    diags{iDATA, iSURF('OMEGA')} = d;
    zs{iDATA, iSURF('OMEGA')} = z;
    
    fprintf(fileID, '(%d,%d,%.4fm): OMEGA done in time %.2f\n', ...
      i0, j0, z0, mytime / nrep);
  end
  
  
  %% Calculate topobaric surface
  if 1    
    mytime = 0;
    if REPEAT
      nrep = max(1, ceil(log2(1024^2 ./ NWC(iDATA)))); % so N x N data gets done once where N=1024.  Others get about same CPU time, if they scale as N^1.
    end
    for rep = 1 : nrep
      [z, ~, ~, RG, ~, ~, ~, d] = topobaric_surface(S, T, Z, z_delta, OPTS);
      mytime = mytime + sum(d.clocktime(1:OPTS.ITER_MIN));
    end
    timer(iDATA, iSURF('TOPOB')) = mytime / nrep;
    diags{iDATA, iSURF('TOPOB')} = d;
    zs{iDATA, iSURF('TOPOB')} = z;
    
    % Show map of regions:
    %arcvec = nan(RG.n_casts,1);  for a = 1:RG.nArcs; arcvec(RG.arc_segment{a}) = a; end; arcmap = nan(ni,nj); arcmap(RG.wet) = arcvec; figi(arcmap)
    
    fprintf(fileID, '(%d,%d,%.4fm): TOPOBARIC done in time %.2f\n', ...
      i0, j0, z0, mytime / nrep);
  end
  
  %end % iZ
end % iDATA

% Extract # of Surface Points per surface
NSP = nan(nDATA, nS);
for iDATA = 1 : nDATA
  for iS = 1 : nS
    NSP(iDATA, iS) = sum(isfinite(flat(zs{iDATA, iS})));
  end
end

cm = [ ... % Colormap for lines
  0.6       0.6       0.6 % grey
  0.0       0.0       0.0 % black
  0.0       0.0       1.0 % blue
  0.5       0.2       0.7 % purple
  1.0       0.0       0.0 % red
  0.0       0.5       0.0 % dark green
  0.85      0.6       0.0 % gold
  ];
markers = 'o*^v+d';

%return

%% Save data
% fn = sprintf('results %s - %s.mat', datestr(now, 'yy-mm-dd hh-MM-ss'), num2str(RES-1));
% save(fn); % save all variables


%% Diagnostics plot: CPU time vs. Resolution

hf = figure('Position', [1920-800, 1080-480, 800, 600]);
ax = axes('Position', [.1, .11, .74, .8]); hold on; grid on; box on;
fn = @log2;
%fn = @(z) z;

% Extract # of iterations per
iters = nan(nDATA, nS);
for iDATA = 1 : nDATA
  for iS = 1 : nS
    d = diags{iDATA, iS};
    if isempty(d); continue; end
    if iS == iSURF('OMEGAGRAD') || iS == iSURF('OMEGA')|| iS == iSURF('TOPOB')
      i = find(d.x_change_L2 <= 1e-3, 1, 'first');
      if isempty(i)
        iters(iDATA, iS) = 99; % flag.  shouldn't happen.
      else
        iters(iDATA, iS) = i;
      end
    else
      % count total # iterations used
      iters(iDATA, iS) = length(d.clocktime);
    end
  end
end

for iS = 1 : nS
  plot(ax, fn(NSP(1:nRES,iS)), fn(timer(1:nRES,iS)), '-', 'color', cm(iS,:), 'MarkerFace', cm(iS,:), 'marker', markers(iS));
end

% MODEL data
model_marker = 'ph';
for iS = 1 : nS
  for iDATA = nRES + 1 : nDATA
    scatter(ax, fn(NSP(iDATA,iS)), fn(timer(iDATA,iS)), 80, cm(iS,:), model_marker(iDATA-nRES), 'MarkerFaceColor', cm(iS,:));
  end
end
xl = ax.XLim;
yl = ax.YLim;

% Add lines for n^slope
x1 = repmat(fn(max(NWC)), nS, 1); x1(iSURF('KMJ')) = 14;
slopes = nan(nS,1);
slopes(iSURF('SIGMA')) = 1;
slopes(iSURF('DELTA')) = 1;
slopes(iSURF('KMJ')) = 1.6;
slopes(iSURF('OMEGAGRAD')) = 1.6;
slopes(iSURF('OMEGA')) = 1.2;
for iS = [ iSURF('DELTA'), iSURF('KMJ'), iSURF('OMEGAGRAD'), iSURF('OMEGA')]
  %for iS = 1 : nS
  i = find(isfinite(timer(1:nRES,iS)), 1, 'last');
  if isempty(i); continue; end
  slope = slopes(iS);
  xx = [fn(NSP(1,iS)), fn(NSP(i,iS)) + .4];
  x2 = fn(NSP(i,iS));
  y2 = fn(timer(i,iS));
  ln = @(x) y2 + slope * (x - x2);
  plot(ax, xx, ln(xx), '--', 'Color', cm(iS,:));  % constant slope line
  
  x = x1(iS) + 0.4;
  y = ln(x) - 0.3;
  txt = ['N^{' num2str(slope) '}'];
  text(ax, x, y, txt, 'Color', cm(iS,:));
end

% Add lines for n log n
for iS = iSURF('TOPOB')
  
  i = find(isfinite(timer(1:nRES,iS)), 1, 'last');
  if isempty(i); continue; end
  
  xx = linspace(log2(NSP(1,iS)), log2(NSP(i,iS)), 100);
  x2 = xx(end);
  y2 = log2(timer(i,iS));
  C = y2 - x2 - log2(x2);
  ln = @(x) C + x + log2(x);
  plot(ax, xx, ln(xx), '--', 'Color', cm(iS,:));  % constant slope line
  
  txt = 'N logN';
  text(ax, x2 + .2, y2 + .2 , txt, 'Color', cm(iS,:));
end


ax.YLim = [-11, 12];
ax.XLim = [log2(min(NWC)) - 1, log2(max(NWC)) + 1.5];
ax.YTick = ax.YLim(1) : ax.YLim(2); % all integers, only
ax.XTick = log2(NWC(1:nRES));
ax.XLabel.String = 'log_2({\itN})';
ax.YLabel.String = 'log_2(CPU time / 1s)';
ax.XLabel.Interpreter = 'tex';
ax.YLabel.Interpreter = 'tex';
leg = [list_surf, MODELS];
legend(leg, 'Location', 'northwest');

% Write # of iterations.  skip SIGMA and DELTA
for iDATA = 1 : nDATA
  for iS = [iSURF('KMJ'), iSURF('OMEGAGRAD'), iSURF('OMEGA'), iSURF('TOPOB')]
    x = fn(NSP(iDATA,iS)) - 0.005 * diff(ax.XLim);
    if iS == iSURF('OMEGA')
      y = fn(timer(iDATA,iS)) - 0.02 * diff(ax.YLim); % below
    else
      y = fn(timer(iDATA,iS)) + 0.02 * diff(ax.YLim); % above
    end
    txt = num2str(iters(iDATA,iS));
    text(ax, x, y, txt, 'Color', cm(iS,:), 'HorizontalAlignment', 'right');
  end
end

%% Add second axes that show non-logged values
ax2 = axes('Position', [.1, .11, .74, .8]);
ax2.Color = 'none';
ax2.XTick = [];
ax2.XAxisLocation = 'top';
ax2.XTick = ax.XTick;
ax2.XLim = ax.XLim;
ax2.XTickLabel = arrayfun(@(c) sprintf('%d^2', 2.^(c/2)), ax.XTick, 'UniformOutput', false);
ax2.XLabel.String = '{\it N} = # surface points';
ax2.XLabel.Interpreter = 'tex';

ax2.YAxisLocation = 'right';
ax2.YTick = ax.YTick;
ax2.YLim = ax.YLim;
for i = 1 : length(ax2.YTickLabel)
  if 2 ^ ax.YTick(i) >= 1
    ax2.YTickLabel{i} = [sprintf('%.1d', 2.^ax.YTick(i)) ' s'];
  else
    ax2.YTickLabel{i} = [sprintf('%.1f', 100 * 2.^ax.YTick(i)) ' ms'];
  end
end
ax2.YLabel.String = 'CPU time';
ax2.YLabel.Interpreter = 'tex';

%%  Save Figure
fn = sprintf('cpu-vs-gridres__log2__%s', datestr(now, 'yy-mm-dd hh-MM-ss'));
export_fig(hf, [PATH_FIGS fn], '-pdf')
%close(hf);

%% Diagnostics plot:  eps L2 norm vs CPU time

iDATA = find(RES == 128+1);  DATA_SOURCE = 'SYNTHRAND';  % Select data source for this plot

hf = figure;
fs = 12;
clear ax;
ax = axes('Position', [.1 .14 .95-.1 .8]);

timefac = ones(nS,1);
timefac(iSURF('KMJ')) = 60;  % measure KMJ in minutes not seconds
timefacstr = repmat('sec', nS, 1);
timefacstr(iSURF('KMJ'),:) = 'min';

% Precompute clocktimes and maximum clocktime of any algorithm
clocktime = cell(nS,1);
maxtime = 0;
for iS = 1:nS
  d = diags{iDATA, iS};
  if ~isempty(d)
    if isscalar(d.clocktime)
      clocktime{iS} = d.clocktime / timefac(iS);
    else
      clocktime{iS} = [0; cumsum(d.clocktime)] / timefac(iS);
    end
    maxtime = max(maxtime, clocktime{iS}(end));
  end
end

epsL2 = nan(nDATA, nS);

iS = 1;
d = diags{iDATA, iS};
if ~isempty(d)
  semilogy(ax, [0 maxtime], [1 1] * d.epsL2(1), '--', 'Color', cm(iS,:), 'LineWidth', 1);
  hold(ax, 'on');  % must be after semilogy
  grid(ax, 'on')
  ax.FontSize = 13;
  epsL2(iDATA, 1) = d.epsL2;
end

for iS = [iSURF('SIGMA'), iSURF('DELTA'), iSURF('TOPOB'), iSURF('KMJ'), iSURF('OMEGAGRAD'), iSURF('OMEGA')]
  d = diags{iDATA, iS};
  if ~isempty(d)
    if iS==2
      foo = semilogy(ax, clocktime{iS}, d.epsL2, '-', 'Color', cm(iS,:), 'LineWidth', 2, 'Marker', markers(iS), 'MarkerFace', cm(iS,:));
    else
      semilogy(ax, clocktime{iS}, d.epsL2, '-', 'Color', cm(iS,:), 'LineWidth', 2, 'Marker', markers(iS), 'MarkerFace', cm(iS,:));
    end
    semilogy(ax, [0 maxtime], [1 1] * d.epsL2(end), '--', 'Color', cm(iS,:), 'LineWidth', 1);
    epsL2(iDATA, iS) = d.epsL2(end);
    
    ha = 'right';
    vertalign = 'top';
    x = maxtime * .99;
    y = 10^(log10(epsL2(iDATA,iS)) + 0.002 * log10(diff(ax.YLim)));
    txt = sprintf('$%s$: %.2e kg m$^{-3}$, %.2f %s', list_surf{iS}, epsL2(iDATA,iS), clocktime{iS}(end), timefacstr(iS,:));
    
    if iS == iSURF('KMJ')
      y = 10^(log10(epsL2(iDATA,iS)) + 0.008 * log10(diff(ax.YLim)));
    elseif iS == iSURF('OMEGAGRAD')
      x = clocktime{iS}(end) / timefac(iS);
      y = 10^(log10(epsL2(iDATA,iS)) - 0.008 * log10(diff(ax.YLim)));
      vertalign = 'bottom';
      ha = 'left';
    elseif iS == iSURF('OMEGA')
      x = clocktime{iS}(end) / timefac(iS);
      ha = 'left';
      vertalign = 'top';
    end
    text(ax, x, y, txt, 'Interpreter', 'latex', 'FontSize', fs, 'Color', cm(iS,:), 'VerticalAlignment', vertalign, 'HorizontalAlignment', ha);
    
  end
end

yl = [1e-11 3e-9];
ax.YLim = yl;

ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.String = 'CPU time [sec or min]';
ax.YLabel.String = '$|\!| \epsilon |\!|_2 \quad [\mathrm{kg} \ \mathrm{m}^{-4}]$';


ax.XLim = [0 maxtime];


fn = sprintf('%srms_dens_error_vs_cputime_%s_(%d,%d)', PATH_FIGS, DATA_SOURCE, ni, nj);
export_fig(hf, fn, '-pdf');
% close(hf);

return

%% Topobaric Surface Subtimers
hf = figure('Position', [1920-800, 1080-480, 800, 480]);
ax = axes('Position', [.1, .11, .77, .8]); hold on; grid on; box on;
fn = @log2;
%fn = @(z) z;

subtimer = nan(nDATA, 5);
iS = iSURF('TOPOB');
for iDATA = 1 : nDATA
  d = diags{iDATA,iS};
  subtimer(iDATA,:) = sum([d.timer_wetting, d.timer_recon, d.timer_reebgraph, d.timer_fitting, d.timer_update], 1);
end
for i = 1 : size(subtimer,2)
  plot(ax, fn(RES), fn(subtimer(:,i)), '-o', 'Color', cm(i,:), 'MarkerFace', cm(i,:));
end
legend({'BFS Wetting & Conncomp', 'Recon', 'Reeb Graph', 'Fitting & Graph Int', 'Vert Update'}, 'Location', 'northwest');

% add slopes
yl = ax.YLim;
slopes = [2, 3, 3, 2, 2];
for i = 1 : size(subtimer,2)
  slope = slopes(i);
  xx = log2(RES([1 nRES]));
  y2 = log2(subtimer(nRES,i));
  plot(xx, [y2 - slope*diff(xx),  y2], '--', 'Color', cm(i,:), 'DisplayName', sprintf('slope %.1f', slope));  % constant slope line
end
ax.YLim = yl;

ax.YTick = ax.YTick(1) : ax.YTick(end);
ax.XTick = log2(RES(1:nRES));
ax.XLabel.String = 'log_2(Grid Resolution)';
ax.YLabel.String = 'log_2(CPU time / [s])';
ax.YLabel.Interpreter = 'tex';


% Add second axes that show non-logged values
ax2 = axes('Position', [.1, .11, .77, .8]);
ax2.Color = 'none';
ax2.XTick = [];
ax2.XAxisLocation = 'top';
ax2.XTick = ax.XTick;
ax2.XLim = ax.XLim;
ax2.XTickLabel = arrayfun(@(c) sprintf('%d', 2.^c), ax.XTick, 'UniformOutput', false);
ax2.XLabel.String = 'Grid Resolution';
ax2.XLabel.Interpreter = 'tex';

ax2.YAxisLocation = 'right';
ax2.YTick = ax.YTick;
ax2.YLim = ax.YLim;
ax2.YTickLabel = arrayfun(@(c) sprintf('%.3g', 2.^c), ax.YTick, 'UniformOutput', false);
ax2.YLabel.String = 'CPU time [s]';
ax2.YLabel.Interpreter = 'tex';

fn = sprintf('cpu-vs-gridres__log2__topobaric-surface-subtimers__%s', datestr(now, 'yy-mm-dd hh-MM-ss'));
%export_fig(hf, [PATH_FIGS fn], '-pdf')
%close(hf);

%% Omega Surface Subtimers
hf = figure('Position', [1920-800, 1080-480, 800, 480]);
ax = axes('Position', [.1, .11, .77, .8]); hold on; grid on; box on;
fn = @log2;
%fn = @(z) z;
fn = @(x) log2(x .* log2(x));

subtimer = nan(nDATA, 3);
iS = iSURF('OMEGA');
for iDATA = 1 : nDATA
  d = diags{iDATA,iS};
  subtimer(iDATA,1) = sum(d.timer_solver);
  subtimer(iDATA,2) = sum(d.timer_update);
  subtimer(iDATA,3) = sum(d.clocktime) - sum(d.timer_solver) - sum(d.timer_update);
end

for i = 1 : size(subtimer,2)
  plot(ax, fn(NSP(:,iS)), fn(subtimer(:,i)), '-o', 'Color', cm(i,:), 'MarkerFace', cm(i,:));
end
legend({'Matrix Solution', 'Vert Update', 'Other'}, 'Location', 'northwest');

% add slopes
yl = ax.YLim;
slopes = [1.18, 1, 1.1, 1];
for i = 1 : size(subtimer,2)
  slope = slopes(i);
  xx = log2(NSP([1 nRES],iS));
  y2 = log2(subtimer(nRES,i));
  plot(xx, [y2 - slope*diff(xx),  y2], '--', 'Color', cm(i,:), 'DisplayName', sprintf('slope %.2f', slope));  % constant slope line
end
ax.YLim = yl;

ax.YTick = ax.YTick(1) : ax.YTick(end);
ax.XTick = log2(NSP(1:nRES,iS).^2);
ax.XLabel.String = 'log_2(# Ocean Casts)';
ax.YLabel.String = 'log_2(CPU time / 1s)';
ax.YLabel.Interpreter = 'tex';

% Add second axes that show non-logged values
ax2 = axes('Position', [.1, .11, .77, .8]);
ax2.Color = 'none';
ax2.XTick = [];
ax2.XAxisLocation = 'top';
ax2.XTick = ax.XTick;
ax2.XLim = ax.XLim;
ax2.XTickLabel = arrayfun(@(c) sprintf('%d^2', 2.^(c/2)), ax.XTick, 'UniformOutput', false);
ax2.XLabel.String = '# Ocean Casts';
ax2.XLabel.Interpreter = 'tex';

ax2.YAxisLocation = 'right';
ax2.YTick = ax.YTick;
ax2.YLim = ax.YLim;
for i = 1 : length(ax2.YTickLabel)
  if 2 ^ ax.YTick(i) >= 1
    ax2.YTickLabel{i} = [sprintf('%.1d', 2.^ax.YTick(i)) ' s'];
  else
    ax2.YTickLabel{i} = [sprintf('%.1f', 100 * 2.^ax.YTick(i)) ' ms'];
  end
end
ax2.YLabel.String = 'CPU time [s]';
ax2.YLabel.Interpreter = 'tex';

fn = sprintf('cpu-vs-gridres__log2__omega-surface-subtimers__%s', datestr(now, 'yy-mm-dd hh-MM-ss'));
%export_fig(hf, [PATH_FIGS fn], '-pdf')
%close(hf);
