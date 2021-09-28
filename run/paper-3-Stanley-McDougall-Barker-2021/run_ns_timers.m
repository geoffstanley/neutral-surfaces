% Comparison of the various surfaces in the Neutral-Surfaces toolbox

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


%% --- BEGIN SETUP --------------------------------------------------------

warning('off', 'MATLAB:nargchk:deprecated')
set(0, 'defaultfigurecolor', [1 1 1]); % white figure background

% Add Neutral Surfaces to MATLAB's path
PATH_NS = '~/work/projects-gfd/neutral-surfaces/'; % << ADJUST AS NEEDED >>
run([PATH_NS 'ns_add_to_path.m']);

% Note: The non-Boussinesq functions densjmd95 and densjmd95_s_t are
% necessary for gamma_n neutrality error calculations and are included in
% PATH_LOCAL alongside this file. 

% Add gamma_n() and nsfcs() MATLAB functions for Jackett and McDougall (1997) neutral density software.
% Obtain from: http://www.teos-10.org/preteos10_software/neutral_density.html
PATH_GAMMA_N = '~/work/MATLAB_notloaded/oceanography_tools/neutral_density/gamma1/'; % << ADJUST AS NEEDED >>
addpath(PATH_GAMMA_N);

%PATH_LOCAL = [fileparts(mfilename('fullpath')) V]; % Get path to this file.
PATH_LOCAL = [PATH_NS 'run/paper-3-Stanley-McDougall-Barker-2021/']; % Manually set path to this file.
cd(PATH_LOCAL);

PATH_OMEGA_KMJ = [PATH_LOCAL 'omega-kmj/'];
addpath(genpath(PATH_OMEGA_KMJ));

% Folder containing the functions eos.m, eos_p.m, and eos_s_t.m
PATH_EOS = '~/work/MATLAB/eos/eos/'; % << ADJUST AS NEEDED >>
addpath(PATH_EOS);  % Doing this last to ensure at top of MATLAB's path, above eos fcns in other dirs

% Folder containing .nc data files
PATH_OCCA = '~/work/data/OCCA/';
PATH_ECCO2 = '~/work/data/ECCO2/latlon/';

% Make a directory for figures and other output
PATH_OUT = [PATH_LOCAL 'output/'];
PATH_OUT = [PATH_OUT datestr(now, 'dd-mm-yyyy') '/'];
if ~exist(PATH_OUT, 'dir'); mkdir(PATH_OUT); end

fileID = 1; % For standard output to the screen
%fileID = fopen([PATH_OUT 'run ' datestr(now, 'yyyy mm dd hh MM ss') '.txt'], 'wt'); % For output to a file

maxNumCompThreads(1);   % Just use one CPU, for apples-to-apples comparison

REPEAT = true; % for paper
% REPEAT = false; nrep = 1; % for quicker tests

%  RES = 2.^(4:7); MODELS = {};  % for quicker tests
% RES = 2.^(4:8); MODELS = {'OCCA'};  % for quicker tests
RES = 2.^(4:11); MODELS = {'OCCA', 'ECCO2'};  % Everything. For paper

nSYNTH = length(RES);
nMODELS = length(MODELS);
nDATA = nSYNTH + nMODELS;
NWC = nan(nDATA,1); % Number of Water Columns

rho_c = 1027.5; % Boussinesq reference density [kg m-3]
grav = 9.81;    % gravitational acceleration [m s-2]

db2Pa = 1e4; % dbar to Pa conversion
Pa2db = 1e-4; % Pa to dbar conversion

lead1 = @(x) reshape(x, [1 size(x)]);

% Choose vertical interpolation method
INTERPFN = @ppc_linterp;


% Selecting surfaces and depths
list_surf = {'\sigma', '\delta', '\gamma^n', '\omega_{KMJ}', '\omega_{\nabla}', '\omega_+', '\tau'};
LIST_SURF = {'SIGMA', 'DELTA', 'GAMMA', 'KMJ', 'OMEGAGRAD', 'OMEGA', 'TOPOB'};

nS = length(list_surf);
iSURF = containers.Map(LIST_SURF, 1:nS);

% Save depths and diagnostic outputs from various surfaces
zs    = cell(nDATA, nS);
diags = cell(nDATA, nS);
timer = nan( nDATA, nS);


%% Set alias functions.  << ENSURE THIS GETS DONE >>
% Choose the Boussinesq densjmd95 and set grav and rho_c in eos.m and eos_p.m
if false % only need to do this once!  Doing every time makes codegen run every time
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq.m'   ] , [PATH_EOS 'eos.m'  ], grav, rho_c); %#ok<UNRCH>
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq_dz.m'] , [PATH_EOS 'eos_p.m'], grav, rho_c);
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq_s_t.m'], [PATH_EOS 'eos_s_t.m'], grav, rho_c);
end
clear eos eos_p eos_s_t % Make sure the copied files get used

%% Begin looping over data
for iDATA = 1 : nDATA
  %% Load grid and hydrographic data
  fprintf(1, '--- Starting grid # %d\n', iDATA);
  
  OPTS = struct();
  
  if iDATA <= nSYNTH
    
    DATA_SOURCE = 'SYNTHRAND';
    WRAP = [false; false]; % non-periodic in both dimensions
    
    
    ni = RES(iDATA);
    nj = ni;
    nk = 50;
    
    [S, T, Z, g] = synthocean_rand(ni, nj, nk, WRAP);
    g.DXC = repmat(g.DXCvec, ni, 1);
    g.DYC = repmat(g.DYCsc, ni, nj);
    g.DXG = repmat(g.DXGvec, ni, 1);
    g.DYG = repmat(g.DYGsc, ni, nj);
        
    z0 = 1500;
    
    i0 = ceil(ni/2); % middle i
    j0 = ceil(nj/2); % middle j
    
  else
    
    DATA_SOURCE = MODELS{iDATA - nSYNTH};
    WRAP = [true; false]; % periodic in longitude, non-periodic in latitude
    
    if DATA_SOURCE(1) == 'O'
      [g, S, T, P, ETAN, ATMP] = load_OCCA(PATH_OCCA);
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
    S = permute(S, [3 1 2]); % [nz,nx,ny]. depth  by  long  by  lat
    T = permute(T, [3 1 2]);
    if exist('P', 'var')
      P = permute(P, [3 1 2]);
    end
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
    x0 = g.XCvec(i0);
    y0 = g.YCvec(j0);
    
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
  
  Sppc = INTERPFN(Z, S);
  Tppc = INTERPFN(Z, T);
  
  Z2P = Pa2db * rho_c * grav ; % Note > 0
  
  Szyx = permute(S, [1 3 2]); % KMJ omega needs dimensions to be [z,y,x].
  Tzyx = permute(T, [1 3 2]);
  Zzyx = repmat(Z(:), [1, nj, ni]);
  
  NWC(iDATA) = sum(flat(isfinite(S(1,:,:))));
  
  
  % --- Setup options neutral surfaces
  OPTS.TOL_P_CHANGE_L2 = 1e-3; % Stop iterations when the L2 change of pressure on the surface is less than this
  OPTS.ITER_MIN = 5;      % .. but do at least 5 iterations
  OPTS.ITER_MAX = 30;     % The maximum number of iterations
  %OPTS.ITER_START_WETTING = 1; % Start wetting on first iteration
  OPTS.ITER_START_WETTING = inf; % No Wetting
  OPTS.REEB = true;       % Do topobaric, not orthobaric, surfaces
  OPTS.GEOSTRF = false;  % Don't add extra criteria that make geostrophic streamfunction well-defined
  OPTS.SIMPLIFY_ARC_REMAIN = Inf; % No simplification
  OPTS.FILL_IJ = [];     % No filling
  OPTS.FILL_PIX = 0;     % No filling.
  OPTS.ML = [];         % Leave the Mixed Layer in.
  OPTS.TOL_P_UPDATE = 1e-4; % error tolerance when updating the surface
  OPTS.P_EXPN = 500;     % expansion of domain to search for solutions in each water column
  OPTS.TOL_LRPD_L1 = 0; % Don't use the LRPD stopping criterion.
  OPTS.VERBOSE = 1;      % Some output while executing functions
  OPTS.FILE_ID = fileID; % Write output to this file
  OPTS.FIGS_SHOW = false; % Don't show figures (for omega surface)
  OPTS.INTERPFN = INTERPFN; % set interpolation function to match this script
  OPTS.Sppc = Sppc;      % Use pre-computed piecewise polynomial for S in terms of X=Z
  OPTS.Tppc = Tppc;      % Use pre-computed piecewise polynomial for T in terms of X=Z
  
  OPTS.DIST1_iJ = g.DXC;   % Distance [m] in 1st dimension centred at (I-1/2, J)
  OPTS.DIST2_Ij = g.DYC;   % Distance [m] in 2nd dimension centred at (I, J-1/2)
  OPTS.DIST2_iJ = g.DYG;   % Distance [m] in 2nd dimension centred at (I-1/2, J)
  OPTS.DIST1_Ij = g.DXG;   % Distance [m] in 1st dimension centred at (I, J-1/2)
  
  % Testing sensitivity to non-uniform grid distances
  %   OPTS.DIST1_iJ = 1;   % Distance [m] in 1st dimension centred at (I-1/2, J)
  %   OPTS.DIST2_Ij = 1;   % Distance [m] in 2nd dimension centred at (I, J-1/2)
  %   OPTS.DIST2_iJ = 1;   % Distance [m] in 2nd dimension centred at (I-1/2, J)
  %   OPTS.DIST1_Ij = 1;   % Distance [m] in 1st dimension centred at (I, J-1/2)
    
  % Call eos() functions to do Just-in-Time Compilation, so that it's not done in one or another function later
  s = squeeze(S(1,:,:));
  t = squeeze(T(1,:,:));
  z = Z(1);
  eos(s, t, z);
  eos_p(s, t, z);
  eos_s_t(s, t, z);
  
  fprintf(1, 'Loaded %s (%d,%d)\n', DATA_SOURCE, ni, nj);
  
  if iDATA == nDATA
    g_ECCO = g;
  end
  
  
  %% Labels for the surface:
  % (a) the Veronis Density, by integrating down a reference cast at (188,-4)
  % (b) from 1997 neutral density value of the bottle at (188,-4) and z0 depth.
  if iDATA > nSYNTH  % either OCCA or ECCO2
    
    i1 = x_to_i(188);  % index closest to 188 deg E
    j1 = y_to_j(-4);   % index closest to  4 deg S
    x1 = g.XCvec(i1);  % longitude at i1
    y1 = g.YCvec(j1);  % latitude at j1
    z1 = z0;           % alias.
    
    % Prepare pressure on cast at (i1,j1)
    if strcmp(DATA_SOURCE, 'OCCA')
      % use model output pressure variable
      P1 = P(:,i1,j1);
    elseif strcmp(DATA_SOURCE, 'ECCO2')
      % trapezoidal integration of hydrostatic balance to get in-situ pressure
      R1 = eos(S(:,i1,j1), T(:,i1,j1), Z);
      P1 = Pa2db * grav * (Z(1) * R1(1) + cumtrapz(Z, R1));
    end
    
    % Calculate Jackett and McDougall (1997) gamma^n neutral density
    s1 = ppc_linterp(Z, S(:,i1,j1), z1); % Evaluate (S,T,P) at (i1,j1,z1)
    t1 = ppc_linterp(Z, T(:,i1,j1), z1);
    p1 = ppc_linterp(Z, P1, z1);
    t1 = eos80_legacy_pt(s1,t1,0,p1);  % convert t1 from potential temperature to in-situ temperature at pressure p1.  Note eos80_legacy_pt is same as sw_ptmp, to machine precision.
    label_gn = eos80_legacy_gamma_n(s1, t1, p1, x1, y1);
    fprintf(fileID, '%s gamma^n label at (188E,4S,%dm): %.16f\n', DATA_SOURCE, z1, label_gn);
    
    % Caclulate Veronis Density. Only internally consistent with this dataset
    z_ref = 0;  % reference depth for potential density
    label_veronis1 = veronis_density(z_ref, S(:,i1,j1), T(:,i1,j1), Z, Z(1), z1);
    fprintf(fileID, '%s Veronis Density at (%.4f,%.4f,%dm): %.16f\n', DATA_SOURCE, x1, y1, z1, label_veronis1);
    
    label_veronis0 = veronis_density(z_ref, S(:,i0,j0), T(:,i0,j0), Z, Z(1), z1);
    fprintf(fileID, '%s Veronis Density at (%.4f,%.4f,%dm): %.16f\n', DATA_SOURCE, x0, y0, z0, label_veronis0);
    
    % For Veronis Density, if using a 4D (with time) ocean dataset, then a reference cast must be selected
    % not just in space but also in time.  A Boussinesq example is as follows.
    % Suppose S and T are 4D arrays, ordered [depth, lon, lat, time].
    % Suppose l0 is the time index for the current time, in which the surface lives
    % Suppose l1 is the time index for the reference cast
    %{
    Sc = squeeze(S(:, i1, j1, l : -1 : l1)); % S casts from current time to reference time
    Tc = squeeze(T(:, i1, j1, l : -1 : l1)); % T casts from current time to reference time
    [z1,s1,t1] = neutral_trajectory(Sc, Tc, Z, z0); % temporal neutral trajectory from current time at z0 to reference time
    s1 = s1(end); % select just the values at the end, the reference time
    t1 = t1(end);
    z1 = z1(end);
    z_ref = 0;  % reference depth for potential density
    label_veronis = veronis_density(z_ref, S(:,i1,j1,l1), T(:,i1,j1,l1), Z, Z(1), z1); % Evaluate Veronis density on cast at the reference time
    fprintf(fileID, '%s Veronis Density at (%.4f, %.4f, %dm, time index %d): %.16f\n', DATA_SOURCE, x1, y1, z1, l1, label_veronis);
    %}
  end
  
  %% gamma^n Neutral Density surface
  
  % Step 1: Convert the model's potential temperature to in-situ temperature,
  % Step 2: Calculate 3D gamma^n
  % Step 3: Calculate gamma^n surface.
  
  % Note: sw_ptmp was updated to now output in ITS-90.
  %   Units for sw_ptmp are as follows:
  %   S  = salinity    [psu      (PSS-78) ]
  %   T  = temperature [degree C (ITS-90)]
  %   P  = pressure    [db]
  %   PR = Reference pressure  [db]
  % OUTPUT:
  %   ptmp = Potential temperature relative to PR [degree C (ITS-90)]
  % THETAj [ITS-90] * 1.00024 -> THETAj [IPTS-68]
  
  % gamma_n call:
  %%%	UNITS:		salinity	psu (IPSS-78)
  %%%			temperature	degrees C (IPTS-68)
  %%%			pressure	db
  %%%			gamma_n		kg m-3
  
  if iDATA > nSYNTH  % either OCCA or ECCO2
    
    d = struct();
    
    % Prepare pressure on cast at (i1,j1)
    if ~strcmp(DATA_SOURCE, 'OCCA')
      % trapezoidal integration of hydrostatic balance to get in-situ pressure
      R = eos(S, T, Z);
      P = (Pa2db * grav) * (Z(1) * R(1,:,:) + cumtrapz(Z, R));
    end
    
    % Evaluate S,T,P at (i0, j0, z0)
    s0 = ppc_linterp(Z, S(:,i0,j0), z0);
    t0 = ppc_linterp(Z, T(:,i0,j0), z0);
    p0 = ppc_linterp(Z, P(:,i0,j0), z0);
    t0 = eos80_legacy_pt(s0, t0, 0, p0) * 1.00024 ; % [deg C IPTS-68]
    
    % Evaluate gamma_n at (i0, j0, p0)
    g0 = gamma_n(s0, t0, p0, g.XCvec(i0), g.YCvec(j0));
    
    mytic = tic;
    s_gamma = nan(ni,nj);
    t_gamma = nan(ni,nj);
    p_gamma = nan(ni,nj);
    for j = 1 : find(g.YCvec <= 64, 1, 'last') % prevent going beyond 64N, or gamma_n() errors
      Sj = S(:,:,j);
      Tj = T(:,:,j);
      Pj = P(:,:,j);
      Tj = eos80_legacy_pt(Sj, Tj, zeros(size(Sj)), Pj) * 1.00024 ; % Convert to in-situ temperature. [deg C IPTS-68]
      XCj = g.XCvec(:);
      YCj = repmat(g.YCvec(j), ni, 1);
      if all(isnan(Sj(:))) % otherwise gamma_n() errors
        continue
      end
      Gj = gamma_n_gjs(Sj, Tj, Pj, XCj, YCj);
      
      % Set any -99 gamma_n values to NaN.  Normally not necessary, but sometimes nsfcs() crashes without this.
      bad = (Gj < 0);
      Sj(bad) = nan;
      Tj(bad) = nan;
      Pj(bad) = nan;
      Gj(bad) = nan;
      
      [s_gamma(:,j), t_gamma(:,j), p_gamma(:,j)] = nsfcs(Sj, Tj, Pj, Gj, g0);
      fprintf(1, 'j=%03d done, time %.2fsec\n', j, toc(mytic));
    end
    s_gamma(s_gamma == -99) = nan;
    t_gamma(t_gamma == -99) = nan;
    p_gamma(p_gamma == -99) = nan;
    t_gamma = eos80_legacy_pt(s_gamma, t_gamma, p_gamma, zeros(size(s_gamma))) / 1.00024 ; % Convert to potential temperature. [deg C ITS-90]
    d.clocktime = toc(mytic);
    
    % Calculate neutrality errors using non-Boussinesq equation of state, and pressure as input. 
    [d.epsL2, d.epsL1] = eps_norms(s_gamma, t_gamma, p_gamma, true, [1 0], {@densjmd95, @densjmd95_s_t}, OPTS.DIST1_iJ, OPTS.DIST2_Ij, OPTS.DIST2_iJ, OPTS.DIST1_Ij);
    
    z_gamma = ppc_linterp(P, Z, lead1(p_gamma)); % Interpolate Pressure to Depth. 
    
    zs{iDATA, iSURF('GAMMA')} = z_gamma;
    diags{iDATA, iSURF('GAMMA')} = d;
    timer(iDATA, iSURF('GAMMA')) = d.clocktime;
  end
  
  %% Potential Density Surface
  if 1
    zref = z0;
    
    mytime = 0;
    if REPEAT
      nrep = max(1, ceil(log2(1024^2 ./ NWC(iDATA)))); % so N x N data gets done once where N=1024.  Others get about same CPU time, if they scale as N^1.
    end
    for rep = 1 : nrep
      [z, s, t, isoval, d] = pot_dens_surf(S, T, Z, zref, [i0, j0, z0], WRAP, OPTS);
      
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
    
    [z, s, t, isoval, s0, t0, d] = delta_surf(S, T, Z, [], [], [i0, j0, z0], WRAP, OPTS);
    
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
  
  
  %% KMJ's Omega surface
  
  if 1 && RES(iDATA) < 256  % Do up to 128^2 and do OCCA.
    if REPEAT
      nrep = max(1, ceil(log2(128^2 ./ NWC(iDATA)))); % so N x N data gets done once where N=128.  Others get about same CPU time, if they scale as N^1.
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
      
      if iDATA <= nSYNTH
        nit = 200; % max iterations
      else
        nit = 50;  % premature stop to KMJ algorithm, about as many iters as it actually takes on SYNTHETIC data
      end
      
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
      [z, ~, ~, d] = omega_surface(S, T, Z, z_delta, I0, WRAP, OPTS);
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
      [z, ~, ~, d] = omega_surface(S, T, Z, z_delta, I0, WRAP, OPTS);
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
      [z, ~, ~, RG, ~, ~, ~, d] = topobaric_surface(S, T, Z, z_delta, I0, WRAP, OPTS);
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
  
end % iDATA

%%
% Extract # of Surface Points per surface
NSP = nan(nDATA, nS);
for iDATA = 1 : nDATA
  for iS = 1 : nS
    NSP(iDATA, iS) = sum(isfinite(flat(zs{iDATA, iS})));
  end
end

%% Save data

% fn = sprintf('results ZS %s - %s.mat', datestr(now, 'yy-mm-dd hh-MM-ss'), num2str(RES-1));
% save(fn, 'zs');
fn = sprintf('results DIAGS %s - %s.mat', datestr(now, 'yy-mm-dd hh-MM-ss'), num2str(RES-1));
save(fn, 'diags', 'timer', 'RES', 'nSYNTH', 'nDATA', 'nS', 'iSURF', 'NSP', 'NWC', 'list_surf', 'MODELS');
return


%% Diagnostics plot: CPU time vs. Resolution

cm = [ ... % Colormap for lines
  0.0       0.0       0.0 % black
  0.85      0.6       0.0 % orange
  0.15      0.85      0.9 % cyan
  0.0       0.0       0.9 % blue
  0.5       0.2       0.7 % purple
  1.0       0.0       0.0 % red
  0.0       0.5       0.0 % dark green
  ]; %#ok<UNRCH>
markers = 'ooooooo';

hf = figure('Position', [1920-800, 1080-700, 800, 600]);
ax = axes('Position', [.1, .11, .74, .8]);
hold on; grid on; box on;
fcn = @log2;

%markersize = [8, 11, 8, 8, 8, 11, 8];
markersize = repmat(6, nS, 1);

surf_order = [ iSURF('SIGMA'), iSURF('DELTA'), iSURF('GAMMA'), iSURF('KMJ'), iSURF('OMEGAGRAD'), iSURF('TOPOB'), iSURF('OMEGA')];  % to make small omega + markers on top of big topobaric diamond markers


% Make markers outside the bounds of the plot, just to set of the first legend
for iS = surf_order
  plot(ax, 0, 0, 's', 'color', cm(iS,:), 'MarkerFace', cm(iS,:), 'MarkerSize', 8);
end


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

% Plot lines
for iS = surf_order
  plot(ax, fcn(NSP(1:nSYNTH,iS)), fcn(timer(1:nSYNTH,iS)), ['-' markers(iS)], 'color', cm(iS,:), 'MarkerFace', cm(iS,:), 'MarkerSize', markersize(iS), 'LineWidth', 1.5);
end

% Add lines for n^slope
slope_color = [1 1 1] * .7;
x1 = repmat(fcn(max(NWC)), nS, 1); x1(iSURF('KMJ')) = 14;
slopes = nan(nS,1);
slopes(iSURF('SIGMA')) = 1;
slopes(iSURF('DELTA')) = 1;
slopes(iSURF('KMJ')) = 1.6;
slopes(iSURF('OMEGAGRAD')) = 1.6;
slopes(iSURF('OMEGA')) = 1.2;
for iS = [ iSURF('DELTA'), iSURF('KMJ'), iSURF('OMEGAGRAD'), iSURF('OMEGA')]
  i = find(isfinite(timer(1:nSYNTH,iS)), 1, 'last');
  if isempty(i); continue; end
  slope = slopes(iS);
  if iS == iSURF('KMJ')
    xx = [fcn(NSP(1,iS)), fcn(NSP(i,iS)) + .4];
  else
    xx = [fcn(NSP(5,iS)), fcn(NSP(i,iS)) + .4];
  end
  x2 = fcn(NSP(i,iS));
  y2 = fcn(timer(i,iS));
  ln = @(x) y2 + slope * (x - x2);
  plot(ax, xx, ln(xx), '--', 'Color', slope_color, 'LineWidth', 1 + .5 * (iS==iSURF('DELTA')));  % constant slope line, thicker for DELTA
  
  x = x1(iS) + 0.4;
  y = ln(x) - 0.3;
  txt = ['N^{' num2str(slope) '}'];
  text(ax, x, y, txt, 'Color', cm(iS,:));
end

% Add lines for n log n
for iS = iSURF('TOPOB')
  
  i = find(isfinite(timer(1:nSYNTH,iS)), 1, 'last');
  if isempty(i); continue; end
  
  xx = linspace(fcn(NSP(5,iS)), fcn(NSP(i,iS)) + .4, 100);
  x2 = fcn(NSP(i,iS));
  y2 = fcn(timer(i,iS));
  C = y2 - x2 - log2(x2);
  ln = @(x) C + x + log2(x);
  plot(ax, xx, ln(xx), '--', 'Color', slope_color);  % constant slope line
  
  txt = 'N logN';
  text(ax, x2 + .2, y2 + .2 , txt, 'Color', cm(iS,:));
end


ax.YLim = [-11, 12];
ax.XLim = [log2(min(NWC)) - 1, log2(max(NWC)) + 1.5];
ax.YTick = ax.YLim(1) : ax.YLim(2); % all integers, only
ax.XTick = log2(NWC(1:nSYNTH));
ax.XLabel.String = 'log_2({\itN})';
ax.YLabel.String = 'log_2(CPU time / 1s)';
ax.XLabel.Interpreter = 'tex';
ax.YLabel.Interpreter = 'tex';
leg = legend(ax, list_surf(surf_order), 'Position', [.115 .71 .12 .18]);

% Write # of iterations.  skip SIGMA and DELTA
for iDATA = 1 : nDATA
  for iS = [iSURF('KMJ'), iSURF('OMEGAGRAD'), iSURF('OMEGA'), iSURF('TOPOB')]
    x = fcn(NSP(iDATA,iS)) - 0.007 * diff(ax.XLim);
    if iS == iSURF('OMEGA')
      y = fcn(timer(iDATA,iS)) - 0.03 * diff(ax.YLim); % below
    else
      y = fcn(timer(iDATA,iS)) + 0.03 * diff(ax.YLim); % above
    end
    txt = num2str(iters(iDATA,iS));
    text(ax, x, y, txt, 'Color', cm(iS,:), 'HorizontalAlignment', 'right');
  end
end


% Add second axes that show non-logged values,  and MODEL data with separate legend
ax2 = axes('Position', [.1, .11, .74, .8]);
hold on;
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
  elseif 2 ^ ax.YTick(i) >= 1e-3
    ax2.YTickLabel{i} = [sprintf('%.1f', 1000 * 2.^ax.YTick(i)) ' ms'];
  end
end

ax2.YTickLabel = {...
  '488 μs', ...
  '977 μs', ...
  '1.95 ms', ...
  '3.9 ms', ...
  '7.8 ms', ...
  '15.6 ms', ...
  '31.3 ms', ...
  '62.5 ms', ...
  '125 ms', ...
  '250 ms', ...
  '500 ms', ...
  '   1 s', ...
  '   2 s', ...
  '   4 s', ...
  '   8 s', ...
  '  16 s', ...
  '  32 s', ...
  '  64 s', ...
  '2.1 min', ...
  '4.3 min', ...
  '8.5 min', ...
  '17 min', ...
  '34 min', ...
  '68 min', ...
  };

ax2.YLabel.String = '\approx CPU time';
ax2.YLabel.Interpreter = 'tex';

% Make markers outside the bounds of the plot, just to set up the second legend
model_marker = 'o^p';
for i = 1 : length(model_marker)
  plot(ax2, 0, 0, model_marker(i), 'color', 'k', 'MarkerFace', 'k', 'MarkerSize', 8);
end

% MODEL data
model_marker = '^p';
model_marker_area = [60, 100];
for iS = [ iSURF('SIGMA'), iSURF('DELTA'), iSURF('GAMMA'), iSURF('KMJ'), iSURF('OMEGAGRAD'), iSURF('TOPOB'), iSURF('OMEGA')]
  for iDATA = nSYNTH + 1 : nDATA
    iMODEL = iDATA - nSYNTH;
    scatter(ax2, fcn(NSP(iDATA,iS)), fcn(timer(iDATA,iS)), model_marker_area(iMODEL), cm(iS,:)*.8, model_marker(iMODEL), 'MarkerFaceColor', cm(iS,:)*.8);
  end
end

leg = legend(ax2, {'Aquaplanet', 'OCCA', 'ECCO2'}, 'Position', [.28, .824, .08, .07]);
leg.String = leg.String(1:3);

%%  Save Figure
fn = sprintf('cpu-vs-gridres__log2__%s', datestr(now, 'yy-mm-dd hh-MM-ss'));
export_fig(hf, [PATH_OUT fn], '-pdf')
savefig(hf, [PATH_OUT fn])
%close(hf);

%% Diagnostics plot:  eps L2 norm vs CPU time

cm = [ ... % Colormap for lines
  0.0       0.0       0.0 % black
  0.85      0.6       0.0 % orange
  0.15      0.85      0.9 % cyan
  0.0       0.0       0.9 % blue
  0.5       0.2       0.7 % purple
  1.0       0.0       0.0 % red
  0.0       0.5       0.0 % dark green
  ];
markers = 'ooooooo';
horzlinestyle = repmat({'--'}, nS, 1);
horzlinestyle{iSURF('GAMMA')} = '-';

hf = figure('Position', [0 0 800 800]);
fs = 12;
clear ax;
timefac = ones(nS,1);
timefacstr = repmat('sec', nS, 1);
markersize = repmat(6, nS, 1);

for i = 1 : 2
  
  %ax = axes('Position', [.1 .14 .95-.1 .8]);
  ax = subaxis(2,1,i, 'margin', 0.07, 'marginright', 0.03, 'marginleft', 0.1, 'spacing', 0.1);
  ax.YAxis.Scale = 'log';
  ax.FontSize = fs;
  hold(ax, 'on');
  grid(ax, 'on');
  box(ax, 'on');
  
  if i == 1
    iDATA = find(RES == 128+1);  DATA_SOURCE = 'SYNTHRAND';  % Select data source for this plot
    timefac(iSURF('GAMMA')) = 60; % measure GAMMA in minutes not seconds
    timefacstr(iSURF('GAMMA'),:) = 'min';
    timefac(iSURF('KMJ')) = 60; % measure KMJ in minutes not seconds
    timefacstr(iSURF('KMJ'),:) = 'min';
    yl = [1e-11 6e-9];
    xtick = 0 : 0.25 : 2.5;
    text(ax, -.05, 1.05, '(a)    Aquaplanet [128 x 128]', 'fontsize', fs+2, 'units', 'normalized');
  else
    iDATA = find(RES == 240);  DATA_SOURCE = 'OCCA';  % Select data source for this plot
    timefac(iSURF('GAMMA')) = 60; % measure GAMMA in minutes not seconds
    timefacstr(iSURF('GAMMA'),:) = 'min';
    timefac(iSURF('KMJ')) = 60;  % measure KMJ in minutes not seconds
    timefacstr(iSURF('KMJ'),:) = 'min';
    yl = [1e-10 5e-8];
    xtick = 0 : 1 : 15;
    text(ax, -.05, 1.05, '(b)    OCCA [360 x 160]', 'fontsize', fs+2, 'units', 'normalized');
  end
  
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
for iS = [iSURF('SIGMA'), iSURF('DELTA'), iSURF('GAMMA'), iSURF('KMJ'), iSURF('TOPOB'), iSURF('OMEGAGRAD'), iSURF('OMEGA')]
  d = diags{iDATA, iS};
  if ~isempty(d)
    j = 1; % for SIGMA, DELTA, GAMMA
    if iS == iSURF('KMJ')
      j = length(d.epsL2);
    elseif iS == iSURF('OMEGAGRAD') || iS == iSURF('OMEGA') || iS == iSURF('TOPOB')
      % Find iteration at which algorithm actually deemed converged, which may be less than OPTS.ITER_MIN
      % +1 epsL2 has one extra element, giving epsL2 on the initial surface (corresponding with clocktime(1) == 0). 
      j = 1 + find(d.x_change_L2 <= OPTS.TOL_P_CHANGE_L2, 1, 'first');
    end
    
    semilogy(ax, clocktime{iS}(1:j), d.epsL2(1:j), '-', 'Color', cm(iS,:), 'LineWidth', 2, 'Marker', markers(iS), 'MarkerSize', markersize(iS), 'MarkerFace', cm(iS,:));
    semilogy(ax, [0 maxtime], [1 1] * d.epsL2(j), horzlinestyle{iS}, 'Color', cm(iS,:), 'LineWidth', 1);
    epsL2(iDATA, iS) = d.epsL2(j);
    
    ha = 'right';
    vertalign = 'top';
    x = maxtime * .99;
    y = 10^(log10(epsL2(iDATA,iS)) + 0.003 * log10(diff(ax.YLim)));
    txt = sprintf('$%s$: %.2e kg m$^{-4}$, %.2f %s', list_surf{iS}, epsL2(iDATA,iS), clocktime{iS}(j), timefacstr(iS,:));
    
    if iS == iSURF('GAMMA')
      x = clocktime{iS}(j);
      ha = 'left';
      y = 10^(log10(epsL2(iDATA,iS)) + 0.005 * log10(diff(ax.YLim)));
    elseif iS == iSURF('KMJ')
      y = 10^(log10(epsL2(iDATA,iS)) + 0.005 * log10(diff(ax.YLim)));
    elseif iS == iSURF('OMEGAGRAD')
      y = 10^(log10(epsL2(iDATA,iS)) - 0.008 * log10(diff(ax.YLim)));
      vertalign = 'bottom';
      if i == 1
        x = clocktime{iS}(end) / timefac(iS);
        ha = 'left';
      end
    elseif iS == iSURF('OMEGA')
      y = 10^(log10(epsL2(iDATA,iS)) + 0.01 * log10(diff(ax.YLim)));
      if i == 1
        x = clocktime{iS}(end) / timefac(iS);
        ha = 'left';
      end
    end
    text(ax, x, y, txt, 'Interpreter', 'latex', 'FontSize', fs, 'Color', cm(iS,:), 'VerticalAlignment', vertalign, 'HorizontalAlignment', ha);
    
  end
end

ax.YLim = yl;

ax.XTick = xtick;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.String = 'CPU time [sec or min]';
ax.YLabel.String = '$|\!| \epsilon |\!|_2 \quad [\mathrm{kg} \ \mathrm{m}^{-4}]$';

ax.XLim = [0 maxtime];

end

%% Save figure
fn = sprintf('rms_dens_error_vs_cputime_SYNTHRAND128x128__OCCA360x160_%s', datestr(now, 'yy-mm-dd hh-MM-ss'));
export_fig(hf, [PATH_OUT fn], '-pdf');
savefig(hf, [PATH_OUT fn])
% close(hf);

%% Maps showing depth of omega surface, depth difference with other surfaces, and gamma^n on on the omega surface

iECCO = nDATA;
assert(g.nx == ni && ni == 1440 && all(size(S) == [50, 1440, 720]), 'The currently loaded data is not for ECCO2');
land = squeeze(isnan(S(1,:,:)));

% --- Compute gamma^n on the omega surface
DO_GAMMAN_ON_OMEGA = false;
if DO_GAMMAN_ON_OMEGA
  
  % Prepare pressure on casts
  % trapezoidal integration of hydrostatic balance to get in-situ pressure
  R = eos(S, T, Z);
  P = (Pa2db * grav) * (Z(1) * R(1,:,:) + cumtrapz(Z, R));
  
  % evaluate S, T, P on the omega surface
  z_omega = zs{iECCO, iSURF('OMEGA')};
  s = ppc_linterp(Z, S, lead1(z_omega));
  t = ppc_linterp(Z, T, lead1(z_omega));
  p = ppc_linterp(Z, P, lead1(z_omega));
  
  % convert to in-situ temperature
  t = eos80_legacy_pt(s, t, 0, p) * 1.00024 ; % Convert to in-situ temperature. [deg C IPTS-68]
  
  x = repmat(g.XCvec, 1, nj);
  y = repmat(g.YCvec, ni, 1);
  good = (y >= -80 & y <= 64) & isfinite(s);
  gamma_on_omega = nan(ni,nj);
  gamma_on_omega(good) = gamma_n_gjs(s(good).', t(good).', p(good).', x(good).', y(good).');  % 35 seconds
end

% --- Make figure: 
fs = 12;

hf = figure('Position', [0 0 800 800]);

margin_v = .05;
margin_l = .07;
margin_r = .03;
spacing_v = 0.06;
spacing_h = 0.02;
subaxopts = {'MarginLeft',margin_l,'MarginRight', margin_r, 'MarginTop', margin_v, 'MarginBottom', margin_v, ...
  'SpacingVertical', spacing_v, 'SpacingHorizontal', spacing_h};

nR = 3;
nC = 2;

for iAX = 1 : (5 + DO_GAMMAN_ON_OMEGA)
  
  ax = subaxis(nR,nC,iAX,subaxopts{:});
  
  if iAX == 1
    dat = zs{iECCO,iSURF('OMEGA')};
    txt1 = 'z[\omega \! = \! \rho_V]';
    clim = [0, max(dat(:))];
    cm = flipud(parula(128));
    
  elseif iAX == 2
    dat = zs{iECCO,iSURF('SIGMA')} - zs{iECCO,iSURF('OMEGA')};
    txt1 = 'z[\sigma \! = \! \rho_V] - z[\omega \! = \! \rho_V]';
    clim = [-200, 200] ;
    cm = bluewhitered(128, clim);
    
  elseif iAX == 3
    dat = zs{iECCO,iSURF('DELTA')} - zs{iECCO,iSURF('OMEGA')};
    txt1 = 'z[\delta \! = \! \rho_V] - z[\omega \! = \! \rho_V]';
    clim = [-300, 100] ;
    cm = bluewhitered(128, clim);
    
  elseif iAX == 4
    dat = zs{iECCO,iSURF('TOPOB')} - zs{iECCO,iSURF('OMEGA')};
    txt1 = 'z[\tau \! = \! \rho_V] - z[\omega \! = \! \rho_V]';
    clim = [-1 1]*10 ;
    cm = bluewhitered(128, clim);
    
  elseif iAX == 5
    dat = zs{iECCO,iSURF('GAMMA')} - zs{iECCO,iSURF('OMEGA')};
    txt1 = 'z[\gamma \! = \! \rho_V] - z[\omega \! = \! \rho_V]';
    clim = [-80, 80] ;
    cm = bluewhitered(128, clim);
    
  elseif iAX == 6 && DO_GAMMAN_ON_OMEGA
    dat = gamma_on_omega;
    txt1 = '\gamma^n [\omega \! = \! \rho_V]';
    clim = [min(dat(:)), max(dat(:))];
    cm = parula(128);
    
  end
  
  fig_map(ax,g.XCvec, g.YCvec, dat, land);
  if diff(clim) > 0
    ax.CLim = clim;
  end
  hcb = colorbar(ax);
  colormap(ax,cm);
  
  if iAX <= nR*nC - 2
    ax.XTickLabel = {};
  end
  if mod(iAX,2) == 0
    ax.YTickLabel = {};
  end
  
  
  text(ax, 0, 1.06, ...
    sprintf('$(%s) \\; %s$', char('a' + iAX-1), txt1), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', fs, 'Interpreter', 'latex');
  
  if iAX == 1
    
    text(ax, 1.1, 1.06, ...
      sprintf('$\\rho_V = %.4f$ kg m$^{-3}$', label_veronis0), ...
      'Units', 'norm', 'HorizontalAlignment', 'right', 'FontSize', fs, 'Interpreter', 'latex');
    
  elseif iAX >= 2 && iAX <= 5
    
    datL2 = sqrt(nanmean(dat(:) .* dat(:)));
    
    text(ax, 1.1, 1.06, ...
      sprintf('$||\\cdot||_2=$ %.1f m', datL2), ...
      'Units', 'norm', 'HorizontalAlignment', 'right', 'FontSize', fs, 'Interpreter', 'latex');
    
  elseif iAX == 6 && DO_GAMMAN_ON_OMEGA
    
    resid = dat(:) - nanmean(dat(:));
    dat_std = sqrt(nanmean(resid .* resid));
    
    text(ax, 1.1, 1.06, ...
      sprintf('SD($\\cdot$)= %.2f g m$^{-3}$', 1000 * dat_std), ...
      'Units', 'norm', 'HorizontalAlignment', 'right', 'FontSize', fs, 'Interpreter', 'latex');
    
    hcb.Ticks = 27.78 : 0.01 : 27.83;
    
  end
  
  
end

%% Save figure
fn = sprintf('omega_surface_comparisons_ECCO2_(%d,%d,%.4fm)_%s', i0, j0, z0, datestr(now, 'yy-mm-dd hh-MM-ss'));
export_fig(hf, [PATH_OUT fn], '-jpg', '-m3');
close(hf);




