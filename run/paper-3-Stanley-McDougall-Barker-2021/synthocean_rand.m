function [S, T, P, g] = synthocean_rand(ni, nj, nk, wrap, pbot, std)
% SYNTHOCAEN_RAND   Synthetic ocean of eddies
%
%
% [S, T, P, g] = synthocean(ni, nj, nk)
% constructs a synthetic ocean of practical / Absolute salinity S,
% potential / Conservative temperature T, pressure P, having ni, nj, nk
% data points in the longitude, latitude depth directions, respectively.
% The grid data is contained in the output struct g.  
%
% [S, T, P, g] = synthocean(ni, nj, nk, wrap, pbot, std)
% also provides additional options determining the periodic nature of the domain, the bottom pressure, 
% and the scale of smoothing.  See code for further details.


% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com

% Set default parameters, and override as specified
if nargin < 4 || isempty(wrap)
  wrap = [false; false]; % non-periodic in both horizontal dimensions, by default
end
if nargin < 5 || isempty(pbot)
  pbot = 4000; % pressure at bottom [dbar]
end
if nargin < 6 || isempty(std)
  std = 6;  % roughly 200km for the 1024x1024 grid, zonally anyways. 
end


X = linspace(-1, 1, ni)';   % scaled longitude
Y = linspace(-80, 80, nj);  % latitude

% Set surface T
Ts = smooth2a(2 * randn(ni,nj), 8*std, std, [], wrap); % 8*std is large enough that the smoothing window is effectively the whole domain
Tmin = 10;
Tmax = 20;
Ts = (Ts - min(Ts(:))) / (max(Ts(:)) - min(Ts(:))) * (Tmax - Tmin) + Tmin;

% Set bottom T: same structure as surface T but different range
Tb = Ts;
Tmin = -1.8;
Tmax = 4;
Tb = (Tb - min(Tb(:))) / (max(Tb(:)) - min(Tb(:))) * (Tmax - Tmin) + Tmin;

% Set surface S: same structure as surface T but different range...
Smin = 34;
Smax = 36;
Ss = Ts;
Ss = (Ss - min(Ss(:))) / (max(Ss(:)) - min(Ss(:))) * (Smax - Smin) + Smin;
Ss = Ss + smooth2a(.1 * randn(ni,nj), 8*std, std, [], wrap); % ... + extra random field to surface S, to produce helicity

% Set bottom S: same as T but different range
Sb = Tb;
Smin = 36;
Smax = 36.5;
Sb = (Sb - min(Sb(:))) / (max(Sb(:)) - min(Sb(:))) * (Smax - Smin) + Smin;

% Add depth structure by linear interpolation of surface->bottom T,S data
T = reshape(linspace(0,1,nk), [1,1,nk]) .* (Tb - Ts) + Ts;
S = reshape(linspace(0,1,nk), [1,1,nk]) .* (Sb - Ss) + Ss;


% Nonlinear spacing of pressure means dTdp and dSdp will be non-uniform
P = linspace(0,1,nk).'.^3 * pbot;


% Re-order data so water-columns are contiguous data:
S = permute(S, [3 1 2]); % [nk,ni,nj]. depth  by  long  by  lat
T = permute(T, [3 1 2]);

if nargout >= 4
  X = (X+1) * 180;  % change longitude, to be [0, 360].
  
  XGorigin = 0;
  YGorigin = Y(1);
  
  g = struct();
  g.nx = ni;
  g.ny = nj;
  g.nz = nk;
  g.XCvec = X;
  g.YCvec = Y;
  g.rSphere = 6370000;
  
  g.resx = ni / 360;
  g.resy = nj / (Y(end) - Y(1));
  
  
  % Build grid variables with accordance with MITgcm's spherical polar
  % grid: See INI_SPHERICAL_POLAR_GRID.F Create our own spherical-polar
  % grid. Short-cuts have been taken in the formulas, appropriate for a
  % uniform latitude-longitude grid. See INI_SPHERICAL_POLAR_GRID.F if
  % adjustments must be made for non-uniform grid.
  deg2rad = pi / 180;
  
  g.XGvec = XGorigin + (0:g.nx-1)' / g.resx ;
  g.YGvec = YGorigin + (0:g.ny-1) / g.resy ;
  g.XCvec = g.XGvec + 1 / (2 * g.resx) ;
  g.YCvec = g.YGvec + 1 / (2 * g.resy) ;
  
  g.DXGvec = g.rSphere * cos(g.YGvec * deg2rad) / g.resx * deg2rad ;
  g.DYGsc = g.rSphere * deg2rad / g.resy;
  g.DXCvec = g.rSphere * cos(g.YCvec * deg2rad) / g.resx * deg2rad ;
  g.DYCsc = g.DYGsc;
  
  g.RACvec = g.rSphere * g.rSphere / g.resx * deg2rad * abs( sin((g.YGvec + 1/g.resy)*deg2rad) - sin(g.YGvec*deg2rad) ) ;
  g.RAWvec = g.RACvec;
  g.RASvec = g.rSphere * g.rSphere / g.resx * deg2rad * abs( sin(g.YCvec * deg2rad) - sin((g.YCvec - 1/g.resy)*deg2rad) );
  g.RAZvec = g.rSphere * g.rSphere / g.resx * deg2rad * abs( sin(g.YCvec*deg2rad) - sin((g.YCvec - 1/g.resy)*deg2rad) );
  
  
end