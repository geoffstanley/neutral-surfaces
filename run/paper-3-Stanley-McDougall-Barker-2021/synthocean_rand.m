function [S, T, P, g] = synthocean_rand(ni, nj, nk, args)
% SYNTHOCAEN_RAND   Synthetic ocean of eddies
%
%
% [S, T, P, g] = synthocean(ni, nj, nk, args)
% constructs a synthetic ocean of practical / Absolute salinity S,
% potential / Conservative temperature T, pressure P, having ni, nj, nk
% data points in the longitude, latitude depth directions, respectively.
% The grid data is contained in the struct g.Extra tuneable parameters are
% provided by args.  See code for further details.

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


defaults = struct();
defaults.pbot = 4000; % pressure at bottom [dbar]
defaults.std = 6;  % roughly 200km for the 1024x1024 grid, zonally anyways. 
if nargin < 4 || ~isstruct(args)
  args = struct();
end
args = catstruct(defaults, args);  % load input arguments, resorting to defaults for any not provided

std = args.std;

X = linspace(-1, 1, ni)';   % scaled longitude
Y = linspace(-80, 80, nj);  % latitude


% Set surface T
Ts = smooth2a(2 * randn(ni,nj), 8*std, std); % 8*std is large enough that the smoothing window is effectively the whole domain
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
Ss = Ss + smooth2a(.1 * randn(ni,nj), 8*std, std); % ... + extra random field to surface S, to produce helicity

% Set bottom S: same as T but different range
Sb = Tb;
Smin = 36;
Smax = 36.5;
Sb = (Sb - min(Sb(:))) / (max(Sb(:)) - min(Sb(:))) * (Smax - Smin) + Smin;

% Add depth structure by linear interpolation of surface->bottom T,S data
T = reshape(linspace(0,1,nk), [1,1,nk]) .* (Tb - Ts) + Ts;
S = reshape(linspace(0,1,nk), [1,1,nk]) .* (Sb - Ss) + Ss;


% Nonlinear spacing of pressure means dTdp and dSdp will be non-uniform
P = linspace(0,1,nk).'.^3 * args.pbot;


% Re-order data so water-columns are contiguous data:
S = permute(S, [3 1 2]); % [nk,ni,nj]. depth  by  long  by  lat
T = permute(T, [3 1 2]);

% Add walls
S(:,1,:) = nan;
S(:,:,1) = nan;
T(isnan(S)) = nan;

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