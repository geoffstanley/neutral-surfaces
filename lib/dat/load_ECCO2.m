function varargout = load_ECCO2(datafolder, vars, timestep, zref)
%LOAD_ECCO2  Load ECCO2 data and/or create its grid.
%
%
% g = load_ECCO2(datafolder, 'grid')
% loads the grid struct g for ECCO2 model data located the directory
% datafolder.
%
% [S, T, ETAN, ATMP, SAP] = load_ECCO2(datafolder, 'casts', timestep)
% loads the practical salinity S, the potential temperature T, the
% sea-surface height ETAN, the atmospheric pressure ATMP, and the Standard
% Atmospheric Pressure SAP at time timestep for ECCO2 model data located
% the directory datafolder.
%
% G = load_ECCO2(datafolder, 'gamma', timestep)
% loads the neutral density G at time timestep for ECCO2 model data
% located the directory datafolder.
%
% U = load_ECCO2(datafolder, 'u', timestep)
% loads the zonal velocity U at time timestep for ECCO2 model data
% located the directory datafolder.
%
% V = load_ECCO2(datafolder, 'v', timestep)
% loads the meridional velocity V at time timestep for ECCO2 model data
% located the directory datafolder.
%
% z = load_ECCO2(datafolder, 'omega', timestep, zref)
% loads the depth z of an omega surface computed from ECCO2 model data
% located the directory datafolder at time timestep, initiated from
% potential density referenced to zref [m] and intersecting a (180 E, 0 N)
% at depth zref [m].
%
%
% --- Input:
% datafolder [string]: path to the folder containing folders containing
%                      ECCO2 netCDF files
% vars [string]: one of 'grid', 'casts', 'gamma', 'u', or 'v'
% timestep [string]: the date in yyyymmdd format, namely that part of the
%                    netcdf file name that specifies the timestep.
% zref [1, 1]: reference depth that initialized an omega-surface.
%
%
% --- Output:
% g [struct]: the model grid, containing standard MITgcm fields, although
%   some fields are replaced by scalar or vector quantities when they are
%   constant in some dimensions.
% S [nx,ny,nz]: practical salinity [psu (PSS-78)]
% T [nx,ny,nz]: potential temperature [degree C (IPTS-68)]
% ETAN [nx,ny]: sea-surface height [m]
% ATMP [1, 1]: zero. ECCO2 does not use atmospheric loading. [dbar]
% SAP  [1, 1]: zero. ECCO2 does not use additional Standard Atmospheric Pressure. [dbar]
% G [nx,ny,nz]: neutral density
% U [nx,ny,nz]: zonal velocity [m s^-1]
% V [nx,ny,nz]: meridional velocity [m s^-1]
% z [nx,ny]   : depth [m, positive] of a pre-computed omega-surface

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

V = filesep(); % /  or  \  depending on OS.
if datafolder(end) ~= V
    datafolder(end+1) = V;
end

switch lower(vars)
    
    case 'grid'
        
        % delR is the ECCO2 readme at
        % ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/readme.txt
        delR   = [10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01, ...
            10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04 , 19.82, 24.85, ...
            31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18, ...
            93.96, 96.58, 98.25, 99.25,100.01,101.33,104.56,111.33,122.83, ...
            139.09,158.94,180.83,203.55,226.50,249.50,272.50,295.50,318.50, ...
            341.50,364.50,387.50,410.50,433.50,456.50];
        
        g = struct();
        g.DRF = delR(:);                          % Vertical Distance (m) between W's
        g.RF = -[0; cumsum(g.DRF(:))];            % Depth (m, negative) of W
        g.RC = (g.RF(1:end-1) + g.RF(2:end)) / 2; % Depth (m, negative) of Tracer at cell centre
        g.DRC = -diff([0; g.RC]);                 % Vertical Distance (m) between T's
        
        g.nx = 1440;           % number of cells in longitude
        g.ny = 720;            % number of cells in latitude
        g.nz = length(g.RC);   % number of cells vertically
        g.grav = 9.81;         % Gravity [m s-2]
        g.rho_c = 1027.5;      % Reference density [kg m-3]
        g.rSphere = 6.3700e6;  % Radius of Earth [m]
        g.sidereal_day = 23.9344696 * 60 * 60; % Number of seconds per sidereal day on Earth.
        g.resx = 4;   % degrees per cell in longitude
        g.resy = 4;   % degrees per cell in latitude
        g.WRAP = [true false]; % Wrapped in longitude, not in latitude
        
        % Calculate the grid in accordance with MITgcm's spherical polar
        % grid: See INI_SPHERICAL_POLAR_GRID.F (This creates the grid
        % correctly in double precision.) Short-cuts taken below,
        % appropriate for a uniform latitude-longitude grid.
        deg2rad = pi / 180;
        
        XGorigin = 0;
        YGorigin = -90;
        g.XGvec = XGorigin + (0:g.nx-1)' / g.resx ; % Longitude (1D vector) at centre of West cells (Longitude of U)
        g.YGvec = YGorigin + (0:g.ny-1) / g.resy ;  % Latitude (1D vector) at centre of South cells (Latitude of V)
        g.XCvec = g.XGvec + 1 / (2 * g.resx) ;      % Longitude (1D vector) at centre of Tracer cells (Longitude of T)
        g.YCvec = g.YGvec + 1 / (2 * g.resy) ;      % Latitude (1D vector) at centre of Tracer cells (Latitude of T)
        
        g.DXGvec = g.rSphere * cos(g.YGvec * deg2rad) / g.resx * deg2rad ;  % X Distance (m, 1D vector) between West Faces (distance between U's)
        g.DYGsc = g.rSphere * deg2rad / g.resy;                             % Y Distance (m, scalar) between South Faces (distance between V's)
        g.DXCvec = g.rSphere * cos(g.YCvec * deg2rad) / g.resx * deg2rad ;  % X Distance (m, 1D vector) between Tracer cell centres
        g.DYCsc = g.DYGsc;                                                  % Y Distance (m, scalar) between Tracer cell centres
        
        g.RACvec = g.rSphere^2 / g.resx * deg2rad * abs( sin((g.YGvec + 1/g.resy)*deg2rad) - sin(g.YGvec*deg2rad) ) ;   % Vertical area of the tracer cells [m^2]
        g.RAWvec = g.RACvec;                                                                                            % Vertical area of the U cells [m^2]
        g.RASvec = g.rSphere^2 / g.resx * deg2rad * abs( sin(g.YCvec * deg2rad) - sin((g.YCvec - 1/g.resy)*deg2rad) );  % Vertical area of the V cells [m^2]
        g.RAZvec = g.rSphere^2 / g.resx * deg2rad * abs( sin(g.YCvec*deg2rad) - sin((g.YCvec - 1/g.resy)*deg2rad) );    % Vertical area of the vorticity cells [m^2]
        
        g.cori    = 2 * (2*pi/g.sidereal_day) * sind(g.YCvec); % Coriolis parameter (s^-1) at Tracer cell centres
        g.cori_YG = 2 * (2*pi/g.sidereal_day) * sind(g.YGvec); % Coriolis parameter (s^-1) at South Faces
        
        varargout{1} = g;
        
    case 'casts'
        S     = ncread([datafolder 'SALT.nc'  V 'SALT.1440x720x50.'  timestep '.nc'], 'SALT');
        T     = ncread([datafolder 'THETA.nc' V 'THETA.1440x720x50.' timestep '.nc'], 'THETA');
        
        n = datenum(timestep, 'yyyymmdd');
        SSH1  = ncread([datafolder 'SSH.nc' V 'SSH.1440x720.' datestr(n-1, 'yyyymmdd')  '.nc'], 'SSH');
        SSH2  = ncread([datafolder 'SSH.nc' V 'SSH.1440x720.' timestep                  '.nc'], 'SSH');
        SSH3  = ncread([datafolder 'SSH.nc' V 'SSH.1440x720.' datestr(n+1, 'yyyymmdd')  '.nc'], 'SSH');
        SSH = double(SSH1 + SSH2 + SSH3) / 3;
        clear SSH1 SSH2 SSH3 n
        
        % Remove really low salinity data, which messes things up.
        % Second check SALT > -1000 is because land data has SALT = -9.999999800000000e+22
        S(S < -1000) = nan;
        nz = size(S,3);
        bad = repmat(any(S < 15,3), [1 1 nz]);
        S(bad) = NaN;
        bad = isnan(S);
        T(bad) = NaN;
        SSH(bad(:,:,1)) = NaN;
        ATMP = 0; % no atmospheric loading
        SAP = 0; % Standard Atmospheric Pressure
        
        S = double(S);
        T = double(T);
        
        varargout{1} = S;
        varargout{2} = T;
        varargout{3} = ATMP;
        varargout{4} = SSH;
        varargout{5} = SAP;
        
    case 'gamma'
        GAMMA = load([datafolder 'GAMMA' V 'GAMMA.1440x720x50.' timestep '.mat'], 'GAMMA');
        GAMMA = double(GAMMA.GAMMA);
        GAMMA = sort(GAMMA,3); % ensure monotonic!
        varargout{1} = GAMMA;
        
    case 'u'
        varargout{1} = double(ncread([datafolder 'UVEL.nc' V 'UVEL.1440x720x50.' timestep '.nc'], 'UVEL'));
        
    case 'v'
        varargout{1} = double(ncread([datafolder 'VVEL.nc' V 'VVEL.1440x720x50.' timestep '.nc'], 'VVEL'));
        
    case 'omega'
        LOAD = load(sprintf('%somega.1440x720.%s.from_SIGMA%04d_through_(180,0,%04d).Boussinesq.mat', [datafolder 'omega_v1.1gjs' V], timestep, zref, zref), 'zns_i');
        varargout{1} = squeeze(LOAD.zns_i).';
end