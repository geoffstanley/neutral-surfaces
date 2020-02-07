function [g, S, T, P, ETAN, ATMP, SAP, GAMMA] = load_OCCA(datafolder)
%LOAD_OCCA  Load OCCA grid and data.
%
%
% [g, S, T, P, ETAN, ATMP, SAP] = load_OCCA(datafolder)
% loads the grid struct g, the practical salinity S, the potential
% temperature T, the pressure P, the sea-surface height ETAN, the
% atmospheric pressure ATMP, and the Standard Atmospheric Pressure SAP, for
% the OCCA 2004 - 2006 climatology output located the directory datafolder.
%
%
% --- Input:
% datafolder [string]: path to the folder containing OCCA netCDF files
%
%
% --- Output:
% g [struct]: the model grid, containing standard MITgcm fields, although
%   some fields are replaced by scalar or vector quantities when they are
%   constant in some dimensions.
% S [nx,ny,nz]: practical salinity [psu (PSS-78)]
% T [nx,ny,nz]: potential temperature [degree C (IPTS-68)]
% P [nx,ny,nz]: pressure [dbar]
% ETAN [nx,ny]: sea-surface height [m]
% ATMP [1, 1]: zero. OCCA does not use atmospheric loading. [dbar]
% SAP  [1, 1]: zero. OCCA does not use additional Standard Atmospheric Pressure. [dbar]

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

V = filesep(); % /  or  \  depending on OS.
if datafolder(end) ~= V
    datafolder(end+1) = V;
end

Pa2db = 1e-4;
deg2rad = pi / 180;

g.rho_c = 1027.5; % A guess. Same as ECCO2
g.grav = 9.81; % A guess. Same as ECCO2
g.rSphere = 6370000; % A guess. Same as ECCO2
g.WRAP = [true false];

varname = 'theta';
ncid = netcdf.open(sprintf('%sDD%s.0406annclim.nc', datafolder, varname));
g.XCvec = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'Longitude_t')));
g.XGvec = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'Longitude_u')));
g.YCvec = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'Latitude_t')));
g.YGvec = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'Latitude_v')));
g.XCvec = g.XCvec(:);
g.YCvec = g.YCvec(:).';
g.XGvec = g.XGvec(:);
g.YGvec = g.YGvec(:).';
Depth_c = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'Depth_c'));
%Depth_l = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'Depth_l'));
%Depth_u = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'Depth_u'));
%Depth_w = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'Depth_w'));
g.RC = double(-Depth_c);
g.DRC = -diff([0; g.RC]);

% Build our own grid, as the MITgcm does it
g.resx = 1; % 1 grid cell per zonal degree
g.resy = 1; % 1 grid cell per meridional degree
g.DXGvec = g.rSphere * cos(g.YGvec * deg2rad) / g.resx * deg2rad ;
g.DYGsc = g.rSphere * deg2rad / g.resy;
g.DXCvec = g.rSphere * cos(g.YCvec * deg2rad) / g.resx * deg2rad ;
g.DYCsc = g.DYGsc;

g.RACvec = g.rSphere^2 / g.resx * deg2rad * abs( sin((g.YGvec + 1/g.resy)*deg2rad) - sin(g.YGvec*deg2rad) ) ;   % Vertical area of the tracer cells [m^2]
g.RAWvec = g.RACvec;                                                                                            % Vertical area of the U cells [m^2]
g.RASvec = g.rSphere^2 / g.resx * deg2rad * abs( sin(g.YCvec * deg2rad) - sin((g.YCvec - 1/g.resy)*deg2rad) );  % Vertical area of the V cells [m^2]
g.RAZvec = g.rSphere^2 / g.resx * deg2rad * abs( sin(g.YCvec*deg2rad) - sin((g.YCvec - 1/g.resy)*deg2rad) );    % Vertical area of the vorticity cells [m^2]


g.nx = length(g.XCvec);
g.ny = length(g.YCvec);
g.nz = length(g.RC);

% Continue reading from theta file
varid = netcdf.inqVarID(ncid, varname);
T = netcdf.getVar(ncid, varid);
mv = netcdf.getAtt(ncid,varid,'missing_value');
T(T == mv) = NaN;
T = double(T);
netcdf.close(ncid);

varname = 'salt';
ncid = netcdf.open(sprintf('%sDD%s.0406annclim.nc', datafolder, varname));
varid = netcdf.inqVarID(ncid, varname);
S = netcdf.getVar(ncid, varid);
mv = netcdf.getAtt(ncid,varid,'missing_value');
S(S == mv) = NaN;
S = double(S);
netcdf.close(ncid);

% phihyd = Pres / rho_c +  grav * z
varname = 'phihyd';
ncid = netcdf.open(sprintf('%sDD%s.0406annclim.nc', datafolder, varname));
varid = netcdf.inqVarID(ncid, varname);
P = netcdf.getVar(ncid, varid);
mv = netcdf.getAtt(ncid,varid,'missing_value');
P(P == mv) = NaN;
P = double(P);
netcdf.close(ncid);
P = (P - g.grav * permute(g.RC(:), [3 2 1])) * g.rho_c * Pa2db;
clear mv ncid varname varid

varname = 'etan';
ncid = netcdf.open(sprintf('%sDD%s.0406annclim.nc', datafolder, varname));
varid = netcdf.inqVarID(ncid, varname);
ETAN = netcdf.getVar(ncid, varid);
mv = netcdf.getAtt(ncid,varid,'missing_value');
ETAN(ETAN == mv) = NaN;
ETAN = double(ETAN);
netcdf.close(ncid);

ATMP = 0;
SAP = 0;

if nargout >= 8
    LOADED = load(sprintf('%sGAMMA.0406annclim.mat', datafolder));
    GAMMA = LOADED.GAMMA;
    clear LOADED
    GAMMA = permute(GAMMA, [2 3 1]);  % change back from [z,x,y] to [x,y,z] for consistency with above. 
    GAMMA(GAMMA < 0) = NaN; % set bad values (should be -99) to NaN
end