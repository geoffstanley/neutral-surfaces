function MLP = mixed_layer_pressure(S, T, P, OPTS)
%MIXED_LAYER_PRESSURE  The pressure at the bottom of the mixed layer.
%
%
% MLP = mixed_layer_pressure(S, T, P) 
% calculates the pressure of the mixed layer, MLP, in an ocean of practical
% / absolute salinity S, potential / Conservative temperature T, and
% pressure P. MLP is pressure at which the potential density (referenced to
% 100 dbar) exceeds the potential density near the surface (the second
% bottle on each cast) by 0.03 kg m^-3.
%
% For a Boussinesq ocean, pass P as the Boussinesq pressure, namely 
% P = Z * Z2P
% where Z is the depth [m] of the S and T data (positive and increasing
% downwards) and Z2P = (1e-4 * rho_c * grav) where grav is the
% gravitational acceleration and rho_c is the Boussinesq reference density.
% Should the mixed layer depth be desired instead of the mixed layer
% pressure, simply divide MLP by Z2P.
%
% mlz = mixed_layer_pressure(S, T, P, Z2P) 
% calculates the depth of the mixed layer, mlz, in a Boussinesq ocean with
% practical / absolute salinity S and potential / Conservative temperature
% T at depth Z, given a depth to pressure conversion factor of Z2P, which
% typically is given by Z2P = 1e-4 * rho_c * grav where grav is the
% gravitational acceleration and rho_c the Boussinesq reference density.
% mlz is the depth at which the potential density (referenced to 100 m)
% exceeds the potential density near the surface (the second bottle on each
% cast) by 0.03 kg m^-3.
% 
% ... = mixed_layer_pressure(..., OPTS) overwrites the default parameters
% according to the struct OPTS (see below).
% 
%
% --- Input:
% S [nz, nx, ny]: practical / Absolute salinity
% T [nz, nx, ny]: potential / Conservative temperature
% P [nz, nx, ny]: pressure [dbar]
% OPTS [struct]: options
%   OPTS.POT_DENS_REF [1, 1]: the reference pressure or depth for potential
%    density [dbar or m, positive]
%   OPTS.POT_DENS_DIFF [1, 1]: the potential density difference that
%    determines the mixed layer [kg m^-3]
%   OPTS.BOTTLE_NUM [1, 1]: the bottle number on each cast where the "near
%    surface" potential density is calculated [integer]
%   OPTS.EOS_INSITUDENS [1, 1]: true when eos.m calculates the in-situ
%    density, and false when eos.m calculates the specific volume [logical]
%
% Note: nz is the maximum number of data points per water column,
%       nx is the number of data points in longitude,
%       ny is the number of data points in latitude.
%
%
% --- Output:
% mlp [nx, ny]: the pressure of the mixed layer [dbar]
% mlz [nx, ny]: the depth of the mixed layer [m, positive]
%
% 
% --- Requirements:
% eos
% interp1qn - https://www.mathworks.com/matlabcentral/fileexchange/69713

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


% Set defaults:
POT_DENS_DIFF = 0.03;   % [kg m^-3]
POT_DENS_REF = 100;     % [dbar]
BOTTLE_NUM = 2;
EOS_INSITUDENS = true;

if nargin == 4 && isstruct(OPTS)
    % Load parameters from struct
    if isfield(OPTS, 'POT_DENS_DIFF')
        POT_DENS_DIFF = OPTS.POT_DENS_DIFF;
    end
    if isfield(OPTS, 'POT_DENS_REF')
        POT_DENS_REF = OPTS.POT_DENS_REF;
    end
    if isfield(OPTS, 'BOTTLE_NUM')
        BOTTLE_NUM = OPTS.BOTTLE_NUM;
    end
    if isfield(OPTS, 'EOS_INSITUDENS')
        EOS_INSITUDENS = OPTS.EOS_INSITUDENS;
    end
end


% If not for these three lines, some eos functions would try to transpose a
% 3D matrix.
lead1 = @(x) reshape(x, [1 size(x)]); 
SB = squeeze(S(BOTTLE_NUM,:,:));
TB = squeeze(T(BOTTLE_NUM,:,:));

% Calculate potential density difference between each data point and the
% near-surface bottle
if EOS_INSITUDENS
    DD =      eos(S, T, POT_DENS_REF) - lead1(eos(SB, TB, POT_DENS_REF)) ;
else
    DD = 1 ./ eos(S, T, POT_DENS_REF) - lead1(1 ./ eos(SB, TB, POT_DENS_REF)) ;
end

% Find the pressure at which the potential density difference exceeds the
% threshold POT_DENS_DIFF
MLP = interp1qn(POT_DENS_DIFF, DD, P);

