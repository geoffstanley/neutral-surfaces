function MLX = mixed_layer(S, T, X, OPTS)
%MIXED_LAYER  The pressure or depth at the bottom of the mixed layer.
%
%
% MLX = mixed_layer(S, T, X)
% calculates the pressure or depth of the mixed layer, MLX, in an ocean
% with practical / absolute salinity S and potential / Conservative
% temperature T located at datasites where the pressure or depth is X.  The
% equation of state for either the in-situ density or the specific volume
% is given by eos.m in the path, which accepts S, T, X as its 3 inputs.
% MLX is the pressure or depth at which the potential density (referenced
% to X = 100 dbar or X = 100 m) exceeds the potential density near the
% surface (the second bottle on each cast) by 0.03 kg m^-3.
%
% ... = mixed_layer(..., OPTS) overwrites the default parameters
% according to the struct OPTS (see below).
%
%
% --- Input:
% S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% X [nk, ni, nj]: pressure [dbar] or depth [m, positive]
% OPTS [struct]: options
%   OPTS.POT_DENS_REF [1, 1]: the reference pressure or depth for potential
%    density [dbar or m, positive]
%   OPTS.POT_DENS_DIFF [1, 1]: the potential density difference that
%    determines the mixed layer [kg m^-3]
%   OPTS.BOTTLE_NUM [1, 1]: the bottle number on each cast where the "near
%    surface" potential density is calculated [integer]
%   OPTS.INTERPFN [function handle]: the vertical interpolation function
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: X must increase monotonically along the first dimension.
%
%
% --- Output:
% MLX [ni, nj]: the pressure [dbar] or depth [m, positive] of the mixed layer

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
% Version   : 2.0.0
%
% Modified by : --
% Date        : --
% Changes     : --


% Set defaults:
POT_DENS_DIFF = 0.03;   % [kg m^-3]
POT_DENS_REF = 100;     % [dbar] or [m], depending on eos()
BOTTLE_NUM = 2;         % reference to the second bottle in the vertical
INTERPFN = @ppc_linterp;   % use linear interpolation

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
    if isfield(OPTS, 'INTERPFN')
        INTERPFN = OPTS.INTERPFN;
    end
end

% These three lines prepare data for a silly case where eos.m transposes
% its inputs.
lead1 = @(x) reshape(x, [1 size(x)]);
SB = squeeze(S(BOTTLE_NUM,:,:));
TB = squeeze(T(BOTTLE_NUM,:,:));

% Calculate potential density difference between each data point and the
% near-surface bottle
if eos(34.5, 3, 1000) > 1  % eos computes in-situ density
    DD =      eos(S, T, POT_DENS_REF) - lead1(eos(SB, TB, POT_DENS_REF)) ;
else                % eos computes specific volume
    DD = 1 ./ eos(S, T, POT_DENS_REF) - lead1(1 ./ eos(SB, TB, POT_DENS_REF)) ;
end

% Find the pressure or depth at which the potential density difference
% exceeds the threshold POT_DENS_DIFF
MLX = INTERPFN(DD, X, POT_DENS_DIFF);

