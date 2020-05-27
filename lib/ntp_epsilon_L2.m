function epsL2 = ntp_epsilon_L2(x, s, t, wrap)
% NTP_EPSILON_L2 Calculate the L2 norm of epsilon neutrality errors
%
%
% epsL2 = ntp_epsilon_L2(x, s, t)
% calculates the L2 norm of the epsilon neutrality error on a surface whose
% depth or pressure is x, whose practical / Absolute salinity is s, and whose
% potential / Conservative temperature is t. 
%
% epsL2 = ntp_epsilon_L2(x, s, t, wrap)
% also specifies periodicity.  The i'th dimension is treated as periodic
% iff wrap(i) is true.
%
% See also: ntp_errors

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
% Version   : 2.2.0
%
% Modified by : --
% Date        : --
% Changes     : --


if nargin < 4 || isempty(wrap)
  % doubly periodic domain, by default
  wrap = [1 1]; 
end

% Calculate epsilon neutrality errors
[eps_i, eps_j] = ntp_errors(s, t, x, 1, 1, true, false, wrap);

% Calculate number of valid neutral links
neq = sum(isfinite(eps_i(:))) + sum(isfinite(eps_j(:)));

% Calculate L2 norm of the epsilon vector field
epsL2 = sqrt((nansum(eps_i(:).^2) + nansum(eps_j(:).^2)) / neq);

