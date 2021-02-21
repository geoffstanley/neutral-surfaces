function epsL2 = ntp_epsilon_L2(x, s, t, WRAP, DX, DY)
% NTP_EPSILON_L2 Calculate the L2 norm of epsilon neutrality errors
%
%
% epsL2 = ntp_epsilon_L2(x, s, t)
% calculates the L2 norm of the epsilon neutrality error on a surface whose
% depth or pressure is x, whose practical / Absolute salinity is s, and whose
% potential / Conservative temperature is t. 
%
% epsL2 = ntp_epsilon_L2(x, s, t, WRAP)
% also specifies periodicity.  The i'th dimension is treated as periodic
% iff WRAP(i) is true.
%
% See also: ntp_errors

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



if nargin < 4 || isempty(WRAP)
  % doubly periodic domain, by default
  WRAP = [1 1]; 
end

if nargin < 5 || isempty(DX)
  DX = 1;
end
if nargin < 6 || isempty(DY)
  DY = 1;
end

% Calculate epsilon neutrality errors
[eps_i, eps_j] = ntp_errors(s, t, x, DX, DY, true, false, WRAP);

% Calculate L2 norm of the epsilon vector field
epsL2 = nanrms([eps_i(:); eps_j(:)]);

