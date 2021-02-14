function [epsL2, epsL1] = eps_norms(s, t, x, use_s_t, wrap, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij)
%EPS_NORMS  The L1 and L2 norms of the epsilon neutrality error
%
%
% [epsL2, epsL1] = eps_norms(s, t, x, use_s_t, wrap, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij)
% computes the L2 and L1 norms of the neutrality error vector, which itself
% is calculated by ntp_errors.m, on approximately neutral surface on which
% the practical / Absolute salinity is s, the potential / Conservative
% temperature is t, and the pressure or depth is x. If use_s_t is true,
% errors are calculated using s and t differences multiplied by the partial
% derivatives of in-situ density w.r.t s and t, respectively; if use_s_t is
% false, errors are calculated using in-situ density differences and x
% differences multiplied by the partial derivative of in-situ density
% w.r.t. x. Data is treated periodic in the i'th (i=1 or 2) dimension of x
% if and only if wrap(i) is true. The distances and areas of the grid are
% used to weight the neutrality errors when calculating its norm.  See
% Input section for details.
%
% The equation of state for the in-situ density is given by eos.m in the
% path.  If use_s_t is true, eos_s_t.m must also exist in the path and
% return the partial derivatives of in-situ density w.r.t s and t;
% otherwise, eos_x.m must also exist in the path and return the partial
% derivative of in-situ density w.r.t x. The inputs to these eos*.m
% functions must be (s, t, x).  Note, eos*.m can involve specific volume
% instead of in-situ density, which merely changes the units of [ex,ey].
%
%
% For a non-Boussinesq ocean, x is pressure [dbar].  
% For a Boussinesq ocean, x is depth [m, positive].  
%
%
% --- Input:
% s     [ni, nj]: practical / Absolute salinity on the surface
% t     [ni, nj]: potential / Conservative temperature on the surface
% x     [ni, nj]: pressure [dbar] or depth [m, positive] of the surface
% use_s_t [1, 1]: true to compute ex,ey using s and t differences,
%                 false to use eos(s,t,x) and x differences.
% wrap    [2, 1]: wrap(i) is true when the data is periodic in i'th
%                 dimension of x
% DIST1_iJ [ni, nj]: Distance [m] in 1st dimension centred at (I-1/2, J)
% DIST2_Ij [ni, nj]: Distance [m] in 2nd dimension centred at (I, J-1/2)
% DIST2_iJ [ni, nj]: Distance [m] in 2nd dimension centred at (I-1/2, J)
% DIST1_Ij [ni, nj]: Distance [m] in 1st dimension centred at (I, J-1/2)
% AREA_iJ [ni, nj]: Area [m^2] centred at (I-1/2, J). Optional. Should = DIST1_iJ .* DIST2_iJ
% AREA_Ij [ni, nj]: Area [m^2] centred at (I, J-1/2). Optional. Should = DIST1_Ij .* DIST2_Ij
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in 1st horizontal dimension,
%       nj is the number of data points in 2nd horizontal dimension.
%
% Note: physical units for s,t,x are determined by eos*.m.
%
%
% --- Output:
% epsL2 [ni, nj]: L2 norm of epsilon neutrality error [(kg m^-3) / m or (m^3 / kg) / m]
% epsL1 [ni, nj]: L1 norm of epsilon neutrality error [(kg m^-3) / m or (m^3 / kg) / m]


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

flat = @(F) F(:);

% Compute epsilon neutrality errors without handling grid distances
[eps_iJ, eps_Ij] = ntp_errors(s, t, x, 1, 1, use_s_t, false, wrap);

if nargin < 11
  AREA_iJ = DIST1_iJ .* DIST2_iJ;   % Area [m^2] centred at (I-1/2, J)
  AREA_Ij = DIST1_Ij .* DIST2_Ij;   % Area [m^2] centred at (I, J-1/2)
end

% L2 norm of vector [a_i], weighted by vector [w_i], is sqrt( sum( w_i * a_i^2 ) / sum( w_i ) )
% Here, weights are AREA_iJ and AREA_Ij.
% But also need to divide epsilon by grid distances DIST1_iJ and DIST2_Ij.
% Thus, the numerator of L2 norm needs to multiply epsilon^2 by
%     AREA_iJ ./ DIST1_iJ.^2 = DIST2_iJ ./ DIST1_iJ , 
% and AREA_Ij ./ DIST2_Ij.^2 = DIST1_Ij ./ DIST2_Ij .
epsL2 = sqrt( ...
  (sum(flat( DIST2_iJ ./ DIST1_iJ .* eps_iJ.^2), 'omitnan') + sum(flat( DIST1_Ij ./ DIST2_Ij .* eps_Ij.^2 ), 'omitnan')) ./ ...
  (sum(flat(  AREA_iJ .* isfinite(eps_iJ) ))                + sum(flat(  AREA_Ij .* isfinite(eps_Ij))     ) ) ...
  );

if nargout < 2; return; end

% L1 norm of vector [a_i], weighted by vector [w_i], is sum( w_i * |a_i| ) / sum( w_i )
% Here, weights are AREA_iJ and AREA_Ij.
% But also need to divide epsilon by grid distances DIST1_iJ and DIST2_Ij.
% Thus, the numerator of L1 norm needs to multiply epsilon by
%     AREA_iJ ./ DIST1_iJ = DIST2_iJ , 
% and AREA_Ij ./ DIST2_Ij = DIST1_Ij .
epsL1 =  ...
  (sum(flat( DIST2_iJ .* abs(eps_iJ) ), 'omitnan' ) + sum(flat( DIST1_Ij .* abs(eps_Ij) ), 'omitnan' )) ./ ...
  (sum(flat(  AREA_iJ .* isfinite(eps_iJ) ))        + sum(flat(  AREA_Ij .* isfinite(eps_Ij) ) ) ) ;
