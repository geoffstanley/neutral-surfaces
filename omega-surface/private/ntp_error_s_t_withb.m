function [eps_i, eps_j, beps_i, beps_j] = ntp_error_s_t_withb(s, t, x, WRAP, b)
%NTP_ERROR_S_T_WITHB  Calculate zonal and meridional errors from the 
%                     neutral tangent plane with the integrating factor b

% [eps_i, eps_j] = ntp_error_s_t_withb(s, t, x, WRAP)
% calculates the neutrality error epsilon for a surface with practical /
% Absolute salinity s, potential / Conservative temperature t, and pressure
% or depth x. At water column (i,j),
% eps_i(i,j) = eos_s * (s(i,j) - s(i-1,j)) + eos_t * (t(i,j) - t(i-1,j))
% where eos_s and eos_t are the partial derivatives of eos with respect to
% s and t, calculated from the function eos_s_t.m in the same folder as
% this function, and evaluated at
% (s(i,j) - s(i-1,j))/2, (t(i,j) - t(i-1,j))/2, (x(i,j) - x(i-1,j))/2.
% Likewise for eps_j but swapping (i-1,j) for (i,j-1).
% Data is treated periodic in the n'th dimension (n = 1 or n = 2) if and
% only if WRAP(n) is true.
%
% [eps_i, eps_j, beps_i, beps_j] = ntp_error_s_t_withb(s, t, x, WRAP, b)
% also calculates b * eps_i and b * eps_j, where b is a non-dimensional
% scalar field the same size as s, t, and x, namely the integrating factor.
%
%
% --- Input
% s [ni,nj]: practical / Absolute salinity on the surface
% t [ni,nj]: potential / Conservative temperature on the surface
% x [ni,nj]: pressure or depth on the surface
% WRAP [2 element vector]: determines which dimensions are periodic.
% b [ni,nj]: integrating factor on the surface [dimensionless]
%
%
% --- Output
% eps_i [ni,nj]: neutrality error in 1st dimension
% eps_j [ni,nj]: neutrality error in 2nd dimension
% beps_i [ni,nj]: b times neutrality error in 1st dimension
% beps_j [ni,nj]: b times neutrality error in 2nd dimension
%
%
% --- Units
% The units of s, t, and x are determined by the functions eos.m and
% eos_s_t.m. The units of eps_i and eps_j are the same as the units output
% from the function eos.m.  Note, eps_i and eps_j incorporate finite
% differences of s and t, not finite difference approximations of the
% gradients of s and t.
%
%
% --- Requirements:
% eos.m, eos_s_t.m

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

USE_INTEGRATING_FACTOR = (nargin == 5);

shift_im1 = @(F) circshift(F, [+1, 0]); % if x is [lat x lon], this shifts south
shift_jm1 = @(F) circshift(F, [0, +1]); % if x is [lat x lon], this shifts west

% --- Errors in 1st dimension (meridional, if x is [lat x lon])
s_sh = shift_im1(s);
t_sh = shift_im1(t);
x_sh = shift_im1(x);

s_dif = s - s_sh;        % backward finite difference
t_dif = t - t_sh;        % backward finite difference
s_avg = (s + s_sh) / 2;  % backward average
t_avg = (t + t_sh) / 2;  % backward average
x_avg = (x + x_sh) / 2;  % backward average

if ~WRAP(1) % Handle non-periodic boundary in first dimension
    % Could also make SAns_dif, CTns_dif, SAns_avg, and CTns_avg nan here,
    % but this nan will infect them all.
    x_avg(1,:) = nan;
end

[eos_s, eos_t] = eos_s_t(s_avg, t_avg, x_avg);

eps_i = (eos_s .* s_dif) + (eos_t .* t_dif);

if USE_INTEGRATING_FACTOR
    b_sh = shift_im1(b);
    b_avg = (b + b_sh) / 2;  % backward average
    beps_i = b_avg .* eps_i;
end

% --- Errors in 2nd dimension (zonal, if x is [lat x lon])
s_sh = shift_jm1(s);
t_sh = shift_jm1(t);
x_sh = shift_jm1(x);

s_dif = s - s_sh;        % backward finite difference
t_dif = t - t_sh;        % backward finite difference
s_avg = (s + s_sh) / 2;  % backward average
t_avg = (t + t_sh) / 2;  % backward average
x_avg = (x + x_sh) / 2;  % backward average

if ~WRAP(2) % Handle non-periodic boundary in second dimension
    % Could also make SAns_dif, CTns_dif, SAns_avg, and CTns_avg nan here,
    % but this nan will infect them all.
    x_avg(:,1) = nan;
end

[eos_s, eos_t] = eos_s_t(s_avg, t_avg, x_avg);

eps_j = (eos_s .* s_dif) + (eos_t .* t_dif);

if USE_INTEGRATING_FACTOR
    b_sh = shift_im1(b);
    b_avg = (b + b_sh) / 2;  % backward average
    beps_j = b_avg .* eps_j;
end