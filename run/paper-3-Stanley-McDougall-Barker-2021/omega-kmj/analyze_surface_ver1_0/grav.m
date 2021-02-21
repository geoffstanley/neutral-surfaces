function [g] = grav(lats)

%           Calculate gravitational acceleration
%
% Usage:    [g] = grav(lats)
%
%           Calculate gravitational acceleration according to 'Geodetic 
%           reference system 1980, Journal of Geodesy, 74, 128-133'. 
%
%
% Input:    lats        latitude
%
% Output:   g           gravitational acceleration    
%
% Calls:    none
%
% Units:    g           m/s^2     
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%


%% check input arguments

if ~(nargin ==1)
  error('grav.m: requires 1 input argument')
end 

%% calculate conversion factors

deg2rad = pi/180;
lats_abs = abs(lats);
lats_rad = lats_abs * deg2rad; % convert to radians
x = sin(lats_rad);  
xx = sin(2 .* lats_rad);
sin2 = x .* x;
sin2x = xx .* xx;

%% calculate gravity

g = 9.780327 .* (1 + 5.3024e-3 .* sin2 - 5.8e-6 .* sin2x);
