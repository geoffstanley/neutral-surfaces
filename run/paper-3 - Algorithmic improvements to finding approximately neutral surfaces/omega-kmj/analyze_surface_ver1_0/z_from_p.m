function [z] = z_from_p(p,g)

%           Calculate depth from pressure
%
% Usage:    [z] = z_from_p(p,g)
%
%           Calculate depth from pressure according to 'Depth-pressure 
%           relationships in the oceans and seas, Leroy and Parthiot, 
%           J. Acoust. Soc. Am., 1998'
%
% Input:    p           pressure
%           g           gravitational acceleration
%
% Output:   z           depth    
%
% Calls:    none
%
% Units:    pressure                    dbar
%           depth                       m
%           gravitational acceleration  m/s^2
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments
if ~(nargin ==2)
  error('z_from_p.m: requires 2 input arguments')
end

%% change pressure from dbar into MPa

p = p ./ 100;  

p2 = p .* p; p3 = p2 .* p; p4 = p3 .* p;

%% change g-matrix to be the same size as p-matrix

if (length(size(g)) == 2)
    [k,dummy,dummy] = size(p); %#ok
    g_rep = repmat(g,[1 1 k]);
    g = permute(g_rep,[3 1 2]);
end

%% calculate numerator

znum = p  .* 9.72659e002 + ...
       p2 .* 2.512e-001 + ...
       p3 .* 2.279e-004 + ...
       p4 .* 1.82e-007;
   
%% calculate denumerator

zden = g + ...
       p .* 1.092e-004;

%% calculate depth

z = znum ./ zden;