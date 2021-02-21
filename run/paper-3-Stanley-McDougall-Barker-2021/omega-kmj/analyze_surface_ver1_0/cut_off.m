function [var_new] = cut_off(var,pns,press)

%           Cut off data on a surface above a certain depth. 
%
% Usage:    [var_new] = cut_off(var,pns,depth)
%
%           Cut off data on a surface above a certain depth (e.g. to 
%           exclude data above the mixed-layer depth). 
%
% Input:    var         variable to cut off
%           pns         pressure on surface
%           press       cut-off pressure
%
% Output:   var_new     cut-off variable
%
% Units:    pressure    dbar
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin==3)
  error('cut_off.m: requires 3 input arguments')
end 

%% find points shallower than 'depth' and replace with nan

mask = change(pns,'<=',press,nan);
mask = change(mask,'>=',press,1);

var_new = var .* mask;