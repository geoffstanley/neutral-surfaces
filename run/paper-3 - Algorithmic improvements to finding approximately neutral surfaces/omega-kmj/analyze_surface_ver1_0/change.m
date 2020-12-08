function [matrix_new] = change(matrix,relation,old,new)

%           Change values in matrix
%
% Usage:    [matrix_new] = change(matrix,relation,old,new)
%
%           Changes the 'old' values in the 'matrix' to the 'new' values 
%           according to the 'relation'.
%
%
% Input:    matrix        matrix in which values are changed
%           relation      string relation e.g. '<', '>', '=='
%           old           old values
%           new           new values
%
% Output:   matrix_new    matrix with changed values
%
% Calls:    none
%
% Units:    
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 4)
    error('change.m: requires 4 input arguments')
end

%% check if relation given is valid

if (strcmp(relation,'==') || strcmp(relation,'>=') || strcmp(relation,'<=') || strcmp(relation,'>') || strcmp(relation,'<'))
  else
    error(['change.m: Relation {' relation '} not valid'])
end

%% find 'old' values in 'matrix'

if isnan(old)
   replace = find(isnan(matrix));
else
   eval(['replace = find(matrix',relation,'old);']);
end

nreplace = length(replace);
matrix_new = matrix;

%% replace values

if nreplace>0
   matrix_new(replace) = new*ones(1,nreplace);
end 