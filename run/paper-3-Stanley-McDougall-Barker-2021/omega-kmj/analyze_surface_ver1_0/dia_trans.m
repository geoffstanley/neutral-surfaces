function [transport,transport_sum] = dia_trans(dia_vel,e1t,e2t)

%           Calculate dia-surface transport and integrate over the domain
%
% Usage:    [transport,transport_sum] = dia_trans(dia_vel,e1t,e2t)
%
%           Calculate dia-surface transport and integrate over the domain
%
% Input:    dia_vel         dia-surface velocity
%           e1t             distance between grid points in latitudinal
%                           direction
%           e2t             distance between grid points in longitudinal
%                           direction
%
% Output:   transport       dia-surface transport
%           transport_sum   integral of dia-surface transport
%
% Calls:    change.m
%
% Units:    dia_vel         m/s
%           e1t             m
%           e2t             m
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin==3)
  error('dia_trans.m: requires 3 input arguments')
end 

%% initialize

[yi,xi] = size(e1t);
transport = nan(1,yi,xi);

%% calculate diapycnal transports

transport(1,:,:) = e1t .* e2t .* squeeze(dia_vel(1,:,:));

%% get rid of Infs
transport(1,:,:) = change(transport(1,:,:),'==',Inf,NaN);
transport(1,:,:) = change(transport(1,:,:),'==',-Inf,NaN);

%% integral of diapycnal transport

transport_sum = nansum(transport(:));