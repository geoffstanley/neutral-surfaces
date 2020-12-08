function [alpha,beta,aonb,T_b] = ab_from_ct(s,ct,p,Z2P)  % GJS added Z2P

%           Calculate alpha, beta, alpha/beta and T_b
%
% Usage:    [alpha,beta,aonb,T_b] = ab_from_ct(s,ct,p)
%
%           Calculate the thermal expansion coefficient (alpha), the saline
%           contraction coefficient (beta), alpha over beta and the
%           thermobaric coefficient (T_b = (alpha/beta)_p * beta)
%
%
% Input:    s           salinity
%           ct          conservative temperature
%           p           pressure
%
% Output:   alpha       thermal expansion coefficient
%           beta        saline contraction coefficient
%           aonb        alpha over beta
%           T_b         thermobaric parameter
%                 
% Calls:    eosall_from_ct.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

%if ~(nargin == 3) % GJS removed
%    error('ab_from_ct.m: requires 3 input arguments') % GJS removed
%end  % GJS removed

if nargin < 4; Z2P = 1; end  % GJS
fac = 0.01 * Z2P; % GJS
p = p * Z2P; % GJS

%% calculate density and its derivatives

[rho,rho_s,rho_ct,rho_p] = eosall_from_ct(s,ct,p); %#ok
[rho1,rho_s1,rho_ct1,rho_p1] = eosall_from_ct(s,ct,p-fac); %#ok
[rho2,rho_s2,rho_ct2,rho_p2] = eosall_from_ct(s,ct,p+fac); %#ok

%% calculate alpha, beta and aonb

alpha = - rho_ct ./ rho;
beta = rho_s ./rho;
aonb = alpha ./beta;

alpha1 = - rho_ct1 ./ rho1;
beta1 = rho_s1 ./rho1;
aonb1 = alpha1 ./beta1;

alpha2 = - rho_ct2 ./ rho2;
beta2 = rho_s2 ./rho2;
aonb2 = alpha2 ./beta2;

aonb_p = (aonb2 - aonb1) ./ (2*fac);
T_b = aonb_p .* beta;