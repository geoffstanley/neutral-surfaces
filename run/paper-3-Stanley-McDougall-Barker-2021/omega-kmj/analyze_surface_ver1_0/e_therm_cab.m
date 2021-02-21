function [e_cab,e_therm] = e_therm_cab(sns,ctns,pns,n2_ns,g,e1t,e2t,keyword,wrap)

%           Calculate cabbeling and thermobaricity
%
% Usage:    [cabbeling,thermobaricity] = e_therm_cab(sns,ctns,pns,n2_ns,g,e1t,e2t,keyword,wrap)
%
%           Calculate the diapycnal velocity caused by cabbeling and
%           thermobaricity
%
%
% Input:    sns         salinity on density surface
%           ctns        conservative temperature on density surface
%           pns         pressure on density surface
%           n2_ns       buoyancy frequency on density surface
%           g           gravitational acceleration
%           e1t         zonal scale factor
%           e2t         meridional scale factor
%           keyword     'op' for values on tracer points
%                       'bp' for values between tracer points
%           wrap        'none'
%                       'long' 
%
% Output:   cabbeling            diapycnal velocity caused by cabbeling
%           thermobaricity       diapycnal velocity caused by thermobaricity
%                 
% Calls:    eosall_from_ct.m, grad_surf.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           n2_ns                       s^-1
%           g                           m/s^2
%           scale factors               m
%           e_cab,e_therm               m/s
%           
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin == 9)
    error('e_therm_cab.m: requires 9 input arguments')
end 

%% set parameters

fac = 0.01;
K = 1000; % isoneutral eddy diffusivity

%% make g 3-dim.

if (length(size(g)) == 2)
    [k,dummy,dummy] = size(pns); %#ok
    g_rep = repmat(g,[1 1 k]);
    g = permute(g_rep,[3 1 2]);
end

%% calculate density, alpha and beta

[rho,rho_s,rho_ct,rho_p] = eosall_from_ct(sns,ctns,pns); %#ok

alpha = - rho_ct ./ rho;
beta = rho_s ./rho;
aonb = alpha ./beta;

%% calculate dadT

[rho1,dummy,rho_ct1,dummy] = eosall_from_ct(sns,ctns-fac,pns); %#ok
[rho2,dummy,rho_ct2,dummy] = eosall_from_ct(sns,ctns+fac,pns); %#ok

alpha1 = - rho_ct1 ./ rho1;
alpha2 = - rho_ct2 ./ rho2;

dadT = (alpha2 - alpha1) ./ (2.*fac);

%% calculate dadS, dbdS

[rho1,rho_s1,rho_ct1,dummy] = eosall_from_ct(sns-fac,ctns,pns); %#ok
[rho2,rho_s2,rho_ct2,dummy] = eosall_from_ct(sns+fac,ctns,pns); %#ok

beta1 = rho_s1 ./ rho1;
beta2 = rho_s2 ./ rho2;

alpha1 = - rho_ct1 ./ rho1;
alpha2 = - rho_ct2 ./ rho2;

dbdS = (beta2 - beta1) ./ (2.*fac);
dadS = (alpha2 - alpha1) ./ (2.*fac);

%% calculate dadp, dbdp

[rho1,rho_s1,rho_ct1,dummy] = eosall_from_ct(sns,ctns,pns-fac); %#ok
[rho2,rho_s2,rho_ct2,dummy] = eosall_from_ct(sns,ctns,pns+fac); %#ok

beta1 = rho_s1 ./ rho1;
beta2 = rho_s2 ./ rho2;

alpha1 = - rho_ct1 ./ rho1;
alpha2 = - rho_ct2 ./ rho2;

dbdp = (beta2 - beta1) ./ (2.*fac);
dadp = (alpha2 - alpha1) ./ (2.*fac);

%% calculate gradients of ct and p on density surface

[ctns_x,ctns_y] = grad_surf(ctns,e1t,e2t,keyword,wrap);
[pns_x,pns_y] = grad_surf(pns,e1t,e2t,keyword,wrap);

%% calculate cabbeling and thermobaricity

cab_para = dadT + 2.* aonb .* dadS - ((alpha.^2)./(beta.^2)) .* dbdS; % cabbeling paramter
e_cab = - ((K .* g) ./ n2_ns) .* cab_para .* ((abs(ctns_x) .* abs(ctns_x)) + abs(ctns_y) .* abs(ctns_y));

thermo_para = dadp - aonb .* dbdp; % thermobaric parameter
e_therm = - ((K .* g) ./ n2_ns) .* thermo_para .* ((ctns_x .* pns_x) + (ctns_y .* pns_y));