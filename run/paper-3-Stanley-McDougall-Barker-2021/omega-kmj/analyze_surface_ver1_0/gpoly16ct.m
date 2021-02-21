function [gamma,gnum,gden] = gpoly16ct(s,ct)

%           Calculate neutral density (gamma^rf) through a rational function
%
% Usage:    [gamma,gnum,gden] = gpoly16ct(s,ct)
%
%           Calculate neutral density (gamma^rf) through a rational function 
%           in terms of salinity and conservative temperature according to 
%           'The material derivative of neutral density, McDougall and 
%           Jackett, Journal of Marine Research, 2005'.
%
% Input:    s           salinity
%           ct          conservative temperature
%
% Output:   gamma       gamma^rf       
%           gnum        numerator  
%           gden        denominator
%
% Calls:          
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           gamma                       kg/m^3
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% check input arguments

if ~(nargin ==2)
  error('gpoly16ct.m: requires 2 input arguments')
end 

zcheck(s,ct)

%% calculate multiples and square roots of ct and s

ct2 = ct.*ct; s2 = s.*s; sqrts = sqrt(s);

%% calculate the numerator

gnum =        1.0022048243661291e003 + ...
        ct.*( 2.0634684367767725e-001 + ...
        ct.*( 8.0483030880783291e-002 + ...
        ct.*(-3.6670094757260206e-004))) + ...
         s.*(-1.4602011474139313e-003 + ...
        ct.* -2.5860953752447594e-003) + ...
        s2.* -3.0498135030851449e-007;
     
%% calculate the denumerator
gden =         1.0 + ...
        ct.*( 4.4946117492521496e-005 + ...
        ct.*( 7.9275128750339643e-005 + ...
        ct.*(-1.2358702241599250e-007 + ...
        ct.*(-4.1775515358142458e-009)))) + ...
         s.*(-4.3024523119324234e-004 + ...
        ct.*( 6.3377762448794933e-006 + ...
       ct2.* -7.2640466666916413e-010) + ...
     sqrts.*(-5.1075068249838284e-005 + ...
       ct2.* -5.8104725917890170e-009));
   
%% calculate gamma^rf

gamma = (gnum./gden) - 1000;