 function [sns,tns,pns,dsns,dtns,dpns] = neutral_surfaces(s,t,p,gamma,glevels)
 
%         [sns,tns,pns,dsns,dtns,dpns) = neutral_surfaces(s,t,p,gamma,glevels)   
%
%         For a section of hydrographic data {s,t,p} that has been labelled with gamma,
%         find the salinities, conservative temperatures and pressures of neutral 
%         density surfaces specified by glevels.
%
%
% Input:          s           array of cast salinities
%                 t           array of cast conservative temperatures
%                 p           array of cast pressures
%                 gamma       array of cast gamma values
%                 glevels     array of neutral density values
%
% Output:         sns         salinity on the neutral density surfaces
%                 tns         conservative temperature on the surfaces
%                 pns         pressure on the surfaces
%                 dsns        surface salinity errors
%                 dtns        surface temperature errors
%                 dpns        surface pressure errors
%
%                 NOTE:       sns, tns and pns values of nan
%                             denotes under or outcropping
%
%                             non-zero dsns, dtns and dpns values
%                             indicates multiply defined surfaces,
%                             and file 'ns-multi.dat' contains
%                             information on the multiple solutions
%
% Units:          salinity    psu (IPSS-78)
%                 temperature degrees C (IPS-90)
%                 pressure    db
%                 gamma       kg m-3
%
%
% DRJ on 16/06/03
%
%

[m,n] = size(s); ng = length(glevels);

sns = nan*ones(ng,n); tns = sns; pns = sns;
dsns = zeros(ng,n); dtns = dsns; dpns = dsns;

for k = 1:n
    inds = find(isfinite(gamma(:,k)));
    if length(inds)>0
      [sns(:,k),tns(:,k),pns(:,k),dsns(:,k),dtns(:,k),dpns(:,k)] = ...
                        neutral_surfaces0(s(inds,k),t(inds,k),p(inds,k),gamma(inds,k),glevels);
    end
end

return