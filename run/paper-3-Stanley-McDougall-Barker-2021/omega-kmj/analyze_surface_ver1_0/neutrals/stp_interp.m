function [s0,t0] = stp_interp(s,t,p,p0)

%%  Linearly interpolate salinity and conservative temperature on a cast 
%%  to a specified pressure
%%
%%  Usage :         [s0,t0] = stp_interp(s,t,p,p0)
%%
%%  Input :         s             cast salinities
%%                  t             cast conservative temperatures
%%                  p             cast pressures
%%                  p0            pressure level
%%
%%  Output :        s0            interpolated salinity
%%                  t0            interpolated conservative temperature
%%
%%  Units :         salinity      psu (IPSS-78)
%%                  temperature   degrees C (IPS-68)
%%                  pressure      db
%%                  density       kg m-3

%%  DRJ on 17/06/03


k = indx(p,p0);

r = (p0-p(k))/(p(k+1)-p(k));

s0 = s(k) + r*(s(k+1)-s(k));

t0 = t(k) + r*(t(k+1)-t(k));


return