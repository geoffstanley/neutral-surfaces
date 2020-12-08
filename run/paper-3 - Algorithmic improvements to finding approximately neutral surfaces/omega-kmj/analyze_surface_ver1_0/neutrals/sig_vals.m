function [sig1,sig2] =  sig_vals(s1,t1,p1,s2,t2,p2)

%%	Computes the sigma values of two neighbouring bottles w.r.t. the mid pressure
%%
%%  Usage :         [sig1,sig2] =  sig_vals(s1,t1,p1,s2,t2,p2)
%%
%%	Input :			s1,s2         bottle salinities
%%					t1,t2         bottle conservative temperatures
%%					p1,p2         bottle pressures
%%
%%	Output :		sig1,sig2     bottle potential density values
%%
%%	Units :			salinity      psu (IPSS-78)
%%					temperature   degrees C (IPTS-68)
%%					pressure      db
%%                  density       kg m-3

%%  DRJ on 17/06/03


pmid = 0.5*(p1+p2);

sig1 = rho_from_ct(s1,t1,pmid)-1000;

sig2 = rho_from_ct(s2,t2,pmid)-1000;


return