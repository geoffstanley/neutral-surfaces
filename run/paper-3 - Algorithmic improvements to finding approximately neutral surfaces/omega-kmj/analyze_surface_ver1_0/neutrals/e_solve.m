function [sns,tns,pns,iter] = e_solve(s,t,p,e,k,s0,t0,p0)

%%	Find the zero of the e function using a bisection method
%%
%%  Usage :         [sns,tns,pns,iter] = e_solve(s,t,p,e,k,s0,t0,p0)
%%
%%	Input :			s		array of cast salinities
%%				t		array of cast conservative temperatures
%%				p		array of cast pressures
%%				e		array of cast e values
%%				k		interval (k-1,k) contains the zero
%%				s0		the bottle salinity
%%				t0		the bottle conservative temperature
%%				p0		the bottle pressure
%%
%%	Output :		sns		salinity of e zero
%%				tns		conservative temperature of e zero
%%				pns		pressure of e zero
%%
%%	Units :			salinities	psu (IPSS-78)
%%				temperatures	degrees C (ITS-90)
%%				pressures	db

%%  DRJ on 17/06/03


pl = p(k-1); el = e(k-1); pu = p(k); eu = e(k);

iter = 0; success = 0; 

while success==0
  iter = iter+1; pm = 0.5*(pl+pu);
  [sm,tm] = stp_interp([s(k-1),s(k)],[t(k-1),t(k)],[p(k-1),p(k)],pm);
  [sigl,sigu] = sig_vals(s0,t0,p0,sm,tm,pm);
  em = sigu-sigl;
  if el*em<0
    pu = pm; eu = em;
  elseif em*eu<0
    pl = pm; el = em;
  elseif em==0
    sns = sm; tns = tm; pns = pm; success = 1;
  end
  if success==0
    if abs(em)<=5.d-5&abs(pu-pl)<=5.d-3
      sns = sm; tns = tm; pns = pm; success = 1;
    elseif iter<=20
      success = 0;
    else
      disp('WARNING in e-solve.f')
	  disp(['iter: ', int2str(iter), '  em: ', num2str(em), '  dp: ', num2str(abs(pu-pl))])
      sns = -99; tns = -99; pns = -99; success = 0; return
    end
   end
end


return