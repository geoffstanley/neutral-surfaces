function [sns,tns,pns] = depth_ntp(s0,t0,p0,s,t,p)

%%  Find the position where the neutral tangent plane passing through a bottle 
%%  intersects a neighbouring cast 
%%
%%  Usage :        [sns,tns,pns] = depth_ntp(s0,t0,p0,s,t,p)
%%
%%  Input :        s0    the bottle salinity
%%                 t0    the bottle conservative temperature
%%                 p0    the bottle pressure
%%                 s     vector of cast salinities
%%                 t     vector of cast conservative temperatures
%%                 p     vector of cast pressures
%%
%%  Output :       sns   salinity of the ntp intersection
%%                 tns   conservative temperature of the intersection
%%                 pns   pressure of the intersection
%%
%%  Units :        salinities	  psu (IPSS-78)
%%                 temperatures   degrees C (IPS-90)
%%                 pressures      db

%%  DRJ on 17/06/03


n = length(s); e = zeros(n,1);

%		find the bottle pairs containing a crossing

ncr = 0;
for k = 1:n
  [sigl,sigu] = sig_vals(s0,t0,p0,s(k),t(k),p(k));
  e(k) = sigu-sigl;
  if k>1    
    if e(k-1)==0                   %  an exact crossing at the k-1 bottle
      ncr = ncr+1; sns = s(k-1); tns = t(k-1); pns = p(k-1);
    elseif e(k)*e(k-1)<0           %  a crossing between k-1 and k bottles
      ncr = ncr+1;
                                   %  some Newton-Raphson iterations
      pc0 = p(k-1)-e(k-1)*(p(k)-p(k-1))/(e(k)-e(k-1));
      iter = 0; success = 0;
      while success==0
        iter = iter+1;
        [sc0,tc0] = stp_interp([s(k-1),s(k)],[t(k-1),t(k)],[p(k-1),p(k)],pc0);
        [sigl,sigu] = sig_vals(s0,t0,p0,sc0,tc0,pc0);
		ec0 = sigu-sigl;
        p1 = 0.5*(p(k-1)+pc0);
	    ez1 = (e(k-1)-ec0)/(pc0-p(k-1));
		p2 = 0.5*(pc0+p(k));
		ez2 = (ec0-e(k))/(p(k)-pc0);
		r = (pc0-p1)/(p2-p1);
		ecz_0 = ez1+r*(ez2-ez1);
        if iter==1
          ecz0 = ecz_0;
	    else
          ecz0 = -(ec0-ec_0)/(pc0-pc_0);
          if ecz0==0, ecz0 = ecz_0; end
	    end
        pc1 = pc0+ec0/ecz0;
                                   %  strategy when iteration jumps out of inteval
        if pc1<p(k-1) | pc1>p(k)
          indsp = find(finite(s));
          [sns,tns,pns,niter] = e_solve(s(indsp),t(indsp),p(indsp),e(indsp),k,s0,t0,p0);
		  if pns<p(k-1) | pns>p(k)
            disp('****ERROR**** in depth-ntp.m')
            return
          else
            success = 1;
          end
        else
                                   %  test accuracy of the iterate
          eps = abs(pc1-pc0);
          if abs(ec0)<=5.d-5&eps<=5.d-3
            sns = sc0; tns = tc0; pns = pc0;
            success = 1; niter = iter;
          elseif iter>10
            [sns,tns,pns,niter] = e_solve(s,t,p,e,k,s0,t0,p0);
		    success = 1;
	      else
            pc_0 = pc0; ec_0 = ec0; pc0 = pc1;
		    success = 0;
          end
        end
      end
    end
  end
  if k==n&e(k)==0                  %  the last bottle
    ncr = ncr+1;
	sns = s(k); tns = t(k); pns = p(k);
  end
end

                                   %  multiple and no crossings

if ncr==0
  if e(1)>0                                 %  outcropping
    sns = -99.1; tns = -99.1; pns = -99.1;
  else                                      %  undercropping
    sns = -99.2; tns = -99.2; pns = -99.2;
  end
elseif ncr>=2                               %  multiple crossings
	sns = -99.3; tns = -99.3; pns = -99.3;
end

    
return