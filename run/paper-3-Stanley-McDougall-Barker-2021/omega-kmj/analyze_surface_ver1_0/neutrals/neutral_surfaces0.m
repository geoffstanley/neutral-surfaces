 function [sns,tns,pns,dsns,dtns,dpns] = neutral_surfaces0(s,t,p,gamma,glevels)
 
%         [sns,tns,pns,dsns,dtns,dpns) = neutral_surfaces0(s,t,p,gamma,glevels)   
%
%         For a cast of hydrographic data {s,t,p} that has been labelled with gamma,
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



n2= 2;  ptol = 1d-3; n = length(s); ng = length(glevels);


%         detect error condition & adjust cast to avoid an exact crossing

in_error = 0;
for k = 1:n
   if gamma(k)<0.d0, in_error = 1; end
      for ig = 1:ng
         if gamma(k)==glevels(ig), gamma(k) = gamma(k)+1.d-10; end
   end 
end
   
if in_error==1
    gamma,glevels
   ERROR = ' in neutral-surfaces.f : missing gamma value'
   return
end

        
%         loop over the surfaces

ierr = 0;

        for ig = 1:ng
%         find the intervals of intersection
          nint = 0;
          for k = 1:n-1
            gmin = min(gamma(k),gamma(k+1));
            gmax = max(gamma(k),gamma(k+1));
            if gmin<=glevels(ig)&glevels(ig)<=gmax
              nint = nint+1; intz(nint) = k;
            end
          end

        
%         find point(s) of intersection

          if nint==0
            sns(ig) = nan; tns(ig) = nan; pns(ig) = nan;
            dsns(ig) = 0; dtns(ig) = 0; dpns(ig) = 0;
          else

%         choose the central interval

            if mod(nint,2)==0&intz(1)>=n/2
              int_middle = (nint+2)/2;
            else
              int_middle = floor((nint+1)/2);
            end

%         loop over all intersections

            for i_int = 1:nint
              k = intz(i_int);
%         coefficients of a quadratic for gamma
              [rho,rho_s,rho_t,rho_p] =  eosall_from_ct(s(k),t(k),p(k));
              alfa_l = -rho_t/rho; beta_l = rho_s/rho;

              [rho,rho_s,rho_t,rho_p] = eosall_from_ct(s(k+1),t(k+1),p(k+1));
              alfa_u = -rho_t/rho; beta_u = rho_s/rho;
 
              smid = (s(k)+s(k+1))/2; tmid = (t(k)+t(k+1))/2; pmid = (p(k)+p(k+1))/2;
              [rhomid,rho_s,rho_t,rho_p] = eosall_from_ct(smid,tmid,pmid);
              alfa_mid = -rho_t/rhomid; beta_mid = rho_s/rho;

              dels = s(k+1)-s(k); delth = t(k+1)-t(k);

              pl = p(k); pu = p(k+1); delp = pu-pl; delp2 = delp*delp;

              bden = rhomid*(beta_mid*dels-alfa_mid*delth);

              if abs(bden)<=1.d-6, bden = 1.d-6; end

              bmid = (gamma(k+1)-gamma(k))/bden;

%          coefficients

              a = dels*(beta_u-beta_l)-delth*(alfa_u-alfa_l);
              a = (a*bmid*rhomid)/(2*delp2);

              b = dels*(pu*beta_l-pl*beta_u) - delth*(pu*alfa_l-pl*alfa_u);
              b = (b*bmid*rhomid)/delp2;

              c = dels*(beta_l*(pl-2*pu)+beta_u*pl) -  ...  
                                    delth*(alfa_l*(pl-2*pu)+alfa_u*pl); 
              c = gamma(k) + (bmid*rhomid*pl*c)/(2*delp2);
              c = c - glevels(ig);

              delta =  b*b-4*a*c;

%         solve the quadratic

              if a~=0&bden~=1d-6&delta>=0
                q = -(b+sign(b)*sqrt(delta))/2;
                pns1 = q/a; pns2 = c/q;
                if pns1>=p(k)-ptol&pns1<=p(k+1)+ptol
                  pns(ig) = min(p(k+1),max(pns1,p(k)));
                elseif pns2>=p(k)-ptol&pns2<=p(k+1)+ptol
                  pns(ig) = min(p(k+1),max(pns2,p(k)));
                else
                  rg = (glevels(ig)-gamma(k))/(gamma(k+1)-gamma(k));
                  pns(ig) = p(k)+rg*(p(k+1)-p(k));
                end
              else
                rg = (glevels(ig)-gamma(k))/(gamma(k+1)-gamma(k));
                pns(ig) = p(k)+rg*(p(k+1)-p(k));
              end

              [sns(ig),tns(ig)] = stp_interp(s,t,p,pns(ig));

              
%         write multiple values to file

            if nint>1

%              if(ierr.eq.0) then
%                ierr = 1
%                call system('rm -f ns-multi.dat')
%                lun = 21
%                open(lun,file='ns-multi.dat',status='unknown')
%              end if
%              if(i_int.eq.1) write(lun,*) ig,nint
%              write(lun,*) sns(ig),tns(ig),pns(ig)


%
%         find median values and errors
%

              if i_int==1
                sns_top = sns(ig); tns_top = tns(ig); pns_top = pns(ig);
              end

              if i_int==int_middle
                sns_middle = sns(ig); tns_middle = tns(ig); pns_middle = pns(ig);
              end

              if i_int==nint
                if (pns_middle-pns_top)>(pns(ig)-pns_middle)
                  dsns(ig) = sns_middle-sns_top;
                  dtns(ig) = tns_middle-tns_top;
                  dpns(ig) = pns_middle-pns_top;
                else
                  dsns(ig) = sns(ig)-sns_middle;
                  dtns(ig) = tns(ig)-tns_middle;
                  dpns(ig) = pns(ig)-pns_middle;
                end
                sns(ig) = sns_middle;
                tns(ig) = tns_middle;
                pns(ig) = pns_middle;
              end

            else
              dsns(ig) = 0; dtns(ig) = 0; dpns(ig) = 0;
            end
        end
    end
end

sns = sns(:); tns = tns(:); pns = pns(:);
dsns = dsns(:); dtns = dtns(:); dpns = dpns(:);


return