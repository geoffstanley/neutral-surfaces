function k = indx(x,x0)
	  
%%  Find the index of a scalar in a monotonically increasing array
%%  
%%  Usage :      k = indx(x,x0)
%%
%%  Input :      x        array of increasing values
%%               x0       scalar
%%  
%%  Output :     k        when x(k) <= x0 < x(k+1), or
%%               n-1      when x0 = x(n)
%%               error    otherwise

%%  DRJ on 17/06/03
  

n = length(x);

if x(1)<x0&x0<x(n)
  inds_x = find(x>=x0);
  k = inds_x(1)-1;
elseif x0==x(1)
  k = 1;
elseif x0==x(n)
  k = n-1;
else
  disp('ERROR in indx.m : out of range')
  x, x0
end


return