function zcheck(x1,x2,x3,x4,x5)

%%   checks arrays have the same dimensions
%%
%%   usage        : zcheck(x1,x2,x3,x4,x5)
%%
%%   x1-5         : up to 5 arrays

%%   DRJ on 31/5/02


n1 = prod(size(x1)); n2 = prod(size(x2));

message = 'ERROR: input array dimensions are not identical';

if n1~=n2, error(message), end

if nargin>=3
  n3 = prod(size(x3));
  if n1~=n3, error(message), end	
  if nargin>=4
    n4 = prod(size(x4));
    if n1~=n4, error(message), end
    if nargin>=5
      n5 = prod(size(x5));
      if n1~=n5, error(message), end	
    end
  end
end


return