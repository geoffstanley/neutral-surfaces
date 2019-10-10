function [x, perturbed] = perturb_until_unique(x)
%PERTURB_UNTIL_UNIQUE  Perturb data to remove duplicate values.
%
%
% [x, perturbed] = perturb_until_unique(x)
% adds small numbers to duplicate elements of x until all non-NaN elements
% of x are unique. If all non-NaN elements of x are unique to begin with,
% perturbed is false; otherwise, perturbed is true. 
%
%
% --- Input:
% x: an array of any size and dimension
% 
%
% --- Output:
% x: the same size as input x, but perturbed to have unique values
% perturbed [1, 1]: true if output x differs from input x, else false. 

% --- Copyright:
% Copyright 2019 Geoff Stanley
%
% This file is part of Topobaric Surface.
% 
% Topobaric Surface is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
% 
% Topobaric Surface is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with Topobaric Surface.  If not, see
% <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au 
% Email     : geoffstanley@gmail.com
% Version   : 1.0
%
% Modified by : --
% Date        : --
% Changes     : --

perturbed = false;
for rep = 0:100
    [datsort, I] = sort(x(:));
    duplicates = find(diff(datsort) == 0);
    if isempty(duplicates)
        break
    else
        for i = duplicates(:)'
            q = x(I(i));
            x(I(i)) = q + eps(q) ;
        end
        perturbed = true;
    end
end
end