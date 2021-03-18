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

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


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