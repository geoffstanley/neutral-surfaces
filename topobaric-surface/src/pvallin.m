function y = pvallin(f,x)
%PVALLIN  Evaluate polynomial, linearly extrapolating outside domain.
%
%
% y = pvallin(f,x)
% evaluates the polynomial f at location site x if x is in the domain of f,
% or linearly extrapolates f to x if x is outside the domain of f. The
% domain of f is f(1) to f(2). If f(1) <= x <= f(2), then
% y = f(3)*(x-f(1))^(N-1) + ... + f(N+1)*(x-f(1)) + f(N+2,l).
% If x < f(1), y is linearly extrapolated from the left end of f.
% If f(2) < x , y is linearly extrapolated from the right end of f.
%
%
% --- Input:
% f [N+2, 1]: polynomial of order N, the valid in a domain f(1) to f(2) and
%             with coefficients f(3:N+2) starting with the highest order
%             monomial
% x [1  , 1]: evaluation site
%
%
% --- Output:
% y [1  , 1]: f evaluated at x

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2019  Geoff Stanley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com
% Version   : 2.1.0
%
% Modified by : --
% Date        : --
% Changes     : --

% Np2-2 is the polynomial's order. (e.g. Np2==4 for linear)
Np2 = size(f,1);

% Evaluate local coordinate:
distL = x - f(1);

if distL < 0
    % Interpolate linearly from left.
    %   val@left + slope@left * (dist from left)
    y = f(end) + f(end-1) * distL;
    
else
    
    xp = 1;
    y = f(end);
    distR = x - f(2);
    
    if distR > 0
        % Interpolate linearly from right.
        
        % Change to evaluate at right, and change to local coordinates
        x = f(2) - f(1);
        
        % Evaluate polynomial AND slope at right:
        slope = 0;
        for k = 1:Np2-3
            slope = slope + k * f(Np2-k) * xp; % xp == (f(2)-f(1)) ^ (k-1)
            xp = xp .* x;
            y = y + f(Np2-k) * xp; % xp == (f(2)-f(1)) ^ k
        end
        
        %   val@right + slope@right * (dist from right)
        y = y + slope * distR;
        
    else
        
        % Evaluate polynomial:
        for k = 1:Np2-3
            xp = xp .* distL; % xp == (x-f(1)) ^ k
            y = y + f(Np2-k) * xp;
        end
        
    end
    
end

