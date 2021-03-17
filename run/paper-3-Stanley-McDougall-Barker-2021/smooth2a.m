function [matrixOut,nrmlize] = smooth2a(matrixIn,N,std,preserveNaN,wrap)
% Smooths 2D array data.  Ignores NaN's when averaging.
%
%function matrixOut = smooth2a(matrixIn,Nr,Nc,preserveNaN,wrapr,wrapc)
% 
% This function smooths the data in matrixIn using a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
% element (i,j) by the mean of the rectange centered on (i,j).  Any NaN
% elements are ignored in the averaging. If element (i,j) of matrixIn is
% NaN, it will be NaN in matrixOut if preserveNaN is true; otherwise it
% takes on the average centred at (i,j) and will be NaN only if the entire
% (2*Nr-1)-by-(2*Nc-1)  neighbourhood is NaN. At the edges of the matrix,
% where you cannot build a full rectangle, as much of the rectangle that
% fits on your matrix is used (similar to the default on Matlab's builtin
% function "smooth"). This behaviour is modified by the logical inputs
% "wrapr" and "wrapc" which allow circular wrapping in the rows and
% columns, respectively.
% 
% "matrixIn": original matrix 
% "Nr": number of points used to smooth rows
% "Nc": number of points to smooth columns.  Default: Nr
% "preserveNaN" : boolean indicating whether NaN's in input remain NaN's in
% output. Default: true. 
% "wrapr" : boolean indicating whether to wrap data
% in rows. Default: false. 
% "wrapc" : boolean indicating whether to wrap
% data in columns. Default: false.
% 
% "matrixOut": smoothed version of original matrix
% 
% 
% 	Written by Greg Reeves, March 2009. Division of Biology Caltech
% 
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004 Applied
% 	Research Laboratory Penn State University
% 
% 	Developed from code written by Olof Liungman, 1997 Dept. of
% 	Oceanography, Earth Sciences Centre G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se
%
%   Modifications by Geoff Stanley, September 2016
%   Added std, preserveNaN, and wrap options.
%   E-mail: geoff.stanley@physics.ox.ac.uk
%


%
% Initial error statements and definitions
%
narginchk(2,8);

Nr = N(1);
if isscalar(N)
    Nc = Nr;
else
    Nc = N(2);
end

if nargin < 3 || isempty(std)
    std = 0; % box smoothing.
end
stdr = std(1);
if isscalar(std)
    stdc = stdr;
else
    stdc = std(2);
end

if nargin < 4 || isempty(preserveNaN)
    preserveNaN = true;
end

if nargin < 5 || isempty(wrap)
    wrapr = false;
    wrapc = false;
elseif isscalar(wrap)
    wrapr = wrap;
    wrapc = wrap;
else
    wrapr = wrap(1);
    wrapc = wrap(2);
end

%
% Add data to either end of the matrix if wrapping is requested
%
if wrapr
    matrixIn = cat(1, matrixIn(end-Nr+1:end,:), matrixIn, matrixIn(1:Nr,:));
end
if wrapc
    matrixIn = cat(2, matrixIn(:,end-Nc+1:end), matrixIn, matrixIn(:,1:Nc));
end

%
% Building matrices that will compute running sums.  The left-matrix, eL,
% smooths along the rows.  The right-matrix, eR, smooths along the columns.
% You end up replacing element "i" by the mean of a (2*Nr+1)-by- (2*Nc+1)
% rectangle centered on element "i".
%
[row,col] = size(matrixIn);

if stdr > 0
    eL = spdiags(repmat(gaussfilter(2*Nr+1, stdr), row, 1), (-Nr:Nr), row, row);
else
    eL = spdiags(ones(row,2*Nr+1),(-Nr:Nr),row,row);
end

if stdc > 0
    eR = spdiags(repmat(gaussfilter(2*Nc+1, stdc), col, 1), (-Nc:Nc), col, col);
else
    eR = spdiags(ones(col,2*Nc+1),(-Nc:Nc),col,col);
end

%
% Setting all "NaN" elements of "matrixIn" to zero so that these will not
% affect the summation.  (If this isn't done, any sum that includes a NaN
% will also become NaN.)
%
A = isnan(matrixIn);
matrixIn(A) = 0;

%
% For each element, we have to count how many non-NaN elements went into
% the sums.  This is so we can divide by that number to get a mean.  We use
% the same matrices to do this (ie, "eL" and "eR").
%
nrmlize = eL*(~A)*eR;
if preserveNaN
    nrmlize(A) = NaN;
end

%
% Actually taking the mean.
%
matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;

%
% Throw out the edges that were added for wrapping
%
if wrapr && wrapc
    matrixOut = matrixOut(Nr+1:end-Nr, Nc+1:end-Nc);
elseif wrapr
    matrixOut = matrixOut(Nr+1:end-Nr,:);
elseif wrapc
    matrixOut = matrixOut(:,Nc+1:end-Nc);
end


function h = gaussfilter(n, std)
% Construct one-dimensional gaussian filter, of size n and standard deviation std
% Create the underlying coordinate
x = linspace(-(n-1)/2, (n-1)/2, n);

% Create the Gaussian mask
h = exp(-(x.*x) / (2*std*std));

% Normalise the mask
h = h/sum(h);
