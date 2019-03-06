function bw = select_conn_cpt(bw, WRAP, REF_IJ)
%SELECT_CONN_CPT  Select only one connected region of the ocean.
%
%
% bw = select_conn_cpt(bw, WRAP) 
% selects the largest connected component of a binary image bw that is
% periodic in its first dimension iff WRAP(1) is true, and periodic in its
% second dimension if WRAP(2) is true.
%
% bw = select_conn_cpt(bw, WRAP, REF_IJ) 
% as above if REF_IJ is empty. If REF_IJ is a 2 element vector, the
% connected component of bw containing the pixel at row REF_IJ(1) and
% column REF_IJ(2) is selected.
%
%
% --- Input:
% bw [nx, ny]: binary image [logical]
% WRAP [1, 2]: periodicity of bw [logical]
% REF_IJ [] or [1, 2]: select the largest component or the component
%                      containing a given pixel
%
%
% --- Output:
% bw [nx, ny]: binary image containing one connected component [logical]
%
%
% --- Requirements:
% bwconncomp, labelmatrix - Image Processing toolbox
% CC2periodic - https://www.mathworks.com/matlabcentral/fileexchange/66079

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

CC = CC2periodic(bwconncomp(bw, 4), WRAP, 'CC');
bw(:) = false;
if nargin < 3 || isempty(REF_IJ)
    % Select the largest region
    [~,id] = max(cellfun(@(idx) length(idx), CC.PixelIdxList));
else
    % Select the chosen region
    i = REF_IJ(1);
    j = REF_IJ(2);
    L = labelmatrix(CC);
    id = L(i,j);
end
bw(CC.PixelIdxList{id}) = true;