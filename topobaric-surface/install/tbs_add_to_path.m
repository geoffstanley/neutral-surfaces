function topobaric_surface_add_to_path()
%TOPOBARIC_SURFACE_ADD_TO_PATH  Add Topobaric Surface and its libraries to
%MATLAB's search path.
%
% This file should exist in the 'run' sub-folder of Topobaric Surface.

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

% Variable names are intentially chosen to be obtuse, hoping they won't
% conflict with existing variables in the base workspace. 

V = filesep;

% Get path to this file. 
PATH_LOCAL = fileparts(mfilename('fullpath'));

% Get path to one directory up, containing all of Topobaric Surface
PATH_PROJECT = PATH_LOCAL(1 : find(PATH_LOCAL == V, 1, 'last'));

addpath([PATH_PROJECT 'run']);
addpath([PATH_PROJECT 'src']);
addpath([PATH_PROJECT 'lib']);
FEX_DIR = dir([PATH_PROJECT 'fex']); % rather than genpath, only go one level deep: avoid codegen folders
for i = 1 : length(FEX_DIR)
    if FEX_DIR(i).isdir && FEX_DIR(i).name(1) ~= '.' % Ignore . and .. and hidden files
        addpath([PATH_PROJECT 'fex' V FEX_DIR(i).name]);
    end
end