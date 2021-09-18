function ns_add_to_path()
%NS_ADD_TO_PATH  Add Neutral Surfaces and its libraries to MATLAB's search path
%
% This file should exist in the root directory of Neutral Surfaces

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


V = filesep(); % /  or  \  depending on OS

% Get path to this file.
PATH_PROJECT = [fileparts(mfilename('fullpath')) V];

% Add Neutral Surfaces libraries and scripts to MATLAB's path
addpath(PATH_PROJECT);

addpath([PATH_PROJECT 'lib']);
addpath([PATH_PROJECT 'lib' V 'bfs']);
addpath([PATH_PROJECT 'lib' V 'dat']);
addpath([PATH_PROJECT 'lib' V 'file-ops']);
addpath([PATH_PROJECT 'lib' V 'fzero_cg']);
addpath([PATH_PROJECT 'lib' V 'gsf']);
addpath([PATH_PROJECT 'lib' V 'ntp-bottle-to-cast']);

FEX_DIR = dir([PATH_PROJECT 'fex']); % rather than genpath, only go one level deep: avoid codegen folders
for i = 1 : length(FEX_DIR)
    if FEX_DIR(i).isdir && FEX_DIR(i).name(1) ~= '.' % Ignore . and .. and hidden files
        addpath([PATH_PROJECT 'fex' V FEX_DIR(i).name]);
    end
end

addpath([PATH_PROJECT 'potential-density-surface']);
addpath([PATH_PROJECT 'specific-volume-anomaly-surface']);

addpath([PATH_PROJECT 'omega-surface']);

addpath([PATH_PROJECT 'topobaric-surface']);
addpath([PATH_PROJECT 'topobaric-surface' V 'src']);
