function bfs_wet_codegen(nk, ni, nj, Xvec, r, OPTS)
%BFS_WET_CODEGEN  Create MEX function for bfs_wet
%
%
% bfs_wet_codegen(nk, ni, nj, false, r)
% runs codegen on bfs_wet.m, appropriate for a grid of ni by nj points in
% the horizontal and nk points in the vertical.  bfs_wet is compiled with
% its final argument (also called r) provided iff r is true.
%
% bfs_wet_codegen(nk, ni, nj, true, r)
% specifies that X in bfs_wet.m is just a vector: X(k) specifies the
% pressure or depth of all grid points having vertical index k. Use this
% for simple Z-level models (not hybrid coordinate models).
%
% bfs_wet_codegen(..., OPTS)
% overrides default verbosity by OPTS.VERBOSE and file output by
% OPTS.FILE_ID.

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

% Set defaults
VERBOSE = 1; % verbose mode
FILE_ID = 1; % output to MATLAB terminal

% Override defaults
if nargin == 6 && isstruct(OPTS)
    if isfield(OPTS, 'VERBOSE')
        VERBOSE = OPTS.VERBOSE;
    end
    if isfield(OPTS, 'FILE_ID')
        FILE_ID = OPTS.FILE_ID;
    end
end

V = filesep();

% Get info about which functions are on the path
name = 'bfs_wet';
if r
    name_mex = [name '_one_mex'];
else
    name_mex = [name '_all_mex'];
end
file_mex = dir(which(name_mex));
file_mat = dir(which(name));
assert(~isempty(file_mat), ['Cannot locate ' name '.m']);
which_eos   = which('eos');
file_eos    = dir(which_eos);

% Test values
s = 34.5;
t = 3;
x = 1000;
m = eos(s, t, x);

% Create textual identifier for this build of the MEX function.
build_text = sprintf('%s_k%d_i%d_j%d_%dD_m=%.59e', name_mex, nk, ni, nj, (1-Xvec)*2+1, m);
fileName_build_text = [file_mat.folder V name '_info.txt'];
fileID_build_text = fopen(fileName_build_text, 'rt');


if isempty(file_mex) ... % No mex file yet
        || file_mex.datenum < max([file_mat.datenum, file_eos.datenum]) ... % MEX is too old
        || (fileID_build_text < 0) ... % No text file specification
        || (fileID_build_text >= 0 && ~strcmp(fgetl(fileID_build_text), build_text) && ~fclose(fileID_build_text)) % Requested parameters do not match those recorded in text file. Also close the text file.
    
    if VERBOSE
        mytic = tic;
        fprintf(FILE_ID, 'Compiling MEX for %s, with\n', name);
        fprintf(FILE_ID, ' %s in %s\n', read_function_name(which_eos), which_eos);
        fprintf(FILE_ID, ' eos(%g,%g,%g) = %e\n', s, t, x, m);
    end
    
    % Note: resulting MEX is actually faster with variable size arrays
    % (using vs = true, below)
    vs = true;
    t_SppX   = coder.typeof(0, [2, nk-1, ni, nj], [false, vs, vs, vs]);
    if Xvec
        t_X = coder.typeof(0, [nk, 1], [vs, vs]);
    else
        t_X = t_S;
    end
    nij = ni * nj;
    t_x    = coder.typeof(0, [ni, nj], [vs, vs]);
    t_A    = coder.typeof(0, [nij, 4], [vs, false]);
    t_BotK = coder.typeof(uint16(0), [ni, nj], [vs, vs]);
    t_q    = coder.typeof(0, [nij, 1], [vs, vs]);
    
    if r        %(SppX,   TppX,   X,   s,   t,  x,X_TOL, A,   BotK, r,  qu)
        args = {t_SppX, t_SppX, t_X, t_x, t_x, t_x, 0, t_A, t_BotK, 0, t_q};
    else
        args = {t_SppX, t_SppX, t_X, t_x, t_x, t_x, 0, t_A, t_BotK, [], t_q};
    end
    
    % Configure MEX for speed.
    mexconfig = coder.config('mex');
    mexconfig.ExtrinsicCalls = false;
    mexconfig.ResponsivenessChecks = false;
    mexconfig.IntegrityChecks = false;
    
    % Compile the MEX function
    codegen(name, '-args', args, '-config', mexconfig, '-o', [file_mat.folder V name_mex], '-d', [file_mat.folder V 'codegen' V 'mex' V name]);
    clear(name_mex)  % Ensure new mex file gets used
    
    % Save textual identifier for this MEX function
    fileID_build_text = fopen(fileName_build_text, 'wt'); % overwrite
    fprintf(fileID_build_text, build_text);
    fclose(fileID_build_text);
    
    if VERBOSE
        fprintf(FILE_ID, ' ...done compiling, time %.2f\n', toc(mytic));
    end
end


