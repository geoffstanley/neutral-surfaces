function omega_vertsolve_codegen(nk, ni, nj, Xvec, OPTS)
%OMEGA_VERTSOLVE_CODEGEN  Create MEX function for omega_vertsolve
%
%
% omega_vertsolve_codegen(nk, ni, nj, false)
% runs codegen on omega_vertsolve.m, appropriate for a grid
% of ni by nj points in the horizontal and nk points in the vertical.
%
% omega_vertsolve_codegen(nk, ni, nj, true)
% specifies that X in omega_vertsolve.m is just a vector: X(k)
% specifies the pressure or depth of all grid points having vertical index
% k. Use this for simple Z-level models (not hybrid coordinate models).
%
% omega_vertsolve_codegen(..., OPTS)
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
% Version   : 2.0.0
%
% Modified by : --
% Date        : --
% Changes     : --

% Set defaults
VERBOSE = 1; % verbose mode
FILE_ID = 1; % output to MATLAB terminal

% Override defaults
if nargin == 5 && isstruct(OPTS)
    if isfield(OPTS, 'VERBOSE')
        VERBOSE = OPTS.VERBOSE;
    end
    if isfield(OPTS, 'FILE_ID')
        FILE_ID = OPTS.FILE_ID;
    end
end

V = filesep();
folder_start = pwd();

try
    
    % Get info about which functions are on the path
    name = 'omega_vertsolve';
    name_mex = [name '_mex'];
    file_mex = dir(which(name_mex));
    file_mat = dir(which(name));
    assert(~isempty(file_mat), ['Cannot locate ' name '.m']);
    which_eos   = which('eos');
    file_eos    = dir(which_eos);
    which_interp = which('interp_firstdim_twovars');
    file_interp = dir(which_interp);
    
    % Test values
    x = 2;
    X = [0; 1; 3; 6];
    Y = [0; 1; 2; 3];
    y = interp_firstdim_twovars(x,X,Y,Y);
    s = 34.5;
    t = 3;
    p = 1000; % or z
    m = eos(s, t, p);
    
    % Create textual identifier for this build of the MEX function.
    build_text = sprintf('%s_k%d_i%d_j%d_%dD_m=%.59e_y=%.59e', name, nk, ni, nj, (1-Xvec)*2+1, m, y);
    fileName_build_text = [file_mat.folder V name '_info.txt'];
    fileID_build_text = fopen(fileName_build_text, 'rt');
    
    
    if isempty(file_mex) ... % No mex file yet
            || file_mex.datenum < max([file_mat.datenum, file_eos.datenum, file_interp.datenum]) ... % MEX is too old
            || (fileID_build_text < 0) ... % No text file specification
            || (fileID_build_text >= 0 && ~strcmp(fgetl(fileID_build_text), build_text) && ~fclose(fileID_build_text)) % Requested parameters do not match those recorded in text file. Also close the text file.
        
        if VERBOSE
            mytic = tic;
            fprintf(FILE_ID, 'Compiling MEX for %s, with\n', name);
            fprintf(FILE_ID, ' %s in %s\n', read_function_name(which_eos), which_eos);
            fprintf(FILE_ID, ' eos(%g,%g,%g) = %e\n', s, t, p, m);
            fprintf(FILE_ID, ' %s in %s\n', read_function_name(which_interp), which_interp);
            fprintf(FILE_ID, ' interp_firstdim_twovars(%g,%s,%s) = %g\n', x, mat2str(X), mat2str(Y), y);
        end
        t_S    = coder.typeof(0, [nk, ni, nj], [false, false, false]);
        t_x    = coder.typeof(0, [ni, nj], [false, false]);
        if Xvec
            t_X = coder.typeof(0, [nk, 1], [false, false]);
        else
            t_X = t_S;
        end
        args = {t_S, t_S, t_X, t_x, t_x, t_x, t_x, 0, t_x};
        
        % Configure MEX for speed.
        mexconfig = coder.config('mex');
        mexconfig.ExtrinsicCalls = false;
        mexconfig.ResponsivenessChecks = false;
        mexconfig.IntegrityChecks = false;
        
        % Compile the MEX function
        cd(file_mat.folder);
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
    
    cd(folder_start);
    
catch err
    cd(folder_start);
    rethrow(err);
end


