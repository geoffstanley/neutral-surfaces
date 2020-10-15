function ntp_bottle_to_cast_codegen(nk, OPTS)
%NTP_BOTTLE_TO_CAST_CODEGEN  Create MEX function for ntp_bottle_to_cast
%
%
% ntp_bottle_to_cast_codegen(nk,false)
% runs codegen on ntp_bottle_to_cast.m, appropriate for a grid
% of nk points in the vertical.
%
% ntp_bottle_to_cast_codegen(nk,true)
% specifies that X in ntp_bottle_to_cast.m is just a vector: X(k) specifies
% the pressure or depth of all grid points having vertical index k. Use
% this for simple Z-level models (not hybrid coordinate models).
%
% ntp_bottle_to_cast_codegen(..., OPTS)
% overrides default verbosity by OPTS.VERBOSE and file output by
% OPTS.FILE_ID.

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2020  Geoff Stanley
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


% Set defaults
VERBOSE = 1; % verbose mode
FILE_ID = 1; % output to MATLAB terminal

% Override defaults
if nargin == 2 && isstruct(OPTS)
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
    name = 'ntp_bottle_to_cast';
    name_mex = [name '_mex'];
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
    build_text = sprintf('%s_k%d_m=%.59e', name, nk, m);
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

        t_SppX = coder.typeof(0, [8, nk-1], [true, true]);
        t_X    = coder.typeof(0, [nk  , 1], [true, false]);
        
        args = {t_SppX, t_SppX, t_X, 0, 0, 0, 0, 0};
        
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




