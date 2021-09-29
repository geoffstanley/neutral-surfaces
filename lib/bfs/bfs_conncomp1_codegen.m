function bfs_conncomp1_codegen(nk, ni, nj, OPTS)
%BFS_CONNCOMP1_CODEGEN  Create MEX function for bfs_conncomp1
%
%
% bfs_conncomp1_codegen(nk, ni, nj)
% runs codegen on bfs_conncomp1.m, appropriate for a grid of ni by nj points
% in the horizontal and nk points in the vertical.
%
% bfs_conncomp1_codegen(..., OPTS)
% overrides default verbosity by OPTS.VERBOSE and file output by
% OPTS.FILE_ID.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


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
folder_start = pwd();

try
    
    % Get info about which functions are on the path
    name = 'bfs_conncomp1';
    name_mex = [name '_mex'];
    file_mex = dir(which(name_mex));
    file_mat = dir(which(name));
    assert(~isempty(file_mat), ['Cannot locate ' name '.m']);
    
    % Create textual identifier for this build of the MEX function.
    build_text = sprintf('%s_k%d_i%d_j%d', name_mex, nk, ni, nj);
    fileName_build_text = [file_mat.folder V name_mex '_info.txt'];
    fileID_build_text = fopen(fileName_build_text, 'rt');

    
    if isempty(file_mex) ... % No mex file yet
            || file_mex.datenum < file_mat.datenum ... % MEX is too old
            || (fileID_build_text < 0) ... % No text file specification
            || (fileID_build_text >= 0 && ~strcmp(fgetl(fileID_build_text), build_text) && ~fclose(fileID_build_text)) % Requested parameters do not match those recorded in text file. Also close the text file.
        
        if VERBOSE
            mytic = tic;
            fprintf(FILE_ID, 'Compiling MEX for %s\n', name);
        end
        
        % Note: resulting MEX is actually faster with variable size arrays
        % (using vs = true, below)
        vs = true;
        t_G    = coder.typeof(true, [ni, nj], [vs, vs]);
        nij = ni * nj;
        t_A    = coder.typeof(0, [9, nij, nj], [true, vs, vs]); % handle [? x nij] or [? x ni x nj]
        t_q    = coder.typeof(0, [nij, 1], [vs, vs]);
        
        args = {t_G, t_A, 0, t_q};

        % Configure MEX for speed.
        mexconfig = coder.config('mex');
        mexconfig.ExtrinsicCalls = false;
        mexconfig.ResponsivenessChecks = false;
        mexconfig.IntegrityChecks = false;
        
        % Compile the MEX function
        cd(folder_start);
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

