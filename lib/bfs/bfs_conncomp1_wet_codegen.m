function bfs_conncomp1_wet_codegen(nk, ni, nj, OPTS)
%BFS_CONNCOMP1_WET_CODEGEN  Create MEX function for bfs_conncomp1_wet
%
%
% bfs_conncomp1_wet_codegen(nk, ni, nj)
% runs codegen on bfs_conncomp1_wet.m, appropriate for a grid of ni by nj
% points in the horizontal and nk points in the vertical.
%
% bfs_wet_codegen(..., OPTS)
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

% Get info about which functions are on the path
name = 'bfs_conncomp1_wet';
name_mex = [name '_mex'];
file_mex = dir(which(name_mex));
file_mat = dir(which(name));
assert(~isempty(file_mat), ['Cannot locate ' name '.m']);
which_eos   = which('eos');
file_eos    = dir(which_eos);

% Test values
s = 34.5;
t = 3;
p = 1000;
m = eos(s, t, p);

% Create textual identifier for this build of the MEX function.
build_text = sprintf('%s_k%d_i%d_j%d_m=%.59e', name_mex, nk, ni, nj, m);
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
        fprintf(FILE_ID, ' eos(%g,%g,%g) = %e\n', s, t, p, m);
    end
    
    % Note: resulting MEX is actually faster with variable size arrays
    % (using vs = true, below)
    vs = true;
    t_SppX   = coder.typeof(0, [8, nk-1, ni, nj], [true, vs, vs, vs]);
    t_X = coder.typeof(0, [nk, ni, nj], [vs, vs, vs]);
    nij = ni * nj;
    t_x    = coder.typeof(0, [ni, nj], [vs, vs]);
    t_A    = coder.typeof(0, [9, nij, nj], [true, vs, vs]); % [? p nij] or [? p ni p nj]
    t_BotK = coder.typeof(0, [ni, nj], [vs, vs]);
    t_q    = coder.typeof(0, [nij, 1], [vs, vs]);
    
    args = {t_SppX, t_SppX, t_X, t_x, t_x, t_x, t_x, 0, t_A, t_BotK, 0, t_q};
    
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


