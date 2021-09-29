function omega_vertsolve_codegen(nk, ni, nj, OPTS)
%OMEGA_VERTSOLVE_CODEGEN  Create MEX function for omega_vertsolve
%
%
% omega_vertsolve_codegen(nk, ni, nj)
% runs codegen on omega_vertsolve.m, appropriate for a grid
% of ni by nj points in the horizontal and nk points in the vertical.
%
% omega_vertsolve_codegen(..., OPTS)
% overrides default verbosity by OPTS.VERBOSE and file output by
% OPTS.FILE_ID.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


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
    
    if strcmpi(file_eos.folder, pwd())
        % If eos() is found in the current working directory, make sure
        % that it is also on the top of MATLAB's path. This is to ensure
        % that eos() used in the MEX function compiled with MATLAB's code
        % generation (codegen) is the intended one, in the current working
        % directory. Interpreted MATLAB uses the current working directory
        % at the top of its search path, but MEX functions that are
        % compiled in a separate directory do not know about this current
        % working directory, so they could inadvertantly call the wrong
        % functions.
        addpath(pwd)
    end
    
    % Test values
    s = 34.5;
    t = 3;
    p = 1000;
    m = eos(s, t, p);
    
    % Create textual identifier for this build of the MEX function.
    build_text = sprintf('%s_k%d_i%d_j%d_m=%.59e', name, nk, ni, nj, m);
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
        
        % Note: for omega surfaces, resulting MEX is the same speed whether
        % arrays are variable size or fixed size (vs = true, or vs = false)
        vs = true;
        t_Sppc   = coder.typeof(0, [8, nk-1, ni, nj], [true, vs, vs, vs]);
        t_P  = coder.typeof(0, [nk, ni, nj], [vs, vs, vs]);
        t_BotK   = coder.typeof(0, [ni, nj], [vs, vs]);
        t_p      = coder.typeof(0, [ni, nj], [vs, vs]);
        
        args = {t_Sppc, t_Sppc, t_P, t_BotK, t_p, t_p, t_p, 0, t_p};
        
        % Configure MEX for speed.
        mexconfig = coder.config('mex');
        mexconfig.ExtrinsicCalls = false;
        mexconfig.ResponsivenessChecks = false;
        mexconfig.IntegrityChecks = false;
        
        % Compile the MEX function
        cd(file_mat.folder);
        codegen(name, '-report', '-args', args, '-config', mexconfig, '-o', [file_mat.folder V name_mex], '-d', [file_mat.folder V 'codegen' V 'mex' V name]);
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


