function ns_install()
%NS_INSTALL  Handle Java, MEX, path, and data installation steps for Neutral Surfaces.
%
% This file should exist in the root directory of Neutral Surfaces.

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
% but WITHOUT ANY WARRANTY; withouasdft even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% Step 0: establish path things
V = filesep(); % /  or  \  depending on OS

% Get path to this file, in ./
PATH_PROJECT = [fileparts(mfilename('fullpath')) V];

% Get path to ./lib/dat/
PATH_DAT = [PATH_PROJECT 'lib' V 'dat' V];

% Get path to ./topobaric-surface
PATH_TOPOBARIC_SURFACE = [PATH_PROJECT 'topobaric-surface' V];

% Step 1: Check MATLAB version
VER_MATLAB_STR = version();
j = find(VER_MATLAB_STR == '.', 2, 'first');
VER_MATLAB = str2double(VER_MATLAB_STR(1 : j(2)-1));
assert(VER_MATLAB >= 9.1, 'MATLAB 9.1 (2016b) or later is required (for automatic expansion).');


% Step 2: Install recon.jar
VER_JAVA_STR = version('-java');
i = find(VER_JAVA_STR == ' ', 1, 'first');
j = find(VER_JAVA_STR == '.', 2, 'first');
VER_JAVA = VER_JAVA_STR(i+1 : j(2)-1);
fprintf('* Found MATLAB uses version %s of the Java Runtime Environment (JRE)\n', VER_JAVA);
if ~strcmp(VER_JAVA, '1.7') && ~strcmp(VER_JAVA, '1.8')
    % Must be that VER_JAVA > 1.8, since we already checked MATLAB >= 2016b
    % which had Java 1.7.
    disp('ReCon has only been pre-compiled (making a .jar file) for JRE versions 1.7 and 1.8.')
    disp('Either compile ReCon from the included source code (see README.md for instructions),');
    disp('or contact the author (Geoff Stanley, g.stanley@unsw.edu.au or geoffstanley@gmail.com)');
    disp('kindly requesting an update.');
    
    WARN_JAVA = true;
else
    
    PATH_RECON_VER = [PATH_TOPOBARIC_SURFACE 'src' V 'recon' V 'build' V 'recon_' VER_JAVA '.jar'];
    PATH_RECON = [PATH_TOPOBARIC_SURFACE 'src' V 'recon' V 'build' V 'recon.jar'];
    
    fprintf('Attempting to link or copy \n %s \nto \n %s\n', PATH_RECON_VER, PATH_RECON);
    if ispc
        status = system(['mklink ' PATH_RECON_VER ' ' PATH_RECON]);
    else
        status = system(['ln -sf ' PATH_RECON_VER ' ' PATH_RECON]);
    end
    if status ~= 0
        if exist(PATH_RECON, 'file')
            delete(PATH_RECON)
        end
        copyfile(PATH_RECON_VER, PATH_RECON, 'f');
    end
    
    WARN_JAVA = false;
end


% Step 3: Add ReCon to MATLAB's javaclasspath
PATH_PREF = prefdir(); % Get path to preferences directory where javaclasspath.txt must be
fprintf('* Shall I add \n %s\nto MATLAB''s javaclasspath in\n %s%sjavaclasspath.txt\n', PATH_RECON, PATH_PREF, V);
reply = input('? Answer Y unless you have done this already. [Y/n]: ', 's');
if isempty(reply) || lower(reply(1)) == 'y'
    fileID = fopen([PATH_PREF V 'javaclasspath.txt'], 'at');
    fprintf(fileID, '\n%s\n', [PATH_TOPOBARIC_SURFACE 'src' V 'recon' V 'build' V 'recon.jar']);
    fclose(fileID);
end


% Step 4: Increase memory that MATLAB allocates to Java
disp('* Increase the memory that MATLAB allocates to Java:');
disp('Go to MATLAB Preferences -> General -> Java HEAP Memory');
disp('Move the slider far to the right.');
disp('e.g. the default value is 512 MB, and you increase it to 2048 MB.');
input('Press any key to continue.')


% Step 5: Setup MEX with C
disp('* About to setup MEX for use with C. You should have a working C compiler');
disp('installed on this system, such as gcc. On MacOS, you can install Xcode to get clang.');
disp('Note: Topobaric Surface can run without MEX compiled code, but it will be far slower.');
input('Press any key to continue.', 's');
mex('-setup', 'C');


% Step 6: Add Neutral Surfaces to MATLAB's search path
fprintf('* Adding\n %s\nand appropriate sub-folders to MATLAB''s search path.\n', PATH_PROJECT)
run([PATH_PROJECT 'ns_add_to_path.m']);

reply = input('* Shall I edit your startup.m file to add these folders to MATLAB''s search path every time MATLAB starts? [Y/n]: ', 's');
if isempty(reply) || lower(reply(1)) == 'y'
    PATH_USER = userpath();
    file_startup = fopen([PATH_USER V 'startup.m'], 'at');
    fprintf(file_startup, '\n%% Add Topobaric Surface\n');
    fprintf(file_startup, 'run(''%s'');\n', [PATH_PROJECT 'ns_add_to_path.m']);
    fclose(file_startup);
end


% Step 7: Download data
reply = input('* Shall I download data used by the scripts in this folder? [Y/n]: ', 's');
if isempty(reply) || lower(reply(1)) == 'y'
    
    % Determine if wget or curl is available
    while true
        disp('Checking for system command: wget');
        [status, ~] = system('wget --version');
        if ~status
            download = @(url, target) system(sprintf('wget -P %s %s', target, url));
            break
        else
            disp('Checking for system command: curl');
            [status, ~] = system('curl --version');
            if ~status
                download = @(url, target) system(sprintf('curl -o %s %s', target, url));
                break
            else
                disp('Cannot find system commands "curl" or "wget".')
                input('Please install curl or install wget then return to MATLAB and press any key.')
                reply = input('* Hello again. Is curl or wget installed now? [Y/n]: ', 's');
                if ~(isempty(reply) || lower(reply(1)) == 'y') % no
                    disp('Okay. Data can be manually downloaded at a later time if you choose.');
                    disp('See README.md for more information.');
                    download = false;
                    break
                end
            end
        end
    end
    
    if isa(download, 'function_handle')
        % Download OCCA data
        URL = 'ftp://mit.ecco-group.org/ecco_for_las/OCCA_1x1_v2/2004-6/annual/';
        reply = input('* Do you want to download OCCA data (used by illustrative_surface.m)? [Y/n]: ', 's');
        if isempty(reply) || lower(reply(1)) == 'y'
            disp(['OCCA data will be downloaded from ' URL]);
            while true
                PATH = input('* Enter path to store OCCA data: ', 's');
                if PATH(end) ~= V % Add trailing / or \ if not already present
                    PATH(end+1) = V; %#ok<AGROW>
                end
                if exist(PATH, 'dir') || mkdir(PATH)
                    vars = {'salt', 'theta', 'phihyd', 'etan'};
                    for var = vars
                        download(sprintf('%s/DD%s.0406annclim.nc', URL, var{1}), PATH)
                    end
                    
                    if ~exist([PATH 'omega_v1.1gjs' V], 'dir')
                        mkdir([PATH 'omega_v1.1gjs' V]);
                    end
                    download('https://ndownloader.figshare.com/files/14536133', [PATH 'omega_v1.1gjs' V]); % Omega 1000
                    movefile([PATH 'omega_v1.1gjs' V '14536133'], [PATH 'omega_v1.1gjs' V 'omega.0406annclim.from_SIGMA1000_through_(180,0,1000).nonBoussinesq.mat']);
                    
                    % Record path to OCCA data in a text file
                    fprintf('Recording\n %s\nas the path to OCCA data in\n %s\n', PATH, PATH_PATHS);
                    fprintf('You can modify this later if you move the data.\n')
                    file_id = fopen([PATH_DAT V 'PATH_OCCA.txt'], 'wt'); % overwrite contents
                    fprintf(file_id, '%s', PATH);
                    fclose(file_id);
                    
                    break
                end
                reply = input('* Invalid directory. Try again? [Y/n]: ', 's');
                if ~(isempty(reply) || lower(reply(1)) == 'y') % no
                    disp('Skipping OCCA download.');
                    break
                end
            end
        end
        
        % Download ECCO data
        URL = 'ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/';
        reply = input('Do you want to download ECCO2 data (used by run.m, pitch.m, example.m)? [Y/n]: ', 's');
        if isempty(reply) || lower(reply(1)) == 'y'
            disp(['ECCO2 data will be downloaded from ' URL]);
            while true
                PATH = input('* Enter path to store ECCO2 data: ', 's');
                if PATH(end) ~= V % Add trailing / or \ if not already present
                    PATH(end+1) = V; %#ok<AGROW>
                end
                if exist(PATH, 'dir') || mkdir(PATH)
                    TIMESTEP = '20021223';
                    vars = {'SALT', 'THETA', 'UVEL', 'VVEL'};
                    for var = vars
                        download(sprintf('%s.nc/%s.1440x720x50.%s.nc', var{1}, var{1}, TIMESTEP), PATH);
                    end
                    var = {'SSH'};
                    n = datenum(TIMESTEP, 'yyyymmdd');
                    download(sprintf('%s.nc/%s.1440x720x50.%s.nc', var{1}, var{1}, datestr(n-1, 'yyyymmdd')), PATH);
                    download(sprintf('%s.nc/%s.1440x720x50.%s.nc', var{1}, var{1}, TIMESTEP                ), PATH);
                    download(sprintf('%s.nc/%s.1440x720x50.%s.nc', var{1}, var{1}, datestr(n+1, 'yyyymmdd')), PATH);
                    
                    if ~exist([PATH 'GAMMA' V], 'dir')
                        mkdir([PATH 'GAMMA' V]);
                    end
                    if ~exist([PATH 'omega_v1.1gjs' V], 'dir')
                        mkdir([PATH 'omega_v1.1gjs' V]);
                    end
                    download('https://ndownloader.figshare.com/files/14536058', [PATH 'GAMMA' V]); % gamma^n
                    movefile([PATH 'GAMMA' V '14536058'], [PATH 'GAMMA' V 'GAMMA.1440x720x50.20021223.mat']);
                    download('https://ndownloader.figshare.com/files/14536061', [PATH 'omega_v1.1gjs' V]); % Omega 1000
                    movefile([PATH 'omega_v1.1gjs' V '14536061'], [PATH 'omega_v1.1gjs' V 'omega.1440x720.20021223.from_SIGMA1000_through_(180,0,1000).Boussinesq.mat']);
                    download('https://ndownloader.figshare.com/files/14536064', [PATH 'omega_v1.1gjs' V]); % Omega 2000
                    movefile([PATH 'omega_v1.1gjs' V '14536064'], [PATH 'omega_v1.1gjs' V 'omega.1440x720.20021223.from_SIGMA2000_through_(180,0,2000).Boussinesq.mat']);
                    
                    % Record path to ECCO2 data in a text file
                    fprintf('Recording\n %s\nas the path to ECCO2 data in\n %s\n', PATH, PATH_PATHS);
                    fprintf('You can modify this later if you move the data.\n')
                    file_id = fopen([PATH_DAT V 'PATH_ECCO2.txt'], 'wt'); % overwrite contents
                    fprintf(file_id, '%s', PATH);
                    fclose(file_id);
                    
                    break
                end
                reply = input('* Invalid directory. Try again? [Y/n]: ', 's');
                if ~(isempty(reply) || lower(reply(1)) == 'y') % no
                    disp('Skipping ECCO2 download.');
                    break
                end
            end
        end
    end
end

% Step 8: Restart MATLAB
disp('----------------------------------------------------------------')
disp('To complete installation, MATLAB must be restarted.')
disp('Upon restarting, execute')
disp('>> javaclasspath')
disp('in MATLAB and confirm the path to recon.jar appears in the final line of the output.')
if WARN_JAVA
    disp('Reminder: topobaric_surface cannot run without recon.jar!')
end
reply = input('* Exit MATLAB now? [Y/n]: ', 's');
if isempty(reply) || lower(reply(1)) == 'y'
    exit;
end

