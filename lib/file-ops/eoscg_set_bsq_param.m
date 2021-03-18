function eoscg_set_bsq_param(eos_in, eos_out, grav, rhob)
%EOSCG_SET_BSQ_PARAM  Copy a .m file but modify 'grav =' and 'rhob =' lines.
%
%
% eoscg_set_bsq_param(eos_in, eos_out, grav, rhob)
% creates a new file whose full path is given by eos_out that is identical
% to the file whose full path is eos_in, except that lines that start with
% 'grav =' or 'rhob =' are modified to assign the input values grav and
% rhob.  These are the gravitational acceleration and Boussinesq reference
% density, both of which determine the common modification of an equation
% of state for seawater into a Boussinesq equation of state for seawater,
% wherein the pressure [dbar] is taken to be 1e-4 * grav * rhob * z, where
% z [m, positive] is the depth.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% Ensure the input file exists
assert(exist(eos_in, 'file') > 0, ['Cannot locate ' eos_in]);

[~, fcn_name] = fileparts(eos_in); % Get the name of the file without path or extension
folder_out = fileparts(eos_out);   % Get the path to the output file 

if ~exist(folder_out, 'dir')
    mkdir(folder_out);     % Make directory for the output file, if necessary
end

fi = fopen(eos_in, 'rt');  % Open input file for textual reading
fo = fopen(eos_out, 'wt'); % Open output file for textual writing
textline = fgetl(fi); % Read first line
while ischar(textline)
    % If the line specifes grav or rhob, use the new values in the output
    % file.  Otherwise, copy the line.
    if strncmp(textline, 'grav =', 6)
        fprintf(fo, 'grav = %.15f; %% gravitational acceleration [m / s^2]\n', grav);
    elseif strncmp(textline, 'rhob =', 6)
        fprintf(fo, 'rhob = %.15f; %% Boussinesq reference density [kg / m^3]\n', rhob);
    else
        fprintf(fo, '%s\n', textline);
    end
    textline = fgetl(fi); % Continue reading the next line
end
fclose(fi); % Close the files
fclose(fo);
clear(fcn_name); % Make sure the new file gets used