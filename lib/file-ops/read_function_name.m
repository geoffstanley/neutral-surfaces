function fcn_name = read_function_name(file_name)
% Read the name of the function in the function declaration line of the
% specified file, or throw an error.
%
% fcn_name = read_function_name(fullpath)
% returns the name of the function as written into the file specified by
% file_name.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


fcn_name = '';
fid = fopen(file_name);
tline = fgetl(fid);
while ischar(tline)
    token = regexp(tline,'.*function.*=\s*(\w+).*', 'tokens');
    if ~isempty(token)
        fcn_name = token{1}{1};
        break
    end
    tline = fgetl(fid);
end
fclose(fid);
if isempty(fcn_name)
    error('Failed to read function declaration line in %s', file_name);
end