function fcn_name = read_function_name(file_name)
% Read the name of the function in the function declaration line of the
% specified file, or throw an error.
%
% fcn_name = read_function_name(fullpath)
% returns the name of the function as written into the file specified by
% file_name.

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