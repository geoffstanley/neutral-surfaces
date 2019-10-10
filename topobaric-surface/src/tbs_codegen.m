function run_codegen(nz, nx, ny, Pvec, fcn_list)
%RUN_CODEGEN  Create MEX functions for Topobaric Surface.
%
%
% run_codegen(nz,nx,ny)
% runs codegen on all the functions topobaric_update, orthobaric_update,
% interp_casts, pchi1d1, isopycnal, and deltasurf, appropriate for a grid
% of nz points in the vertical and nx by ny points in the horizontal.
%
% run_codegen(nz,nx,ny,true)
% specifies that the pressure data is actually just a vector specifying
% data at each vertical grid point, and taking the same value at all
% horizontal locations. Use this for simple Z-level models (not hybrid
% coordinate models).
%
% run_codegen(...,fcn_list)
% only runs codegen on the functions specified in the cell array of
% strings given by fcn_list.

% --- Copyright:
% Copyright 2019 Geoff Stanley
%
% This file is part of Topobaric Surface.
%
% Topobaric Surface is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% Topobaric Surface is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with Topobaric Surface.  If not, see
% <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com
% Version   : 1.0
%
% Modified by : --
% Date        : --
% Changes     : --

start = pwd;

if nargin < 4
    Pvec = false;
end

if nargin < 5
    fcn_list = {...
        'topobaric_update', ...
        'orthobaric_update', ...
        'interp1qn2', ...
        'pchi1d1', ...
        'isopycnal', ...
        'deltasurf' };
elseif ischar(fcn_list)
    fcn_list = {fcn_list};
end

try
    
    % file_eos = dir(which('eos.m'));
    
    mexconfig = coder.config('mex');
    mexconfig.ExtrinsicCalls = false;
    mexconfig.ResponsivenessChecks = false;
    mexconfig.IntegrityChecks = false;
    
    
    for i = 1 : length(fcn_list)
        
        fcn = fcn_list{i};
        
        switch lower(fcn)
            
            case 'topobaric_update'
                file_mat = dir(which('topobaric_update'));
                assert(~isempty(file_mat), 'Cannot locate topobaric_update');
                % file_mex = dir(which('topobaric_update_mex'));
                % if isempty(file_mex) || file_mex.datenum < file_mat.datenum || file_mex.datenum < file_eos.datenum
                cd(file_mat.folder)
                if Pvec
                    type_P  = coder.typeof(0, [nz, 1], [false, false]);
                else
                    type_P  = coder.typeof(0, [nz, nx, ny], [false, true, true]);
                end
                type_S      = coder.typeof(0, [nz, nx, ny], [false, false, false]);
                type_K      = coder.typeof(0, [nx, ny], [false, false]);
                type_p      = coder.typeof(0, [nx, ny], [false, false]);
                type_branch = coder.typeof(0, [nx, ny], [false, false]);
                type_vafnp  = coder.typeof(0, [5, nx*ny], [false true]);
                codegen('topobaric_update', '-args', {type_S, type_S, type_P, type_K, type_p, type_branch, type_vafnp, 0, 0, 0}, '-config', mexconfig, '-o', 'topobaric_update_mex');
                clear type_P type_S type_K type_p type_branch type_vafnp
                % end
                
            case 'orthobaric_update'
                file_mat = dir(which('orthobaric_update'));
                assert(~isempty(file_mat), 'Cannot locate orthobaric_update');
                % file_mex = dir(which('orthobaric_update_mex'));
                % if isempty(file_mex) || file_mex.datenum < file_mat.datenum || file_mex.datenum < file_eos.datenum
                cd(file_mat.folder)
                if Pvec
                    type_P  = coder.typeof(0, [nz, 1], [false, false]);
                else
                    type_P  = coder.typeof(0, [nz, nx, ny], [false, true, true]);
                end
                type_S      = coder.typeof(0, [nz, nx, ny], [false, false, false]);
                type_K      = coder.typeof(0, [nx, ny], [false, false]);
                type_p      = coder.typeof(0, [nx, ny], [false, false]);
                type_breaks = coder.typeof(0, [1, 64], [false, true]);
                type_coefs  = coder.typeof(0, [64, 8], [true, true]);
                codegen('orthobaric_update', '-args', {type_S, type_S, type_P, type_K, type_p, type_breaks, type_coefs, 0, 0, 0}, '-config', mexconfig, '-o', 'orthobaric_update_mex');
                clear type_P type_S type_K type_p type_breaks type_coefs
                % end
                
            case 'interp1qn2'
                file_mat = dir(which('interp1qn2'));
                assert(~isempty(file_mat), 'Cannot locate interp1qn2');
                % file_mex = dir(which('interp1qn2_mex'));
                % if isempty(file_mex) || file_mex.datenum < file_mat.datenum || file_mex.datenum < file_eos.datenum
                cd(file_mat.folder)
                if Pvec
                    type_P  = coder.typeof(0, [nz, 1], [false, false]);
                else
                    type_P  = coder.typeof(0, [nz, nx, ny], [false, true, true]);
                end
                type_S  = coder.typeof(0, [nz, nx, ny], [false, false, false]);
                type_p = coder.typeof(0, [1, nx, ny], [false, false, false]);
                codegen('interp1qn2', '-args', {type_p, type_P, type_S, type_S}, '-config', mexconfig, '-o', 'interp1qn2_mex');
                clear type_P type_S type_p
                % end
                
            case 'pchi1d1'
                file_mat = dir(which('pchi1d1'));
                assert(~isempty(file_mat), 'Cannot locate pchi1d1');
                % file_mex = dir(which('pchi1d1_mex'));
                % if isempty(file_mex) || file_mex.datenum < file_mat.datenum
                cd(file_mat.folder)
                type_U = coder.typeof(0, [256, nx, ny], [true  true true]);
                type_X = coder.typeof(0, [nz,  nx, ny], [false true true]);
                codegen('pchi1d1', '-args', {type_U, type_X, type_X}, '-o', 'pchi1d1_mex');
                clear type_U type_X
                % end
                
            case 'isopycnal'
                file_mat = dir(which('isopycnal'));
                assert(~isempty(file_mat), 'Cannot locate isopycnal');
                % file_mex = dir(which('isopycnal_ijp_mex'));
                % if isempty(file_mex) || file_mex.datenum < file_mat.datenum || file_mex.datenum < file_eos.datenum
                cd(file_mat.folder)
                type_S = coder.typeof(0, [nz, nx, ny], [false, false, false]);
                if Pvec
                    type_P  = coder.typeof(0, [nz, 1], [false, false]);
                else
                    type_P  = coder.typeof(0, [nz, nx, ny], [false, true, true]);
                end
                codegen('isopycnal', '-args', {type_S, type_S, type_P, 0, [0, 0, 0], 0, 0}, '-config', mexconfig, '-o', 'isopycnal_ijp_mex');
                clear type_S type_P
                % end
                
            case 'deltasurf'
                file_mat = dir(which('deltasurf'));
                assert(~isempty(file_mat), 'Cannot locate deltasurf');
                % file_mex = dir(which('deltasurf_mex'));
                % if isempty(file_mex) || file_mex.datenum < file_mat.datenum || file_mex.datenum < file_eos.datenum
                cd(file_mat.folder)
                type_S = coder.typeof(0, [nz, nx, ny], [false, false, false]);
                if Pvec
                    type_P  = coder.typeof(0, [nz, 1], [false, false]);
                else
                    type_P  = coder.typeof(0, [nz, nx, ny], [false, true, true]);
                end
                codegen('deltasurf', '-args', {type_S, type_S, type_P, 0, 0, [0, 0, 0], 0, 0}, '-config', mexconfig, '-o', 'deltasurf_mex');
                clear type_S type_P
                % end
                
            otherwise
                error('run_codegen:unknownFunction', 'Unknown function name %s', lower(fcn));
        end
    end
    
    cd(start);
catch err
    cd(start);
    rethrow(err);
end