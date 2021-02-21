function varargout = analyze_surface(varargin)

%           analyze_surface.m -- The main GUI for the analyze_surface toolbox.
%
% Usage:    analyze_surface.m
%
%           A graphical user interface (GUI) used to analyze water masses
%           on density surfaces. This GUI uses m-files from the analyze_surface toolbox
%           to label hydrographies with a density variable, calculate surfaces of
%           constant density, optimize these surfaces according to 'Klocker et al.
%           A new method of forming approximately neutral surfaces.
%           submitted.' and calculate several properties on these surfaces. 
%           Detailed instructions on how to use this GUI can be found in 
%           analyze_surface.pdf (which is part of the toolbox).
%
% Input:    s             salinity
%           t or ct       potential temperature or conservative
%                         temperature
%           p             pressure
%           lats          latitude
%           longs         longitude
%
% optional 
% Input:    u             lateral velocity in x-direction
%           v             lateral velocity in y-direction
%           pre           prelabelled density
%           n2            buoyancy frequency
%           g             gravitational acceleration
%
% Output:   ocean               structure with all variables
%           ocean.settings      structure with all settings used to
%                               calculate variables/properties
% Output 
% in detail:ocean.s             salinity
%           ocean.t             potential temperature
%           ocean.ct            conservative temperature
%           ocean.lats          latitude
%           ocean.longs         longitude
%           ocean.u             lateral velocity in x-direction
%           ocean.v             lateral velocity in y-direction
%           ocean.e1t           scale factor in x-direction
%           ocean.e2t           scale factor in y-direction
%           ocean.n2            buoyancy frequency
%           ocean.g             gravitational acceleration
%           ocean.p_mid         pressure at mid-points
%           ocean.pre           pre-labelled density
%           ocean.sigma         potential density
%           ocean.pr            reference pressure
%           ocean.gamma         neutral density
%           ocean.glevels       levels of density surfaces
%           ocean.mld           mixed-layer depth
%           ocean.sns           salinity on density surface
%           ocean.ctns          conservative temperature on
%                               density surface
%           ocean.pns           pressure on density surface
%           ocean.n2_ns         buoyancy frequency on density surface
%           ocean.hel           neutral helicity on density
%                               surface
%           ocean.ss            slope error on density surface
%           ocean.sx            x-component of slope error on
%                               density surface
%           ocean.sy            y-component of slope error on
%                               density surface
%           ocean.curl_s        curl of slope error on density
%                               surface
%           ocean.ee            density gradient error on density
%                               surface
%           ocean.ex            x-component of density gradient
%                               error on density surface
%           ocean.ey            y-component of density gradient
%                               error on density surface
%           ocean.curl_e        curl of density gradient error
%                               on density surface
%           ocean.e_hel         diapycnal velocity due to
%                               neutral helicity
%           ocean.e_hel_x       x-component of diapycnal velocity
%                               due to neutral helicity
%           ocean.e_hel_y       y-component of diapycnal velocity
%                               due to neutral helicity
%           ocean.e_therm       diapycnal velocity due to
%                               thermobaricity
%           ocean.e_cab         diapycnal velocity due to
%                               cabbeling
%           ocean.*_trans       transport caused by the
%                               respective diapycnal velocity
%           ocean.*_trans_sum   integral of the transport
%                               caused by the respective
%                               diapycnal velocity
%           ocean.streamfunc    geostrophic streamfunction
%           ocean.geo_vel_x     zonal component of geostrophic
%                               velocity
%           ocean.geo_vel_y     meridional component of geostrophic
%                               velocity
%           ocean.ref_level     reference level for calculation
%                               of geostrophic velocities
%           ocean.*_i           same variables as above but for
%                               the optimized surface
%
%
% Calls:    bfrq_ct.m, e_hel.m, dia_trans.m, gamma_n_3d.m, gpoly16ct.m,
%           grav.m, helicity_surface.m, input_settings.m, mld.m, ns_3d.m,
%           optimize_surface.m, scale_fac.m, slope_error.m, e_therm_cab.m,
%           var_on_surf.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

%% initialization

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @analyze_surface_OpeningFcn, ...
    'gui_OutputFcn',  @analyze_surface_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT

%% set defaults

function analyze_surface_OpeningFcn(hObject, eventdata, handles, varargin) %#ok

handles.output = hObject;
handles.density = 'potdens';
handles.helicity = 0;
handles.slope = 0;
handles.e_hel = 0;
handles.velocity = 0;
handles.bfrq = 0;
handles.transport = 0;
handles.e_therm_cab = 0;
handles.rho = [];
handles.nit = 150;
handles.which_surface = 'initial';
handles.streamfunc = 0;
handles.geo_vel = 0;

global ocean_options %#ok

% Update handles structure

guidata(hObject, handles);

%% ouput function

function varargout = analyze_surface_OutputFcn(hObject, eventdata, handles) %#ok

varargout{1} = handles.output;

%% GUI menu

function file_Callback(hObject, eventdata, handles) %#ok

% --------------------------------------------------------------------

function settings_Callback(hObject, eventdata, handles) %#ok
% --------------------------------------------------------------------

function help_menu_Callback(hObject, eventdata, handles) %#ok

% --------------------------------------------------------------------

function input_settings_Callback(hObject, eventdata, handles) %#ok

input_settings

global ocean_options

handles.ocean.settings.grid = ocean_options.grid;
handles.ocean.settings.wrap = ocean_options.wrap;

guidata(hObject, handles);

% --------------------------------------------------------------------

function contents_Callback(hObject, eventdata, handles) %#ok

type Contents

% --------------------------------------------------------------------

function license_Callback(hObject, eventdata, handles) %#ok

type analyze_surface_license

% --------------------------------------------------------------------

function version_Callback(hObject, eventdata, handles) %#ok

helpdlg('analyze_surface by A. Klocker version 1.00','Version information')

% --------------------------------------------------------------------

function help_Callback(hObject, eventdata, handles) %#ok

open analyze_surface.pdf

%% open file

function open_Callback(hObject, eventdata, handles) %#ok

clear handles.ocean
clear ocean

[filename,pathname] = uigetfile('*.mat','Select the M-file');
eval(['load ',pathname,filename])

global ocean_options

if exist('ocean') %#ok if structure ocean exists then use it

    handles.ocean = ocean %#ok

    % open input_settings GUI if settings are not defined
    
    if ~isfield(ocean,'settings') || ~isfield(ocean,'settings') %#ok
        
        input_settings % open GUI
        
        handles.ocean.settings.grid = ocean_options.grid;
        handles.ocean.settings.wrap = ocean_options.wrap;
        
    end

    % convert potential temperature to conservative temperature if input 
    % is potential temperature 
    
    if ~isfield(ocean,'ct') && isfield(ocean,'t') %#ok
        
        tic
        disp('calculating conservative temperature');
        [handles.ocean.ct] = ct_from_t(handles.ocean.s,handles.ocean.t,handles.ocean.p); %#ok
        disp(['calculating conservative temperature took ',int2str(toc),' seconds']);
        
    end

    % calculate scale factors from longitudes/latitudes if they don't exist
    
    if ~isfield(ocean,'e1t') %#ok
        
        tic
        disp('calculating scale factors');
        [handles.ocean.e2t,handles.ocean.e1t] = scale_fac(handles.ocean.lats,handles.ocean.longs,handles.ocean.settings.wrap); %#ok
        disp(['calculating scale factors took ',int2str(toc),' seconds']);
        
    end

    % calculate gravitational acceleration if not an input
        
        tic
        display('calculating gravitational acceleration');
        handles.ocean.g = grav(handles.ocean.lats);
        display(['calculating gravitational acceleration took ',int2str(toc),' seconds']);
        

    % calculate buoyancy frequency if not an input
    
    if ~isfield(ocean,'n2') %#ok

        tic
        display('calculating buoyancy frequency');
        [handles.ocean.n2,handles.ocean.p_mid] = bfrq_ct(handles.ocean.s,handles.ocean.ct, ...
            handles.ocean.p,handles.ocean.g);
        display(['calculating buoyancy frequency took ',int2str(toc),' seconds']);
        
    end

else % if structure ocean does not exist then create it

    % choose input settings

    input_settings

    handles.ocean.settings.grid = ocean_options.grid;
    handles.ocean.settings.wrap = ocean_options.wrap;

    % read 3-dim hydrography

    handles.ocean.s = s;

    if exist('t','var')
        
        handles.ocean.t = t;
        
    end

    % make pressure matrix correct size

    if (length(size(p)) == 2)
        
        [zi,yi,xi] = size(handles.ocean.s);
        handles.ocean.p = repmat(p, [1 yi xi]);
        
    else
        
        handles.ocean.p = p;
        
    end

    % convert potential temperature into conservative temperature

    if exist('ct','var') && ~isempty(ct)
        
        handles.ocean.ct = ct;
        
    else
        
        disp('calculating conservative temperature');
        tic;
        [handles.ocean.ct] = ct_from_t(handles.ocean.s,handles.ocean.t,handles.ocean.p);
        disp(['calculating conservative temperature took ',int2str(toc),' seconds']);
        
    end

    % read lateral velocities if they exist

    if exist('u','var')
        
        handles.ocean.u = u;
        handles.ocean.v = v;
        
    end

    % read prelabelled density if it exists

    if exist('pre','var')
        
        handles.ocean.pre = pre;
        
    end

    % produce 2-dimensional longitude/latitude grid

    [yi,xi] = size(lats);
    
    if (yi == 1) || (xi == 1)
        
        [handles.ocean.longs,handles.ocean.lats] = meshgrid(longs,lats);
        
    else
        
        handles.ocean.lats = lats;
        handles.ocean.longs = longs;
        
    end

    % calculate scale factors

    tic
    disp('calculating scale factors');
    [handles.ocean.e2t,handles.ocean.e1t] = scale_fac(handles.ocean.lats,handles.ocean.longs,handles.ocean.settings.wrap); %#ok
    disp(['calculating scale factors took ',int2str(toc),' seconds']);

    % calculate gravitational acceleration

    if ~exist('g','var') %#ok
        
        tic
        display('calculating gravitational acceleration');
        handles.ocean.g = grav(handles.ocean.lats);
        display(['calculating gravitational acceleration took ',int2str(toc),' seconds']);
        
    else
        
        handles.ocean.g = g;
        
    end

    % calculate buoyancy frequency

    if ~exist('n2','var') %#ok
        
        tic
        display('calculating buoyancy frequency');
        [handles.ocean.n2,handles.ocean.p_mid] = bfrq_ct(handles.ocean.s,handles.ocean.ct, ...
            handles.ocean.p,handles.ocean.g);
        display(['calculating buoyancy frequency took ',int2str(toc),' seconds']);
        
    else
        
        handles.ocean.n2 = n2;
        
    end

end

% set more defaults

handles.ocean.ref_level = [];
handles.ocean.settings.slope = 'epsilon';

% display ocean structure

ocean = handles.ocean %#ok

% update handles

guidata(hObject, handles);

%% save file

function save_Callback(hObject, eventdata, handles) %#ok

switch handles.density

    case 'potdens' % if initial density surface is a potential density surface

        handles.ocean.sigma = handles.rho;

    case 'neutral' % if initial density surface is a neutral density surface

        handles.ocean.gamma = handles.rho;

    case 'poly' % if initial density surface is a gamma^rf density surface

        handles.ocean.poly = handles.rho;

end

ocean = handles.ocean; %#ok

% choose location and file name under which file is saved

[filename,pathname] = uiputfile('*.mat','Save to');
eval(['save ',pathname,filename,' ocean'])

% display ocean structure which is saved

display(['Saving in ',filename,':']);
disp(ocean);

%% choose density variable for initial surface

function potential_density_Callback(hObject, eventdata, handles) %#ok

handles.density = 'potdens'; % potential density

guidata(hObject, handles);

%--------------------------------------------------------------------------
function neutral_density_Callback(hObject, eventdata, handles) %#ok

handles.density = 'neutral'; % neutral density

guidata(hObject, handles);

%--------------------------------------------------------------------------
function gamma_polynomial_Callback(hObject, eventdata, handles) %#ok

handles.density = 'poly'; % gamma^rf

guidata(hObject, handles);

%--------------------------------------------------------------------------
function prelabelled_density_Callback(hObject, eventdata, handles) %#ok

handles.density = 'prelabelled'; % prelabelled density

guidata(hObject, handles);

%% choose reference pressure for potential density

function edit_pr_Callback(hObject, eventdata, handles) %#ok

prtext = get(hObject,'String');
handles.ocean.pr = eval(['[',prtext,']']);

guidata(hObject, handles);

%--------------------------------------------------------------------------
function edit_pr_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% labe data with density variable

function go_density_Callback(hObject, eventdata, handles) %#ok

switch handles.density

    case 'potdens' % potential density

        if (isempty(handles.ocean.pr))
            
            errordlg('Potential Density needs a reference pressure!','Error Message','modal');
            
        else
            
            tic
            display('calculating potential density');
            pr_rep = handles.ocean.pr * ones(size(handles.ocean.s));
            handles.rho = rho_from_ct(handles.ocean.s,handles.ocean.ct,pr_rep) - 1000;
            display(['calculating potential density took ',int2str(toc),' seconds']);
            
        end

    case 'neutral' % neutral density

        tic
        display('calculating neutral density');
        handles.rho = gamma_n_3d(handles.ocean.s,handles.ocean.ct,...
            handles.ocean.p,handles.ocean.lats,handles.ocean.longs);
        display(['calculating neutral density took ',int2str(toc),' seconds']);

    case 'poly' % gamma^rf

        tic
        display('calculating gamma polynomial');
        handles.rho = gpoly16ct(handles.ocean.s,handles.ocean.ct);
        display(['calculating gamma polynomial took ',int2str(toc),' seconds']);

    case 'prelabelled' % prelabelled density

        if (isempty(handles.ocean.pre))
            
            errordlg('Prelabelled Density does not exist!','Error Message','modal');
            
        else
            
            handles.rho = handles.ocean.pre;
            display('using prelabelled density');
            
        end

    otherwise

        errordlg('You should never get here','Error Message','modal');

end

% display density distribution

axes(handles.density_dist);
cla;
hist(handles.rho(:),100)
set(gca,'YTick',[])

guidata(hObject, handles);

%% edit levels for calculation of surfaces on which the chosen density variable is constant

function edit_levels_Callback(hObject, eventdata, handles) %#ok

gtext = get(hObject,'String');
handles.ocean.glevels = eval(['[',gtext,']']);

guidata(hObject, handles);

%--------------------------------------------------------------------------
function edit_levels_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% calculation of surfaces on which the chosen density variable is constant

function go_surfaces_Callback(hObject, eventdata, handles) %#ok

if (handles.ocean.glevels < nanmin(handles.rho(:)))
    
    errordlg('Selected density surface is too light','Error Message','modal');
    
elseif (handles.ocean.glevels > nanmax(handles.rho(:)))
    
    errordlg('Selected density surface is too dense','Error Message','modal');
    
else
    
    tic
    display('calculating s, ct and p of density surface');
    [handles.ocean.sns,handles.ocean.ctns,handles.ocean.pns,handles.ocean.dsns,handles.ocean.dctns,handles.ocean.dpns] = ...
        ns_3d(handles.ocean.s,handles.ocean.ct,handles.ocean.p,handles.rho,handles.ocean.glevels);
    display(['calculating s, ct and p on density surface took ',int2str(toc),' seconds']);

end

% calculate and display average pressure of density surface

p_ave = nanmean(handles.ocean.pns(:));
display(['the average pressure of this surface is ',int2str(p_ave),' meters']);

% exclude data above mixed-layer depth

display('excluding data above mixed-layer depth');

[zi,yi,xi] = size(handles.ocean.s);
handles.ocean.mld(1,1:yi,1:xi) = mld(handles.ocean.s,handles.ocean.ct,handles.ocean.p);
[handles.ocean.sns] = cut_off(handles.ocean.sns,handles.ocean.pns,handles.ocean.mld);
[handles.ocean.ctns] = cut_off(handles.ocean.ctns,handles.ocean.pns,handles.ocean.mld);
[handles.ocean.pns] = cut_off(handles.ocean.pns,handles.ocean.pns,handles.ocean.mld);

guidata(hObject, handles);

%% choose if using density gradient error or slope error to calculate omega-surface

function choose_epsilon_Callback(hObject, eventdata, handles) %#ok

handles.ocean.settings.slope = 'epsilon'; % density gradient error

guidata(hObject, handles);

%--------------------------------------------------------------------------
function choose_s_Callback(hObject, eventdata, handles) %#ok

handles.ocean.settings.slope = 's'; % slope error

guidata(hObject, handles);

%% edit number of iterations

function edit_niter_Callback(hObject, eventdata, handles) %#ok

nittext = get(hObject,'String');
handles.nit = eval(['[',nittext,']']);

guidata(hObject, handles);

function edit_niter_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% optimize surfaces (construct omega-surfaces)

function go_improve_Callback(hObject, eventdata, handles) %#ok

display('optimizing density surface');
tic
[handles.ocean.sns_i,handles.ocean.ctns_i,handles.ocean.pns_i] = ...
    optimize_surface(handles.ocean.s,handles.ocean.ct,handles.ocean.p,...
    handles.ocean.g,handles.ocean.n2,handles.ocean.sns,handles.ocean.ctns,...
    handles.ocean.pns,handles.ocean.e1t,handles.ocean.e2t,handles.nit,...
    handles.ocean.settings.slope,handles.ocean.settings.wrap);
display(['optimizing density surface took ',int2str(toc),' seconds']);

guidata(hObject, handles);

%% choose which properties on density surfaces to calculate

function helicity_switch_Callback(hObject, eventdata, handles) %#ok

handles.helicity = get(hObject,'Value'); % neutral helicity

guidata(hObject, handles);

%--------------------------------------------------------------------------
function slope_switch_Callback(hObject, eventdata, handles) %#ok

handles.slope = get(hObject,'Value'); % slope and density gradient errors

guidata(hObject, handles);

%--------------------------------------------------------------------------
function pdvv_switch_Callback(hObject, eventdata, handles) %#ok

handles.e_hel = get(hObject,'Value'); % diapycnal velocity due to neutral helicity

guidata(hObject, handles);

%--------------------------------------------------------------------------
function velocity_switch_Callback(hObject, eventdata, handles) %#ok

handles.velocity = get(hObject,'Value'); % lateral velocity on density surface

guidata(hObject, handles);

%--------------------------------------------------------------------------
function buoyancy_switch_Callback(hObject, eventdata, handles) %#ok

handles.bfrq = get(hObject,'Value'); % buoyancy frequency

guidata(hObject, handles);

%--------------------------------------------------------------------------
function transport_switch_Callback(hObject, eventdata, handles) %#ok

handles.transport = get(hObject,'Value'); % diapycnal transports

guidata(hObject, handles);

%--------------------------------------------------------------------------
function e_therm_cab_switch_Callback(hObject, eventdata, handles) %#ok

handles.e_therm_cab = get(hObject,'Value'); % cabbeling/thermobaricity 

guidata(hObject, handles);

%--------------------------------------------------------------------------
function streamfunction_switch_Callback(hObject, eventdata, handles) %#ok

handles.streamfunc = get(hObject,'Value'); % geostrophic streamfunction

guidata(hObject, handles);

%--------------------------------------------------------------------------
function geo_vel_switch_Callback(hObject, eventdata, handles) %#ok

handles.geo_vel = get(hObject,'Value'); % geostrophic velocities

guidata(hObject, handles);

%% choose reference level for geostrophic velocity calculation

function edit_ref_level_Callback(hObject, eventdata, handles) %#ok

ref_level_text = get(hObject,'String');
handles.ocean.ref_level = eval(['[',ref_level_text,']']);

guidata(hObject, handles);

function edit_ref_level_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% choose if calculating properties on initial, optimized or both density surfaces

function which_surface_popup_Callback(hObject, eventdata, handles) %#ok

str = get(hObject,'String');
val = get(hObject,'Value');

handles.which_surface = str{val};

guidata(hObject, handles);

function which_surface_popup_CreateFcn(hObject, eventdata, handles) %#ok

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% calculate properties on density surfaces

function go_properties_Callback(hObject, eventdata, handles) %#ok

switch handles.which_surface

    case 'initial' % on initial density surface

        tic
        display('calculating buoyancy frequency on initial density surface');
        [handles.ocean.n2_ns] = var_on_surf(handles.ocean.pns,handles.ocean.p_mid,handles.ocean.n2);
        display(['calculating buoyancy frequency on initial density surface took ',int2str(toc),' seconds']);

        if (handles.helicity == 1)

            tic
            display('calculating neutral helicity on initial density surface');
            [handles.ocean.hel] = ...
                helicity_surface(handles.ocean.sns,handles.ocean.ctns, ...
                handles.ocean.pns,handles.ocean.n2_ns,handles.ocean.g, ...
                handles.ocean.e1t,handles.ocean.e2t,handles.ocean.settings.wrap);
            display(['calculating neutral helicity on initial density surface took ',int2str(toc),' seconds']);

        end

        if (handles.slope == 1)

            tic
            display('calculating slope errors on initial density surface');
            [handles.ocean.ss,handles.ocean.sx,handles.ocean.sy,handles.ocean.curl_s, ...
                handles.ocean.ee,handles.ocean.ex,handles.ocean.ey,handles.ocean.curl_e,handles.ocean.fdd] = ...
                slope_error(handles.ocean.p,handles.ocean.g,handles.ocean.n2,handles.ocean.sns, ...
                handles.ocean.ctns,handles.ocean.pns,handles.ocean.e1t,handles.ocean.e2t, ...
                'op',handles.ocean.settings.wrap);
            display(['calculating slope errors on initial density surface took ',int2str(toc),' seconds']);

        end

        if (handles.velocity == 1)

            if isfield(handles.ocean,'u')

                tic
                display('calculating velocity on initial density surface');
                [handles.ocean.uns,handles.ocean.vns] = var_on_surf(handles.ocean.pns,handles.ocean.p,handles.ocean.u,handles.ocean.v);
                display(['calculating velocity on initial density surface took ',int2str(toc),' seconds']);

            else

                errordlg('You can not calculate velocities on a density surface without velocity data!','Error Message','modal');

            end

        end

        if (handles.e_hel == 1)

            if isfield(handles.ocean,'u')

                if (handles.velocity ~= 1)

                    tic
                    display('calculating velocity on initial density surface');
                    [handles.ocean.uns,handles.ocean.vns] = var_on_surf(handles.ocean.pns,handles.ocean.p,handles.ocean.u,handles.ocean.v);
                    display(['calculating velocity on initial density surface took ',int2str(toc),' seconds']);

                end

                if (handles.slope ~=1)

                    tic
                    display('calculating slope errors on initial density surface');
                    [handles.ocean.ss,handles.ocean.sx,handles.ocean.sy,handles.ocean.curl_s, ...
                        handles.ocean.ee,handles.ocean.ex,handles.ocean.ey,handles.ocean.curl_e,handles.ocean.fdd] = ...
                        slope_error(handles.ocean.p,handles.ocean.g,handles.ocean.n2,handles.ocean.sns, ...
                        handles.ocean.ctns,handles.ocean.pns,handles.ocean.e1t,handles.ocean.e2t, ...
                        'op',handles.ocean.settings.wrap);
                    display(['calculating slope errors on initial density surface took ',int2str(toc),' seconds']);

                end

                tic
                display('calculating v*s on initial density surface');
                [handles.ocean.e_hel,handles.ocean.e_hel_x,handles.ocean.e_hel_y] = ...
                    e_hel(handles.ocean.sx,handles.ocean.sy,handles.ocean.uns, ...
                    handles.ocean.vns,handles.ocean.settings.grid,handles.ocean.settings.wrap);

                display(['calculating v*s on initial density surface took ',int2str(toc),' seconds']);

            else

                errordlg('You can not calculate v*s without velocity data!','Error Message','modal');

            end

        end

        if (handles.e_therm_cab == 1)

            tic
            display('calculating diapycnal velocity due to thermobaricity/cabbeling through initial density surface');
       
            [handles.ocean.e_cab,handles.ocean.e_therm] = e_therm_cab(handles.ocean.sns,handles.ocean.ctns,...
                handles.ocean.pns,handles.ocean.n2_ns,handles.ocean.g,handles.ocean.e1t,handles.ocean.e2t,'op','none');
            display(['calculating diapycnal velocity due to thermobaricity/cabbeling through initial density surface took ',int2str(toc),' seconds']);

        end

        % calculate diapycnal transports for diapycnal velocities which
        % have already been calculated

        if (handles.transport == 1)

            tic
            display('calculating diapycnal transports through initial density surface');

            if ~isempty(handles.ocean.e_cab)

                [handles.ocean.e_cab_trans,handles.ocean.e_cab_trans_sum] = ...
                    dia_trans(handles.ocean.e_cab,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('transport caused cabbeling is not calculated on initial density surface');

            end

            if ~isempty(handles.ocean.e_therm)

                [handles.ocean.e_therm_trans,handles.ocean.e_therm_trans_sum] = ...
                    dia_trans(handles.ocean.e_therm,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('transport caused by thermobaricity is not calculated on initial density surface');

            end

            if isfield(handles.ocean,'e_hel')

                [handles.ocean.e_hel_trans,handles.ocean.e_hel_trans_sum] = ...
                    dia_trans(handles.ocean.e_hel,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('transport caused by e^hel is not calculated on initial density surface');

            end

            display(['calculating diapycnal transports through initial density surface took ',int2str(toc),' seconds']);

        end

        if (handles.streamfunc == 1) && (handles.geo_vel == 0)

            tic
            display('calculating geostrophic streamfunction on initial density surface');
            [handles.ocean.streamfunc] = optimize_streamfunc(handles.ocean.s,handles.ocean.ct,...
                handles.ocean.p,handles.ocean.sns,handles.ocean.ctns,handles.ocean.pns,...
                handles.ocean.settings.wrap);
            display(['calculating geostrophic streamfunction on initial density surface took ',int2str(toc),' seconds']);

        elseif ((handles.streamfunc == 1) && (handles.geo_vel == 1)) || ((handles.streamfunc == 0) && (handles.geo_vel == 1))
            
            tic
            display('calculating geostrophic streamfunction on initial density surface');
            [handles.ocean.streamfunc] = optimize_streamfunc(handles.ocean.s,handles.ocean.ct,...
                handles.ocean.p,handles.ocean.sns,handles.ocean.ctns,handles.ocean.pns,...
                handles.ocean.settings.wrap);
            display(['calculating geostrophic streamfunction on initial density surface took ',int2str(toc),' seconds']);

            if (isempty(handles.ocean.ref_level))

                errordlg('Geostrophic velocity calculation needs a reference level!','Error Message','modal');

            else

                tic
                display('calculating geostrophic velocities on initial density surface');
                [handles.ocean.geo_vel_x,handles.ocean.geo_vel_y] = geo_vel(handles.ocean.streamfunc,...
                    handles.ocean.s,handles.ocean.ct,handles.ocean.p,nanmean(handles.ocean.sns(:)),...
                    nanmean(handles.ocean.ctns(:)),handles.ocean.lats,handles.ocean.e1t,handles.ocean.e2t,...
                    handles.ocean.ref_level,handles.ocean.settings.wrap);
                display(['calculating geostrophic velocities on initial density surface took ',int2str(toc),' seconds']);

            end

        end
        
        handles.ocean %DELETE

    case 'optimized' % on omega surface


        tic
        display('calculating buoyancy frequency on omega surface');
        [handles.ocean.n2_ns_i] = var_on_surf(handles.ocean.pns_i,handles.ocean.p_mid,handles.ocean.n2);
        display(['calculating buoyancy frequency on omega surface took ',int2str(toc),' seconds']);

        if (handles.helicity == 1)

            tic
            display('calculating neutral helicity on omega surface');
            [handles.ocean.hel_i] = ...
                helicity_surface(handles.ocean.sns_i,handles.ocean.ctns_i, ...
                handles.ocean.pns_i,handles.ocean.n2_ns_i,handles.ocean.g, ...
                handles.ocean.e1t,handles.ocean.e2t,handles.ocean.settings.wrap);
            display(['calculating neutral helicity on omega surface took ',int2str(toc),' seconds']);

        end

        if (handles.slope == 1)

            tic
            display('calculating slope errors on omega surface');
            [handles.ocean.ss_i,handles.ocean.sx_i,handles.ocean.sy_i,handles.ocean.curl_s_i, ...
                handles.ocean.ee_i,handles.ocean.ex_i,handles.ocean.ey_i,handles.ocean.curl_e_i,handles.ocean.fdd_i] = ...
                slope_error(handles.ocean.p,handles.ocean.g,handles.ocean.n2,handles.ocean.sns_i, ...
                handles.ocean.ctns_i,handles.ocean.pns_i,handles.ocean.e1t,handles.ocean.e2t, ...
                'op',handles.ocean.settings.wrap);
            display(['calculating slope errors on omega surface took ',int2str(toc),' seconds']);

        end

        if (handles.velocity == 1)

            if isfield(handles.ocean,'u')

                tic
                display('calculating velocity on omega surface');
                [handles.ocean.uns_i,handles.ocean.vns_i] = var_on_surf(handles.ocean.pns_i,handles.ocean.p,handles.ocean.u,handles.ocean.v);
                display(['calculating velocity on omega surface took ',int2str(toc),' seconds']);

            else

                errordlg('You can not calculate velocities on a density surface without velocity data!','Error Message','modal');

            end

        end

        if (handles.e_hel == 1)

            if isfield(handles.ocean,'u')

                if (handles.velocity ~= 1)
                    tic
                    display('calculating velocity on omega surface');
                    [handles.ocean.uns_i,handles.ocean.vns_i] = var_on_surf(handles.ocean.pns_i,handles.ocean.p,handles.ocean.u,handles.ocean.v);
                    display(['calculating velocity on omega surface took ',int2str(toc),' seconds']);
                end

                if (handles.slope ~=1)
                    tic
                    display('calculating slope errors on omega surface');
                    [handles.ocean.ss_i,handles.ocean.sx_i,handles.ocean.sy_i,handles.ocean.curl_s_i, ...
                        handles.ocean.ee_i,handles.ocean.ex_i,handles.ocean.ey_i,handles.ocean.curl_e_i,handles.ocean.fdd_i] = ...
                        slope_error(handles.ocean.p,handles.ocean.g,handles.ocean.n2,handles.ocean.sns_i, ...
                        handles.ocean.ctns_i,handles.ocean.pns_i,handles.ocean.e1t,handles.ocean.e2t, ...
                        'op',handles.ocean.settings.wrap);
                    display(['calculating slope errors on omega surface took ',int2str(toc),' seconds']);
                end

                tic
                display('calculating v*s on omega surface');
                [handles.ocean.e_hel_i,handles.ocean.e_hel_x_i,handles.ocean.e_hel_y_i] = ...
                    e_hel(handles.ocean.sx_i,handles.ocean.sy_i,handles.ocean.uns_i, ...
                    handles.ocean.vns_i,handles.ocean.settings.grid,handles.ocean.settings.wrap);

                display(['calculating v*s on omega surface took ',int2str(toc),' seconds']);

            else

                errordlg('You can not calculate v*s without velocity data!','Error Message','modal');

            end

        end

        if (handles.e_therm_cab == 1)

            tic
            display('calculating diapycnal velocity due to thermobaricity/cabbeling through initial density surface');
            [handles.ocean.e_cab_i,handles.ocean.e_therm_i] = e_therm_cab(handles.ocean.sns_i,handles.ocean.ctns_i,...
                handles.ocean.pns_i,handles.ocean.n2_ns_i,handles.ocean.g,handles.ocean.e1t,handles.ocean.e2t,'op',handles.ocean.settings.wrap);

            display(['calculating diapycnal velocity due to thermobaricity/cabbeling through initial density surface took ',int2str(toc),' seconds']);

        end

        % calculate diapycnal transports for diapycnal velocities which
        % have already been calculated

        if (handles.transport == 1)

            tic
            display('calculating diapycnal transports through omega surface');

            if ~isempty(handles.ocean.e_cab_i)

                [handles.ocean.e_cab_trans_i,handles.ocean.e_cab_trans_sum_i] = ...
                    dia_trans(handles.ocean.e_cab_i,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('Transport caused cabbeling is not calculated on omega surface');

            end

            if ~isempty(handles.ocean.e_therm_i)

                [handles.ocean.e_therm_trans_i,handles.ocean.e_therm_trans_sum_i] = ...
                    dia_trans(handles.ocean.e_therm_i,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('Transport caused by thermobaricity is not calculated on omega surface');

            end
            if isfield(handles.ocean,'e_hel_i')

                [handles.ocean.e_hel_trans_i,handles.ocean.e_hel_trans_sum_i] = ...
                    dia_trans(handles.ocean.e_hel_i,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('Transport caused by e^hel is not calculated on omega surface');

            end

            display(['calculating diapycnal transports through omega surface took ',int2str(toc),' seconds']);

        end

        if (handles.streamfunc == 1) && (handles.geo_vel == 0)

            tic
            display('calculating geostrophic streamfunction on omega surface');
            [handles.ocean.streamfunc_i] = optimize_streamfunc(handles.ocean.s,handles.ocean.ct,...
                handles.ocean.p,handles.ocean.sns_i,handles.ocean.ctns_i,handles.ocean.pns_i,...
                handles.ocean.settings.wrap);
            display(['calculating geostrophic streamfunction on omega surface took ',int2str(toc),' seconds']);

        elseif ((handles.streamfunc == 1) && (handles.geo_vel == 1)) || ((handles.streamfunc == 0) && (handles.geo_vel == 1))

            tic
            display('calculating geostrophic streamfunction on omega surface');
            [handles.ocean.streamfunc_i] = optimize_streamfunc(handles.ocean.s,handles.ocean.ct,...
                handles.ocean.p,handles.ocean.sns_i,handles.ocean.ctns_i,handles.ocean.pns_i,...
                handles.ocean.settings.wrap);
            display(['calculating geostrophic streamfunction on omega surface took ',int2str(toc),' seconds']);
            
            if (isempty(handles.ocean.ref_level))

                errordlg('Geostrophic velocity calculation needs a reference level!','Error Message','modal');

            else

                tic
                display('calculating geostrophic velocities on omega surface');
                [handles.ocean.geo_vel_x_i,handles.ocean.geo_vel_y_i] = geo_vel(handles.ocean.streamfunc_i,...
                    handles.ocean.s,handles.ocean.ct,handles.ocean.p,nanmean(handles.ocean.sns_i(:)),...
                    nanmean(handles.ocean.ctns_i(:)),handles.ocean.lats,handles.ocean.e1t,handles.ocean.e2t,...
                    handles.ocean.ref_level,handles.ocean.settings.wrap);
                display(['calculating geostrophic velocities on omega surface took ',int2str(toc),' seconds']);

            end

        end

    case 'both' % on initial and omega surface

        tic
        display('calculating buoyancy frequency on initial and omega surface');
        [handles.ocean.n2_ns] = var_on_surf(handles.ocean.pns,handles.ocean.p_mid,handles.ocean.n2);
        [handles.ocean.n2_ns_i] = var_on_surf(handles.ocean.pns_i,handles.ocean.p_mid,handles.ocean.n2);
        display(['calculating buoyancy frequency on initial and omega surface took ',int2str(toc),' seconds']);

        if (handles.helicity == 1)

            tic
            display('calculating neutral helicity on initial density surface');
            [handles.ocean.hel] = ...
                helicity_surface(handles.ocean.sns,handles.ocean.ctns, ...
                handles.ocean.pns,handles.ocean.n2_ns,handles.ocean.g, ...
                handles.ocean.e1t,handles.ocean.e2t,handles.ocean.settings.wrap);
            display(['calculating neutral helicity on initial density surface took ',int2str(toc),' seconds']);

        end

        if (handles.slope == 1)

            tic
            display('calculating slope errors on initial density surface');
            [handles.ocean.ss,handles.ocean.sx,handles.ocean.sy,handles.ocean.curl_s, ...
                handles.ocean.ee,handles.ocean.ex,handles.ocean.ey,handles.ocean.curl_e,handles.ocean.fdd] = ...
                slope_error(handles.ocean.p,handles.ocean.g,handles.ocean.n2,handles.ocean.sns, ...
                handles.ocean.ctns,handles.ocean.pns,handles.ocean.e1t,handles.ocean.e2t, ...
                'op',handles.ocean.settings.wrap);
            display(['calculating slope errors on initial density surface took ',int2str(toc),' seconds']);

        end

        if (handles.velocity == 1)

            if isfield(handles.ocean,'u')

                tic
                display('calculating velocity on initial density surface');
                [handles.ocean.uns,handles.ocean.vns] = var_on_surf(handles.ocean.pns,handles.ocean.p,handles.ocean.u,handles.ocean.v);
                display(['calculating velocity on initial density surface took ',int2str(toc),' seconds']);

            else

                errordlg('You can not calculate velocities on a density surface without velocity data!','Error Message','modal');

            end

        end

        if (handles.e_hel == 1)

            if isfield(handles.ocean,'u')

                if (handles.velocity ~= 1)

                    tic
                    display('calculating velocity on initial density surface');
                    [handles.ocean.uns,handles.ocean.vns] = var_on_surf(handles.ocean.pns,handles.ocean.p,handles.ocean.u,handles.ocean.v);
                    display(['calculating velocity on initial density surface took ',int2str(toc),' seconds']);

                end

                if (handles.slope ~=1)

                    tic
                    display('calculating slope errors on initial density surface');
                    [handles.ocean.ss,handles.ocean.sx,handles.ocean.sy,handles.ocean.curl_s, ...
                        handles.ocean.ee,handles.ocean.ex,handles.ocean.ey,handles.ocean.curl_e,handles.ocean.fdd] = ...
                        slope_error(handles.ocean.p,handles.ocean.g,handles.ocean.n2,handles.ocean.sns, ...
                        handles.ocean.ctns,handles.ocean.pns,handles.ocean.e1t,handles.ocean.e2t, ...
                        'op',handles.ocean.settings.wrap);
                    display(['calculating slope errors on initial density surface took ',int2str(toc),' seconds']);

                end

                tic
                display('calculating v*s on initial density surface');
                [handles.ocean.e_hel,handles.ocean.e_hel_x,handles.ocean.e_hel_y] = ...
                    e_hel(handles.ocean.sx,handles.ocean.sy,handles.ocean.uns, ...
                    handles.ocean.vns,handles.ocean.settings.grid,handles.ocean.settings.wrap);

                display(['calculating v*s on initial density surface took ',int2str(toc),' seconds']);

            else

                errordlg('You can not calculate v*s without velocity data!','Error Message','modal');

            end

        end

        if (handles.e_therm_cab == 1)

            tic
            display('calculating diapycnal velocity due to thermobaricity/cabbeling through initial density surface');
            [handles.ocean.e_cab,handles.ocean.e_therm] = e_therm_cab(handles.ocean.sns,handles.ocean.ctns,...
                handles.ocean.pns,handles.ocean.n2_ns,handles.ocean.g,handles.ocean.e1t,handles.ocean.e2t,'op',handles.ocean.settings.wrap);
            display(['calculating diapycnal velocity due to thermobaricity/cabbeling through initial density surface took ',int2str(toc),' seconds']);

        end

        if (handles.transport == 1)

            tic
            display('calculating diapycnal transports through initial density surface');

            if ~isempty(handles.ocean.e_cab)

                [handles.ocean.e_cab_trans,handles.ocean.e_cab_trans_sum] = ...
                    dia_trans(handles.ocean.e_cab,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('transport caused cabbeling is not calculated on initial density surface');

            end

            if ~isempty(handles.ocean.e_therm)

                [handles.ocean.e_therm_trans,handles.ocean.e_therm_trans_sum] = ...
                    dia_trans(handles.ocean.e_therm,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('transport caused by thermobaricity is not calculated on initial density surface');

            end

            if isfield(handles.ocean,'e_hel')

                [handles.ocean.e_hel_trans,handles.ocean.e_hel_trans_sum] = ...
                    dia_trans(handles.ocean.e_hel,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('transport caused by e^hel is not calculated on initial density surface');

            end

            display(['calculating diapycnal transports through initial density surface took ',int2str(toc),' seconds']);

        end

        if (handles.streamfunc == 1) && (handles.geo_vel == 0)

            tic
            display('calculating geostrophic streamfunction on initial density surface');
            [handles.ocean.streamfunc] = optimize_streamfunc(handles.ocean.s,handles.ocean.ct,...
                handles.ocean.p,handles.ocean.sns,handles.ocean.ctns,handles.ocean.pns,...
                handles.ocean.settings.wrap);
            display(['calculating geostrophic streamfunction on initial density surface took ',int2str(toc),' seconds']);


        elseif ((handles.streamfunc == 1) && (handles.geo_vel == 1)) || ((handles.streamfunc == 0) && (handles.geo_vel == 1)) 

             tic
            display('calculating geostrophic streamfunction on initial density surface');
            [handles.ocean.streamfunc] = optimize_streamfunc(handles.ocean.s,handles.ocean.ct,...
                handles.ocean.p,handles.ocean.sns,handles.ocean.ctns,handles.ocean.pns,...
                handles.ocean.settings.wrap);
            display(['calculating geostrophic streamfunction on initial density surface took ',int2str(toc),' seconds']);
            
            if (isempty(handles.ocean.ref_level))

                errordlg('Geostrophic velocity calculation needs a reference level!','Error Message','modal');

            else

                tic
                display('calculating geostrophic velocities on initial density surface');
                [handles.ocean.geo_vel_x,handles.ocean.geo_vel_y] = geo_vel(handles.ocean.streamfunc,...
                    handles.ocean.s,handles.ocean.ct,handles.ocean.p,nanmean(handles.ocean.sns(:)),...
                    nanmean(handles.ocean.ctns(:)),handles.ocean.lats,handles.ocean.e1t,handles.ocean.e2t,...
                    handles.ocean.ref_level,handles.ocean.settings.wrap);
                display(['calculating geostrophic velocities on initial density surface took ',int2str(toc),' seconds']);

            end

        end
        
        if (handles.helicity == 1)

            tic
            display('calculating neutral helicity on omega surface');
            [handles.ocean.hel_i] = ...
                helicity_surface(handles.ocean.sns_i,handles.ocean.ctns_i, ...
                handles.ocean.pns_i,handles.ocean.n2_ns,handles.ocean.g, ...
                handles.ocean.e1t,handles.ocean.e2t,handles.ocean.settings.wrap);
            display(['calculating neutral helicity on omega surface took ',int2str(toc),' seconds']);

        end

        if (handles.slope == 1)

            tic
            display('calculating slope errors on omega surface');
            [handles.ocean.ss_i,handles.ocean.sx_i,handles.ocean.sy_i,handles.ocean.curl_s_i, ...
                handles.ocean.ee_i,handles.ocean.ex_i,handles.ocean.ey_i,handles.ocean.curl_e_i,handles.ocean.fdd_i] = ...
                slope_error(handles.ocean.p,handles.ocean.g,handles.ocean.n2,handles.ocean.sns_i, ...
                handles.ocean.ctns_i,handles.ocean.pns_i,handles.ocean.e1t,handles.ocean.e2t, ...
                'op',handles.ocean.settings.wrap);
            display(['calculating slope errors on omega surface took ',int2str(toc),' seconds']);

        end

        if (handles.velocity == 1)

            if isfield(handles.ocean,'u')

                tic
                display('calculating velocity on omega surface');
                [handles.ocean.uns_i,handles.ocean.vns_i] = var_on_surf(handles.ocean.pns_i,handles.ocean.p,handles.ocean.u,handles.ocean.v);
                display(['calculating velocity on omega surface took ',int2str(toc),' seconds']);

            else

                errordlg('You can not calculate velocities on a density surface without velocity data!','Error Message','modal');

            end

        end

        if (handles.e_hel == 1)

            if isfield(handles.ocean,'u')

                if (handles.velocity ~= 1)

                    tic
                    display('calculating velocity on omega surface');
                    [handles.ocean.uns_i,handles.ocean.vns_i] = var_on_surf(handles.ocean.pns_i,handles.ocean.p,handles.ocean.u,handles.ocean.v);
                    display(['calculating velocity on omega surface took ',int2str(toc),' seconds']);

                end

                if (handles.slope ~=1)

                    tic
                    display('calculating slope errors on omega surface');
                    [handles.ocean.ss_i,handles.ocean.sx_i,handles.ocean.sy_i,handles.ocean.curl_s_i, ...
                        handles.ocean.ee_i,handles.ocean.ex_i,handles.ocean.ey_i,handles.ocean.curl_e_i,handles.ocean.fdd_i] = ...
                        slope_error(handles.ocean.s,handles.ocean.ct,handles.ocean.p,handles.ocean.g,handles.ocean.n2,handles.ocean.sns_i, ...
                        handles.ocean.ctns_i,handles.ocean.pns_i,handles.ocean.e1t,handles.ocean.e2t, ...
                        'op',handles.ocean.settings.wrap);
                    display(['calculating slope errors on omega surface took ',int2str(toc),' seconds']);

                end

                tic
                display('calculating v*s on omega surface');
                [handles.ocean.e_hel_i,handles.ocean.e_hel_x_i,handles.ocean.e_hel_y_i] = ...
                    e_hel(handles.ocean.sx_i,handles.ocean.sy_i,handles.ocean.uns_i, ...
                    handles.ocean.vns_i,handles.ocean.settings.grid,handles.ocean.settings.wrap);
                handles.vel_done = 1;
                display(['calculating v*s on omega surface took ',int2str(toc),' seconds']);

            else

                errordlg('You can not calculate v*s without velocity data!','Error Message','modal');

            end

        end

        if (handles.e_therm_cab == 1)

            tic
            display('calculating diapycnal velocity due to thermobaricity/cabbeling through initial density surface');
            [handles.ocean.e_cab_i,handles.ocean.e_therm_i] = e_therm_cab(handles.ocean.sns_i,handles.ocean.ctns_i,...
                handles.ocean.pns_i,handles.ocean.n2_ns_i,handles.ocean.g,handles.ocean.e1t,handles.ocean.e2t,'op',handles.ocean.settings.wrap);
            display(['calculating diapycnal velocity due to thermobaricity/cabbeling through initial density surface took ',int2str(toc),' seconds']);

        end

        % calculate diapycnal transports for diapycnal velocities which
        % have already been calculated

        if (handles.transport == 1)

            tic
            display('calculating diapycnal transports through omega surface');

            if ~isempty(handles.ocean.e_cab_i)

                [handles.ocean.e_cab_trans_i,handles.ocean.e_cab_trans_sum_i] = ...
                    dia_trans(handles.ocean.e_cab_i,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('Transport caused cabbeling is not calculated on omega surface');

            end

            if ~isempty(handles.ocean.e_therm_i)

                [handles.ocean.e_therm_trans_i,handles.ocean.e_therm_trans_sum_i] = ...
                    dia_trans(handles.ocean.e_therm_i,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('Transport caused by thermobaricity is not calculated on omega surface');

            end

            if isfield(handles.ocean,'e_hel_i')

                [handles.ocean.e_hel_trans_i,handles.ocean.e_hel_trans_sum_i] = ...
                    dia_trans(handles.ocean.e_hel_i,handles.ocean.e1t,handles.ocean.e2t);

            else

                display('Transport caused by e^hel is not calculated on omega surface');

            end

            display(['calculating diapycnal transports through omega surface took ',int2str(toc),' seconds']);

        end

        if (handles.streamfunc == 1) && (handles.geo_vel == 0)

            tic
            display('calculating geostrophic streamfunction on omega surface');
            [handles.ocean.streamfunc_i] = optimize_streamfunc(handles.ocean.s,handles.ocean.ct,...
                handles.ocean.p,handles.ocean.sns_i,handles.ocean.ctns_i,handles.ocean.pns_i,...
                handles.ocean.settings.wrap);
            display(['calculating geostrophic streamfunction on omega surface took ',int2str(toc),' seconds']);

        elseif (handles.streamfunc == 1) && (handles.geo_vel == 1) || (handles.streamfunc == 0) && (handles.geo_vel == 1)

            tic
            display('calculating geostrophic streamfunction on omega surface');
            [handles.ocean.streamfunc_i] = optimize_streamfunc(handles.ocean.s,handles.ocean.ct,...
                handles.ocean.p,handles.ocean.sns_i,handles.ocean.ctns_i,handles.ocean.pns_i,...
                handles.ocean.settings.wrap);
            display(['calculating geostrophic streamfunction on omega surface took ',int2str(toc),' seconds']);
            
            if (isempty(handles.ocean.ref_level))

                errordlg('Geostrophic velocity calculation needs a reference level!','Error Message','modal');

            else

                tic
                display('calculating geostrophic velocities on omega surface');
                [handles.ocean.geo_vel_x_i,handles.ocean.geo_vel_y_i] = geo_vel(handles.ocean.streamfunc_i,...
                    handles.ocean.s,handles.ocean.ct,handles.ocean.p,nanmean(handles.ocean.sns_i(:)),...
                    nanmean(handles.ocean.ctns_i(:)),handles.ocean.lats,handles.ocean.e1t,handles.ocean.e2t,...
                    handles.ocean.ref_level,handles.ocean.settings.wrap);
                display(['calculating geostrophic velocities on omega surface took ',int2str(toc),' seconds']);

            end

        end

end

guidata(hObject, handles);

%% others

function density_panel_CreateFcn(hObject, eventdata, handles) %#ok

function density_dist_CreateFcn(hObject, eventdata, handles) %#ok
