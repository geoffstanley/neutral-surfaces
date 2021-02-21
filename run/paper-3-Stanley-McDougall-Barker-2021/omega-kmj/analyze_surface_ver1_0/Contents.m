% analyze_surface toolbox
%
% Contents:
%
%   ab_from_ct              - Calculate alpha, beta, alpha/beta and T_b
%   analyze_surface         - analyze_surface.m -- The main GUI for the 
%                             analyze_surface toolbox
%   analyze_surface_license - License statement and permissions for 
%                             analyze_surface package.
%   bernoulli_streamfunc    - Calculate the Bernoulli streamfunction
%   bfrq_ct                 - Calculate the buoyancy frequency, N^2 
%   change                  - Change values in matrix
%   coriolis                - Calculate the coriolis parameter
%   cut_off                 - Cut off data on a surface above a certain depth. 
%   delta_streamfunc        - Calculate local difference in streamfunctions
%   dia_trans               - Calculate dia-surface transport and integrate 
%                             over the domain
%   dyn_height              - Calculate dynamic height
%   e_hel                   - Calculate the diapycnal velocity due to 
%                             neutral helicity
%   e_therm_cab             - Calculate cabbeling and thermobaricity
%   gamma_ew                - Calculate Eden/Willebrand neutral density
%   gamma_n_3d              - Calculate neutral density (gamma^n)
%   geo_streamfunc          - Calculate the geostrophic streamfunction on
%                             approximately neutral surfaces
%   geo_vel                 - Calculate geostrophic velocities
%   gpoly16ct               - Calculate neutral density (gamma^rf) through 
%                             a rational function
%   grad_surf               - Find the surface gradient of a variable in x 
%                             and y-direction.
%   grav                    - Calculate gravitational acceleration
%   helicity_pressure       - Calculate neutral helicity on pressure levels
%   helicity_surface        - Calculate neutral helicity on a density surface
%   input_settings          - input_settings.m -- GUI for setting input settings.
%   mld                     - Calculate mixed-layer depth
%   montgomery_streamfunc   - Calculate the Montgomery/Zhang & Hogg streamfunction
%   ns_3d                   - Find s, ct and p on a density surface
%   optimize_streamfunc     - Optimize streamfunction on density surfaces 
%   optimize_surface        - Optimize density surfaces to minimize 
%                             fictitious diapycnal diffusivity
%   scale_fac               - Find distances between gridpoints of given 
%                             latitude/longitude in a 2-dim data set 
%   slope_error             - Calculate slope errors
%   var_on_surf             - Vertically interpolate a variable onto a 
%                             density surface
%   z_from_p                - Calculate depth from pressure
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%