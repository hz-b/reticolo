% Single Layer Laminar (Trapezoidal) Grating — Geometry and Meshgrid Visualisation
%
% Builds a laminar (trapezoidal-profile) grating with 15-degree sidewalls and
% a 1 nm carbon contamination layer, then renders the refractive-index
% cross-section at a chosen photon energy.
% No RCWA solve is performed — this script is purely for geometry inspection.

clear; warning('off', 'all');

base    = fileparts(mfilename('fullpath'));
oc_path = fullfile(base, '..', '..', 'Optical_Constants');
addpath(fullfile(base, '..', '..', '..', 'helpers'));


% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_res_nm = 0.5;

grPeriod_lpermm   = 400;
depth_nm          = 14.9;
sidewall_angle_deg = 15;           % from horizontal
groove_width_ratio = 0.67;         % groove floor width / period

period_nm = 1e6 / grPeriod_lpermm;

fprintf('Grating period : %.2f nm  (%.0f l/mm)\n', period_nm, grPeriod_lpermm);
fprintf('Groove depth   : %.2f nm  |  Sidewall angle: %.1f deg\n', depth_nm, sidewall_angle_deg);

grating = build_grating('trapezoidal', period_nm, depth_nm, ...
                         groove_width_ratio, sidewall_angle_deg, x_res_nm);


% Layer stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

substrate_file = fullfile(oc_path, 'n_Si_cxro.txt');

stack = build_stack(grating);
stack = add_layer(stack, fullfile(oc_path, 'n_C_cxro.txt'), 1);   % 1 nm C contamination


% Meshgrid visualisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

photon_eV   = 500;
z_res_nm    = 0.5;
results_dir = fullfile(base, 'Results');
save_png    = fullfile(results_dir, sprintf('meshgrid_%.0feV.png', photon_eV));

plot_meshgrid(stack, substrate_file, photon_eV, z_res_nm, save_png);
