% Single Layer Blazed Grating — Geometry and Meshgrid Visualisation


clear; warning('off', 'all');

base    = fileparts(mfilename('fullpath'));
oc_path = fullfile(base, '..', '..', 'Optical_Constants');
addpath(fullfile(base, '..', '..', '..', 'helpers'));


% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_res_nm = 0.5;

grPeriod_lpermm      = 600;
grBlazeAngle_deg     = 0.73;
grAntiBlazeAngle_deg = 5.60;

period_nm = 1e6 / grPeriod_lpermm;
tan_blaze = tand(grBlazeAngle_deg);
tan_anti  = tand(grAntiBlazeAngle_deg);
w_blaze_nm = period_nm / (1 + tan_blaze / tan_anti);
depth_nm   = w_blaze_nm * tan_blaze;

fprintf('Grating period : %.2f nm  (%.0f l/mm)\n', period_nm, grPeriod_lpermm);
fprintf('Groove depth   : %.4f nm\n', depth_nm);

grating = build_grating('blazed', period_nm, depth_nm, ...
                         grBlazeAngle_deg, grAntiBlazeAngle_deg, x_res_nm);


% Layer stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

substrate_file = fullfile(oc_path, 'n_Si_cxro.txt');

stack = build_stack(grating);
stack = add_layer(stack, fullfile(oc_path, 'n_Au_cxro.txt'), 31);   % 31 nm Au coating


% Meshgrid visualisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

photon_eV   = 500;
z_res_nm    = 0.5;
results_dir = fullfile(base, 'Results');
save_png    = fullfile(results_dir, sprintf('meshgrid_%.0feV.png', photon_eV));

plot_meshgrid(stack, substrate_file, photon_eV, z_res_nm, save_png);
