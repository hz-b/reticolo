% Multilayer Blazed Grating — Geometry and Meshgrid Visualisation

clear; warning('off', 'all');

base    = fileparts(mfilename('fullpath'));
oc_path = fullfile(base, '..', '..', 'Optical_Constants');
addpath(fullfile(base, '..', '..', '..', 'helpers'));


% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_res_nm = 0.1;

grPeriod_lpermm      = 2400;
grBlazeAngle_deg     = 1.37;
grAntiBlazeAngle_deg = 3.25;

period_nm  = 1e6 / grPeriod_lpermm;
tan_blaze  = tand(grBlazeAngle_deg);
tan_anti   = tand(grAntiBlazeAngle_deg);
w_blaze_nm = period_nm / (1 + tan_blaze / tan_anti);
depth_nm   = w_blaze_nm * tan_blaze;

fprintf('Grating period : %.2f nm  (%.0f l/mm)\n', period_nm, grPeriod_lpermm);
fprintf('Groove depth   : %.4f nm\n', depth_nm);

grating = build_grating('blazed', period_nm, depth_nm, ...
                         grBlazeAngle_deg, grAntiBlazeAngle_deg, x_res_nm);


% Layer stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

substrate_file = fullfile(oc_path, 'n_Si_cxro.txt');

layers_per_period  = 2;
num_multilayers    = 60;
layer_a_file       = fullfile(oc_path, 'n_Cr_cxro.txt');
layer_a_nm         = 1.9;
layer_b_file       = fullfile(oc_path, 'n_C_cxro.txt');
layer_b_nm         = 2.9;

stack = build_stack(grating);
for i = 1:(layers_per_period * num_multilayers)
    stack = add_layer(stack, layer_a_file, layer_a_nm);
    stack = add_layer(stack, layer_b_file, layer_b_nm);
end

fprintf('Multilayer stack: %d periods of Cr(%.1fnm) / C(%.1fnm)\n', ...
        num_multilayers, layer_a_nm, layer_b_nm);


% Meshgrid visualisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

photon_eV   = 500;
z_res_nm    = 0.1;
results_dir = fullfile(base, 'Results');
save_png    = fullfile(results_dir, sprintf('meshgrid_%.0feV.png', photon_eV));

plot_meshgrid(stack, substrate_file, photon_eV, z_res_nm, save_png);
