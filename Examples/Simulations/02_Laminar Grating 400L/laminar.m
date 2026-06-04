clear; warning('off', 'all');

base    = fileparts(mfilename('fullpath'));
oc_path = fullfile(base, '..', '..', 'Optical_Constants');
addpath(fullfile(base, '..', '..', 'helpers'));


% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_res_nm = 0.5;
z_res_nm = 0.5;

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
stack = add_layer(stack, fullfile(oc_path, 'n_Pt_cxro.txt'), 29);   % 29 nm Pt coating
% stack = add_layer(stack, fullfile(oc_path, 'n_C_cxro.txt'), 1);   % 1 nm C contamination


% Sweep parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sweep.type      = 'energy';
sweep.values    = 50:5:1000;
% sweep.alpha_deg = 4;        % fixed grazing incidence angle in degrees
sweep.Cff     = 2.25;     % uncomment to use Cff-based angle instead

% Solver options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.FourierOrders = 11;
options.pol           = -1;          % -1 = TM,  +1 = TE
options.GR_Order      = -1;
options.z_res_nm      = z_res_nm;
options.reticolo_path = fullfile(base, '..', '..', '..', 'V9', 'reticolo_allege_v9');
options.output_dir    = fullfile(base, 'Results');
options.oc_path       = oc_path;
options.verbose       = true;


% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = run_rcwa(stack, substrate_file, sweep, options);


% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

photon_eV = sweep.values(round(end/2));   % mid-range energy for cross-section preview
save_png  = fullfile(options.output_dir, 'meshgrid_laminar.png');
plot_meshgrid(stack, substrate_file, photon_eV, z_res_nm, save_png);

if isempty(results.efficiency)
    error('No valid efficiency data — check angle range and optical constant files.');
end

figure(1); clf;
plot(results.sweep_values, results.efficiency * 100, ...
     'b-o', 'LineWidth', 1.2, 'MarkerSize', 2);

xlabel('Photon Energy (eV)', 'FontSize', 12);
ylabel('Diffraction Efficiency (%)', 'FontSize', 12);

pol_str = 'TM'; if options.pol == 1; pol_str = 'TE'; end
if isfield(sweep, 'Cff')
    inc_tag = sprintf('Cff=%.2f', sweep.Cff);
else
    inc_tag = sprintf('α=%.2f°', sweep.alpha_deg);
end
title(sprintf('Energy Sweep | %d l/mm | %s | %s | Order %+d', ...
    grPeriod_lpermm, inc_tag, pol_str, options.GR_Order), 'FontSize', 11);

grid on;
set(gca, 'FontSize', 11);
xlim([min(results.sweep_values), max(results.sweep_values)]);

saveas(gcf, fullfile(options.output_dir, 'efficiency_fixed_angle_energy_sweep.png'));
fprintf('Plot saved.\n');
