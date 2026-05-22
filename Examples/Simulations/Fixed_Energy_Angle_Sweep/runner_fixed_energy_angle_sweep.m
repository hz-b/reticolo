% Fixed Energy Angle Sweep — Single Layer Blazed Grating Simulation
%
% Sweeps grazing incidence angle at a fixed photon energy.
% Use this to find the optimal working angle or to reproduce a rocking-curve
% measurement at a given photon energy.
%
% Grating : 600 l/mm single-layer Au on Si
% Sweep   : grazing angle 0.5–8° at fixed photon energy
% Output  : efficiency CSV + plot saved to Results/

clear; warning('off', 'all');

base    = fileparts(mfilename('fullpath'));
oc_path = fullfile(base, '..', '..', 'Optical_Constants');
addpath(fullfile(base, '..', '..', '..', 'helpers'));


% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_res_nm = 0.5;
z_res_nm = 0.5;

grPeriod_lpermm      = 600;
grBlazeAngle_deg     = 0.73;
grAntiBlazeAngle_deg = 5.60;

period_nm  = 1e6 / grPeriod_lpermm;
tan_blaze  = tand(grBlazeAngle_deg);
tan_anti   = tand(grAntiBlazeAngle_deg);
w_blaze_nm = period_nm / (1 + tan_blaze / tan_anti);
depth_nm   = w_blaze_nm * tan_blaze;

fprintf('Derived groove depth: %.4f nm\n', depth_nm);

grating = build_grating('blazed', period_nm, depth_nm, ...
                         grBlazeAngle_deg, grAntiBlazeAngle_deg, x_res_nm);


% Layer stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

substrate_file = fullfile(oc_path, 'n_Si_cxro.txt');

stack = build_stack(grating);
stack = add_layer(stack, fullfile(oc_path, 'n_Au_cxro.txt'), 31);   % 31 nm Au coating


% Sweep parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fixed photon energy: solve efficiency as a function of grazing angle.

sweep.type      = 'alpha';
sweep.values    = 0.5:0.1:8.0;    % grazing incidence angles in degrees
sweep.energy_eV = 500;            % fixed photon energy in eV


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

if isempty(results.efficiency)
    error('No valid efficiency data — check angle range and optical constant files.');
end

figure(1); clf;
plot(results.sweep_values, results.efficiency * 100, ...
     'r-o', 'LineWidth', 1.2, 'MarkerSize', 3);

xlabel('Grazing Angle (deg)', 'FontSize', 12);
ylabel('Diffraction Efficiency (%)', 'FontSize', 12);

pol_str = 'TM'; if options.pol == 1; pol_str = 'TE'; end
title(sprintf('Fixed Energy Angle Sweep | %d l/mm | %.0f eV | %s | Order %+d', ...
    grPeriod_lpermm, sweep.energy_eV, pol_str, options.GR_Order), 'FontSize', 11);

grid on;
set(gca, 'FontSize', 11);
xlim([min(results.sweep_values), max(results.sweep_values)]);

saveas(gcf, fullfile(options.output_dir, 'efficiency_fixed_energy_angle_sweep.png'));
fprintf('Plot saved.\n');
