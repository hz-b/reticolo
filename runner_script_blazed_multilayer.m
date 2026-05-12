clear;
warning('off', 'all')

% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolution parameters--------------------------------

x_res_nm             = .1;
y_res_nm             = .01;

% Blazed grating parameters--------------------------------

blaze_grPeriod_lpermm      = 2400;
grBlazeAngle_deg     = 1.37;
grAntiBlazeAngle_deg = 3.25;

blaze_period_nm  = 1e6 / blaze_grPeriod_lpermm;
tan_blaze  = tand(grBlazeAngle_deg);
tan_anti   = tand(grAntiBlazeAngle_deg);
w_blaze_nm = blaze_period_nm / (1 + tan_blaze / tan_anti);
depth_nm   = w_blaze_nm * tan_blaze;

fprintf('Derived groove depth : %.4f nm\n', depth_nm);

grating = build_grating('blazed', blaze_period_nm, depth_nm, grBlazeAngle_deg, grAntiBlazeAngle_deg, x_res_nm);

% Layer stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layers_per_period = 2;
num_multilayers = 4;

layer_a = 'n_Cr_cxro.txt';
layer_a_thickness_nm = 1.9;
layer_b = 'n_C_cxro.txt';
layer_b_thickness_nm = 2.9;

substrate_file = 'n_Si_cxro.txt';

stack = build_stack(grating);
stack = add_layer(stack, 'n_SiO2_cxro.txt', 1);

for i = 1:(layers_per_period * num_multilayers)
    stack = add_layer(stack, layer_a, layer_a_thickness_nm);
    stack = add_layer(stack, layer_b, layer_b_thickness_nm);
end

% Sweep parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lookup_file = 'bragg_lookup_table.csv';

% FIX: Octave-compatible CSV reading (no readtable)
if isfile(lookup_file)
    % Skip header row, read numeric data
    data = dlmread(lookup_file, ',', 1, 0);  % Skip header row (row 1)
    lookup_energy = data(:, 1);      % Energy column (first column)
    lookup_alpha = data(:, 4);       % alpha column (4th column)
    
    fprintf('Loaded lookup table: %.1f - %.1f eV (%d points)\n', ...
        min(lookup_energy), max(lookup_energy), length(lookup_energy));
else
    error('Lookup table file not found: %s', lookup_file);
end

% Energy range within lookup table bounds
energy_start = 500;
energy_end = 538;
energy_step = 2;

sweep_energy = energy_start:energy_step:energy_end;
sweep_alpha = interp1(lookup_energy, lookup_alpha, sweep_energy, 'linear', 'extrap');

% Validate interpolation
if any(isnan(sweep_alpha))
    warning('Some alpha values are NaN - check energy range');
end

% Bragg sweep configuration
sweep.type = 'bragg';
sweep.values = sweep_energy;        % Energy array (eV)
sweep.alpha_deg = sweep_alpha;      % Grazing angle array (degrees)

fprintf('Bragg sweep: %d points from %.1f to %.1f eV\n', ...
    length(sweep_energy), min(sweep_energy), max(sweep_energy));

% Solver options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.FourierOrders = 15;
options.pol           = -1;       % -1 = TM,  +1 = TE
options.GR_Order      = -1;       % Verify this matches your lookup table
options.z_res_nm      = y_res_nm;
options.reticolo_path = fullfile(pwd, 'V9', 'reticolo_allege_v9');
options.output_dir    = pwd;
options.verbose       = true;

% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = run_rcwa(stack, substrate_file, sweep, options);

% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(results.efficiency)
    error('No valid efficiency data — check energy range and optical constant files.');
end

figure(1); clf;
plot(results.sweep_values, results.efficiency * 100, ...
     'b-o', 'LineWidth', 0.2, 'MarkerSize', 1);

xlabel('Photon Energy (eV)', 'FontSize', 12);
ylabel('Diffraction Efficiency (%)', 'FontSize', 12);

pol_str = 'TM'; if options.pol == 1; pol_str = 'TE'; end

% Build layer string from stack
layer_str = '';
for li = 1:numel(stack.layers)
    layer_str = [layer_str, sprintf('%s %.1fnm', stack.layers{li}.label, stack.layers{li}.thickness_nm)];
    if li < numel(stack.layers); layer_str = [layer_str, ' / ']; end
end

% Handle bragg sweep type in incidence descriptor
if strcmpi(sweep.type, 'bragg')
    inc_tag = sprintf('Bragg_%.0f-%.0feV', min(sweep_energy), max(sweep_energy));
elseif strcmpi(sweep.type, 'energy')
    if isfield(sweep, 'Cff')
        inc_tag = sprintf('Cff=%.2f', sweep.Cff);
    else
        inc_tag = sprintf('alpha=%.2fdeg', sweep.alpha_deg);
    end
else
    inc_tag = sprintf('E=%.0feV', sweep.energy_eV);
end

% Use blaze_grPeriod_lpermm instead of laminar_grPeriod_lpermm
title_str = sprintf('%s | %d l/mm | %s | %s | Order %+d', ...
    grating.type, blaze_grPeriod_lpermm, layer_str, inc_tag, options.GR_Order);

title(title_str, 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);
xlim([min(results.sweep_values), max(results.sweep_values)]);
ylim([0, max(results.efficiency * 100) * 1.15 + 1]);

png_name = sprintf('rcwa_%s_%dlmm_%s_%s_order%+d.png', ...
    grating.type, blaze_grPeriod_lpermm, inc_tag, pol_str, options.GR_Order);
saveas(gcf, png_name);
fprintf('Plot saved to: %s\n', png_name);