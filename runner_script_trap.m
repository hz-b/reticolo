clear;
warning('off', 'all')


% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolution parameters--------------------------------

x_res_nm             = .1;
y_res_nm             = .1;


% Blazed grating parameters--------------------------------

blaze_grPeriod_lpermm      = 150;
grBlazeAngle_deg     = 0.729;
grAntiBlazeAngle_deg = 5.597;


blaze_period_nm  = 1e6 / blaze_grPeriod_lpermm;
tan_blaze  = tand(grBlazeAngle_deg);
tan_anti   = tand(grAntiBlazeAngle_deg);
w_blaze_nm = blaze_period_nm / (1 + tan_blaze / tan_anti);
depth_nm   = w_blaze_nm * tan_blaze;

fprintf('Derived groove depth : %.4f nm\n', depth_nm);

% grating = build_grating('blazed', blaze_period_nm, depth_nm, grBlazeAngle_deg, grAntiBlazeAngle_deg, x_res_nm);


% trapezoidal grating parameters--------------------------------

laminar_grPeriod_lpermm      = 150;
laminar_period_nm      = 1e6 / laminar_grPeriod_lpermm;
depth_nm             = 60;
trapezoid_angle      = 10;
ratio_groove_width_to_groove_period = 0.65;

grating = build_grating('trapezoidal', laminar_period_nm, depth_nm ,ratio_groove_width_to_groove_period, trapezoid_angle, x_res_nm);


% Layer stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
substrate_file = 'n_Si_cxro.txt';

stack = build_stack(grating); 
% stack = add_layer(stack, 'n_SiO2_cxro.txt', 1);   %adhesion layer
stack = add_layer(stack, 'n_Au_cxro.txt', 30); 
% stack = add_layer(stack, 'n_Pt_cxro.txt', 28.77);
% stack = add_layer(stack, 'n_C_cxro.txt', .7);     %carbon contamination layer on top of the grating



% Sweep parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sweep.type   = 'energy';    % if sweep type is 'energy' fix either sweep.Cff or sweep.angle_deg at desired value. If sweep type is 'angle', fix sweep.energy_eV at desired value.  
sweep.values = 50:50:1000;
sweep.Cff    = 1.45;    % only used if sweep.type is 'energy' and no alpha specified
% sweep.alpha_deg = 4;      % only used if sweep.type is 'energy' and no Cff specified

%  Angle sweep at fixed energy:
%sweep.type      = 'alpha';
%sweep.values    = 0.5:0.1:5.0;
%sweep.energy_eV = 500;



%  4. Solver options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All are optional with defaults, but you can specify any of these to override the defaults.  See run_rcwa.m for details.
options.FourierOrders = 41;
options.pol           = -1;       % -1 = TM,  +1 = TE
options.GR_Order      = -1;
options.z_res_nm      = y_res_nm; %default is 0.1 nm
options.reticolo_path = fullfile(pwd, 'V9', 'reticolo_allege_v9');
options.output_dir    = pwd;
options.verbose       = true;



%  Run  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = run_rcwa(stack, substrate_file, sweep, options);



% Plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(results.efficiency)
    error('No valid efficiency data — check energy range and optical constant files.');
end

plot_meshgrid(stack, 'n_Si_cxro.txt', 500, 0.01, 'meshgrid_500eV.png');


figure(1); clf;
plot(results.sweep_values, results.efficiency * 100, ...
     'b-o', 'LineWidth', 0.2, 'MarkerSize', 1);
 
xlabel(results.sweep_label, 'FontSize', 12);
ylabel('Diffraction Efficiency (%)', 'FontSize', 12);
 
pol_str = 'TM'; if options.pol == 1; pol_str = 'TE'; end
 
% Build layer string from stack  e.g. "Pt 28.8nm / C 0.9nm"
layer_str = '';
for li = 1:numel(stack.layers)
    layer_str = [layer_str, sprintf('%s %.1fnm', stack.layers{li}.label, stack.layers{li}.thickness_nm)];
    if li < numel(stack.layers); layer_str = [layer_str, ' / ']; end
end
 
% Sweep-mode-dependent incidence descriptor
if strcmpi(sweep.type, 'energy')
    if isfield(sweep, 'Cff')
        inc_tag = sprintf('Cff=%.2f', sweep.Cff);
    else
        inc_tag = sprintf('alpha=%.2fdeg', sweep.alpha_deg);
    end
elseif strcmpi(sweep.type, 'bragg')
    inc_tag = sprintf('bragg_%.0f-%.0feV', min(sweep.values), max(sweep.values));
else  % 'alpha' sweep
    inc_tag = sprintf('E=%.0feV', sweep.energy_eV);
end
 
title_str = sprintf('%s | %d l/mm | %s | %s | %s | Order %+d', ...
    grating.type, blaze_period_nm, layer_str, inc_tag, pol_str, options.GR_Order);
 
title(title_str, 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);
xlim([min(results.sweep_values), max(results.sweep_values)]);
ylim([0, max(results.efficiency * 100) * 1.15 + 1]);
 
png_name = sprintf('rcwa_%s_%dlmm_%s_%s_order%+d.png', ...
    grating.type, blaze_period_nm, inc_tag, pol_str, options.GR_Order);
saveas(gcf, png_name);
fprintf('Plot saved to: %s\n', png_name);