clear;

% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Blazed grating parameters

grPeriod_lpermm      = 600;
grBlazeAngle_deg     = 0.729;
grAntiBlazeAngle_deg = 5.597;
x_res_nm             = 1;

period_nm  = 1e6 / grPeriod_lpermm;
tan_blaze  = tand(grBlazeAngle_deg);
tan_anti   = tand(grAntiBlazeAngle_deg);
w_blaze_nm = period_nm / (1 + tan_blaze / tan_anti);
depth_nm   = w_blaze_nm * tan_blaze;

fprintf('Derived groove depth : %.4f nm\n', depth_nm);

% grating = build_grating('blazed', period_nm, depth_nm, grBlazeAngle_deg, grAntiBlazeAngle_deg, x_res_nm);

% trapezoidal grating with same period
grPeriod_lpermm      = 600;
period_nm            = 1e6 / grPeriod_lpermm;
depth_nm             = 14.9;
trapezoid_angle      = 75;
ratio_groove_width_to_groove_period = 0.67;

grating = build_grating('trapezoidal', period_nm, depth_nm, trapezoid_angle, ratio_groove_width_to_groove_period, x_res_nm);


% Layer stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


stack = build_stack(grating);
stack = add_layer(stack, 'n_Pt_cxro.txt', 28.77);
% stack = add_layer(stack, 'n_C_cxro.txt', .9);     %incase of a dielectric coating on top of the grating

substrate_file = 'n_Si_cxro.txt';


% Sweep parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sweep.type   = 'energy';    % if sweep type is 'energy' fix either sweep.Cff or sweep.angle_deg at desired value. If sweep type is 'angle', fix sweep.energy_eV at desired value.  
sweep.values = 50:100:2000;
% sweep.Cff    = 2.25;    % only used if sweep.type is 'energy' and no alpha specified
sweep.alpha_deg = 4;      % only used if sweep.type is 'energy' and no Cff specified

%  Angle sweep at fixed energy:
%sweep.type      = 'alpha';
%sweep.values    = 0.5:0.1:5.0;
%sweep.energy_eV = 500;



%  4. Solver options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All are optional with defaults, but you can specify any of these to override the defaults.  See run_rcwa.m for details.
options.FourierOrders = 11;
options.pol           = -1;       % -1 = TM,  +1 = TE
options.GR_Order      = -1;
options.z_res_nm      = .5;
options.reticolo_path = fullfile(pwd, 'V9', 'reticolo_allege_v9');
options.output_dir    = pwd;
options.verbose       = true;



%  Run  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = run_rcwa(stack, substrate_file, sweep, options);



% Plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(results.efficiency)
    error('No valid efficiency data — check energy range and optical constant files.');
end

figure(1); clf;
plot(results.sweep_values, results.efficiency * 100, ...
     'b-o', 'LineWidth', 0.2, 'MarkerSize', 1);

xlabel(results.sweep_label, 'FontSize', 12);
ylabel('Diffraction Efficiency (%)', 'FontSize', 12);
title(sprintf('Blazed Grating RCWA | %d l/mm | Pt | Cff = %.2f | TM | Order %+d', ...
    grPeriod_lpermm, results.sweep_label, options.GR_Order), 'FontSize', 12);
grid on;
set(gca, 'FontSize', 11);
xlim([min(results.sweep_values), max(results.sweep_values)]);
ylim([0, max(results.efficiency * 100) * 1.15 + 1]);
saveas(gcf, "blazed_grating_results.png");
  

