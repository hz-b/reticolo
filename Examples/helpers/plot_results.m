function plot_results(stack, substrate_file, sweep, options, results)
% PLOT_RESULTS  Save the grating cross-section meshgrid and efficiency curve.
%
% USAGE
%   plot_results(stack, substrate_file, sweep, options, results)
%
% INPUTS
%   stack          struct from build_stack / add_layer
%   substrate_file CXRO file for the substrate
%   sweep          struct with fields: type ('energy' or 'alpha'), values, and
%                    alpha_deg / Cff  for energy sweeps
%                    energy_eV        for angle sweeps
%   options        struct with: output_dir, pol, GR_Order, z_res_nm
%   results        struct from run_rcwa with: efficiency, sweep_values

if isempty(results.efficiency)
    error('plot_results: no valid efficiency data.');
end

% --- photon energy for the cross-section preview ----------------------------
if strcmp(sweep.type, 'energy')
    photon_eV = sweep.values(round(end/2));
else
    photon_eV = sweep.energy_eV;
end

% --- meshgrid ----------------------------------------------------------------
meshgrid_png = fullfile(options.output_dir, ...
    sprintf('meshgrid_%s_sweep.png', sweep.type));
plot_meshgrid(stack, substrate_file, photon_eV, options.z_res_nm, meshgrid_png);

% --- efficiency plot ---------------------------------------------------------
lpermm  = round(1e6 / stack.grating.period_nm);
pol_str = 'TM'; if options.pol == 1; pol_str = 'TE'; end

if strcmp(sweep.type, 'energy')
    line_spec = 'b-o';
    x_label   = 'Photon Energy (eV)';
    if isfield(sweep, 'Cff')
        inc_tag = sprintf('Cff=%.2f', sweep.Cff);
    else
        inc_tag = sprintf('α=%.2f°', sweep.alpha_deg);
    end
    ttl = sprintf('Energy Sweep | %d l/mm | %s | %s | Order %+d', ...
        lpermm, inc_tag, pol_str, options.GR_Order);
else
    line_spec = 'r-o';
    x_label   = 'Grazing Angle (deg)';
    ttl = sprintf('Angle Sweep | %d l/mm | %.0f eV | %s | Order %+d', ...
        lpermm, sweep.energy_eV, pol_str, options.GR_Order);
end

figure('Name', 'Diffraction Efficiency'); clf;
plot(results.sweep_values, results.efficiency * 100, ...
     line_spec, 'LineWidth', 1.2, 'MarkerSize', 2);

xlabel(x_label, 'FontSize', 12);
ylabel('Diffraction Efficiency (%)', 'FontSize', 12);
title(ttl, 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11);
xlim([min(results.sweep_values), max(results.sweep_values)]);

eff_png = fullfile(options.output_dir, ...
    sprintf('efficiency_%s_sweep.png', sweep.type));
saveas(gcf, eff_png);
fprintf('Efficiency plot saved to: %s\n', eff_png);

end
