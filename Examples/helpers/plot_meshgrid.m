function plot_meshgrid(stack, substrate_file, photon_eV, z_res_nm, save_png)
% PLOT_MESHGRID  Visualise the refractive-index cross-section of the stack.
%
% Reproduces the Im(n) visualisation from the original monolithic script.
% Calls the same build_reticolo_input internals used by run_rcwa so what
% you see is exactly what the solver receives.
%
% USAGE
%   plot_meshgrid(stack, substrate_file, photon_eV)
%   plot_meshgrid(stack, substrate_file, photon_eV, z_res_nm)
%   plot_meshgrid(stack, substrate_file, photon_eV, z_res_nm, save_png)
%
% INPUTS
%   stack          struct from build_stack / add_layer
%   substrate_file CXRO file for the substrate
%   photon_eV      photon energy at which to evaluate optical constants
%   z_res_nm       z slice thickness in nm  (default 0.5)
%   save_png       filename to save PNG, e.g. 'meshgrid.png'
%                  omit or pass '' to skip saving

if nargin < 4 || isempty(z_res_nm); z_res_nm = 0.5; end
if nargin < 5; save_png = ''; end

% --- load optical constants --------------------------------------------------
n_sub = load_cxro_local(substrate_file, photon_eV);
if isnan(real(n_sub))
    error('plot_meshgrid: substrate OC out of range at %.1f eV', photon_eV);
end

n_layers = cell(1, numel(stack.layers));
for li = 1:numel(stack.layers)
    n_layers{li} = load_cxro_local(stack.layers{li}.material_file, photon_eV);
    if isnan(real(n_layers{li}))
        error('plot_meshgrid: layer %d OC out of range at %.1f eV', li, photon_eV);
    end
end

n_inc = 1;
grating = stack.grating;

% --- build meshgrid ----------------------------------------------------------
x    = grating.x;
Prf0 = grating.z_surface;

total_coat_nm = 0;
for li = 1:numel(stack.layers)
    total_coat_nm = total_coat_nm + stack.layers{li}.thickness_nm;
end
th_total_nm = grating.depth_nm + total_coat_nm + 5;

z      = linspace(th_total_nm, 0, round(th_total_nm / z_res_nm) + 1);
[X, Z] = meshgrid(x, z);
n_grid = X .* 0;

% Fill: substrate -> vacuum above surface -> coating layers
n_grid(:) = n_sub;
P = find(Z >= Prf0);  n_grid(P) = n_inc;

Prf_bot = Prf0;
for li = 1:numel(stack.layers)
    Prf_top = Prf_bot + stack.layers{li}.thickness_nm;
    P = find(Z >= Prf_bot & Z < Prf_top);
    n_grid(P) = n_layers{li};
    Prf_bot = Prf_top;
end

% --- plot --------------------------------------------------------------------
figure('Name', 'Grating Cross-Section', 'Position', [100, 100, 1000, 800]);

imagesc(x, z, imag(n_grid));
axis xy; axis tight;
colormap(jet(256));
cb = colorbar;
ylabel(cb, 'Im(n)  [absorption]', 'FontSize', 10);

hold on;
% contour(x, z, imag(n_grid), 'k', 'LineWidth', 0.5);

xlabel('x  (nm)', 'FontSize', 12);
ylabel('z  (nm)', 'FontSize', 12);

% Build a descriptive title from stack parameters
layer_str = '';
for li = 1:numel(stack.layers)
    layer_str = [layer_str, sprintf('%s %.0fnm', ...
        stack.layers{li}.label, stack.layers{li}.thickness_nm)];
    if li < numel(stack.layers); layer_str = [layer_str, ' / ']; end
end
if isempty(layer_str); layer_str = 'bare'; end

title(sprintf('%s grating | %.0f l/mm | %s | %.1f eV', ...
    stack.grating.type, round(1e6 / stack.grating.period_nm), ...
    layer_str, photon_eV), 'FontSize', 12);

set(gca, 'FontSize', 11);

if ~isempty(save_png)
    saveas(gcf, save_png);
    fprintf('Meshgrid plot saved to: %s\n', save_png);
end

end


% -------------------------------------------------------------------------
%  Local copy of load_cxro so plot_meshgrid is self-contained
% -------------------------------------------------------------------------
function n_cmpl = load_cxro_local(filepath, photon_eV)
    if exist(filepath, 'file') ~= 2
        error('load_cxro: file not found: %s', filepath);
    end
    raw   = importdata(filepath);
    data  = raw.data;
    delta = interp1(data(:,1), data(:,2), photon_eV, 'linear', NaN);
    beta  = interp1(data(:,1), data(:,3), photon_eV, 'linear', NaN);
    n_cmpl = 1 - delta + 1i * beta;
end