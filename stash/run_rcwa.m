function results = run_rcwa(stack, substrate_file, sweep, options)

%  Defaults

opt = struct( ...
    'FourierOrders', 11,   ...
    'pol',           -1,   ...
    'GR_Order',      -1,   ...
    'z_res_nm',      0.5,  ...
    'reticolo_path', fullfile(pwd, 'V9', 'reticolo_allege_v9'), ...
    'output_dir',    pwd,  ...
    'verbose',       true  ...
);

if nargin >= 4 && ~isempty(options)
    fnames = fieldnames(options);
    for fi = 1:numel(fnames)
        opt.(fnames{fi}) = options.(fnames{fi});
    end
end

%  RETICOLO initialisation

if exist(opt.reticolo_path, 'dir') ~= 7
    error('run_rcwa: RETICOLO path not found:\n  %s', opt.reticolo_path);
end
addpath(genpath(opt.reticolo_path));
retio;
warning('off', 'all');


%  Sweep mode detection

use_energy_sweep = strcmpi(sweep.type, 'energy');
use_angle_sweep  = strcmpi(sweep.type, 'alpha');

if use_energy_sweep
    use_cff = isfield(sweep, 'Cff') && ~isfield(sweep, 'alpha_deg');
    if use_cff
        Cff = sweep.Cff;
    else
        sweep.alpha_deg = 90 - sweep.alpha_deg;  % default alpha for energy sweep if not specified
        fixed_alpha_deg = sweep.alpha_deg; %% convert from grazing angle to alpha convention
    end
elseif use_angle_sweep
    fixed_energy_eV = sweep.energy_eV;
else
    error('run_rcwa: sweep.type must be ''energy'' or ''alpha''');
end


% assign geometry and stack to local variables for convenience
grating   = stack.grating;
p_nm      = grating.period_nm;
nn        = opt.FourierOrders;
pol       = opt.pol;
gr_order  = opt.GR_Order;


%  Output accumulators
out_sweep = [];
out_eff   = [];
out_alpha = [];
out_beta  = [];


%  Main sweep loop

sweep_values = sweep.values;
global meshgrid_saved;  % flag to save meshgrid plot only once
meshgrid_saved = false;


for sv = sweep_values

    if use_energy_sweep
        photon_eV = sv;
        lambda_nm = 1239.8 / photon_eV;
    else  % angle sweep
        photon_eV = fixed_energy_eV;
        lambda_nm = 1239.8 / photon_eV;
    end

    % resolve alpha
    if use_energy_sweep && use_cff
        current_alpha_deg =resolve_alpha_cff(lambda_nm, p_nm, gr_order, Cff, photon_eV);
        if isnan(current_alpha_deg); continue; end
    elseif use_energy_sweep
        current_alpha_deg = fixed_alpha_deg;
    else
        current_alpha_deg = sv;   % angle sweep
    end

    k_parallel = sin(deg2rad(current_alpha_deg));

    % load substrate optical constants 
    n_sub = load_cxro(substrate_file, photon_eV);
    if isnan(real(n_sub))
        if opt.verbose
            fprintf('  Substrate OC out of range at %.1f eV, skipping.\n', photon_eV);
        end
        continue;
    end

    % load layer optical constants 
    n_layers = cell(1, numel(stack.layers));
    skip = false;
    for li = 1:numel(stack.layers)
        n_layers{li} = load_cxro(stack.layers{li}.material_file, photon_eV);
        if isnan(real(n_layers{li}))
            if opt.verbose
                fprintf('  Layer %d OC out of range at %.1f eV, skipping.\n', li, photon_eV);
            end
            skip = true; break;
        end
    end
    if skip; continue; end

    n_inc = 1;   % vacuum above
    
    %  build meshgrid and fill refractive indices 
    [textures, profile, ~] = build_reticolo_input(grating, stack.layers, n_layers, n_sub, n_inc, opt.z_res_nm);

    % RETICOLO call 
    parm            = res0(pol);
    parm.res1.trace = 0;

    aa = res1(lambda_nm, p_nm, textures, nn, k_parallel, parm);
    ef = res2(aa, profile, parm);


    % extract order efficiency
    orders    = ef.inc_top_reflected.order(:, 1);
    idx_order = find(orders == gr_order);

    if isempty(idx_order)
        if opt.verbose
            fprintf('  Order %d not found at %.1f eV, skipping.\n', gr_order, photon_eV);
        end
        continue;
    end
     
    idx_order
    eff_val  = ef.inc_top_reflected.efficiency(idx_order);
    beta_val = 90 - ef.inc_top_reflected.theta(idx_order);
    current_alpha_deg = 90 - current_alpha_deg;  % convert back to grazing angle for output

    if opt.verbose
        if use_energy_sweep
            fprintf('E = %6.1f eV | alpha = %.3f deg | eff(%+d) = %.4f (%.2f%%) | beta = %.3f deg\n', ...
                photon_eV, current_alpha_deg, gr_order, eff_val, eff_val*100, beta_val);
        else
            fprintf('alpha = %.3f deg | eff(%+d) = %.4f (%.2f%%) | beta = %.3f deg\n', ...
                current_alpha_deg, gr_order, eff_val, eff_val*100, beta_val);
        end
    end

    out_sweep = [out_sweep, sv];
    out_eff   = [out_eff,   eff_val];
    out_alpha = [out_alpha, current_alpha_deg];
    out_beta  = [out_beta,  beta_val];
end

%  Assemble results struct

results.efficiency  = out_eff;
results.alpha_deg   = 90 - out_alpha;  % convert back from alpha convention to grazing angle
results.beta_deg    = 90 - out_beta;

if use_energy_sweep
    results.sweep_values = out_sweep;
    results.sweep_label  = 'PhotonEnergy_eV';
else
    results.sweep_values = out_sweep;
    results.sweep_label  = 'GrazingAngle_deg';
end

%  CSV output

csv_name = build_csv_name(stack, substrate_file, sweep, opt);
csv_path = fullfile(opt.output_dir, csv_name);

fid = fopen(csv_path, 'w');
if use_energy_sweep
    fprintf(fid, 'PhotonEnergy_eV,GrazingAlpha_deg,DiffractionEfficiency,ExitAngle_beta_deg\n');
    for k = 1:numel(out_sweep)
        fprintf(fid, '%.4f,%.6f,%.6f,%.6f\n', out_sweep(k),  out_alpha(k), out_eff(k),  out_beta(k));
    end
else
    fprintf(fid, 'GrazingAngle_deg,PhotonEnergy_eV,DiffractionEfficiency,ExitAngle_beta_deg\n');
    for k = 1:numel(out_sweep)
        fprintf(fid, '%.6f,%.4f,%.6f,%.6f\n',  out_sweep(k), fixed_energy_eV, out_eff(k),  out_beta(k));
    end
end
fclose(fid);

results.csv_file = csv_path;
fprintf('\nResults saved to: %s\n', csv_path);

end



%  Internal helper functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [textures, profile, X, Z, n_grid] = build_reticolo_input(grating, layers, n_layers, n_sub, n_inc, z_res_nm)

% Filling order:
%   1. Everything = n_sub  (substrate)
%   2. Everything above grating surface (Prf0) = n_inc  (vacuum)
%   3. For each coating layer i: band between Prf_{i-1} and Prf_i = n_layers{i}

    x        = grating.x;
    Prf0     = grating.z_surface;   % grating surface height at each x [1 x Nx]
    depth_nm = grating.depth_nm;

    % Total z height: groove depth + all coatings + 5 nm vacuum buffer
    total_coat_nm = 0;
    for li = 1:numel(layers)
        total_coat_nm = total_coat_nm + layers{li}.thickness_nm;
    end
    th_total_nm = depth_nm + total_coat_nm + 5;

    z        = linspace(th_total_nm, 0, round(th_total_nm / z_res_nm) + 1);
    N_z      = numel(z);

    [X, Z]  = meshgrid(x, z);

    n_grid  = X .* 0;   % real zero, same size — matches original idiom

    %fill everything with substrate index first
    n_grid(:) = n_sub;

    %vacuum above grating surface 
    P = find(Z >= Prf0);
    n_grid(P) = n_inc;

    % fill each coating 
    % Prf_bot and Prf_top are [1 x Nx] surfaces that follow the grating shape.
    Prf_bot = Prf0;
    for li = 1:numel(layers)
        Prf_top = Prf_bot + layers{li}.thickness_nm;
        P = find(Z >= Prf_bot & Z < Prf_top);
        n_grid(P) = n_layers{li};
        Prf_bot = Prf_top;
    end

    global meshgrid_saved;  % flag to save meshgrid plot only once

    % Mesh plot of groove profile (optional)
    if ~meshgrid_saved

        figure('Name','Grating Meshgrid','Position',[100 100 1000 800]);

        imagesc(x, z, imag(n_grid));
        axis xy;
        axis tight;

        colormap(jet(6));
        cb = colorbar;
        ylabel(cb, 'Im(n)', 'FontSize', 10);

        hold on;
        % contour(x, z, imag(n_grid), 'k', 'LineWidth', 0.2);

        xlabel('x (nm)');
        ylabel('z (nm)');

        set(gca, 'FontSize', 11);

        saveas(gcf, 'grating_meshgrid.png');

        meshgrid_saved = true;
    end


    % compress into RETICOLO texture cells :
    % RETICOLO needs textures defined by jump positions and values.
    % So we find all the jumps in n_grid along z, then for each layer build a texture cell from jump positions and values.
    % We also build a texture cell for the incident medium (vacuum) and substrate, then assemble into a profile list.
    % then for each z-layer build a texture cell from jump positions and values.

    deltan = diff(n_grid, 1, 2);
    [nonZeroRows, nonZeroCols] = find(deltan ~= 0);
    ProfEdge = sortrows([nonZeroRows, nonZeroCols], 1);

    textures      = cell(1, N_z + 2);
    textures{1}   = {n_inc};
    textures{end} = {n_sub};

    for layer = 1:N_z
        p_idx = find(ProfEdge(:,1) == layer);
        if ~isempty(p_idx)
            x_position = X(layer, ProfEdge(p_idx, 2) + 1);   % position after jump
            n_value    = n_grid(layer, ProfEdge(p_idx, 2));   % value before jump
            textures{layer+1} = {x_position, n_value};
        else
            textures{layer+1} = {n_grid(layer, 1)};
        end
    end

    texture_list = 1:(N_z + 2);
    th_list      = [0, ones(1, N_z) .* z_res_nm, 0];
    profile      = {th_list, texture_list};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n_cmpl = load_cxro(filepath, photon_eV)
% Load a CXRO optical constants file and interpolate at photon_eV.
% Returns n = 1 - delta + i*beta.

    if exist(filepath, 'file') ~= 2
        error('load_cxro: file not found: %s', filepath);
    end

    raw  = importdata(filepath);
    data = raw.data;

    delta = interp1(data(:,1), data(:,2), photon_eV, 'linear', NaN);
    beta  = interp1(data(:,1), data(:,3), photon_eV, 'linear', NaN);

    n_cmpl = 1 - delta + 1i * beta;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alpha_deg = resolve_alpha_cff(lambda_nm, p_nm, gr_order, Cff, photon_eV)
% Compute grazing incidence angle from Cff and the grating equation.
% Returns NaN on failure.

    m         = abs(gr_order);
    ml_over_d = m * lambda_nm / p_nm;

    A =  1 - Cff^2;
    B = -2 * ml_over_d;
    C =  ml_over_d^2 + Cff^2 - 1;
    disc = B^2 - 4*A*C;

    alpha_deg = NaN;

    if disc < 0 || A == 0
        warning('resolve_alpha_cff: no real solution at %.1f eV', photon_eV);
        return;
    end

    candidates = [(-B + sqrt(disc))/(2*A), (-B - sqrt(disc))/(2*A)];
    valid = [];

    for c = candidates
        if abs(c) > 1; continue; end
        ag = asin(c);
        if ag <= 0 || ag >= pi/2; continue; end
        bg_sin = ml_over_d - c;
        if bg_sin <= -1 || bg_sin >= 1; continue; end
        bg = asin(bg_sin);
        if bg <= -pi/2 || bg >= pi/2; continue; end
        if abs(cos(bg)/cos(ag) - Cff) > 1e-6; continue; end
        valid = [valid; ag, bg];
    end

    if isempty(valid)
        warning('resolve_alpha_cff: no physical solution at %.1f eV', photon_eV);
        return;
    end

    [~, idx] = min(valid(:,1));
    alpha_deg = rad2deg(valid(idx,1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function name = build_csv_name(stack, substrate_file, sweep, opt)
% Build a descriptive CSV filename from simulation parameters.

    g = stack.grating;

    % grating type tag
    switch lower(g.type)
        case 'blazed'
            grating_tag = sprintf('blazed_b%.2fdeg_ab%.2fdeg', ...
                g.params.blaze_deg, g.params.antiblaze_deg);
        case 'trapezoidal'
            grating_tag = sprintf('trap_dc%.2f', g.params.width_ratio);
    end

    period_tag = sprintf('%.0flmm', round(1e6 / g.period_nm));
    depth_tag  = sprintf('d%.0fnm', g.depth_nm);

    % substrate label
    [~, sfname, ~] = fileparts(substrate_file);
    sparts = strsplit(sfname, '_');
    if numel(sparts) >= 2; sub_tag = sparts{2}; else; sub_tag = sfname; end

    % layer stack tag:  Au30nm_Pt10nm etc.
    layer_tag = '';
    for li = 1:numel(stack.layers)
        layer_tag = [layer_tag, stack.layers{li}.label, ...
            sprintf('%.0fnm', stack.layers{li}.thickness_nm)];
        if li < numel(stack.layers); layer_tag = [layer_tag, '_']; end
    end
    if isempty(layer_tag); layer_tag = 'bare'; end

    % polarisation
    if opt.pol == 1; pol_tag = 'TE'; else; pol_tag = 'TM'; end
    order_tag = sprintf('order%+d', opt.GR_Order);

    % sweep tag
    switch lower(sweep.type)
        case 'energy'
            if isfield(sweep, 'Cff')
                sweep_tag = sprintf('Cff%.2f_%.0f-%.0feV', sweep.Cff, ...
                    min(sweep.values), max(sweep.values));
            else
                sweep_tag = sprintf('alpha%.2fdeg_%.0f-%.0feV', sweep.alpha_deg, ...
                    min(sweep.values), max(sweep.values));
            end
        case 'alpha'
            sweep_tag = sprintf('%.0feV_alpha%.2f-%.2fdeg', sweep.energy_eV, ...
                min(sweep.values), max(sweep.values));
    end

    name = sprintf('rcwa_%s_%s_%s_%s_%s_%s_%s_%s.csv', ...
        grating_tag, period_tag, depth_tag, sub_tag, layer_tag, pol_tag, order_tag, sweep_tag);

    % Replace any characters that are awkward in filenames
    name = strrep(name, ' ', '_');
end