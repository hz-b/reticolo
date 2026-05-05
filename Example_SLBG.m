clear;
warning('off', 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Blazed Grating RCWA Simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Grating geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.grPeriod_lpermm      = 600;
metadata.GR_Order             = -1;
metadata.grBlazeAngle_deg     = 0.729;
metadata.grAntiBlazeAngle_deg = 5.597;
metadata.material_sub         = 'Si';

% Groove depth derived from blaze and anti-blaze angles and the period.
% For a sawtooth: depth = period * tan(blaze).

p_nm_check          = 1e6 / metadata.grPeriod_lpermm;          
tan_blaze           = tand(metadata.grBlazeAngle_deg);
tan_anti            = tand(metadata.grAntiBlazeAngle_deg);
% w_b + w_a = p
% depth = w_b * tan_blaze = w_a * tan_anti  //  w_b = p / (1 + tan_blaze/tan_anti)
w_blaze_nm          = p_nm_check / (1 + tan_blaze / tan_anti);
metadata.grDepth_nm = w_blaze_nm * tan_blaze;

fprintf('Derived groove depth : %.4f nm\n', metadata.grDepth_nm);
fprintf('Blaze face width     : %.4f nm\n', w_blaze_nm);
fprintf('Anti-blaze face width: %.4f nm\n', p_nm_check - w_blaze_nm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Coating layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.material_layer      = 'Au';
metadata.layerThickness_nm   = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation resolution / Fourier orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.z_resolution_nm = 0.5;   
metadata.x_resolution_nm = 1;
metadata.FourierOrders   = 11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Photon energy range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.photonEnergy_eV = 50:5:1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Incidence angle — conditional: alpha (default) or Cff
%
%  Set alpha_deg to the desired grazing incidence angle.
%  If alpha_deg == 0 the script falls back to Cff to compute alpha from
%  the grating equation at each photon energy (energy-dependent alpha).
%
%  Cff = cos(beta) / cos(alpha)   with  sin(alpha) + sin(beta) = m*lambda/d
%  Solving: alpha = acos( sqrt( (m*lambda/d)^2 / (Cff^2 - 1) + ... ) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_deg = 0;    % grazing incidence angle in degrees  (set 0 to use Cff)
Cff       = 2.25; % used only when alpha_deg == 0

use_cff = (alpha_deg == 0);

if use_cff
    fprintf('Mode: Cff = %.4f  (alpha computed per energy from grating eq.)\n', Cff);
else
    fprintf('Mode: fixed alpha = %.4f deg\n', alpha_deg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RETICOLO setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reticolo_path = fullfile(pwd, 'V9', 'reticolo_allege_v9');
if exist(reticolo_path, 'dir') ~= 7
    error('RETICOLO path not found: %s', reticolo_path);
end
addpath(genpath(reticolo_path));
retio;

pol  = -1;   % TM polarisation  (convention: +1 = TE, -1 = TM)
nn   = metadata.FourierOrders;
a    = metadata;

eff  = [];
En   = [];
beta_out = [];   % store exit angles

meshgrid_saved = false;  % flag to save meshgrid plot only once

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Main loop over photon energies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for photonEnergy_eV = metadata.photonEnergy_eV

    a.photonEnergy_eV = photonEnergy_eV;
    lambda_nm         = 1239.8 / photonEnergy_eV;
    p_nm              = 1e6   / a.grPeriod_lpermm;   % period in nm

    % --- incidence angle resolution -------------------------------------------
    m   = abs(a.GR_Order);   % |order| for grating equation
    sgn = sign(a.GR_Order);  % -1 for inside order

    if use_cff
        % Grating equation (grazing angles):
        %   sin(alpha_g) + sin(beta_g) = m*lambda/d      [grazing convention]
        % where alpha_g, beta_g are grazing angles and d = p_nm.
        % Cff = cos(beta_g) / cos(alpha_g)
        % Combined: alpha_g = asin( (m*lambda/d - sqrt((Cff^2-1) +
        %           (m*lambda/d)^2)) / ... )  -- use the standard analytic form.

        ml_over_d = m * lambda_nm / p_nm;

        % Solve: sin(a) + sin(b) = S,  cos(b)/cos(a) = Cff         S = m*lambda/d
        % =>  sin(b) = S - sin(a),  cos(b) = Cff*cos(a)
        %   (S-x)^2 + Cff^2*(1-x^2) = 1                            x = sin(alpha_g)
        %   S^2 - 2Sx + x^2 + Cff^2 - Cff^2*x^2 = 1
        %   (1-Cff^2)*x^2 - 2S*x + (S^2 + Cff^2 - 1) = 0

        A_coeff =  1 - Cff^2;
        B_coeff = -2 * ml_over_d;
        C_coeff =  ml_over_d^2 + Cff^2 - 1;
        discriminant = B_coeff^2 - 4*A_coeff*C_coeff;

        if discriminant < 0 || A_coeff == 0
            % No real solution at this energy — skip
            warning('No real alpha solution at %.1f eV, skipping.', photonEnergy_eV);
            continue;
        end

        % Two roots; pick the smallest positive grazing angle
        x1 = (-B_coeff + sqrt(discriminant)) / (2*A_coeff);
        x2 = (-B_coeff - sqrt(discriminant)) / (2*A_coeff);

        % Valid root: |x| <= 1 and gives beta_g in (0, pi/2)
        candidates = [x1, x2];
        valid = [];

        for c = candidates
            
            if abs(c) > 1
                continue;
            end

            ag = asin(c);


            if ag <= 0 || ag >= pi/2
                continue;
            end

            % beta from grating equation
            bg_sin = ml_over_d - c;


            if bg_sin <= -1 || bg_sin >= 1
                continue;
            end

            bg = asin(bg_sin);

            if bg <= -pi/2 || bg >= pi/2
                continue;
            end

            % Check Cff consistency`
            if abs(cos(bg)/cos(ag) - Cff) > 1e-6
                continue;
            end

            % store valid solution
            valid = [valid; ag, bg];
        end

        if isempty(valid)
            warning('Could not find physical alpha at %.1f eV, skipping.', photonEnergy_eV);
            continue;
        end

        [~, idx] = min(valid(:,1)); % smallest alpha_g
        alpha_g_rad = valid(idx,1);
        current_alpha_deg = rad2deg(alpha_g_rad);
        % disp(['Using alpha = ', num2str(current_alpha_deg), ' deg at ', num2str(photonEnergy_eV), ' eV']);
    else
        current_alpha_deg = alpha_deg;
        disp(['Using fixed alpha: ', num2str(current_alpha_deg), ' deg (from Cff) at ', num2str(photonEnergy_eV), ' eV']);
    end

    % polar angle from normal
    theta0_deg = current_alpha_deg;
    k_parallel = sin(deg2rad(theta0_deg));

    % --- load optical constants -----------------------------------------------
    n_sub = NaN;
    n_HZ  = NaN;

    for i = 1:2
        if i == 1
            nfile = ['n_', a.material_sub, '_cxro.txt'];
        else
            nfile = ['n_', a.material_layer, '_cxro.txt'];
        end

        if exist(nfile, 'file') == 2
            nData = importdata(nfile);
            nData = nData.data;
        else
            error('Index file not found: %s', nfile);
        end

        n_real = interp1(nData(:,1), nData(:,2), photonEnergy_eV, 'linear', NaN);
        n_imag = interp1(nData(:,1), nData(:,3), photonEnergy_eV, 'linear', NaN);
        n_cmpl = 1 - n_real + 1i*n_imag;

        if i == 1
            n_sub = n_cmpl;
        else
            n_HZ = n_cmpl;
        end
    end

    if isnan(real(n_sub)) || isnan(real(n_HZ))
        warning('Optical constants out of range at %.1f eV, skipping.', photonEnergy_eV);
        continue;
    end

    n_inc = 1;

    % --- build blazed sawtooth profile ----------------------------------------
    %
    % Profile convention 
    %   x = 0         -> height = 0           (groove valley / start)
    %   x = x_peak    -> height = grDepth_nm  (groove tip / blaze peak)
    %   x = p_nm      -> height = 0           (next valley, hard step)
    %
    % The gentle slope (0 -> x_peak) is the blaze face.
    % The steep drop  (x_peak -> p_nm) is the anti-blaze face.
    % x_peak is derived from the two angles so both are honoured exactly.

    x_peak_nm = w_blaze_nm;   % computed above from both angles
    edge_step  = max(a.x_resolution_nm, p_nm * 1e-6);
    x_peak_nm  = min(x_peak_nm, p_nm - edge_step);

    % Discretise one period
    x = linspace(0, p_nm, round(p_nm / a.x_resolution_nm) + 1);

    % Interpolate sawtooth: linear rise then hard drop
    Prf_pts  = [0, x_peak_nm, p_nm];
    Prf_hts  = [0, a.grDepth_nm, 0];
    Prf_z    = interp1(Prf_pts, Prf_hts, x, 'linear');

    

    % --- build z stack and texture --------------------------------------------

    th_nm    = a.grDepth_nm + a.layerThickness_nm + 5;
    z        = linspace(th_nm, 0, round(th_nm / a.z_resolution_nm) + 1);
    N_layers = length(z);

    [X, Z]   = meshgrid(x, z);
    n_grid   = X .* 0;

    Prf0 = Prf_z;                            % substrate surface
    Prf1 = Prf_z + a.layerThickness_nm;      % top of coating

    % Assign refractive indices layer by layer
    P = find(Z <  Prf0);            n_grid(P) = n_sub;
    P = find(Z >= Prf0);            n_grid(P) = n_inc;
    P = find(Z >= Prf0 & Z < Prf1); n_grid(P) = n_HZ;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUALIZE MESHGRID (first iteration only)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    if ~meshgrid_saved
        figure('Name', 'Blazed Grating Meshgrid Visualization', ...
            'Position', [100, 100, 1000, 800]);
 
        imagesc(x, z, imag(n_grid));    
        axis xy;
        axis tight;
 
        colormap(jet(6));
        cb = colorbar;
        ylabel(cb, 'Im(n)  [absorption]', 'FontSize', 10);
 
        hold on;
        contour(x, z, imag(n_grid), 'k', 'LineWidth', 0.01);
 
        xlabel('x Position (nm)', 'FontSize', 12);
        ylabel('z Height (nm)',   'FontSize', 12);
        title(sprintf('Blazed Grating Cross-Section | Im(n) | %.1f eV | Pt/Au', ...
            photonEnergy_eV), 'FontSize', 12);
 
        
        set(gca, 'FontSize', 11);
 
        saveas(gcf, 'blazed_grating_meshgrid.png');
 
        disp('Saved meshgrid visualization to: blazed_grating_meshgrid.png');
 
        meshgrid_saved = true;
    end

    % --- compress into RETICOLO texture cells ---------------------------------

    deltan = diff(n_grid, 1, 2);
    [nonZeroRowIndices, nonZeroColIndices] = find(deltan ~= 0);
    ProfEdge = [nonZeroRowIndices, nonZeroColIndices];
    ProfEdge = sortrows(ProfEdge, 1);

    textures      = cell(1, N_layers + 2);
    textures{1}   = {n_inc};
    textures{end} = {n_sub};

    for layer = 1:N_layers
        p_idx = find(ProfEdge(:,1) == layer);
        if ~isempty(p_idx)
            x_position = X(layer, ProfEdge(p_idx, 2) + 1);
            n_value    = n_grid(layer, ProfEdge(p_idx, 2));
            textures{layer+1} = {x_position, n_value};
        else
            textures{layer+1} = {n_grid(layer, 1)};
        end
    end

    texture_list = 1:(N_layers + 2);
    th_list      = [0, ones(1, N_layers) .* a.z_resolution_nm, 0];
    profile      = {th_list, texture_list};

    % --- run RETICOLO ---------------------------------------------------------

    parm            = res0(pol);
    parm.res1.trace = 0;

    aa = res1(lambda_nm, p_nm, textures, nn, k_parallel, parm);
    ef = res2(aa, profile, parm);

    % --- extract -1st order efficiency ----------------------------------------
    %
    % ef.inc_top_reflected.order contains signed diffraction orders.
    % We look for order == GR_Order (i.e. -1).

    orders      = ef.inc_top_reflected.order(:, 1);
    idx_order   = find(orders == a.GR_Order);

    if isempty(idx_order)
        warning('Order %d not found in RETICOLO output at %.1f eV.', a.GR_Order, photonEnergy_eV);
        continue;
    end

    eff_val  = ef.inc_top_reflected.efficiency(idx_order);
    beta_val = 90 - ef.inc_top_reflected.theta(idx_order);  % exit grazing angle

    fprintf('E = %6.1f eV | alpha = %.3f deg | eff(order %d) = %.4f (%.2f%%) | beta = %.3f deg\n', ...
        photonEnergy_eV, current_alpha_deg, a.GR_Order, eff_val, eff_val*100, beta_val);

    eff      = [eff,      eff_val];
    En       = [En,       photonEnergy_eV];
    beta_out = [beta_out, beta_val];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(eff)
    error('No valid efficiency data computed. Check energy range and optical constant files.');
end

figure(1);
clf;

plot(En, eff * 100, 'b-o', 'LineWidth', .2, 'MarkerSize', 1, 'DisplayName', 'Simulation (TM, order -1)');

xlabel('Photon Energy (eV)', 'FontSize', 12);
ylabel('Diffraction Efficiency (%)', 'FontSize', 12);

if use_cff
    title_str = sprintf('Blazed Grating RCWA | 600 l/mm | Pt/Au | Cff = %.2f | TM | Order -1', Cff);
else
    title_str = sprintf('Blazed Grating RCWA | 600 l/mm | Pt/Au | \\alpha = %.1f° | TM | Order -1', alpha_deg);
end

title(title_str, 'FontSize', 12);
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 11);
xlim([min(En), max(En)]);
ylim([0, max(eff*100)*1.15 + 1]);

out_fig = 'blazed_grating_efficiency.png';
saveas(gcf, out_fig);
fprintf('\nPlot saved as: %s\n', out_fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Save results to CSV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_csv = 'blazed_grating_results.csv';
fid = fopen(out_csv, 'w');
fprintf(fid, 'PhotonEnergy_eV,DiffractionEfficiency,ExitAngle_beta_deg\n');
for k = 1:length(En)
    fprintf(fid, '%.4f,%.6f,%.6f\n', En(k), eff(k), beta_out(k));
end
fclose(fid);
fprintf('Results saved as: %s\n', out_csv);