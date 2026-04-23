clear;
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input grating structure
%
%     ideal or measured profile for:
%     lamellar grating
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


metadata.grPeriod_lpermm = 400; %('simModel: input 1) grating period in l/mm: ')
metadata.GR_Order = 1; %('simModel:input 3)diffraction order: ')
metadata.grWidthtoD = 0.67; %('simModel:input 4)ratio of grWidth to grPeriod: ')
metadata.grDepth_nm = 14.9; %('simModel:input 5)grDepth in nm: ')
metadata.grTrapezoidAngL_deg = 75; %('simModel:input 6.1)grTrapezoid left: ')
metadata.grTrapezoidAngR_deg = 75; %('simModel:input 6.2)grTrapezoid right: ')
metadata.material_sub = 'Si';%('simModel:input 7) grating substrate meaterial(type Si for silicon): ');
metadata.material_sub_density = [2.33]; % optional OC_ELISA density suffix for substrate
metadata.reference_dir = ''; % optional override for the refractive index folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input layer structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metadata.material_layer = 'Pt'; %input('layer material, type Au for gold: ');
metadata.layerThickness_nm = 28.77; %input('single layer thickness in nm: ');
metadata.material_layer_density = [20.132]; % optional OC_ELISA density suffix for layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input simulation parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.z_resolution_nm = .1 ;%input('z_slicing_step_in_nm: ');
metadata.x_resolution_nm = 1 ;%input('x_slicing_step_in_nm: ');
metadata.FourierOrders = 15 ;%input('Harmonics: '); % total Fourier harmonics, should be odd
metadata.GR_groove = 1; % number of periods in simulation cell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input incidence angle and photon energy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.photonEnergy_eV = 60:5:2000; %input(' which photon energy(range) wants to comput_in eV: ');:
grazing_angle_deg = 4;% input('desired grazing incidence angle range) in deg: ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    build model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = metadata;
reticolo_path = fullfile(pwd, 'V9', 'reticolo_allege_v9');
if exist(reticolo_path, 'dir') ~= 7
    error('RETICOLO path not found: %s', reticolo_path);
end
addpath(genpath(reticolo_path))

% Store energy range for metadata
metadata.photonEnergy_min = min(metadata.photonEnergy_eV);
metadata.photonEnergy_max = max(metadata.photonEnergy_eV);
metadata.photonEnergy_step = metadata.photonEnergy_eV(2) - metadata.photonEnergy_eV(1);

eff=[];
En=[];
i=1;
for photonEnergy_eV=metadata.photonEnergy_eV
    a.photonEnergy_eV=photonEnergy_eV;

% Fourier truncation used by RETICOLO: orders run from -nn to +nn
fourier_harmonics = metadata.FourierOrders;
nn = floor(fourier_harmonics / 2);

% input refractive index files
[nData_sub, nfile_sub] = load_oc_elisa_refractive_index(a.material_sub, a.reference_dir, a.material_sub_density);
disp(['Loaded substrate refractive index: ', nfile_sub]);
if a.layerThickness_nm > 0
    [nData_layer, nfile_layer] = load_oc_elisa_refractive_index(a.material_layer, a.reference_dir, a.material_layer_density);
    disp(['Loaded layer refractive index: ', nfile_layer]);
else
    nData_layer = [];
    n_HZ = 1;
end

%wavelength
lambda_nm = 1239.8/a.photonEnergy_eV;

%index
n_sub_real = interp1(nData_sub(:,1),nData_sub(:,2),a.photonEnergy_eV);
n_sub_imag = interp1(nData_sub(:,1),nData_sub(:,3),a.photonEnergy_eV);
n_sub = 1-n_sub_real+n_sub_imag*1i;
if a.layerThickness_nm > 0
    n_HZ_real = interp1(nData_layer(:,1),nData_layer(:,2),a.photonEnergy_eV);
    n_HZ_imag = interp1(nData_layer(:,1),nData_layer(:,3),a.photonEnergy_eV);
    n_HZ = 1-n_HZ_real+n_HZ_imag*1i;
end
n_inc = 1;

%grating pitch
p_nm = 1/a.grPeriod_lpermm*1E6;

%angle of incidence
theta0_deg = 90 - grazing_angle_deg;
k_parallel = n_inc*sin(theta0_deg*pi/180);


% thicknes mean grooveHeight
th_nm = a.grDepth_nm+a.layerThickness_nm+5; %a.layerThickness_nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lamellar grating profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- derived geometry ---
W = a.grWidthtoD * p_nm;   % groove width at bottom

if a.grTrapezoidAngL_deg == 0 && a.grTrapezoidAngR_deg == 0
    % ideal rectangular case: walls are vertical, no lateral foot offset         % bottom flat width (split equally on both sides)
    
    A = [W/2,          0           ];
    B = [W/2,          a.grDepth_nm];
    C = [p_nm - W/2,   a.grDepth_nm];
    D = [p_nm - W/2,   0           ];
    
else
    % trapezoidal case: walls have finite slope
    foot_L = a.grDepth_nm / tand(a.grTrapezoidAngL_deg);  % horizontal extent of left wall
    foot_R = a.grDepth_nm / tand(a.grTrapezoidAngR_deg);  % horizontal extent of right wall


    A = [W/2,                       0           ];   % left foot of left wall
    B = [W/2 + foot_L,              a.grDepth_nm];  % left top of ridge
    C = [p_nm - ((W/2)+foot_R),     a.grDepth_nm];  % right top of ridge
    D = [p_nm - (W/2),              0           ];   % right foot of right wall
      
end

% full Prf matrix: columns are [x, z], rows are corner points in order
% unit cell spans exactly [0, p_nm] in x, [0, grDepth_nm] in z
Prf = [0,      0;       % left edge of period (start of left flat)
       A;                % left foot of left wall
       B;                % left top of ridge
       C;                % right top of ridge
       D;                % right foot of right wall
       p_nm,   0];       % right edge of period (end of right flat)

% sanity check: x-coordinates must be strictly increasing for interp1
assert(all(diff(Prf(:,1)) > 0), ...
    'Prf x-coordinates are not strictly increasing. Check geometry parameters.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build spatial grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = linspace(0, p_nm * a.GR_groove, ...
        round(p_nm * a.GR_groove / a.x_resolution_nm) + 1);
z = linspace(th_nm* a.GR_groove, 0 , ...
        round(th_nm * a.GR_groove / a.z_resolution_nm) + 1);

[X, Z] = meshgrid(x, z);
N_layers = length(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate Prf onto x grid -> surface height at each x position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prf_z = interp1(Prf(:,1), Prf(:,2), x, 'linear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign refractive indices on the (Z, X) grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start with everything as substrate, then overwrite upward
n = repmat(n_sub, size(X));

% background medium above the grating surface
n(Z >= repmat(Prf_z, N_layers, 1)) = n_inc;

% coating layer: a band of thickness layerThickness_nm sitting on top of the surface
Prf_z_top = Prf_z + a.layerThickness_nm;
in_layer = (Z >= repmat(Prf_z, N_layers, 1)) & ...
           (Z <  repmat(Prf_z_top, N_layers, 1));
n(in_layer) = n_HZ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build Reticolo textures: one per z-layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find lateral index discontinuities per row
deltan = diff(n, 1, 2);   % differences along x (columns)

textures        = cell(1, N_layers + 2);
textures{1}     = {n_inc};   % superstrate (homogeneous)
textures{end}   = {n_sub};   % substrate   (homogeneous)

for layer = 1:N_layers
    row_diffs = find(deltan(layer, :) ~= 0);   % column indices where n changes
    
    if isempty(row_diffs)
        % homogeneous layer: just store the index value
        textures{layer + 1} = {n(layer, 1)};
    else
        % patterned layer: discontinuity x-positions and left-side indices
        x_disc  = x(row_diffs + 1);     % x-position of each boundary
        n_left  = n(layer, row_diffs);  % index to the LEFT of each boundary
        textures{layer + 1} = {x_disc, n_left};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build Reticolo profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

texture_list = 1 : N_layers + 2;                          % [superstrate, layers..., substrate]
th_list      = [0, ones(1, N_layers) * a.z_resolution_nm, 0];  % thicknesses
profile      = {th_list, texture_list};


pol = -1;  %polarisation -1 -> TM, 1-> TE the parameter is a multipupose argumen for both dimensionality and polarization 1 for 1D and 2 for 2D; + for TE and - for TM
parm = res0(pol);    % -1 → else branch → 1D mode, TM implied by sign convention
parm.not_io = 1;
param.fourier_orders = 3;
% initialisation
parm.res1.trace =  0; % refrax curve for each layer
aa = res1(lambda_nm, p_nm * a.GR_groove, textures, nn, k_parallel, parm);
ef = res2(aa, profile, parm);

% % ========== PLOT GRATING PROFILE ==========
% figure('Name', 'Laminar Grating Profile', 'Position', [100, 100, 800, 600]);
% plot(x, Prf_z, 'b-', 'LineWidth', 2);
% hold on;
% xlabel('Position (nm)', 'FontSize', 12);
% ylabel('Height (nm)', 'FontSize', 12);
% title('Laminar Grating Profile', 'FontSize', 14);
% grid on;
% saveas(gcf, 'SLAG_Laminar_profile.png');




idx_DesignOrder = nn;

if isempty(idx_DesignOrder)
    error('Unable to locate diffraction order %d in the RETICOLO output.', target_diffraction_order);
end

% BOUNDS CHECK - Prevent out-of-bounds error
num_orders = length(ef.inc_top_reflected.efficiency);
if idx_DesignOrder > num_orders
    fprintf('WARNING: Requested order %d, but only %d orders available.\n', ...
        target_diffraction_order, num_orders);
    fprintf('Available efficiency vectors: %s\n', mat2str(ef.inc_top_reflected.efficiency));
    error('Order index out of bounds.');
end

Output(i,1) = ef.inc_top_reflected.efficiency(idx_DesignOrder);
Output(i,2) = 90 - ef.inc_top_reflected.theta(idx_DesignOrder);

% disp(ef.inc_top_reflected.efficiency)
eff_pct = sprintf('%.3f', Output(i,1) * 100);

disp(['En',num2str(a.photonEnergy_eV),'eV, diffraction efficiency at the first order: ', eff_pct, '%,  Diffraction angle : ',num2str(Output(i,2)),'deg'])
eff=[eff,Output(i,1)];
En=[En,a.photonEnergy_eV];

i=i+1;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %     LOAD EXPERIMENTAL DATA
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Load experimental data from CSV file
% exp_filename = 'Re__ELISA,_400l_mm_laminar_grating_from_HORIBA/lG400-HZB-ELISA_ascan-energy_alpha-4deg_1-order.csv';

% % Read the file, skipping header lines
% fid = fopen(exp_filename, 'r');
% if fid == -1
%     error('Cannot open experimental data file: %s', exp_filename);
% end

% fgetl(fid);  % Energy;alpha = 4 deg
% fgetl(fid);  % eV;
% fgetl(fid);  % ;1 order

% % Read data with semicolon separator and comma decimal
% exp_data = [];
% while ~feof(fid)
%     line = fgetl(fid);
%     if ~isempty(line) && line(1) ~= ';'
%         % Parse line: Energy;Efficiency (with comma as decimal separator)
%         parts = strsplit(line, ';');
%         if length(parts) >= 2
%             energy = str2double(strrep(parts{1}, ',', '.'));
%             efficiency = str2double(strrep(parts{2}, ',', '.'));
%             if ~isnan(energy) && ~isnan(efficiency)
%                 exp_data = [exp_data; energy, efficiency];
%             end
%         end
%     end
% end
% fclose(fid);

% exp_energy = exp_data(:,1);
% exp_efficiency = exp_data(:,2);



figure(20);
clf(20);

% Plot simulation data
plot(En, eff, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 2, 'DisplayName', 'Simulation');
% hold on;

% Plot experimental data
% plot(exp_energy, exp_efficiency, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 2, 'DisplayName', 'Experimental Data');

% Add labels and legend
xlabel('Photon Energy (eV)', 'FontSize', 12);
ylabel('Diffraction Efficiency', 'FontSize', 12);
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 11);
drawnow;
saveas(gcf, 'Example_SLAG_efficiency.png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     SAVE DATA WITH METADATA TO CSV FILE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create output filename with metadata and timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
output_filename = sprintf('SLAG_simulation_%dlmm_alpha%.1fdeg_order%d_sub_%s_layer_%s_%s.csv', ...
    metadata.grPeriod_lpermm, grazing_angle_deg, metadata.GR_Order, ...
    metadata.material_sub, metadata.material_layer, timestamp);

% Open file for writing
fid = fopen(output_filename, 'w');
if fid == -1
    error('Cannot create output file: %s', output_filename);
end

% Write metadata as comments (lines starting with #)
fprintf(fid, '# SLAG Simulation Data\n');
fprintf(fid, '# Generated with Example_SLAG_save_data.m\n');
fprintf(fid, '#\n');
fprintf(fid, '# === GRATING PARAMETERS ===\n');
fprintf(fid, '# grating_period_lpermm: %.1f\n', metadata.grPeriod_lpermm);
fprintf(fid, '# diffraction_order: %d\n', metadata.GR_Order);
fprintf(fid, '# grating_width_to_period_ratio: %.2f\n', metadata.grWidthtoD);
fprintf(fid, '# grating_depth_nm: %.1f\n', metadata.grDepth_nm);
fprintf(fid, '# trapezoid_angle_left_deg: %.1f\n', metadata.grTrapezoidAngL_deg);
fprintf(fid, '# trapezoid_angle_right_deg: %.1f\n', metadata.grTrapezoidAngR_deg);
fprintf(fid, '# substrate_material: %s\n', metadata.material_sub);
fprintf(fid, '# substrate_density: %.3f\n', metadata.material_sub_density(1));
fprintf(fid, '#\n');
fprintf(fid, '# === LAYER PARAMETERS ===\n');
fprintf(fid, '# layer_material: %s\n', metadata.material_layer);
fprintf(fid, '# layer_thickness_nm: %.1f\n', metadata.layerThickness_nm);
fprintf(fid, '# layer_density: %.3f\n', metadata.material_layer_density(1));
fprintf(fid, '#\n');
fprintf(fid, '# === SIMULATION PARAMETERS ===\n');
fprintf(fid, '# z_resolution_nm: %.1f\n', metadata.z_resolution_nm);
fprintf(fid, '# x_resolution_nm: %.1f\n', metadata.x_resolution_nm);
fprintf(fid, '# fourier_orders: %d\n', metadata.FourierOrders);
fprintf(fid, '#\n');
fprintf(fid, '# === INCIDENCE CONDITIONS ===\n');
fprintf(fid, '# grazing_angle_deg: %.1f\n', grazing_angle_deg);
fprintf(fid, '# photon_energy_min_eV: %.1f\n', metadata.photonEnergy_min);
fprintf(fid, '# photon_energy_max_eV: %.1f\n', metadata.photonEnergy_max);
fprintf(fid, '# photon_energy_step_eV: %.1f\n', metadata.photonEnergy_step);
fprintf(fid, '# polarization: %.1f\n', pol);
fprintf(fid, '#\n');
fprintf(fid, '# === DATA COLUMNS ===\n');
fprintf(fid, '# Column 1: photon_energy_eV\n');
fprintf(fid, '# Column 2: diffraction_efficiency\n');
fprintf(fid, '#\n');
fprintf(fid, '# photon_energy_eV;diffraction_efficiency\n');

% Write data with semicolon separator
for j = 1:length(En)
    fprintf(fid, '%.4f;%.8f\n', En(j), eff(j));
end

try
    close all;  % Close all open figures
catch
    % Ignore errors if no figures exist
end

disp('');
disp(['Simulation complete!']);
disp(['Data saved to: ', output_filename]);
disp('');
disp('use plot_simulation_data.py to visualize the results.');
disp('');

% Display summary
fprintf('Simulation Summary:\n');
fprintf('  Grating: %d l/mm, %.1f° blaze, %s substrate, %s layer\n', ...
    metadata.grPeriod_lpermm, metadata.grTrapezoidAngL_deg, metadata.material_sub, metadata.material_layer);
fprintf('  Angle: %.1f° grazing incidence\n', grazing_angle_deg);
fprintf('  Energy range: %.1f - %.1f eV (%.1f eV steps)\n', ...
    metadata.photonEnergy_min, metadata.photonEnergy_max, metadata.photonEnergy_step);
fprintf('  Data points: %d\n', length(En));
fprintf('  Max efficiency: %.4f (%.2f%%) at %.1f eV\n', ...
    max(eff), max(eff)*100, En(find(eff == max(eff), 1)));
