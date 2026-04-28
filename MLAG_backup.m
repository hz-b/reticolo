clear;
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     MULTILAYER GRATING SIMULATION (MLAG)
%
%     Ideal or measured profile for:
%     Lamellar grating with up to 3 coating layers
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
%     INPUT LAYER STRUCTURE (UP TO 3 LAYERS)
%     Layers are activated when thickness > 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Layer 1 (bottom-most coating, directly on grating surface)
metadata.layer1_material = 'Pt'; 
metadata.layer1_Thickness_nm = 28.77; % Set > 0 to activate
metadata.layer1_density = [20.132]; 

% Layer 2 (middle coating)
metadata.layer2_material = 'CO'; 
metadata.layer2_Thickness_nm = 10;  % Set > 0 to activate
metadata.layer2_density = [1.38]; 

% Layer 3 (top-most coating)
metadata.layer3_material = 'Cr'; 
metadata.layer3_Thickness_nm = 0;  % Set > 0 to activate
metadata.layer3_density = [7.139]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     INPUT SIMULATION PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.z_resolution_nm = .1 ;%input('z_slicing_step_in_nm: ');
metadata.x_resolution_nm = .1 ;%input('x_slicing_step_in_nm: ');
metadata.FourierOrders = 11 ;%input('Harmonics: '); % total Fourier harmonics, should be odd
metadata.GR_groove = 1; % number of periods in simulation cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     INPUT INCIDENCE ANGLE AND PHOTON ENERGY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.photonEnergy_eV = 100:20:500; %input(' which photon energy(range) wants to comput_in eV: ');:
grazing_angle_deg = 4;% input('desired grazing incidence angle range) in deg: ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    BUILD MODEL
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

% Count active layers
active_layers = 0;
if a.layer1_Thickness_nm > 0, active_layers = active_layers + 1; end
if a.layer2_Thickness_nm > 0, active_layers = active_layers + 1; end
if a.layer3_Thickness_nm > 0, active_layers = active_layers + 1; end

eff=[];
En=[];
i=1;
for photonEnergy_eV=metadata.photonEnergy_eV
    a.photonEnergy_eV=photonEnergy_eV;

% Fourier truncation used by RETICOLO: orders run from -nn to +nn
fourier_harmonics = metadata.FourierOrders;
nn = floor(fourier_harmonics / 2);

% ========== LOAD REFRACTIVE INDICES FOR ALL MATERIALS ==========
[nData_sub, nfile_sub] = load_oc_elisa_refractive_index(a.material_sub, a.reference_dir, a.material_sub_density);
disp(['Loaded substrate refractive index: ', nfile_sub]);

% Load layer 1 refractive index (if active)
if a.layer1_Thickness_nm > 0
    [nData_layer1, nfile_layer1] = load_oc_elisa_refractive_index(a.layer1_material, a.reference_dir, a.layer1_density);
    disp(['Loaded layer 1 refractive index: ', nfile_layer1]);
else
    nData_layer1 = [];
end

% Load layer 2 refractive index (if active)
if a.layer2_Thickness_nm > 0
    [nData_layer2, nfile_layer2] = load_oc_elisa_refractive_index(a.layer2_material, a.reference_dir, a.layer2_density);
    disp(['Loaded layer 2 refractive index: ', nfile_layer2]);
else
    nData_layer2 = [];
end

% Load layer 3 refractive index (if active)
if a.layer3_Thickness_nm > 0
    [nData_layer3, nfile_layer3] = load_oc_elisa_refractive_index(a.layer3_material, a.reference_dir, a.layer3_density);
    disp(['Loaded layer 3 refractive index: ', nfile_layer3]);
else
    nData_layer3 = [];
end

% ========== CALCULATE REFRACTIVE INDICES AT CURRENT ENERGY ==========
lambda_nm = 1239.8/a.photonEnergy_eV;

% Substrate
n_sub_real = interp1(nData_sub(:,1),nData_sub(:,2),a.photonEnergy_eV);
n_sub_imag = interp1(nData_sub(:,1),nData_sub(:,3),a.photonEnergy_eV);
n_sub = 1-n_sub_real+n_sub_imag*1i;

% Layer 1
if a.layer1_Thickness_nm > 0
    n_layer1_real = interp1(nData_layer1(:,1),nData_layer1(:,2),a.photonEnergy_eV);
    n_layer1_imag = interp1(nData_layer1(:,1),nData_layer1(:,3),a.photonEnergy_eV);
    n_layer1 = 1-n_layer1_real+n_layer1_imag*1i;
end

% Layer 2
if a.layer2_Thickness_nm > 0
    n_layer2_real = interp1(nData_layer2(:,1),nData_layer2(:,2),a.photonEnergy_eV);
    n_layer2_imag = interp1(nData_layer2(:,1),nData_layer2(:,3),a.photonEnergy_eV);
    n_layer2 = 1-n_layer2_real+n_layer2_imag*1i;
end

% Layer 3
if a.layer3_Thickness_nm > 0
    n_layer3_real = interp1(nData_layer3(:,1),nData_layer3(:,2),a.photonEnergy_eV);
    n_layer3_imag = interp1(nData_layer3(:,1),nData_layer3(:,3),a.photonEnergy_eV);
    n_layer3 = 1-n_layer3_real+n_layer3_imag*1i;
end

n_inc = 1;

% ========== GRATING PARAMETERS ==========
p_nm = 1/a.grPeriod_lpermm*1E6;
theta0_deg = 90 - grazing_angle_deg;
k_parallel = n_inc*sin(theta0_deg*pi/180);

% Total thickness for z-grid (grating depth + all active layers + margin)
total_layer_thickness = 0;
if a.layer1_Thickness_nm > 0, total_layer_thickness = total_layer_thickness + a.layer1_Thickness_nm; end
if a.layer2_Thickness_nm > 0, total_layer_thickness = total_layer_thickness + a.layer2_Thickness_nm; end
if a.layer3_Thickness_nm > 0, total_layer_thickness = total_layer_thickness + a.layer3_Thickness_nm; end

th_nm = a.grDepth_nm + total_layer_thickness + 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE LAMELLAR GRATING % Layer 1 (bottom-most coating, directly on grating surface)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = a.grWidthtoD * p_nm;   % groove width at bottom

if a.grTrapezoidAngL_deg == 0 && a.grTrapezoidAngR_deg == 0
    % Ideal rectangular case
    A = [W/2,          0           ];
    B = [W/2,          a.grDepth_nm];
    C = [p_nm - W/2,   a.grDepth_nm];
    
    D = [p_nm - W/2,   0           ];
else
    % Trapezoidal case
    foot_L = a.grDepth_nm / tand(a.grTrapezoidAngL_deg);
    foot_R = a.grDepth_nm / tand(a.grTrapezoidAngR_deg);

    A = [W/2,                       0           ];
    B = [W/2 + foot_L,              a.grDepth_nm];
    C = [p_nm - ((W/2)+foot_R),     a.grDepth_nm];
    D = [p_nm - (W/2),              0           ];
end

% Full Prf matrix: columns are [x, z]
Prf = [0,      0;
       A;
       B;
       C;
       D;
       p_nm,   0];

assert(all(diff(Prf(:,1)) > 0), ...
    'Prf x-coordinates are not strictly increasing. Check geometry parameters.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD SPATIAL GRIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = linspace(0, p_nm * a.GR_groove, ...
        round(p_nm * a.GR_groove / a.x_resolution_nm) + 1);
z = linspace(th_nm * a.GR_groove, 0 , ...
        round(th_nm * a.GR_groove / a.z_resolution_nm) + 1);

[X, Z] = meshgrid(x, z);
N_layers = length(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATE SURFACE HEIGHT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Prf_z = interp1(Prf(:,1), Prf(:,2), x, 'linear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN REFRACTIVE INDICES ON (Z,X) GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with everything as substrate
n = repmat(n_sub, size(X));
real_substrate = real(n_sub);

% Incident medium above all layers
n(Z >= repmat(Prf_z + total_layer_thickness, N_layers, 1)) = n_inc;
real_incident = real(n_inc);

% Layer 3 
if a.layer3_Thickness_nm > 0
    layer3_bottom = Prf_z + a.layer1_Thickness_nm;
    if a.layer2_Thickness_nm > 0
        layer3_bottom = layer3_bottom + a.layer2_Thickness_nm;
    end
    layer3_top = layer3_bottom + a.layer3_Thickness_nm;
    
    in_layer3 = (Z >= repmat(layer3_bottom, N_layers, 1)) & ...
                (Z <  repmat(layer3_top, N_layers, 1));
    n(in_layer3) = n_layer3;
    real_layer3  = real(n_layer3);
end

% Layer 2
if a.layer2_Thickness_nm > 0
    layer2_bottom = Prf_z + a.layer1_Thickness_nm;
    layer2_top = layer2_bottom + a.layer2_Thickness_nm;
    
    in_layer2 = (Z >= repmat(layer2_bottom, N_layers, 1)) & ...
                (Z <  repmat(layer2_top, N_layers, 1));
    n(in_layer2) = n_layer2;
    real_layer2  = real(n_layer2);
end

% Layer 1 (bottom-most coating, directly on grating surface)
if a.layer1_Thickness_nm > 0
    layer1_top = Prf_z + a.layer1_Thickness_nm;
    
    in_layer1 = (Z >= repmat(Prf_z, N_layers, 1)) & ...
                (Z <  repmat(layer1_top, N_layers, 1));
    n(in_layer1) = n_layer1;
    real_layer1  = real(n_layer1); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZE MESHGRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if i == 1
    figure('Name', 'MLAG Meshgrid Visualization', 'Position', [100, 100, 1000, 800]);
    

    imagesc(x, z, imag(n));
    set(gca, 'YDir', 'normal');
    set(gca, 'CLim', [min(imag(n(:))), max(imag(n(:)))]);
    colormap(jet);

    cb = colorbar;
    ylabel(cb, 'Im(n)');
    set(cb, 'FontSize', 11);

    axis tight;
    hold on;

    grid on;

    saveas(gcf, 'MLAG_meshgrid_visualization.png');
    disp('Saved meshgrid visualization to: MLAG_meshgrid_visualization.png');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD RETICOLO TEXTURES (ONE PER Z-LAYER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find lateral index discontinuities per row
deltan = diff(n, 1, 2);   % differences along x (columns)

textures        = cell(1, N_layers + 2);
textures{1}     = {n_inc};   % superstrate (homogeneous)
textures{end}   = {n_sub};   % substrate   (homogeneous)

for layer = 1:N_layers
    row_diffs = find(deltan(layer, :) ~= 0);   % column indices where n changes
    
    if isempty(row_diffs)
        % Homogeneous layer: just store the index value
        textures{layer + 1} = {n(layer, 1)};
    else
        % Patterned layer: discontinuity x-positions and left-side indices
        x_disc  = x(row_diffs + 1);     % x-position of each boundary
        n_left  = n(layer, row_diffs);  % index to the LEFT of each boundary
        textures{layer + 1} = {x_disc, n_left};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD RETICOLO PROFILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

texture_list = 1 : N_layers + 2;                          % [superstrate, layers..., substrate]
th_list      = [0, ones(1, N_layers) * a.z_resolution_nm, 0];  % thicknesses
profile      = {th_list, texture_list};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETICOLO SOLVER CALLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pol = -1;  % 1D TM polarization
parm = res0(pol);
parm.not_io = 1;
parm.FourierOrders = metadata.FourierOrders;  % Pass Fourier orders explicitly
parm.res1.trace = 0;

% res1: Calculate grating profile & k-vectors
aa = res1(lambda_nm, p_nm * a.GR_groove, textures, nn, k_parallel, parm);

% res2: Calculate diffraction efficiency
ef = res2(aa, profile, parm);

% ========== EXTRACT DIFFRACTION EFFICIENCY ==========
% Order arrangement: last index is specular (order 0), previous are diffraction orders
num_orders = length(ef.inc_top_reflected.efficiency);

% Calculate target order index based on order mapping
order_mapping = [-(num_orders-1): -1, 0];  % Orders from negative to 0
target_order = -a.GR_Order * a.GR_groove;
idx_DesignOrder = find(order_mapping == target_order);

if isempty(idx_DesignOrder)
    error('Unable to locate diffraction order %d in the RETICOLO output. Available orders: %s', ...
        target_order, mat2str(order_mapping));
end

% BOUNDS CHECK
if idx_DesignOrder > num_orders
    error('Requested order index %d exceeds available orders (%d)', idx_DesignOrder, num_orders);
end

Output(i,1) = ef.inc_top_reflected.efficiency(idx_DesignOrder);
Output(i,2) = 90 - ef.inc_top_reflected.theta(idx_DesignOrder);

eff_pct = sprintf('%.3f', Output(i,1) * 100);
disp(['En',num2str(a.photonEnergy_eV),'eV, diffraction efficiency at the first order: ', eff_pct, '%,  Diffraction angle : ',num2str(Output(i,2)),'deg'])
eff=[eff,Output(i,1)];
En=[En,a.photonEnergy_eV];

i=i+1;
end

% ========== PLOT SIMULATION RESULTS ==========
figure(20);
clf(20);
plot(En, eff, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 2, 'DisplayName', 'Simulation');
xlabel('Photon Energy (eV)', 'FontSize', 12);
ylabel('Diffraction Efficiency', 'FontSize', 12);
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 11);
drawnow;
saveas(gcf, 'Example_MLAG_efficiency.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     SAVE DATA WITH METADATA TO CSV FILE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create output filename with metadata and timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
output_filename = sprintf('MLAG_simulation_%dlmm_alpha%.1fdeg_order%d_sub_%s_layers_', ...
    metadata.grPeriod_lpermm, grazing_angle_deg, metadata.GR_Order, ...
    metadata.material_sub);

% Add layer materials to filename
layer_count = 0;
if a.layer1_Thickness_nm > 0
    layer_count = layer_count + 1;
    output_filename = [output_filename, sprintf('L%d_%s_%.1fnm', layer_count, metadata.layer1_material, metadata.layer1_Thickness_nm)];
end
if a.layer2_Thickness_nm > 0
    layer_count = layer_count + 1;
    output_filename = [output_filename, sprintf('_L%d_%s_%.1fnm', layer_count, metadata.layer2_material, metadata.layer2_Thickness_nm)];
end
if a.layer3_Thickness_nm > 0
    layer_count = layer_count + 1;
    output_filename = [output_filename, sprintf('_L%d_%s_%.1fnm', layer_count, metadata.layer3_material, metadata.layer3_Thickness_nm)];
end

if layer_count == 0
    output_filename = [output_filename, 'none'];
end

output_filename = [output_filename, '_', timestamp, '.csv'];

% Open file for writing
fid = fopen(output_filename, 'w');
if fid == -1
    error('Cannot create output file: %s', output_filename);
end

% Write metadata as comments (lines starting with #)
fprintf(fid, '# SLAG Simulation Data\n');
fprintf(fid, '# Generated with Example_MLAG.m (Multilayer Version)\n');
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
fprintf(fid, '# === LAYER CONFIGURATION ===\n');
fprintf(fid, '# active_layers: %d\n', active_layers);
fprintf(fid, '#\n');
fprintf(fid, '# === LAYER 1 (Bottom-most coating) ===\n');
fprintf(fid, '# layer1_material: %s\n', metadata.layer1_material);
fprintf(fid, '# layer1_Thickness_nm: %.2f\n', metadata.layer1_Thickness_nm);
fprintf(fid, '# layer1_density: %.3f\n', metadata.layer1_density(1));
fprintf(fid, '#\n');
fprintf(fid, '# === LAYER 2 (Middle coating) ===\n');
fprintf(fid, '# layer2_material: %s\n', metadata.layer2_material);
fprintf(fid, '# layer2_Thickness_nm: %.2f\n', metadata.layer2_Thickness_nm);
fprintf(fid, '# layer2_density: %.3f\n', metadata.layer2_density(1));
fprintf(fid, '#\n');
fprintf(fid, '# === LAYER 3 (Top-most coating) ===\n');
fprintf(fid, '# layer3_material: %s\n', metadata.layer3_material);
fprintf(fid, '# layer3_Thickness_nm: %.2f\n', metadata.layer3_Thickness_nm);
fprintf(fid, '# layer3_density: %.3f\n', metadata.layer3_density(1));
fprintf(fid, '#\n');
fprintf(fid, '# === SIMULATION PARAMETERS ===\n');
fprintf(fid, '# z_resolution_nm: %.1f\n', metadata.z_resolution_nm);
fprintf(fid, '# x_resolution_nm: %.1f\n', metadata.x_resolution_nm);
fprintf(fid, '# fourier_orders: %d\n', metadata.FourierOrders);
fprintf(fid, '# polarization: %.1f\n', pol);
fprintf(fid, '#\n');
fprintf(fid, '# === INCIDENCE CONDITIONS ===\n');
fprintf(fid, '# grazing_angle_deg: %.1f\n', grazing_angle_deg);
fprintf(fid, '# photon_energy_min_eV: %.1f\n', metadata.photonEnergy_min);
fprintf(fid, '# photon_energy_max_eV: %.1f\n', metadata.photonEnergy_max);
fprintf(fid, '# photon_energy_step_eV: %.1f\n', metadata.photonEnergy_step);
fprintf(fid, '#\n');
fprintf(fid, '# === REFRACTIVE INDEX FILES ===\n');
fprintf(fid, '# substrate_refractive_index_file: %s\n', nfile_sub);
if a.layer1_Thickness_nm > 0
    fprintf(fid, '# layer1_refractive_index_file: %s\n', nfile_layer1);
end
if a.layer2_Thickness_nm > 0
    fprintf(fid, '# layer2_refractive_index_file: %s\n', nfile_layer2);
end
if a.layer3_Thickness_nm > 0
    fprintf(fid, '# layer3_refractive_index_file: %s\n', nfile_layer3);
end
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
fprintf('  Grating: %d l/mm, %.1f° blaze, %s substrate\n', ...
    metadata.grPeriod_lpermm, metadata.grTrapezoidAngL_deg, metadata.material_sub);
fprintf('  Active Layers: %d\n', active_layers);
if a.layer1_Thickness_nm > 0
    fprintf('    Layer 1: %s, %.2f nm\n', metadata.layer1_material, metadata.layer1_Thickness_nm);
end
if a.layer2_Thickness_nm > 0
    fprintf('    Layer 2: %s, %.2f nm\n', metadata.layer2_material, metadata.layer2_Thickness_nm);
end
if a.layer3_Thickness_nm > 0
    fprintf('    Layer 3: %s, %.2f nm\n', metadata.layer3_material, metadata.layer3_Thickness_nm);
end
fprintf('  Angle: %.1f° grazing incidence\n', grazing_angle_deg);
fprintf('  Energy range: %.1f - %.1f eV (%.1f eV steps)\n', ...
    metadata.photonEnergy_min, metadata.photonEnergy_max, metadata.photonEnergy_step);
fprintf('  Data points: %d\n', length(En));
fprintf('  Max efficiency: %.4f (%.2f%%) at %.1f eV\n', ...
    max(eff), max(eff)*100, En(find(eff == max(eff), 1)));