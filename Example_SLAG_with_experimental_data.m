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
metadata.grTrapezoidAngL_deg = 20; %('simModel:input 6.1)grTrapezoid left: ')
metadata.grTrapezoidAngR_deg = 20; %('simModel:input 6.2)grTrapezoid right: ')
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
metadata.FourierOrders = 5 ;%input('Harmonics: '); %or FourierOrders = GR_Order*GR_groove+3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input incidence angle and photon energy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.photonEnergy_eV = 140:5:160; %input(' which photon energy(range) wants to comput_in eV: ');:
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
retio
eff=[];
En=[];
i=1;
for photonEnergy_eV=metadata.photonEnergy_eV
    a.photonEnergy_eV=photonEnergy_eV;
%number of Fourier orders
nn = a.FourierOrders;
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
theta0_deg = 90-grazing_angle_deg;
k_parallel = n_inc*sin(theta0_deg*pi/180);


% thicknes mean grooveHeight
th_nm = a.grDepth_nm+a.layerThickness_nm+5; %a.layerThickness_nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lamellar grating profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1/tand(a.grTrapezoidAngL_deg) && 1/tand(a.grTrapezoidAngR_deg) == 0 % ideal case
    edge_z = a.grDepth_nm;
    edge_x = a.grWidthtoD*p_nm;
else
    pos1 = [(p_nm-a.grWidthtoD*p_nm)/2-a.grDepth_nm/tand(a.grTrapezoidAngL_deg),0];
    pos2 = [(p_nm-a.grWidthtoD*p_nm)/2,a.grDepth_nm];
    pos3 = [(p_nm+a.grWidthtoD*p_nm)/2,a.grDepth_nm];
    pos4 = [(p_nm+a.grWidthtoD*p_nm)/2+a.grDepth_nm/tand(a.grTrapezoidAngR_deg),0];
    Prf = [0,0;pos1;pos2;pos3;pos4;p_nm,0];
end

a.GR_groove = 1;
x = linspace(0,p_nm*a.GR_groove,round(p_nm*a.GR_groove/a.x_resolution_nm)+1);%
z = linspace(th_nm,0,round(th_nm/a.z_resolution_nm)+1);%

% genergy grid
[X,Z] = meshgrid(x,z);

% number of z layers
N_layers = length(z);

% mat to store n
n = X.*0;

if 1/tand(a.grTrapezoidAngL_deg) && 1/tand(a.grTrapezoidAngR_deg) == 0
    P = find(x<= p_nm/2-edge_x/2 | x>= p_nm/2+edge_x/2);
    Prf_z(P) = 0;
    P = [];
    P = find(x>p_nm/2-edge_x/2 & x<= p_nm/2+edge_x/2);
    Prf_z(P) = edge_z;
else
    Prf_z = interp1(Prf(:,1)',Prf(:,2)',x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate testure for single layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prf0 = Prf_z;
Prf1 = Prf_z+a.layerThickness_nm;
P = find(Z<Prf0);         n(P) = n_sub;%substrate
P = find(Z>= Prf0);         n(P) = n_inc;% background
P = find(Z>= Prf0&Z<Prf1);  n(P) = n_HZ; % single layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reduce unecessaty points in texture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltan = diff(n,1,2);
[nonZeroRowIndices, nonZeroColIndices] = find(deltan~=  0);
ProfEdge = [nonZeroRowIndices,nonZeroColIndices];
ProfEdge = sortrows(ProfEdge, 1);

textures = cell(1,N_layers+2);
textures{1} = {n_inc};
textures{end} = {n_sub};

for layer = 1:N_layers
    if ~isnan(find(ProfEdge(:,1) == layer))
        p = find(ProfEdge(:,1) == layer);
        x_position = X(layer,ProfEdge(p,2)+1);
        n_value = n(layer,ProfEdge(p,2));
        textures{layer+1} = {x_position, n_value};
    else
        n_value = n(layer,1);
        textures{layer+1} = {n_value};
    end
end

texture_list = 1:N_layers+2;
th_list = [0,ones(1,N_layers).*a.z_resolution_nm,0];
profile = {th_list,texture_list};% starting1, from incidence medium


pol = 1 ; % 1:TE   -1:TM
parm = res0(pol);  % parameter init

% initialisation
parm.res1.trace =  0; % refrax curve for each layer
aa = res1(lambda_nm, p_nm*a.GR_groove, textures, nn, k_parallel, parm);
ef = res2(aa, profile, parm);
idx_DesignOrder=find(a.FourierOrders:-1:0==a.GR_Order*a.GR_groove);
Output(i,1)=ef.inc_top_reflected.efficiency(idx_DesignOrder);
Output(i,2)=90-ef.inc_top_reflected.theta(idx_DesignOrder);
eff_pct = sprintf('%.3f', Output(i,1) * 100);
disp(['En',num2str(a.photonEnergy_eV),'eV, diffraction efficiency at the first order: ', eff_pct, '%,  Diffraction angle : ',num2str(Output(i,2)),'deg'])
eff=[eff,Output(i,1)];
En=[En,a.photonEnergy_eV];
i=i+1;
figure(10)
clf(10);
plot(En,eff,'*-');
xlabel('photonEnergy,eV');
ylabel('Diffraction efficiency');
grid on;
drawnow;
saveas(gcf, 'Example_SLAG_efficiency_original.png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     LOAD EXPERIMENTAL DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load experimental data from CSV file
exp_filename = 'Re__ELISA,_400l_mm_laminar_grating_from_HORIBA/lG400-HZB-ELISA_ascan-energy_alpha-4deg_1-order.csv';

% Read the file, skipping header lines
fid = fopen(exp_filename, 'r');
if fid == -1
    error('Cannot open experimental data file: %s', exp_filename);
end

% Skip header lines (first 3 lines are headers)
fgetl(fid);  % Energy;alpha = 4 deg
fgetl(fid);  % eV;
fgetl(fid);  % ;1 order

% Read data with semicolon separator and comma decimal
exp_data = [];
while ~feof(fid)
    line = fgetl(fid);
    if ~isempty(line) && line(1) ~= ';'
        % Parse line: Energy;Efficiency (with comma as decimal separator)
        parts = strsplit(line, ';');
        if length(parts) >= 2
            energy = str2double(strrep(parts{1}, ',', '.'));
            efficiency = str2double(strrep(parts{2}, ',', '.'));
            if ~isnan(energy) && ~isnan(efficiency)
                exp_data = [exp_data; energy, efficiency];
            end
        end
    end
end
fclose(fid);

exp_energy = exp_data(:,1);
exp_efficiency = exp_data(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     PLOT SIMULATION AND EXPERIMENTAL DATA TOGETHER
%     Title is automatically set from grazing_angle_deg variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create title string from grazing_angle_deg variable
plot_title = sprintf('SLAG Simulation vs Experimental Data (400 l/mm Grating, \\alpha = %.1f°)', grazing_angle_deg);

figure(20);
clf(20);

% Plot simulation data
plot(En, eff, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Simulation');
hold on;

% Plot experimental data
plot(exp_energy, exp_efficiency, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Experimental Data');

% Add labels and legend
xlabel('Photon Energy (eV)', 'FontSize', 12);
ylabel('Diffraction Efficiency', 'FontSize', 12);
title(plot_title, 'FontSize', 14);
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% Save the figure with angle in filename
filename = sprintf('Example_SLAG_efficiency_with_experimental_%.1fdeg.png', grazing_angle_deg);
saveas(gcf, filename);

disp(['Plot saved as ', filename]);

% Display comparison statistics
if length(En) > 0 && length(exp_energy) > 0
    % Interpolate simulation data to match experimental energy points
    sim_eff_at_exp = interp1(En, eff, exp_energy, 'linear', 'extrap');
    
    % Calculate mean absolute error
    mae = mean(abs(sim_eff_at_exp - exp_efficiency));
    rmse = sqrt(mean((sim_eff_at_exp - exp_efficiency).^2));
    
    fprintf('\n=== Comparison Statistics ===\n');
    fprintf('Mean Absolute Error (MAE): %.6f\n', mae);
    fprintf('Root Mean Square Error (RMSE): %.6f\n', rmse);
    fprintf('Number of experimental points: %d\n', length(exp_energy));
    fprintf('Number of simulation points: %d\n', length(En));
end

% Also plot on the original figure for comparison
figure(10);
plot(En, eff, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Simulation');
hold on;
plot(exp_energy, exp_efficiency, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Experimental Data');
xlabel('photonEnergy,eV');
ylabel('Diffraction efficiency');
title(plot_title);
legend('Location', 'best');
grid on;
saveas(gcf, 'Example_SLAG_efficiency.png');

disp('Plot saved as Example_SLAG_efficiency.png');
