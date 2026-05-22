function ef = efficiency_bgrML(grPeriod_lpermm, ...
    grBA_deg, grAntiBA_deg, photonEnergy_eV, grazing_angle_deg, ...
    material_sub, material_HZ, material_LZ, ...
    ML_d_nm, ML_d_HZtod, ML_N, ...
    z_resolution_nm, x_resolution_nm, FourierOrders, plSelect)
% efficiency_bgrML - Compute diffraction efficiency using RETICOLO

fprintf('\n[efficiency_bgrML] Starting simulation at %s\n', datestr(now, 'HH:MM:SS'));
fprintf('   Photon energy: %.1f eV | Grating: %d l/mm | Grazing: %.3f°\n', ...
        photonEnergy_eV, grPeriod_lpermm, grazing_angle_deg);

base = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(base, '..', 'V7-reticolo-blazr', 'RETICOLO V7')));
retio

% Number of Fourier orders
nn = FourierOrders;
lambda_nm = 1239.8 / photonEnergy_eV;

fprintf(' -> Step 1/7: Loading refractive indices...\n');

% ------------------ Load refractive indices ------------------
for i = 1:3
    if i == 1
        nfile = ['n_', material_sub, '_cxro.txt'];
    elseif i == 2
        nfile = ['n_', material_HZ, '_cxro.txt'];
    else
        nfile = ['n_', material_LZ, '_cxro.txt'];
    end

    if exist(nfile, 'file') == 2
        nData = importdata(nfile);
        nData = nData.data;
    else
        error(['Index file does not exist: ', nfile]);
    end

    if i == 1
        n_sub_real = interp1(nData(:,1), nData(:,2), photonEnergy_eV);
        n_sub_imag = interp1(nData(:,1), nData(:,3), photonEnergy_eV);
        n_sub = 1 - n_sub_real + n_sub_imag * 1i;
    elseif i == 2
        n_HZ_real = interp1(nData(:,1), nData(:,2), photonEnergy_eV);
        n_HZ_imag = interp1(nData(:,1), nData(:,3), photonEnergy_eV);
        n_HZ = 1 - n_HZ_real + n_HZ_imag * 1i;
    else
        n_LZ_real = interp1(nData(:,1), nData(:,2), photonEnergy_eV);
        n_LZ_imag = interp1(nData(:,1), nData(:,3), photonEnergy_eV);
        n_LZ = 1 - n_LZ_real + n_LZ_imag * 1i;
    end
end
fprintf('    ✓ Index data loaded successfully.\n');

n_inc = 1;

% ------------------ Geometry setup ------------------
fprintf(' -> Step 2/7: Setting up geometry...\n');
p_nm = 1 / grPeriod_lpermm * 1e6;
theta0_deg = 90 - grazing_angle_deg;
k_parallel = n_inc * sind(theta0_deg);

th_nm = p_nm * tand(grBA_deg) * tand(grAntiBA_deg) / ...
         (tand(grBA_deg) + tand(grAntiBA_deg));
th_nm = th_nm + z_resolution_nm * 25 + ML_d_nm * ML_N;

Apex_z = p_nm * tand(grBA_deg) * tand(grAntiBA_deg) / ...
          (tand(grBA_deg) + tand(grAntiBA_deg));
Apex_x = Apex_z / tand(grBA_deg);
Prf = [0, 0; Apex_x, Apex_z; p_nm, 0];
fprintf('    ✓ Geometry defined (period = %.1f nm, thickness = %.1f nm)\n', p_nm, th_nm);

% ------------------ Grid setup ------------------
fprintf(' -> Step 3/7: Building mesh grid...\n');
x = linspace(0, p_nm, round(p_nm / x_resolution_nm) + 1);
z = linspace(th_nm, 0, round(th_nm / z_resolution_nm) + 1);
[X, Z] = meshgrid(x, z);
N_layers = length(z);
n = X .* 0;

Prf_z = interp1(Prf(:,1)', Prf(:,2)', x);
fprintf('    ✓ Mesh grid generated: %d x %d points.\n', size(X,1), size(X,2));

% ------------------ Multilayer construction ------------------
fprintf(' -> Step 4/7: Filling multilayer refractive index map...\n');
for i = 1:ML_N
    Prf0 = Prf_z + 1 + ML_d_nm * (i - 1);
    Prf1 = Prf_z + 1 + ML_d_nm * (i - 1) + ML_d_nm * ML_d_HZtod;
    Prf2 = Prf_z + 1 + ML_d_nm * i;

    if i == 1
        P = find(Z < Prf0);             n(P) = n_sub;
        P = find(Z >= Prf0 & Z < Prf1); n(P) = n_HZ;
        P = find(Z >= Prf1 & Z < Prf2); n(P) = n_LZ;
    elseif i == ML_N
        P = find(Z >= Prf2);            n(P) = n_inc;
        P = find(Z >= Prf0 & Z < Prf1); n(P) = n_HZ;
        P = find(Z >= Prf1 & Z < Prf2); n(P) = n_LZ;
    else
        P = find(Z >= Prf0 & Z < Prf1); n(P) = n_HZ;
        P = find(Z >= Prf1 & Z < Prf2); n(P) = n_LZ;
    end
end
fprintf('    ✓ Multilayer stack (%d periods) built.\n', ML_N);

% ------------------ Texture and profile ------------------
fprintf(' -> Step 5/7: Building RETICOLO texture and profile...\n');
deltan = diff(n, 1, 2);
[nonZeroRowIndices, nonZeroColIndices] = find(deltan ~= 0);
ProfEdge = sortrows([nonZeroRowIndices, nonZeroColIndices], 1);

textures = cell(1, N_layers + 2);
textures{1} = {n_inc};
textures{end} = {n_sub};

for layer = 1:N_layers
    p = find(ProfEdge(:,1) == layer);
    if ~isempty(p)
        x_position = X(layer, ProfEdge(p,2) + 1);
        n_value = n(layer, ProfEdge(p,2));
        textures{layer + 1} = {x_position, n_value};
    else
        n_value = n(layer,1);
        textures{layer + 1} = {n_value};
    end
end

texture_list = 1:N_layers + 2;
th_list = [0, ones(1, N_layers) .* z_resolution_nm, 0];
profile = {th_list, texture_list};
fprintf('    ✓ Texture/profile created.\n');

% ------------------ RETICOLO computation ------------------
fprintf(' -> Step 6/7: Running RETICOLO solver...\n');
pol = 1; % TE
parm = res0(pol);
parm.res1.trace = 0;

aa = res1(lambda_nm, p_nm, textures, nn, k_parallel, parm);
ef = res2(aa, profile, parm);
fprintf('    ✓ RETICOLO main computation completed.\n');

% ------------------ Optional plotting ------------------
% if plSelect == 1
%     fprintf(' -> Step 7/7: Computing field map for visualization...\n');
%     parm.res3.trace = 1;
%     parm.res3.cale = [];
%     parm.res3.npts = [10, 80, 10];

%     if pol == 1
%         einc = ef.inc_top.PlaneWave_E(2);
%     else
%         einc = ef.inc_top.PlaneWave_H(2);
%     end

%     [e, z, o] = res3(x, aa, profile, einc, parm);
%     axis square
%     set(gcf, 'WindowStyle', 'docked')
%     fprintf('    ✓ Field visualization ready.\n');
% end

fprintf('[efficiency_bgrML] Finished at %s\n', datestr(now, 'HH:MM:SS'));
fprintf('------------------------------------------------------------\n');

end
