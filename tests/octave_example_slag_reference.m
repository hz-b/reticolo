warning('off', 'all');
more off;

args = argv();
if numel(args) < 2
    error('Usage: octave_example_slag_reference.m <output_csv> <energy1> [energy2 ...]');
end

output_file = args{1};
energies = zeros(1, numel(args) - 1);
for idx = 2:numel(args)
    energies(idx - 1) = str2double(args{idx});
end

addpath(genpath(fullfile(pwd, 'V9', 'reticolo_allege_v9')));

metadata.grPeriod_lpermm = 400;
metadata.GR_Order = 1;
metadata.grWidthtoD = 0.67;
metadata.grDepth_nm = 14.9;
metadata.grTrapezoidAngL_deg = 15;
metadata.grTrapezoidAngR_deg = 15;
metadata.material_sub = 'Si';
metadata.material_layer = 'Pt';
metadata.layerThickness_nm = 28.77;
metadata.z_resolution_nm = 0.1;
metadata.x_resolution_nm = 0.1;
metadata.FourierOrders = 25;
grazing_angle_deg = 4;

fid = fopen(output_file, 'w');
if fid < 0
    error('Unable to open output file %s', output_file);
end
fprintf(fid, 'energy_ev,efficiency,diffraction_angle_deg\n');

for photonEnergy_eV = energies
    wavelength_nm = 1239.8 / photonEnergy_eV;

    nData = importdata(['n_', metadata.material_sub, '_cxro.txt']);
    nData = nData.data;
    n_sub_real = interp1(nData(:,1), nData(:,2), photonEnergy_eV);
    n_sub_imag = interp1(nData(:,1), nData(:,3), photonEnergy_eV);
    n_sub = 1 - n_sub_real + n_sub_imag * 1i;

    nData = importdata(['n_', metadata.material_layer, '_cxro.txt']);
    nData = nData.data;
    n_layer_real = interp1(nData(:,1), nData(:,2), photonEnergy_eV);
    n_layer_imag = interp1(nData(:,1), nData(:,3), photonEnergy_eV);
    n_layer = 1 - n_layer_real + n_layer_imag * 1i;

    n_inc = 1;
    period_nm = 1 / metadata.grPeriod_lpermm * 1e6;
    theta0_deg = 90 - grazing_angle_deg;
    k_parallel = n_inc * sin(theta0_deg * pi / 180);
    thickness_nm = metadata.grDepth_nm + metadata.layerThickness_nm + 5;

    pos1 = [(period_nm - metadata.grWidthtoD * period_nm) / 2 - metadata.grDepth_nm * tand(metadata.grTrapezoidAngL_deg), 0];
    pos2 = [(period_nm - metadata.grWidthtoD * period_nm) / 2, metadata.grDepth_nm];
    pos3 = [(period_nm + metadata.grWidthtoD * period_nm) / 2, metadata.grDepth_nm];
    pos4 = [(period_nm + metadata.grWidthtoD * period_nm) / 2 + metadata.grDepth_nm * tand(metadata.grTrapezoidAngR_deg), 0];
    prf = [0, 0; pos1; pos2; pos3; pos4; period_nm, 0];

    x = linspace(0, period_nm, round(period_nm / metadata.x_resolution_nm) + 1);
    z = linspace(thickness_nm, 0, round(thickness_nm / metadata.z_resolution_nm) + 1);
    [x_grid, z_grid] = meshgrid(x, z);
    surface = interp1(prf(:,1)', prf(:,2)', x);

    index_grid = x_grid .* 0;
    coating_top = surface + metadata.layerThickness_nm;
    p = find(z_grid < surface);
    index_grid(p) = n_sub;
    p = find(z_grid >= surface);
    index_grid(p) = n_inc;
    p = find(z_grid >= surface & z_grid < coating_top);
    index_grid(p) = n_layer;

    deltan = diff(index_grid, 1, 2);
    [row_idx, col_idx] = find(deltan ~= 0);
    edge_table = sortrows([row_idx, col_idx], 1);

    n_layers = length(z);
    textures = cell(1, n_layers + 2);
    textures{1} = {n_inc};
    textures{end} = {n_sub};

    for layer = 1:n_layers
        matches = find(edge_table(:,1) == layer);
        if isempty(matches)
            textures{layer + 1} = {index_grid(layer, 1)};
        else
            x_position = x_grid(layer, edge_table(matches, 2) + 1);
            n_value = index_grid(layer, edge_table(matches, 2));
            textures{layer + 1} = {x_position, n_value};
        end
    end

    profile = {[0, ones(1, n_layers) .* metadata.z_resolution_nm, 0], 1:(n_layers + 2)};

    parm = res0(1);
    aa = res1(wavelength_nm, period_nm, textures, metadata.FourierOrders, k_parallel, parm);
    ef = res2(aa, profile, parm);

    idx = find(ef.inc_top_reflected.order == -metadata.GR_Order);
    efficiency = ef.inc_top_reflected.efficiency(idx(1));
    diffraction_angle_deg = 90 - ef.inc_top_reflected.theta(idx(1));
    fprintf(fid, '%.6f,%.12f,%.12f\n', photonEnergy_eV, efficiency, diffraction_angle_deg);
end

fclose(fid);
