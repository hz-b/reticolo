clear;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT GRATING STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.grPeriod_lpermm = 2400;      % grooves per mm
metadata.GR_Order = 2;               % diffraction order to analyze
metadata.grBA_deg = 0.9;               % blaze angle
metadata.grAntiBA_deg = 23.6;          % anti-blaze angle
metadata.material_sub = 'Si';        % substrate material

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT MULTILAYER STRUCTURE (TEST: Au/Au)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.material_layers = {'Cr', 'C'};   % alternating materials
metadata.layerThicknesses_nm = [2.3, 2.7];     % thickness per layer
metadata.N_periods = 40;                    % number of bilayers (5×(5+5)=50 nm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.z_resolution_nm = .1;
metadata.x_resolution_nm = .1;
metadata.FourierOrders = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INCIDENCE AND PHOTON ENERGY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.photonEnergy_eV = 500:4500:5000;  % energy sweep for the test
grazing_angle_deg = 1.5;                    % grazing incidence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = metadata;
addpath('MODIFY HERE!!! PATH to \V9\reticolo_allege_v9')
retio
eff = [];
En = [];
i = 1;

for photonEnergy_eV = metadata.photonEnergy_eV
    a.photonEnergy_eV = photonEnergy_eV;
    nn = a.FourierOrders;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD REFRACTIVE INDICES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_inc = 1; % vacuum
    mats = unique([{a.material_sub}, metadata.material_layers]); % list all used materials
    n_map = struct();
    for j = 1:length(mats)
        nfile = ['n_', mats{j}, '_cxro.txt'];
        if exist(nfile, 'file') ~= 2
            error(['Index file ', nfile, ' not found.']);
        end
        nData = importdata(nfile).data;
        n_real = interp1(nData(:,1), nData(:,2), a.photonEnergy_eV);
        n_imag = interp1(nData(:,1), nData(:,3), a.photonEnergy_eV);
        n_map.(mats{j}) = 1 - n_real + 1i * n_imag;
    end

    n_sub = n_map.(a.material_sub);
    n_layers = cellfun(@(m) n_map.(m), metadata.material_layers, 'UniformOutput', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GEOMETRY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_nm = 1/a.grPeriod_lpermm * 1e6;
    theta0_deg = 90 - grazing_angle_deg;
    k_parallel = n_inc * sin(theta0_deg*pi/180);

    % Grating groove height
    th_nm = p_nm * tand(a.grBA_deg) * tand(a.grAntiBA_deg) / (tand(a.grBA_deg) + tand(a.grAntiBA_deg));
    th_nm = th_nm + sum(metadata.layerThicknesses_nm) * metadata.N_periods + 5;

    % Grating profile
    Apex_z = p_nm * tand(a.grBA_deg) * tand(a.grAntiBA_deg) / (tand(a.grBA_deg) + tand(a.grAntiBA_deg));
    Apex_x = Apex_z / tand(a.grBA_deg);
    Prf = [0,0; Apex_x, Apex_z; p_nm,0];

    a.GR_groove = 1;
    x = linspace(0, p_nm * a.GR_groove, round(p_nm * a.GR_groove / a.x_resolution_nm) + 1);
    z = linspace(th_nm, 0, round(th_nm / a.z_resolution_nm) + 1);
    [X, Z] = meshgrid(x, z);
    N_layers = length(z);
    n = X .* 0;

    % Profile height function
    Prf_z = interp1(Prf(:,1)', Prf(:,2)', x);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MULTILAYER ASSIGNMENT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Prf0 = Prf_z;
    z_current = Prf0;

    for p = 1:metadata.N_periods
        for j = 1:length(metadata.material_layers)
            t = metadata.layerThicknesses_nm(j);
            n_val = n_layers{j};
            z_next = z_current + t;
            P = find(Z >= z_current & Z < z_next);
            n(P) = n_val;
            z_current = z_next;
        end
    end

    % Assign substrate and vacuum
    P = find(Z < Prf0);  n(P) = n_sub;
    P = find(Z >= z_current);  n(P) = n_inc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RETICOLO TEXTURE CONSTRUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deltan = diff(n,1,2);
    [nonZeroRowIndices, nonZeroColIndices] = find(deltan ~= 0);
    ProfEdge = sortrows([nonZeroRowIndices, nonZeroColIndices], 1);

    textures = cell(1, N_layers + 2);
    textures{1} = {n_inc};
    textures{end} = {n_sub};

    for layer = 1:N_layers
        if ~isempty(find(ProfEdge(:,1) == layer, 1))
            pidx = find(ProfEdge(:,1) == layer);
            x_position = X(layer, ProfEdge(pidx,2) + 1);
            n_value = n(layer, ProfEdge(pidx,2));
            textures{layer+1} = {x_position, n_value};
        else
            n_value = n(layer,1);
            textures{layer+1} = {n_value};
        end
    end

    texture_list = 1:N_layers + 2;
    th_list = [0, ones(1, N_layers) * a.z_resolution_nm, 0];
    profile = {th_list, texture_list};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RCWA SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pol = 1; % TE
    parm = res0(pol);
    parm.res1.trace = 0;
    aa = res1(1239.8/a.photonEnergy_eV, p_nm*a.GR_groove, textures, nn, k_parallel, parm);
    ef = res2(aa, profile, parm);

    idx_DesignOrder = find(a.FourierOrders:-1:0 == a.GR_Order*a.GR_groove);
    Output(i,1) = ef.inc_top_reflected.efficiency(idx_DesignOrder);
    Output(i,2) = 90 - ef.inc_top_reflected.theta(idx_DesignOrder);

    disp(['En', num2str(a.photonEnergy_eV), 'eV, diffraction efficiency at the first order: ', ...
          num2str(round(Output(i,1)*1000)/1000*100), '%,  Diffraction angle : ', ...
          num2str(Output(i,2)), 'deg']);

    eff = [eff, Output(i,1)];
    En = [En, a.photonEnergy_eV];
    i = i + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
plot(En, eff, '*-');
xlabel('Photon energy (eV)');
ylabel('Diffraction efficiency');
title('Test multilayer grating (Au/Au)');
grid on;
