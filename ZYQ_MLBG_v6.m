clear;

% Add the helper folder to the MATLAB/Octave path
addpath(genpath(fullfile(pwd, 'helpers')))

% ------------------ METADATA SETUP ------------------
metadata.grPeriod_lpermm = 2400;
metadata.grBA_deg = 1;
metadata.grAntiBA_deg = 3;

metadata.material_sub = 'Si';
metadata.material_HZ  = 'Cr';
metadata.material_LZ  = 'C';

metadata.ML_d_nm     = 6.5;
metadata.ML_d_HZtod  = 0.45;
metadata.ML_N        = 60;

metadata.z_resolution_nm = 1 %0.01;
metadata.x_resolution_nm = 1;
metadata.FourierOrders   = 5; % or FourierOrders = GR_order*GR_groove+3
metadata.plSelect        = 1;

% Grating parameters
metadata.ML_order      = 1;
metadata.GR_groove     = 1;
metadata.GR_order      = 1;
metadata.photonEnergy_eV = 2500;

% ------------------ LOOP DETECTION ------------------
fields = fieldnames(metadata);
disp('Parameters with more than one value:');
loopNumber = 1;

for i = 1:numel(fields)
    field = fields{i};
    if ~ischar(metadata.(field)) && numel(metadata.(field)) > 1
        loop.para{loopNumber}  = field;
        loop.Idx{loopNumber}   = i;
        loop.value{loopNumber} = metadata.(field);
        disp(['Field: ' field]);
        disp(['Values: ' num2str(metadata.(field))]);
        loopNumber = loopNumber + 1;
    end
end

loopNumber = loopNumber - 1;

% Define loop bounds
loopBounds = zeros(1, loopNumber);
for i = 1:loopNumber
    IdxCurrent = loop.Idx{loopNumber};
    loopBounds(i) = length(metadata.(fields{IdxCurrent}));
end

% Output data containers
outputdata0 = cell(1, prod(loopBounds, 2));
numLoops    = length(loopBounds);
outputdata1 = zeros(prod(loopBounds, 2), numLoops + 3);
loopIndices = cell(1, numLoops);

for i = 1:numLoops
    loopIndices{i} = 1;
end

k = 1;
idx_outputdata0 = 1;
tic;

% ------------------ MAIN COMPUTATION ------------------
if numLoops == 0
    disp('No loop indices, single-point calculation');
    metadataCurrent = metadata;

    thetaEst = estimateTheta(metadataCurrent.material_HZ, metadataCurrent.material_LZ, ...
        metadataCurrent.ML_order, metadataCurrent.ML_d_nm, metadataCurrent.ML_d_HZtod, ...
        metadataCurrent.GR_order, metadataCurrent.GR_groove, ...
        metadataCurrent.grPeriod_lpermm, metadataCurrent.photonEnergy_eV);

    grazing_angle_deg = thetaEst;

    ef = efficiency_bgrML(metadataCurrent.grPeriod_lpermm, ...
        metadataCurrent.grBA_deg, metadataCurrent.grAntiBA_deg, metadataCurrent.photonEnergy_eV, ...
        grazing_angle_deg, metadataCurrent.material_sub, metadataCurrent.material_HZ, metadataCurrent.material_LZ, ...
        metadataCurrent.ML_d_nm, metadataCurrent.ML_d_HZtod, metadataCurrent.ML_N, ...
        metadataCurrent.z_resolution_nm, metadataCurrent.x_resolution_nm, ...
        metadataCurrent.FourierOrders, metadataCurrent.plSelect);

    idx_DesignOrder = find(metadataCurrent.FourierOrders:-1:0 == metadataCurrent.GR_order);
    eff   = ef.inc_top_reflected.efficiency(idx_DesignOrder);
    theta = grazing_angle_deg;

    plot(theta, eff);
    xlabel('theta (deg)');
    ylabel('-1st order diffraction efficiency');
    title('Diffraction efficiency');
    legend('TE');
    grid on;
    drawnow;
    disp('Plot displayed — press any key to continue.');
    pause;             % Wait for user input before closing
else
    disp('Multi-parameter scan not implemented in this clean version.');
end

elapsedTime = toc;
disp(['Total simulation time: ' datestr(seconds(elapsedTime), 'HH:MM:SS')]);
