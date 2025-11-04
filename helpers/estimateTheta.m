function thetaEst = estimateTheta(material_HZ, material_LZ, ...
    ML_order, ML_d_nm, ML_d_HZtod, ...
    GR_order, GR_groove, grPeriod_lpermm, photonEnergy_eV)
% estimateTheta - Estimate the grazing angle (in degrees) for the multilayer grating

for i = 2:3
    if i == 2
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

    if i == 2
        n_HZ_real = interp1(nData(:,1), nData(:,2), photonEnergy_eV);
        n_HZ_imag = interp1(nData(:,1), nData(:,3), photonEnergy_eV);
    else
        n_LZ_real = interp1(nData(:,1), nData(:,2), photonEnergy_eV);
        n_LZ_imag = interp1(nData(:,1), nData(:,3), photonEnergy_eV);
    end
end

nAvg = n_LZ_real * (1 - ML_d_HZtod) + n_HZ_real * ML_d_HZtod;
ThetaBRAGG = asind(sqrt((ML_order * 1.2398 / (2 * ML_d_nm * photonEnergy_eV / 1000))^2 + 2 * nAvg));
LD = 1000 / grPeriod_lpermm * 1000;

thetaEst = ThetaBRAGG - asind(GR_order * GR_groove * 1.2398 / LD / ...
    (photonEnergy_eV / 1000) / 2 / sind(ThetaBRAGG));
end
