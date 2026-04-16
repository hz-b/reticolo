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
metadata.grTrapezoidAngL_deg = 15; %('simModel:input 6.1)grTrapezoid left: ')
metadata.grTrapezoidAngR_deg = 15; %('simModel:input 6.2)grTrapezoid right: ')
metadata.material_sub = 'Pt';%('simModel:input 7) grating substrate meaterial(type Si for silicon): ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input layer structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.useTopLayer = false; % set true to add a coating layer on top of the Pt grating
metadata.material_layer = 'Au'; % coating material, only used when useTopLayer is true
metadata.layerThickness_nm = .1; % coating thickness in nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input simulation parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.z_resolution_nm = 1 ;%input('z_slicing_step_in_nm: ');
metadata.x_resolution_nm = 1 ;%input('x_slicing_step_in_nm: ');
metadata.FourierOrders = 5 ;%input('Harmonics: '); %or FourierOrders = GR_Order*GR_groove+3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input incidence angle and photon energy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.photonEnergy_eV = 100:10:2000; %input(' which photon energy(range) wants to comput_in eV: ');:
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
energyIdx = 1;
for photonEnergy_eV=metadata.photonEnergy_eV
    a.photonEnergy_eV=photonEnergy_eV;
%number of Fourier orders
nn = a.FourierOrders;

% wavelength
lambda_nm = 1239.8/a.photonEnergy_eV;

% substrate optical constants
subFile = ['n_', a.material_sub, '_cxro.txt'];
if exist(subFile, 'file') ~= 2
    error('Index file does not exist: %s', subFile);
end
subData = importdata(subFile);
subData = subData.data;
n_sub_real = interp1(subData(:,1), subData(:,2), a.photonEnergy_eV);
n_sub_imag = interp1(subData(:,1), subData(:,3), a.photonEnergy_eV);
n_sub = 1 - n_sub_real + n_sub_imag*1i;

% top region is either vacuum or a coating layer
n_inc = 1;
if a.useTopLayer
    layerFile = ['n_', a.material_layer, '_cxro.txt'];
    if exist(layerFile, 'file') ~= 2
        error('Index file does not exist: %s', layerFile);
    end
    layerData = importdata(layerFile);
    layerData = layerData.data;
    n_HZ_real = interp1(layerData(:,1), layerData(:,2), a.photonEnergy_eV);
    n_HZ_imag = interp1(layerData(:,1), layerData(:,3), a.photonEnergy_eV);
    n_HZ = 1 - n_HZ_real + n_HZ_imag*1i;
    layerThickness_nm = a.layerThickness_nm;
else
    n_HZ = n_inc;
    layerThickness_nm = 0;
end

%grating pitch
p_nm = 1/a.grPeriod_lpermm*1E6;

%angle of incidence
theta0_deg = 90-grazing_angle_deg;
k_parallel = n_inc*sin(theta0_deg*pi/180); 


% thicknes mean grooveHeight
th_nm = a.grDepth_nm + layerThickness_nm + 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lamellar grating profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1/tand(a.grTrapezoidAngL_deg) && 1/tand(a.grTrapezoidAngR_deg) == 0 % ideal case
    edge_z = a.grDepth_nm;
    edge_x = a.grWidthtoD*p_nm;
else
    pos1 = [(p_nm-a.grWidthtoD*p_nm)/2-a.grDepth_nm*tand(a.grTrapezoidAngL_deg),0]; %the gomentry was changed to accomodate the vertical wall reference for angles
    pos2 = [(p_nm-a.grWidthtoD*p_nm)/2,a.grDepth_nm];
    pos3 = [(p_nm+a.grWidthtoD*p_nm)/2,a.grDepth_nm];
    pos4 = [(p_nm+a.grWidthtoD*p_nm)/2+a.grDepth_nm*tand(a.grTrapezoidAngR_deg),0];
    Prf = [0,0;pos1;pos2;pos3;pos4;p_nm,0];
end

a.GR_groove = 1;
x = linspace(0,p_nm*a.GR_groove,round(p_nm*a.GR_groove/a.x_resolution_nm)+1);
z = linspace(th_nm,0,round(th_nm/a.z_resolution_nm)+1);%

% genergy grid
[X,Z] = meshgrid(x,z);

% number of z layers
N_layers = length(z);

% mat to store n
n = X.*0;

%I am considering using 
%if a.grTrapezoidAngL_deg == 90 && a.grTrapezoidAngR_deg == 90
%the block below is building a height profile based on the prf previosly built

if a.grTrapezoidAngL_deg == 0 && 1/a.grTrapezoidAngR_deg == 0
    P = find(x<= p_nm/2-edge_x/2 | x>= p_nm/2+edge_x/2);
    Prf_z(P) = 0;
    P = [];
    P = find(x>p_nm/2-edge_x/2 & x<= p_nm/2+edge_x/2);
    Prf_z(P) = edge_z;
else
    Prf_z = interp1(Prf(:,1)',Prf(:,2)',x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate texture for single layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prf0 = Prf_z;
Prf1 = Prf_z + layerThickness_nm;
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
Output(energyIdx,1)=ef.inc_top_reflected.efficiency(idx_DesignOrder);
Output(energyIdx,2)=90-ef.inc_top_reflected.theta(idx_DesignOrder);
eff_pct = sprintf('%.3f', Output(energyIdx,1) * 100);
disp(['En',num2str(a.photonEnergy_eV),'eV, diffraction efficiency at the first order: ', eff_pct, '%,  Diffraction angle : ',num2str(Output(energyIdx,2)),'deg'])
eff=[eff,Output(energyIdx,1)];
En=[En,a.photonEnergy_eV];
energyIdx = energyIdx + 1;
figure(10)
clf(10);
plot(En,eff,'*-');
xlabel('photonEnergy,eV');
ylabel('Diffraction efficiency');
grid on;
drawnow;
saveas(gcf, 'Example_SLAG_efficiency.png');
end


select = 1;
 
if select == 1
    parm.res3.trace = 1;
    parm.res3.cale  = [];
    parm.res3.npts  = [10, 80, 10];
 
    if pol == 1
        einc = ef.inc_top.PlaneWave_E(2);
    else
        einc = ef.inc_top.PlaneWave_H(2);
    end
 
    [e, z_field, o] = res3(x, aa, profile, einc, parm);
end
