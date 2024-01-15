clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input grating structure
%
%     ideal or measured profile for:
%     blazed grating 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.grPeriod_lpermm = 600; %('simModel: input 1) grating period in l/mm: ')
metadata.GR_Order = 1; %('simModel:input 2)diffraction order: ')
metadata.grBA_deg = 1; %('simModel:input 3) grating blaze angle in deg: ')
metadata.grAntiBA_deg = 30; %('simModel:input 4) grating anti-blaze angle in deg: ')
metadata.material_sub = 'Si';%('simModel:input 5) grating substrate meaterial(type Si for silicon): ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     input layer structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata.material_layer = 'Au'; %input('layer material, type Au for gold: ');
metadata.layerThickness_nm = 30; %input('single layer thickness in nm: ');

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

metadata.photonEnergy_eV = 100:10:1000; %input(' which photon energy(range) wants to comput_in eV: ');:
grazing_angle_deg = 2;% input('desired grazing incidence angle range) in deg: ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    build model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = metadata;
addpath('MODIFY HERE!!! PATH to here \V9\reticolo_allege_v9')
retio
eff=[];
En=[];
for photonEnergy_eV=metadata.photonEnergy_eV
    a.photonEnergy_eV=photonEnergy_eV;
%number of Fourier orders
nn = a.FourierOrders;
%input material file
for i = 1:2
    if i == 1
        nfile = ['n_',a.material_sub,'_cxro.txt'];
    elseif i == 2
        nfile = ['n_',a.material_layer,'_cxro.txt'];
    end
    
    if exist(nfile, 'file') ==  2
        nData = importdata(nfile);
        nData = nData.data;
    else
        disp('index File does not exist./type wrong material');
    end
    
    %wavelength
    lambda_nm = 1239.8/a.photonEnergy_eV;
    
    %index
    if i == 1
        n_sub_real = interp1(nData(:,1),nData(:,2),a.photonEnergy_eV);
        n_sub_imag = interp1(nData(:,1),nData(:,3),a.photonEnergy_eV);
        n_sub = 1-n_sub_real+n_sub_imag*1i;
    elseif i == 2
        n_HZ_real = interp1(nData(:,1),nData(:,2),a.photonEnergy_eV);
        n_HZ_imag = interp1(nData(:,1),nData(:,3),a.photonEnergy_eV);
        n_HZ = 1-n_HZ_real+n_HZ_imag*1i;
    end
end
n_inc = 1;

%grating pitch
p_nm = 1/a.grPeriod_lpermm*1E6;

%angle of incidence
theta0_deg = 90-grazing_angle_deg;
k_parallel = n_inc*sin(theta0_deg*pi/180);


% thicknes mean grooveHeight
th_nm = p_nm*tand(a.grBA_deg)*tand(a.grAntiBA_deg)/(tand(a.grBA_deg)+tand(a.grAntiBA_deg));
th_nm = th_nm+a.layerThickness_nm+5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate blazed grating profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Apex_z=p_nm*tand(a.grBA_deg)*tand(a.grAntiBA_deg)/(tand(a.grBA_deg)+tand(a.grAntiBA_deg));
Apex_x=Apex_z/tand(a.grBA_deg);
Prf=[0,0;Apex_x,Apex_z;p_nm,0];
groove=Prf(2:end,:);

a.GR_groove = 1;
x = linspace(0,p_nm*a.GR_groove,round(p_nm*a.GR_groove/a.x_resolution_nm)+1);%
z = linspace(th_nm,0,round(th_nm/a.z_resolution_nm)+1);%

% genergy grid
[X,Z] = meshgrid(x,z);

% number of z layers
N_layers = length(z);

% mat to store n
n = X.*0;

% interplot the z of profilr based on x
Prf_z=interp1(Prf(:,1)',Prf(:,2)',x);

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
disp(['En',num2str(a.photonEnergy_eV),'eV, diffraction efficiency at the first order:', num2str(round(Output(i,1),3)*100),'%,  Diffraction angle : ',num2str(Output(i,2)),'deg'])
eff=[eff,Output(i,1)];
En=[En,a.photonEnergy_eV];
i=i+1;
figure(10)
plot(En,eff,'*'); xlabel('photonEnergy,eV'),ylabel('Diffraction efficiency');
end

select = input('wants to check the model refracx distribution? 1 = yes, 2 = No: ');
if select == 1
    parm.res3.trace = 1 ; % trace automatique
    parm.res3.cale = [];
    parm.res3.npts = [10,80,10];
    
    if pol == 1 % 1:TE   -1:TM
        einc =  ef.inc_top.PlaneWave_E(2);
    else
        einc =  ef.inc_top.PlaneWave_H(2);
    end
    % x x section; aa from res1, profile, 1,
    [e,z,o] = res3(x,aa,profile,einc,parm);
    axis square
    set(gcf,'WindowStyle','docked')
end

