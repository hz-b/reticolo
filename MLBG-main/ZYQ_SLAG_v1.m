clear;
grPeriod_lpermm=150;
grazing_angle_deg=2;

grDepth_nm=60;
grRatio=0.35;
photonEnergy_eV=380;

material='Au';
z_resolution_nm=0.5;
x_resolution_nm=10;
FourierOrders=20;
plSelect=0;

eff=[];
energy_=[];
effAll=[];


% eta = efficiency_bgrSL(grPeriod_lpermm, ...
%     grDepth_nm, grRatio, photonEnergy_eV, grazing_angle_deg, material, z_resolution_nm,x_resolution_nm, FourierOrders,plSelect);

for grPeriod_lpermm=150 %[1000,500,150:-20:20]
    effAll=[];
    figure
    for FourierOrders=5:5:40
        for energy=50:10:1000
            photonEnergy_eV= energy;
            eta = efficiency_bgrSL(grPeriod_lpermm, ...
                grDepth_nm, grRatio, photonEnergy_eV, grazing_angle_deg, material, z_resolution_nm,x_resolution_nm, FourierOrders,plSelect);
            eff=[eff,eta];
            %     energy_=[energy_,energy];
            
        end
        plot(50:10:1000,eff,'-','DisplayName',['number of Harmonic: ',num2str(FourierOrders)]);xlabel('energy');ylabel('1 order efficiency''ordre 0');grid;title('TEST RETICOLO');pause(eps);
        hold on
        effAll=[effAll;eff];
        eff=[];
        energy_=[];
    end
    legend show
    effAll=[50:10:1000;effAll];
    save(['agrCmp_effAll_v1_',num2str(grPeriod_lpermm),'LinePerMm_Depth15.mat'],'effAll');
end

function eta = efficiency_bgrSL(grPeriod_lpermm, ...
    grDepth_nm, grRatio, photonEnergy_eV, grazing_angle_deg, material, z_resolution_nm,x_resolution_nm, FourierOrders,plSelect)


addpath(genpath('RETICOLO V7'))
retio

%number of Fourier orders
nn = FourierOrders;

%material 
nfile=['n_',material,'_cxro.txt'];
if exist(nfile, 'file') == 2
    nData=importdata(nfile);
    nData=nData.data;
else
    disp('index File does not exist./type wrong material');
end

%wavelength
lambda_nm = 1239.8/photonEnergy_eV;

%index
n_inc = 1;
n_gr_real=interp1(nData(:,1),nData(:,2),photonEnergy_eV);
n_gr_imag=interp1(nData(:,1),nData(:,3),photonEnergy_eV);
n_gr = 1-n_gr_real+n_gr_imag*1i;
n_sub = n_gr;

%grating pitch
p_nm = 1/grPeriod_lpermm*1E6;

%angle of incidence
theta0_deg = 90-grazing_angle_deg;
k_parallel = n_inc*sin(theta0_deg*pi/180);

% thicknes mean grooveHeight
th_nm = grDepth_nm;


% polarisation
pol = -1; % 1:TE   -1:TM
parm = res0(pol);  % parameter init


% generate texture and Profile
edge_z=th_nm;
edge_x=grRatio*p_nm;
%Prf=[0,0;edge_x,edge_z;p_nm,0];
%Prf(:,1)=Prf(:,1)-p_nm/2;

x=linspace(0,p_nm,round(p_nm/x_resolution_nm)+1);%
z=linspace(th_nm,0,round(th_nm/z_resolution_nm)+1);%

% genergy grid 
[X,Z]=meshgrid(x,z);

% number of blazed layers
N_layers=length(z);

% mat to store n 
n=X.*0;

% interplot the z of profilr based on x
P=find(x<=p_nm/2-edge_x/2 | x>=p_nm/2+edge_x/2);
Prf_z(P)=0;
P=[];
P=find(x>p_nm/2-edge_x/2 & x<=p_nm/2+edge_x/2);
Prf_z(P)=edge_z;


Idx_P=find(Z<Prf_z);         n(Idx_P)=n_gr; % gr_substrate
Idx_PP = find(Z>=Prf_z);        n(Idx_PP)=n_inc; % inc_medium

deltan = diff(n,1,2);
[nonZeroRowIndices, nonZeroColIndices]=find(deltan~= 0);
ProfEdge=[nonZeroRowIndices,nonZeroColIndices];
ProfEdge=sortrows(ProfEdge, 1);

textures = cell(1,N_layers+2);
textures{1}={n_inc};
textures{end}={n_sub};

for layer=1:N_layers
    if ~isnan(find(ProfEdge(:,1)==layer))
        p=find(ProfEdge(:,1)==layer);
        x_position=X(layer,ProfEdge(p,2)+1);
        n_value=n(layer,ProfEdge(p,2));     
        textures{layer+1} = {x_position, n_value};
    else
        n_value=n(layer,1);
        textures{layer+1} = {n_value};
    end
end


texture_list=1:N_layers+2;
th_list=[0,ones(1,N_layers).*z_resolution_nm,0];
profile={th_list,texture_list};% starting1, from incidence medium

% initialisation
parm.res1.trace= 0;
aa = res1(lambda_nm, p_nm, textures, nn, k_parallel, parm);

ef = res2(aa, profile, parm);
% grating efficiency
eta = ef.inc_top_reflected.efficiency{-1};

if plSelect==1
    %             show if needed
    % % xx=linspace(-p_nm,p_nm,length(x)*2);% on trace 2 periodes
    parm.res3.trace=1 ; % trace automatique
    parm.res3.cale=[];
    parm.res3.npts=[10,80,10];
    
    if pol==1 % 1:TE   -1:TM
        einc= ef.inc_top.PlaneWave_E(2);
    else
        einc= ef.inc_top.PlaneWave_H(2);
    end
    % x x section; aa from res1, profile, 1,
    [e,z,o]=res3(x,aa,profile,einc,parm);
    axis square
    set(gcf,'WindowStyle','docked')
else
end


end