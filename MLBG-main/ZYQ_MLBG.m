clear;
grPeriod_lpermm=150;
grazing_angle_deg=2;

grBA_deg=1;
grAntiBA_deg=89;
photonEnergy_eV=380;

material='Au';
z_resolution_nm=0.5;
x_resolution_nm=10;
FourierOrders=10;
plSelect=0;

eff=[];
energy_=[];
effAll=[];
for FourierOrders=5:5:40
for energy=50:10:250
    photonEnergy_eV= energy;
    eta = efficiency_bgrSL(grPeriod_lpermm, ...
        grBA_deg, grAntiBA_deg, photonEnergy_eV, grazing_angle_deg, material, z_resolution_nm,x_resolution_nm, FourierOrders,plSelect);
    eff=[eff,eta];
%     energy_=[energy_,energy];
    
end
    plot(50:10:250,eff,'-','DisplayName',['number of Harmonic: ',num2str(FourierOrders)]);xlabel('energy');ylabel('1 order efficiency''ordre 0');grid;title('TEST RETICOLO');pause(eps);   
    hold on
    effAll=[effAll;eff];
    eff=[];
    energy_=[];
end
legend show
effAll=[50:10:250;effAll];
save('bgrCmp_effAll.mat','effAll');
function eta = efficiency_bgrSL(grPeriod_lpermm, ...
    grBA_deg, grAntiBA_deg, photonEnergy_eV, grazing_angle_deg, material, z_resolution_nm,x_resolution_nm, FourierOrders,plSelect)


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
th_nm = p_nm*tand(grBA_deg)*tand(grAntiBA_deg)/(tand(grBA_deg)+tand(grAntiBA_deg));


% polarisation
pol = -1; % 1:TE   -1:TM
parm = res0(pol);  % parameter init


% generate texture and Profile
Apex_z=p_nm*tand(grBA_deg)*tand(grAntiBA_deg)/(tand(grBA_deg)+tand(grAntiBA_deg));
Apex_x=Apex_z/tand(grBA_deg);
Prf=[0,0;Apex_x,Apex_z;p_nm,0];

x=linspace(0,p_nm,round(p_nm/x_resolution_nm)+1);%
z=linspace(th_nm,0,round(th_nm/z_resolution_nm)+1);%

% genergy grid 
[X,Z]=meshgrid(x,z);

% number of blazed layers
N_layers=length(z);

% mat to store n 
n=X.*0;

% interplot the z of profilr based on x
Prf_z=interp1(Prf(:,1)',Prf(:,2)',x);

P=find(Z<Prf_z);         n(P)=n_gr; % gr_substrate
P = find(Z>=Prf_z);        n(P)=n_inc; % inc_medium

deltan = diff(n);
P=find(deltan~= 0);
rowdist=ones(1,N_layers);
textures_tmp= [mat2cell(X,rowdist) mat2cell(n,rowdist)];

textures = cell(1,N_layers+2);
textures{1}={n_inc};
textures{end}={n_sub};
for i = 1:length(z)
    textures{i+1} = {textures_tmp{i,1}, textures_tmp{i,2}};
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